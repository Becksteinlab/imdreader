from MDAnalysis.coordinates.base import (
    ReaderBase,
    FrameIteratorIndices,
    FrameIteratorAll,
    FrameIteratorSliced,
)
from MDAnalysis.lib.util import store_init_arguments
from .IMDProtocol import *
from .util import *
import socket
import threading
import numpy as np
import signal
import logging
import time


class IMDReader(ReaderBase):
    """
    Reader for IMD protocol packets.

    Default buffer size is set to 8MB for testing
    Buffer_size kwarg is in bytes

    We are assuming the header will never be sent without the body as in the sample code.
    If this assumption is violated, the producer thread can cause a deadlock.
    """

    format = "IMD"

    @store_init_arguments
    def __init__(
        self,
        filename,
        convert_units=True,
        num_atoms=None,
        buffer_size=2**26,
        socket_bufsize=None,
        **kwargs,
    ):
        super(IMDReader, self).__init__(filename, **kwargs)

        # NOTE: Replace me with header packet which contains this information
        # OR get this information from the topology?
        if not num_atoms:
            raise ValueError("num_atoms must be specified")
        self.n_atoms = num_atoms
        self.ts = self._Timestep(
            self.n_atoms, positions=False, **self._ts_kwargs
        )

        self.units = {
            "time": "ps",
            "length": "nm",
            "force": "kJ/(mol*nm)",
        }

        self._buffer = CircularByteBuf(buffer_size, self.n_atoms, self.ts)

        self._producer = IMDProducer(
            filename, self._buffer, self._n_frames, self.n_atoms, socket_bufsize
        )

        self._frame = -1

    def _read_next_timestep(self):
        if self._frame == -1:
            self._producer.start()

        return self._read_frame(self._frame + 1)

    def _read_frame(self, frame):
        # When we've read all frames,
        # we wait for the producer to finish
        # _full is safe to access here since producer can no longer change it
        if self._buffer.producer_finished and self._buffer._full == 0:
            self._producer.join()
            raise IOError from None

        self._frame = frame

        # loads the timestep with step, positions, and energy
        self._buffer.consume_next_timestep()
        self.ts.frame = self._frame
        self.ts.dimensions = None

        return self.ts

    @property
    def n_frames(self):
        """number of frames in trajectory"""
        raise RuntimeError("IMDReader: n_frames is not unknown for a stream")

    def rewind(self):
        pass

    @staticmethod
    def _format_hint(thing):
        try:
            parse_host_port(thing)
        except:
            return False
        return True

    def close(self):
        """Gracefully shut down the reader."""
        self._producer.stop()
        print("IMDReader shut down gracefully.")

    # Incompatible methods
    def copy(self):
        raise NotImplementedError("IMDReader does not support copying")

    def _reopen(self):
        pass

    def __getitem__(self, frame):
        """This method from ProtoReader must be overridden
        to prevent slicing that doesn't make sense in a stream.
        """
        raise RuntimeError("IMDReader: Trajectory can only be read in for loop")


class IMDProducer(threading.Thread):
    """
    Producer thread for IMDReader. Reads packets from the socket
    and places them into the shared buffer.
    """

    def __init__(
        self,
        filename,
        buffer,
        n_atoms,
        socket_bufsize=None,
        pausable=True,
    ):
        super(IMDProducer, self).__init__(daemon=True)
        self._host, self._port = parse_host_port(filename)
        self._conn = None
        self.running = False

        self._buffer = buffer
        self.n_atoms = n_atoms
        self._expected_data_bytes = 12 * n_atoms
        self._socket_bufsize = socket_bufsize
        self.pausable = pausable
        # < represents little endian and > represents big endian
        # we assume big by default and use handshake to check
        self.parsed_frames = 0
        self._full_frames = 0
        self._parse_frame_time = None
        # Saving memory by preallocating space for the frame
        # we're loading into the buffer
        self._energy_byte_buf = bytearray(IMDENERGYPACKETLENGTH)
        self._energy_byte_view = memoryview(self._energy_byte_buf)
        self._body_byte_buf = bytearray(self._expected_data_bytes)
        self._body_byte_view = memoryview(self._body_byte_buf)

        signal.signal(signal.SIGINT, self._handle_signal)
        signal.signal(signal.SIGTERM, self._handle_signal)
        self._is_disconnected = False

        # The body of a force or position packet should contain
        # (4 bytes per float * 3 atoms * n_atoms) bytes
        self._expected_data_bytes = 12 * self.n_atoms

    def _await_IMD_handshake(self):
        """
        Wait for the server to send a handshake packet, set endianness,
        and check IMD Protocol version.
        """
        print("waiting for handshake...")
        handshake = self._expect_header(expected_type=IMDType.IMD_HANDSHAKE)
        if handshake.length != IMDVERSION:
            # Try swapping endianness
            swapped = struct.unpack("<i", struct.pack(">i", handshake.length))[
                0
            ]
            if swapped != IMDVERSION:
                err_version = min(swapped, handshake.length)
                raise ValueError(
                    f"Incompatible IMD version. Expected {IMDVERSION}, got {err_version}"
                )
            else:
                self._buffer.inform_endianness("<")
        print("handshake received")

    def _send_go_packet(self):
        """
        Send a go packet to the client to start the simulation
        and begin receiving data.
        """
        print("sending go packet...")
        go = create_header_bytes(IMDType.IMD_GO, 0)
        self._conn.sendall(go)

    def _toggle_pause_simulation(self):
        """
        Block the simulation until the buffer has more space.
        """
        pause = create_header_bytes(IMDType.IMD_PAUSE, 0)
        self._conn.sendall(pause)
        self._buffer.wait_almost_empty()
        pause = create_header_bytes(IMDType.IMD_PAUSE, 0)
        self._conn.sendall(pause)

    def _connection_sequence(self):
        """
        Establish connection with the server and perform
        the handshake and go packet exchange.
        """
        self._conn = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        if self._socket_bufsize:
            self._conn.setsockopt(
                socket.SOL_SOCKET, socket.SO_RCVBUF, self._socket_bufsize
            )
        self._conn.connect((self._host, self._port))
        self._await_IMD_handshake()
        self._send_go_packet()
        self.running = True

    def run(self):
        """
        Producer thread method. Reads from the socket and
        sends a 'pause' signal if needed.
        """
        self._connection_sequence()

        while (self.parsed_frames) and self.running:
            try:
                with timeit() as parse_frame:
                    self._expect_header(
                        expected_type=IMDType.IMD_ENERGIES, expected_value=1
                    )
                    energies = self._recv_energy_bytes()
                    self._expect_header(
                        expected_type=IMDType.IMD_FCOORDS,
                        expected_value=self.n_atoms,
                    )
                    pos = self._recv_body_bytes()
                    self._full_frames = self._buffer.insert(energies, pos)
                    self.parsed_frames += 1

                self._parse_frame_time = parse_frame.elapsed
                logging.debug(
                    f"IMDProducer: Added frame #{self.parsed_frames - 1} to buffer in {parse_frame.elapsed} seconds"
                )
            except Exception as e:
                logging.debug(
                    f"IMDProducer: Error parsing frame #{self.parsed_frames - 1}: {e}"
                )
                self.running = False
                self._buffer.producer_finished = True
            # If buffer is more than 90% full, pause the simulation
            # until we have more space
            if self.pausable:
                if self._full_frames >= self._buffer.capacity // (1.1):
                    self._toggle_pause_simulation()

            # Will attempt to reconnect if the socket is closed
            # in case of a network error
            self._check_end_of_simulation()

            # Reset timeout if it was changed during the loop
            self._conn.settimeout(None)

        self._ensure_disconnect()
        return

    def _check_end_of_simulation(self):
        if not self._is_socket_connected():
            # Continue reading if socket contains bytes
            try:
                self._conn.settimeout(10 * self._parse_frame_time)
                self._conn.recv(1, socket.MSG_PEEK)
            except socket.timeout:
                try:
                    # Attempt to reconnect
                    self._connection_sequence()
                except Exception as e:
                    logging.debug(
                        "IMDProducer: Assuming simulation is over at frame "
                        + f"#{self.parsed_frames - 1} due to error: {e}"
                    )
                    self.running = False
                    self._buffer.producer_finished = True

    def stop(self):
        """Method to stop the thread and disconnect safely."""
        self.running = False
        self._ensure_disconnect()

    def _handle_signal(self, signum, frame):
        """Handle SIGINT and SIGTERM signals."""
        self.running = False
        self._ensure_disconnect()

    def _ensure_disconnect(self):
        """Ensure the connection is closed only once."""
        if not self._is_disconnected:
            self._disconnect()
            self._is_disconnected = True

    def _recv_n_bytes(self, num_bytes):
        """Receive an arbitrary number of bytes from the socket."""
        data = bytearray(num_bytes)
        view = memoryview(data)
        total_received = 0
        while total_received < num_bytes:
            chunk = self._conn.recv(num_bytes - total_received)
            if not chunk:
                raise ConnectionError("Socket connection was closed")
            view[total_received : total_received + len(chunk)] = chunk
            total_received += len(chunk)
        return data

    def _recv_body_bytes(self):
        """Efficiently receives the correct number of bytes
        for position data
        """
        total_received = 0
        num_bytes = self._expected_data_bytes
        data = self._body_byte_buf
        view = self._body_byte_view
        while total_received < num_bytes:
            chunk = self._conn.recv(num_bytes - total_received)
            if not chunk:
                raise ConnectionError("Socket connection was closed")
            view[total_received : total_received + len(chunk)] = chunk
            total_received += len(chunk)
        return data

    def _recv_energy_bytes(self):
        """Efficiently receives the correct number of bytes
        for energy data
        """
        total_received = 0
        num_bytes = IMDENERGYPACKETLENGTH
        data = self._energy_byte_buf
        view = self._energy_byte_view
        while total_received < num_bytes:
            chunk = self._conn.recv(num_bytes - total_received)
            if not chunk:
                raise ConnectionError("Socket connection was closed")
            view[total_received : total_received + len(chunk)] = chunk
            total_received += len(chunk)
        return data

    def _expect_header(self, expected_type=None, expected_value=None):
        """
        Read a header packet from the socket.
        """
        header = parse_header_bytes(self._recv_n_bytes(IMDHEADERSIZE))
        if expected_type is not None and header.type != expected_type:
            raise ValueError(
                f"Expected packet type {expected_type}, got {header.type}"
            )
        elif expected_value is not None and header.length != expected_value:
            raise ValueError(
                f"Expected packet length {expected_value}, got {header.length}"
            )
        return header

    def _disconnect(self):
        if self._conn:
            try:
                disconnect = create_header_bytes(IMDType.IMD_DISCONNECT, 0)
                self._conn.sendall(disconnect)
                self._conn.close()
                print("Disconnected from server")
            except Exception as e:
                print(f"Error during disconnect: {e}")

    def _is_socket_connected(self):
        try:
            self._conn.getpeername()
            return True
        except socket.error:
            return False


class CircularByteBuf:
    """
    Acts as interface between producer and consumer threads
    """

    # NOTE: Use 1 buffer for pos, vel, force rather than 3
    def __init__(self, buffer_size, n_atoms, ts):
        # a frame is the number of bytes needed to hold
        # energies + positions
        self._frame_size = 40 + (n_atoms * 12)
        self._n_atoms = n_atoms
        self._ts = ts
        self._ts_buf = bytearray(self._frame_size)
        self._body_bytes = n_atoms * 12

        if self._frame_size > buffer_size:
            raise MemoryError(
                f"Requested buffer size of {buffer_size} "
                + f"doesn't meet memory requirement of {self._frame_size} "
                + f"(energy and position data for {n_atoms})"
            )

        self._buf_frames = (buffer_size) // self._frame_size
        self._buf = bytearray(self._buf_frames * self._frame_size)
        self._memview = memoryview(self._buf)

        self._mutex = threading.Lock()
        self._not_empty = threading.Condition(self._mutex)
        self._not_full = threading.Condition(self._mutex)
        self._has_space = threading.Condition()

        # _empty refers to the number of empty
        # "frames" in the buffer
        # where a frame refers to space in a numpy array
        # for a full ts frame of positions
        self._empty = self._buf_frames
        self._full = 0
        self._fill = 0
        self._use = 0

        # Timing for analysis
        self._t1 = None
        self._t2 = None
        self._start = True
        self._analyze_frame_time = None

        # Syncing reader and producer
        self._producer_finished = False

    def inform_endianness(self, endianness):
        """Producer thread must determine
        endianness before sending data to the buffer"""
        self._end = endianness

    def insert(self, energy, pos):
        with self._not_full:
            while not self._empty:
                self._not_full.wait()

            self._memview[self._fill : self._fill + IMDENERGYPACKETLENGTH] = (
                energy
            )
            self._fill += IMDENERGYPACKETLENGTH

            self._memview[self._fill : self._fill + self._body_bytes] = pos
            self._fill += self._body_bytes

            self._fill %= len(self._buf)

            self._full += 1
            self._empty -= 1
            self._not_empty.notify()

            logging.debug(
                f"IMDProducer: Buffer at ({self._full}/{self._buf_frames}) frames"
            )
        # This doesn't need to be a perfect count
        return self._full

    def consume_next_timestep(self):
        # Start timer- one frame of analysis is starting (including removal
        # from buffer)
        if self._start:
            self._t1 = time.time()

        # Stop the timer- one frame of analysis has finished
        if not self._start:
            self._t2 = time.time()
            logging.debug(
                f"IMDReader: Frame #{self._frame + 1} took {self._t2 - self._t1} seconds to analyze"
            )
            self._analyze_frame_time = self._t2 - self._t1

        self._start = not self._start

        with self._not_empty:
            while not self._full:
                self._not_empty.wait()

            self._ts_buf[:] = self._memview[
                self._use : self._use + self._frame_size
            ]
            self._use = (self._use + self._frame_size) % len(self._buf)

            self._full -= 1
            self._empty += 1
            self._not_full.notify()

        self._ts.data["step"] = np.frombuffer(
            self._ts_buf, dtype=f"{self._end}i4", offset=0, count=1
        )[0]
        # absolute temperature
        self._ts.data["temperature"] = np.frombuffer(
            self._ts_buf, dtype=f"{self._end}f4", offset=4, count=1
        )[0]
        self._ts.data["total_energy"] = np.frombuffer(
            self._ts_buf, dtype=f"{self._end}f4", offset=8, count=1
        )[0]
        self._ts.data["potential_energy"] = np.frombuffer(
            self._ts_buf, dtype=f"{self._end}f4", offset=12, count=1
        )[0]
        self._ts.data["van_der_walls_energy"] = np.frombuffer(
            self._ts_buf, dtype=f"{self._end}f4", offset=16, count=1
        )[0]
        self._ts.data["coulomb_energy"] = np.frombuffer(
            self._ts_buf, dtype=f"{self._end}f4", offset=20, count=1
        )[0]
        self._ts.data["bonds_energy"] = np.frombuffer(
            self._ts_buf, dtype=f"{self._end}f4", offset=24, count=1
        )[0]
        self._ts.data["angles_energy"] = np.frombuffer(
            self._ts_buf, dtype=f"{self._end}f4", offset=28, count=1
        )[0]
        self._ts.data["dihedrals_energy"] = np.frombuffer(
            self._ts_buf, dtype=f"{self._end}f4", offset=32, count=1
        )[0]
        self._ts.data["improper_dihedrals_energy"] = np.frombuffer(
            self._ts_buf, dtype=f"{self._end}f4", offset=36, count=1
        )[0]

        self._ts.positions = np.frombuffer(
            self._ts_buf, dtype=f"{self._end}f4", offset=40
        ).reshape((self._n_atoms, 3))

        self.print_buffer()

    def wait_almost_empty(self):
        with self._not_full:
            while self._full >= 1:
                self._not_full.wait()

            logging.debug(
                "IMDProducer: Buffer almost empty, resuming simulation"
            )

    @property
    def capacity(self):
        return self._buf_frames

    @property
    def analyze_frame_time(self):
        if self._analyze_frame_time is not None:
            return self._analyze_frame_time
        else:
            return None

    @property
    def producer_finished(self):
        return self._producer_finished

    @producer_finished.setter
    def producer_finished(self, value):
        if value:
            # Only producer should call this
            self._producer_finished = True
        else:
            raise ValueError("Cannot unset producer_finished")
