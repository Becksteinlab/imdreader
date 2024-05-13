from MDAnalysis.coordinates.base import (
    ReaderBase,
    FrameIteratorIndices,
    FrameIteratorAll,
    FrameIteratorSliced,
)
from MDAnalysis.lib.util import store_init_arguments
from .IMDProtocol import *
import socket
import threading
import numpy as np
import numbers
import signal


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
        n_frames=None,
        buffer_size=2**26,
        socket_bufsize=None,
        **kwargs,
    ):
        super(IMDReader, self).__init__(filename, **kwargs)

        # NOTE: Replace me with header packet which contains this information
        # OR get this information from the topology?
        if not n_frames:
            raise ValueError("n_frames must be specified")
        self._n_frames = n_frames
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
        # When we've read all frames,
        # we wait for the producer to finish
        if self._frame == self._n_frames - 1:
            self._producer.join()
        return self._read_frame(self._frame + 1)

    def _read_frame(self, frame):
        if frame >= self._n_frames:
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
        return self._n_frames

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
        raise ValueError("IMDReader: Trajectory can only be read in for loop")


def parse_host_port(filename):
    # Check if the format is correct
    parts = filename.split(":")
    if len(parts) == 2:
        host = parts[0]  # Hostname part
        try:
            port = int(parts[1])  # Convert the port part to an integer
            return (host, port)
        except ValueError:
            # Handle the case where the port is not a valid integer
            raise ValueError("Port must be an integer")
    else:
        # Handle the case where the format does not match "host:port"
        raise ValueError("Filename must be in the format 'host:port'")


class IMDProducer(threading.Thread):
    """
    Producer thread for IMDReader. Reads packets from the socket
    and places them into the shared buffer.
    """

    def __init__(
        self,
        filename,
        buffer,
        n_frames,
        n_atoms,
        socket_bufsize=None,
        pausable=True,
    ):
        super(IMDProducer, self).__init__(daemon=True)
        self._host, self._port = parse_host_port(filename)
        self._conn = None
        self.running = False

        self._buffer = buffer
        self._n_frames = n_frames
        self.n_atoms = n_atoms
        self._expected_data_bytes = 12 * n_atoms
        self._socket_bufsize = socket_bufsize
        self.pausable = pausable
        # < represents little endian and > represents big endian
        # we assume big by default and use handshake to check
        self.parsed_frames = 0
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
        Pause/unpause the simulation for buffer filling purposes.
        Put the thread to sleep until the buffer has more space.
        """
        pause = create_header_bytes(IMDType.IMD_PAUSE, 0)
        self._conn.sendall(pause)

    def run(self):
        """
        Producer thread method. Reads from the socket and
        sends a 'pause' signal if needed.
        """
        self._conn = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        if self._socket_bufsize:
            self._conn.setsockopt(
                socket.SOL_SOCKET, socket.SO_RCVBUF, self._socket_bufsize
            )
        self._conn.connect((self._host, self._port))
        self._await_IMD_handshake()
        # NOTE: Buffer allocation ideally goes here
        self._send_go_packet()
        self.running = True
        percent_full = 0

        while (self.parsed_frames < self._n_frames) and self.running:
            # if self.pausable and percent_full > 0.5:
            #   self._toggle_pause_simulation()

            header1 = self._expect_header(
                expected_type=IMDType.IMD_ENERGIES, expected_value=1
            )
            energies = self._recv_energy_bytes()
            header2 = self._expect_header(
                IMDType.IMD_FCOORDS, expected_value=self.n_atoms
            )
            pos = self._recv_body_bytes()
            self._buffer.insert(energies, pos)
            self.parsed_frames += 1

        self._ensure_disconnect()
        return

    def stop(self):
        """Method to stop the thread and disconnect safely."""
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

    def _handle_signal(self, signum, frame):
        self._ensure_disconnect()

    def _disconnect(self):
        if self._conn:
            try:
                disconnect = create_header_bytes(IMDType.IMD_DISCONNECT, 0)
                self._conn.sendall(disconnect)
                self._conn.close()
                print("Disconnected from server")
            except Exception as e:
                print(f"Error during disconnect: {e}")


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

    def consume_next_timestep(self):
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

    def print_buffer(self):
        """
        Debugging method to show buffer fullness
        """
        print(f"Full: {self._full}/{self._buf_frames}")
