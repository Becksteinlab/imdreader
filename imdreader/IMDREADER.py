from MDAnalysis.coordinates.base import (
    ReaderBase,
    FrameIteratorIndices,
    FrameIteratorAll,
    FrameIteratorSliced,
)
from .IMDProtocol import *
import socket
import threading
import numpy as np
import numbers
from MDAnalysis.lib.util import store_init_arguments


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
        buffer_size=2**23,
        **kwargs,
    ):
        super(IMDReader, self).__init__(filename, **kwargs)
        self._host, self._port = parse_host_port(filename)
        self._socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        # NOTE: Replace me with header packet which contains this information
        # OR get this information from the topology?
        if not n_frames:
            raise ValueError("n_frames must be specified")
        self.n_frames = n_frames
        if not num_atoms:
            raise ValueError("num_atoms must be specified")
        self.n_atoms = num_atoms
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)

        # The body of a force or position packet should contain
        # (4 bytes per float * 3 atoms * n_atoms) bytes
        self._expected_data_bytes = 12 * self.n_atoms
        # Preallocate body byte holder
        self._body_byte_buf = bytearray(self._expected_data_bytes)
        self._body_byte_view = memoryview(self._body_byte_buf)
        # < represents little endian and > represents big endian
        # we assume big by default and use handshake to check
        self._endianness = ">"
        self._connect_to_server()
        self._await_IMD_handshake()
        self._start_producer_thread()

        # NOTE: We need to know if we are getting forces or not to construct the correct buffers
        self._buffer = CircularNPBuf(
            buffer_size, self.n_atoms, self._endianness
        )
        self._frame = -1

    def _connect_to_server(self):
        """
        Listen for a connection on the specified port.
        """
        self._socket.connect((self._host, self._port))
        self._conn = self._socket

    def _await_IMD_handshake(self):
        """
        Wait for the server to send a handshake packet, set endianness,
        and check IMD Protocol version.
        """
        print("waiting for handshake...")
        handshake = self._expect_header(IMDType.IMD_HANDSHAKE)
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
                self._endianness = "<"
        print("handshake received")

    def _pause_simulation(self):
        """
        Pause the simulation for buffer filling purposes
        or to exit gracefully.
        """
        pause = create_header_bytes(IMDType.IMD_PAUSE, 0)
        self._conn.sendall(pause)

    def _send_go_packet(self):
        """
        Send a go packet to the client to start the simulation
        and begin receiving data.
        """
        print("sending go packet...")
        go = create_header_bytes(IMDType.IMD_GO, 0)
        self._conn.sendall(go)

    def _start_producer_thread(self):
        """
        Starts the thread which reads packets from the socket and
        places them into a buffer for the parent thread to read.
        """
        self._producer_thread = threading.Thread(target=self._control_flow)
        self._producer_thread.daemon = True
        self._producer_thread.start()

    def _control_flow(self):
        """
        Producer thread method. Reads from the socket and
        sends a 'pause' signal if needed.
        """
        self._parsed_frames = 0

        while self._parsed_frames < self.n_frames:
            header = self._expect_header()
            if header.type == IMDType.IMD_ENERGIES and header.length == 1:
                self._recv_energies()
                header2 = self._expect_header(IMDType.IMD_FCOORDS)
                if header2.length == self.n_atoms:
                    self._recv_fcoords()
                    self._parsed_frames += 1
                else:
                    raise ValueError(
                        f"Unexpected coordinate packet length of "
                        + f"{header2.length} bytes, expected "
                        + f"{self._expected_data_bytes} bytes"
                    )
            else:
                raise ValueError(
                    f"Unexpected packet type {header.type} "
                    + f"of length {header.length}"
                )
        print(self._parsed_frames)

    def _expect_header(self, expected_type=None):
        """
        Read a header packet from the socket.
        """
        header = parse_header_bytes(self._recv_n_bytes(IMDHEADERSIZE))
        if expected_type is not None and header.type != expected_type:
            raise ValueError(
                "Expected packet type {}, got {}".format(
                    expected_type, header.type
                )
            )
        return header

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

    def _recv_fcoords(self):
        self._buffer.insert(self._recv_body_bytes())

    def _recv_energies(self):
        self._recv_n_bytes(IMDENERGYPACKETLENGTH)

    def _read_next_timestep(self):
        if self._frame == -1:
            self._send_go_packet()
        return self._read_frame(self._frame + 1)

    def _read_frame(self, frame):
        # NOTE: read everything provided, don't throw out forces
        self._frame = frame
        pos = self._buffer.get_frame()
        self.ts.positions = pos
        self.ts.frame = self._frame
        self.ts.dimensions = None
        return self.ts

    def _reopen(self):
        # NOTE: try setting self._frame to -1
        pass

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
        self._socket.close()
        print("IMDReader shut down gracefully.")

    # Incompatible methods
    def copy(self):
        raise NotImplementedError("IMDReader does not support copying")

    def __getitem__(self, frame):
        """This method from ProtoReader must be overridden
        to prevent slicing that doesn't make sense in a stream.

        We want the user to only be able to call this once.

        Original docstring
        ##################

        Return the Timestep corresponding to *frame*.

        If `frame` is a integer then the corresponding frame is
        returned. Negative numbers are counted from the end.

        If frame is a :class:`slice` then an iterator is returned that
        allows iteration over that part of the trajectory.

        Note
        ----
        *frame* is a 0-based frame index.
        """
        if isinstance(frame, numbers.Integral):
            frame = self._apply_limits(frame)
            return self._read_frame_with_aux(frame)
        elif isinstance(frame, (list, np.ndarray)):
            if len(frame) != 0 and isinstance(frame[0], (bool, np.bool_)):
                # Avoid having list of bools
                frame = np.asarray(frame, dtype=bool)
                # Convert bool array to int array
                frame = np.arange(len(self))[frame]
            return FrameIteratorIndices(self, frame)
        elif isinstance(frame, slice):
            start, stop, step = self.check_slice_indices(
                frame.start, frame.stop, frame.step
            )
            if start == 0 and stop == len(self) and step == 1:
                return FrameIteratorAll(self)
            else:
                return FrameIteratorSliced(self, frame)
        else:
            raise TypeError(
                "Trajectories must be an indexed using an integer,"
                " slice or list of indices"
            )


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


class CircularNPBuf:
    """
    Thread-safe circular numpy buffer
    """

    # NOTE: Use 1 buffer for pos, vel, force rather than 3
    def __init__(self, buffer_size, n_atoms, endianness):
        # a frame is the number of 32bit floats
        # needed to hold postition data for n_atoms
        self._frame_size = n_atoms * 3
        self._shape = (n_atoms, 3)
        self._n_atoms = n_atoms
        self._end = endianness

        if self._frame_size * 4 > buffer_size:
            raise MemoryError(
                f"Requested buffer size of {buffer_size} "
                + f"doesn't meet memory requirement of {self._frame_size * 4} "
                + f"(position data for {n_atoms})"
            )

        self._buf_frames = (buffer_size // 4) // self._frame_size
        self._buf = np.empty(
            self._buf_frames * self._frame_size, dtype=np.float32
        )

        self._mutex = threading.Lock()
        self._not_empty = threading.Condition(self._mutex)
        self._not_full = threading.Condition(self._mutex)

        # _empty refers to the number of empty
        # "frames" in the buffer
        # where a frame refers to space in a numpy array
        # for a full ts frame of positions
        self._empty = self._buf_frames
        self._full = 0
        self._fill = 0
        self._use = 0

        self._curr_frame = np.empty(self._shape, dtype=np.float32)

    def insert(self, data_bytes):
        with self._not_full:
            while not self._empty:
                self._not_full.wait()

            self._buf[self._fill : self._fill + self._frame_size] = (
                np.frombuffer(data_bytes, dtype=f"{self._end}f4").copy()
            )
            self._fill = (self._fill + self._frame_size) % len(self._buf)

            self._full += 1
            self._empty -= 1
            self._not_empty.notify()

    def get_frame(self):
        # NOTE: Init buffer in CircularNPBuf init()
        with self._not_empty:
            while not self._full:
                self._not_empty.wait()
            self._curr_frame[:] = self._buf[
                self._use : self._use + self._frame_size
            ].reshape(self._shape)
            self._use = (self._use + self._n_atoms) % len(self._buf)

            self._full -= 1
            self._empty += 1
            self._not_full.notify()
            return self._curr_frame
