"""

Example: Streaming an IMD v2 trajectory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To stream a trajectory from GROMACS or another simulation engine that supports 
IMD v2, ensure that the simulation engine is running and waiting for an IMD connection.

For example, in GROMACS, you can use ``gmx mdrun`` with the ``-imdwait`` flag
to ensure that GROMACS will wait for a client before starting the simulation.
In GROMACS, you will know that the simulation is ready and waiting for the
MDAnalysis IMDReader client when this line is printed to the terminal:

.. code-block:: none

    IMD: Will wait until I have a connection and IMD_GO orders.

Once the simulation is ready for a client connection, setup your :class:`Universe`
like this: ::

    import MDAnalysis as mda
    # Pass host and port of the listening GROMACACS simulation
    # server as the trajectory argument
    # Number of atoms must be provided for internal buffer allocation
    u = mda.Universe("topology.tpr", "localhost:8888", num_atoms=100)

For more information on IMD v2 as it is implemented in GROMACS, see the imd command line
arguments described `here <https://manual.gromacs.org/documentation/5.1/onlinehelp/gmx-mdrun.html>`_,
the :ref:`v2_spec` writeup, and this `example simulation  <https://www.mpinat.mpg.de/grubmueller/interactivemd>`_. Note that 
setting the following line in your ``.mdp`` file allows you to change which subset of the
simulation is sent to the client via IMD v2:

.. code-block:: none

    ; Group to display and/or manipulate in interactive MD session
    IMD-group                = System

Classes
^^^^^^^

.. autoclass:: IMDReader
   :members:
   :inherited-members:

"""

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
        """
        Parameters
        ----------
        filename : a string of the form "host:port" where host is the hostname
            or IP address of the listening GROMACS server and port
            is the port number.

        """
        super(IMDReader, self).__init__(filename, **kwargs)

        # NOTE: Replace me with header packet which contains this information
        # OR get this information from the topology?
        if not num_atoms:
            raise ValueError("num_atoms must be specified")
        self.n_atoms = num_atoms
        self.ts = self._Timestep(
            self.n_atoms, positions=True, **self._ts_kwargs
        )

        self.units = {
            "time": "ps",
            "length": "nm",
            "force": "kJ/(mol*nm)",
        }

        self._buffer = CircularByteBuf(buffer_size, self.n_atoms, self.ts)

        self._producer = IMDProducer(
            filename, self._buffer, self.n_atoms, socket_bufsize
        )

        self._frame = -1

    def _read_next_timestep(self):
        if self._frame == -1:
            self._producer.start()

        return self._read_frame(self._frame + 1)

    def _read_frame(self, frame):

        # loads the timestep with step, positions, and energy
        self._buffer.consume_next_timestep()
        self.ts.frame = self._frame + 1
        self.ts.dimensions = None

        # Must set frame after read occurs successfully
        # Since buffer raises IO error
        # after producer is finished and there are no more frames
        self._frame = frame

        return self.ts

    @property
    def n_frames(self):
        """Changes as stream is processed unlike other readers"""
        return self._frame + 1

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
        """Gracefully shut down the reader. Stops the producer thread."""
        self._producer.stop()
        self._producer.join()
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
        super(IMDProducer, self).__init__()
        self._host, self._port = parse_host_port(filename)
        self._conn = None
        self.running = False

        self._buffer = buffer
        self.n_atoms = n_atoms
        self._expected_data_bytes = 12 * n_atoms
        self._socket_bufsize = socket_bufsize
        self.pausable = pausable
        self.paused = False
        # < represents little endian and > represents big endian
        # we assume big by default and use handshake to check
        self.parsed_frames = 0
        self._full_frames = 0
        self._parse_frame_time = 0
        # Saving memory by preallocating space for the frame
        # we're loading into the buffer
        self._energy_byte_buf = bytearray(IMDENERGYPACKETLENGTH)
        self._energy_byte_view = memoryview(self._energy_byte_buf)
        self._body_byte_buf = bytearray(self._expected_data_bytes)
        self._body_byte_view = memoryview(self._body_byte_buf)

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
                # Don't call stop, simulation hasn't started yet
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

    def _pause_simulation(self):
        """
        Block the simulation until the buffer has more space.
        """
        try:
            pause = create_header_bytes(IMDType.IMD_PAUSE, 0)
            self._conn.sendall(pause)
            logger.debug(
                "IMDProducer: Pausing simulation because buffer is almost full"
            )
            self.paused = True
        except ConnectionResetError as e:
            self._is_disconnected = True
            self.stop()
            logger.warning(
                "IMDProducer: Connection reset during pausing, "
                + "data likely lost in frame "
                + f"{self.parsed_frames}"
            )

    def _unpause_simulation(self):
        try:
            self._buffer.wait_almost_empty()
            logger.debug("IMDProducer: Unpausing simulation, buffer has space")
            pause = create_header_bytes(IMDType.IMD_PAUSE, 0)
            self._conn.sendall(pause)
            self.paused = False
            self._conn.settimeout(None)
        except ConnectionResetError as e:
            self._is_disconnected = True
            self.stop()
            logger.warning(
                "IMDProducer: Connection reset during unpausing, "
                + "data likely lost in frame "
                + f"{self.parsed_frames}"
            )

    def _connection_sequence(self):
        """
        Establish connection with the server and perform
        the handshake and go packet exchange.
        """
        self._conn = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        if self._socket_bufsize is not None:
            self._conn.setsockopt(
                socket.SOL_SOCKET, socket.SO_RCVBUF, self._socket_bufsize
            )
        self._conn.settimeout(5)
        try:
            self._conn.connect((self._host, self._port))
        except ConnectionRefusedError:
            self._buffer.producer_finished = True
            raise ConnectionRefusedError(
                f"IMDProducer: Connection to {self._host}:{self._port} refused"
            )
        self._conn.settimeout(None)
        self._await_IMD_handshake()
        self._send_go_packet()
        self.running = True
        self._is_disconnected = False

    def run(self):
        """
        Producer thread method. Reads from the socket and
        sends a 'pause' signal if needed.
        """
        self._connection_sequence()

        while self.running:
            self._parse_and_insert_frame()

            # If buffer is more than 50% full, pause the simulation
            if self.pausable:
                if (
                    not self.paused
                    and self._full_frames >= self._buffer.capacity // 10
                ):
                    self._pause_simulation()

            # Simulation will unpause in _if
            # 1. buffer is almost empty
            # 2. socket is emptied

            # Only check for simulation end if we're not paused
            # otherwise we can get a false positive
            if not self.paused:
                self._check_end_of_simulation()

            # Reset timeout if it was changed during the loop
            self._conn.settimeout(None)

        return

    def _parse_and_insert_frame(self):
        with timeit() as parse_frame:
            self._expect_header(
                expected_type=IMDType.IMD_ENERGIES, expected_value=1
            )
            energies = self._recv_data("energy")
            self._expect_header(
                expected_type=IMDType.IMD_FCOORDS,
                expected_value=self.n_atoms,
            )
            pos = self._recv_data("body")
            self._full_frames = self._buffer.insert(energies, pos)

        self.parsed_frames += 1
        # Use the longest parse frame time to calculate
        # the timeout for the next frame
        self._parse_frame_time = max(
            parse_frame.elapsed, self._parse_frame_time
        )
        logger.debug(
            f"IMDProducer: Added frame #{self.parsed_frames - 1} to buffer in {parse_frame.elapsed} seconds"
        )

    def stop(self):
        # Tell reader not to expect more frames to be added
        self._buffer.producer_finished = True
        # MUST disconnect before stopping run loop
        self._ensure_disconnect()
        self.running = False

    def _check_end_of_simulation(self):
        # It is the server's reponsibility
        # to close the connection, but this may not happen in time
        # for us to check during the last frame
        # Therefore, use a timeout on a peek to check if the server has closed
        try:
            logger.debug(
                f"IMDProducer: Checking for frame #{self.parsed_frames}"
            )
            # Continue reading if socket contains bytes
            self._conn.settimeout(5 * self._parse_frame_time)
            b = self._conn.recv(1, socket.MSG_PEEK)
            # Two cases for ending producer reads from socket

            # case 1: server has closed the connection and we have no more data to read
            if not b:
                logger.debug(
                    "IMDProducer: Assuming simulation is over at frame "
                    + f"#{self.parsed_frames - 1} due to closed connection"
                )
                self.running = False
                self._buffer.producer_finished = True
            # case 2: server has not closed the connection and we have no more data to read
        except socket.timeout:
            logger.debug(
                "IMDProducer: Assuming simulation is over at frame "
                + f"#{self.parsed_frames - 1} due to read timeout"
            )
            self.running = False
            self._buffer.producer_finished = True

    def _handle_signal(self, signum, frame):
        """Handle SIGINT and SIGTERM signals."""
        self.running = False
        self._ensure_disconnect()

    def _recv_data(self, type):
        """Used to receive headers and data packets from the socket.
        For energies and positions, the data is stored in a preallocated
        buffer to avoid memory allocation during the run loop.

        This method will behave differently is self.paused is True,
        since having no more data available in the socket doesn't
        indicate the connection is closed, but rather that
        the simulation is now ready to be unpaused

        ``type`` can be one of: "header", "body", "energy"
        """

        if type == "header":
            data = bytearray(IMDHEADERSIZE)
            view = memoryview(data)

        elif type == "body":
            data = self._body_byte_buf
            view = self._body_byte_view

        elif type == "energy":
            data = self._energy_byte_buf
            view = self._energy_byte_view

        if self.paused:
            self._conn.settimeout(0)

        total_received = 0
        while total_received < len(data):
            try:
                chunk = self._conn.recv(len(data) - total_received)
                if not chunk:
                    self._is_disconnected = True
                    self.stop()
                    logger.warning(
                        f"IMDProducer: Data likely lost in frame {self.parsed_frames}"
                    )
                    raise ConnectionError("Socket connection was closed")
            except BlockingIOError:
                # If we're paused, we can assume the simulation is ready to unpause
                if self.paused:
                    self._unpause_simulation()
                    continue
            view[total_received : total_received + len(chunk)] = chunk
            total_received += len(chunk)

        return data

    def _expect_header(self, expected_type=None, expected_value=None):
        """
        Read a header packet from the socket.
        """
        header = parse_header_bytes(self._recv_data("header"))
        if expected_type is not None and header.type != expected_type:
            self.stop()
            raise ValueError(
                f"Expected packet type {expected_type}, got {header.type}"
            )
        elif expected_value is not None and header.length != expected_value:
            self.stop()
            raise ValueError(
                f"Expected packet length {expected_value}, got {header.length}"
            )
        return header

    def _ensure_disconnect(self):
        """Ensure the connection is closed only once."""
        # We rely on the server to close the connection
        # so only disconnect in cases where it isn't the
        # server's responsibility (stopping mid-simulation)
        if not self._is_disconnected and self.running:
            self._disconnect()
            self._is_disconnected = True

    def _disconnect(self):
        if self._conn:
            try:
                disconnect = create_header_bytes(IMDType.IMD_DISCONNECT, 0)
                self._conn.sendall(disconnect)
                self._conn.close()
                logger.debug("IMDProducer: Disconnected from server")
            except Exception as e:
                print(f"IMDProducer: Error during disconnect: {e}")


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
        self._frame = 0

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

            logger.debug(
                f"IMDProducer: Buffer at ({self._full}/{self._buf_frames}) frames"
            )
        # This doesn't need to be a perfect count
        return self._full

    def consume_next_timestep(self):
        # Start timer- one frame of analysis is starting (including removal
        # from buffer)

        self._t1 = self._t2
        self._t2 = time.time()
        if self._t1 is not None:
            logger.debug(
                f"IMDReader: Frame #{self._frame - 1} analyzed in {self._t2 - self._t1} seconds"
            )
            self._analyze_frame_time = self._t2 - self._t1

        with self._not_empty:
            while not self._full and not self.producer_finished:
                self._not_empty.wait()

            # Buffer is responsible for stopping iteration
            if self.producer_finished and not self._full:
                raise IOError from None

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

        self._frame += 1

    def wait_almost_empty(self):
        with self._not_full:
            while self._full >= 1:
                self._not_full.wait()

            logger.debug(
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
            # Wake up the reader if it's waiting for more frames
            with self._not_full:
                self._producer_finished = True
                self._not_empty.notify()
        else:
            raise ValueError("Cannot unset producer_finished")
