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
    u = mda.Universe("topology.tpr", "localhost:8888")

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
from MDAnalysis.coordinates import core
from MDAnalysis.lib.util import store_init_arguments
from .IMDProtocol import *
from .util import *
import socket
import threading
import numpy as np
import signal
import logging
import time

logger = logging.getLogger(__name__)


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
        n_atoms=None,
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
        self._producer = None
        super(IMDReader, self).__init__(filename, **kwargs)
        if not n_atoms:
            raise ValueError("`n_atoms` kwarg must be specified")
        self.n_atoms = n_atoms

        self.units = {
            "time": "ps",
            "length": "nm",
            "force": "kJ/(mol*nm)",
        }

        self._host, self._port = parse_host_port(filename)
        self._buffer_size = buffer_size
        self._socket_bufsize = socket_bufsize

        self._frame = -1

    def _read_next_timestep(self):
        if self._frame == -1:
            # Reader is responsible for performing handshake
            # and parsing the configuration before
            # passing the connection off to the appropriate producer
            # and allocating an appropriate buffer
            conn = self._connect_to_server()
            imdsinfo = self._await_IMD_handshake(conn)

            if imdsinfo.version == 2:
                self._ts = self._Timestep(self.n_atoms, positions=True)
                self._buffer = CircularByteBuf(
                    self._buffer_size, self.n_atoms, self._ts, imdsinfo
                )
                self._producer = IMDv2Producer(
                    conn,
                    self._buffer,
                    imdsinfo,
                    self.n_atoms,
                )
            # Producer responsible for sending go packet
            self._producer.start()

        return self._read_frame(self._frame + 1)

    def _read_frame(self, frame):

        # loads the timestep with step, positions, and energy
        self._buffer.consume_next_timestep()
        self._ts.frame = self._frame + 1

        # Must set frame after read occurs successfully
        # Since buffer raises IO error
        # after producer is finished and there are no more frames
        self._frame = frame

        return self._ts

    def _connect_to_server(self):
        """
        Establish connection with the server, failing out if this
        does not occur within 5 seconds.
        """
        conn = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        if self._socket_bufsize is not None:
            conn.setsockopt(
                socket.SOL_SOCKET, socket.SO_RCVBUF, self._socket_bufsize
            )
        conn.settimeout(5)
        try:
            conn.connect((self._host, self._port))
        except ConnectionRefusedError:
            logger.error(
                f"IMDReader: Connection to {self._host}:{self._port} refused"
            )
            raise ConnectionRefusedError(
                f"IMDReader: Connection to {self._host}:{self._port} refused"
            )
        conn.settimeout(None)
        return conn

    def _await_IMD_handshake(self, conn) -> IMDSessionInfo:
        """
        Wait for the server to send a handshake packet, then parse
        endianness and version information and IMD session configuration.
        """
        end = ">"
        ver = None

        handshake = parse_header_bytes(conn.recv(IMDHEADERSIZE))
        if handshake.type != IMDType.IMD_HANDSHAKE:
            raise ValueError(
                f"Expected packet type {IMDType.IMD_HANDSHAKE}, got {handshake.type}"
            )

        if handshake.length not in IMDVERSIONS:
            # Try swapping endianness
            swapped = struct.unpack("<i", struct.pack(">i", handshake.length))[
                0
            ]
            if swapped not in IMDVERSIONS:
                err_version = min(swapped, handshake.length)
                # Don't call stop, simulation hasn't started yet
                raise ValueError(
                    f"Incompatible IMD version. Expected version in {IMDVERSIONS}, got {err_version}"
                )
            else:
                end = "<"
                ver = swapped
        else:
            ver = handshake.length

        sinfo = None
        if ver == 2:
            # IMD v2 does not send a configuration packet
            sinfo = IMDSessionInfo(
                version=ver,
                endianness=end,
                imdterm=None,
                imdwait=None,
                imdpull=None,
                wrapped_coords=False,
                energies=1,
                dimensions=0,
                positions=1,
                velocities=0,
                forces=0,
            )
        return sinfo

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
        if self._producer is not None:
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


class AbstractIMDProducer(abc.ABC, threading.Thread):

    def __init__(self, conn, buffer, imdsinfo, n_atoms):
        # call threading.Thread init
        super(AbstractIMDProducer, self).__init__()
        self._conn = conn
        self.imdsinfo = imdsinfo
        self.running = False
        self.paused = False
        self._is_disconnected = False
        self._buffer = buffer

        self.parsed_frames = 0
        self._full_frames = 0
        self._parse_frame_time = 0

        # The body of a force or position packet should contain
        # (4 bytes per float * 3 atoms * n_atoms) bytes
        self.n_atoms = n_atoms
        self._data_bytes = 12 * n_atoms

        self._byte_dict = {}
        if self.imdsinfo.energies > 0:
            self._byte_dict["energies"] = bytearray(40)
        if self.imdsinfo.dimensions > 0:
            self._byte_dict["dimensions"] = bytearray(36)
        if self.imdsinfo.positions > 0:
            self._byte_dict["positions"] = bytearray(self._data_bytes)
        if self.imdsinfo.velocities > 0:
            self._byte_dict["velocities"] = bytearray(self._data_bytes)
        if self.imdsinfo.forces > 0:
            self._byte_dict["forces"] = bytearray(self._data_bytes)

    @abc.abstractmethod
    def _parse_and_insert_frame(self):
        pass

    @abc.abstractmethod
    def _check_end_of_simulation(self):
        pass

    def _send_go_packet(self):
        """
        Send a go packet to the client to start the simulation
        and begin receiving data.
        """
        print("sending go packet...")
        go = create_header_bytes(IMDType.IMD_GO, 0)
        self._conn.sendall(go)
        logger.debug("IMDProducer: Sent go packet to server")

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
            # make the socket non-blocking so we can use BlockingIOError to check
            # if we can unpause
            self._conn.settimeout(0)
            self.paused = True
        except ConnectionResetError as e:
            self._is_disconnected = True
            self.stop()
            logger.warning(
                "IMDProducer: Connection reset during pausing, "
                "data likely lost in frame {}".format(self.parsed_frames)
            )

    def _unpause_simulation(self):
        try:
            logger.debug(
                "IMDProducer: Waiting to unpause until buffer almost empty"
            )
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
                "data likely lost in frame {}".format(self.parsed_frames)
            )

    def run(self):
        """
        Producer thread method. Reads from the socket and
        sends a 'pause' signal if needed.
        """
        self._send_go_packet()
        self.running = True
        while self.running:
            self._parse_and_insert_frame()

            # If buffer is more than 50% full, pause the simulation
            if (
                not self.paused
                and self._full_frames >= self._buffer.capacity // 2
            ):
                self._pause_simulation()

            # Simulation will unpause during frame parsing if
            # 1. buffer is empty
            # 2. socket is empty

            # Only check for simulation end if we're not paused
            # otherwise we can get a false positive
            if not self.paused:
                self._check_end_of_simulation()

        return

    def stop(self):
        # Tell reader not to expect more frames to be added
        self._buffer.producer_finished = True
        # MUST disconnect before stopping run loop
        self._ensure_disconnect()
        self.running = False

    def _handle_signal(self, signum, frame):
        """Handle SIGINT and SIGTERM signals."""
        self.running = False
        self._ensure_disconnect()

    def _recv_data(self, bytearray):
        """Used to receive headers and data packets from the socket
        into self._frame_buf at a specified offset.

        This method will behave differently is self.paused is True,
        since having no more data available in the socket doesn't
        indicate the connection is closed, but rather that
        the simulation is now ready to be unpaused
        """
        logger.debug(f"IMDProducer: Receiving {len(bytearray)} bytes")

        memview = memoryview(bytearray)
        total_received = 0
        while total_received < len(bytearray):
            try:
                chunk = self._conn.recv(len(bytearray) - total_received)
                if not chunk:
                    self._is_disconnected = True
                    self.stop()
                    logger.warning(
                        "IMDProducer: Data likely lost in frame {}".format(
                            self.parsed_frames
                        )
                    )
                    raise ConnectionError("Socket connection was closed")
            except BlockingIOError:
                # If we're paused, we can assume the simulation is ready to unpause
                if self.paused:
                    self._unpause_simulation()
                    continue
            except Exception as e:
                logger.error(
                    f"IMDProducer: Error receiving data: {e}. Stopping producer"
                )
            memview[total_received : total_received + len(chunk)] = chunk
            total_received += len(chunk)
            logger.debug(
                f"IMDProducer: Receiving data. Total received: {total_received}"
            )

    def _expect_header(self, expected_type=None, expected_value=None):
        """
        Read a header packet from the socket.
        """
        header_bytes = bytearray(IMDHEADERSIZE)
        self._recv_data(header_bytes)
        header = parse_header_bytes(header_bytes)
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
            except ConnectionResetError:
                logger.debug(
                    f"IMDProducer: Server already terminated the connection"
                )


class IMDv2Producer(AbstractIMDProducer):
    """
    Producer thread for IMDReader. Reads packets from the socket
    and places them into the shared buffer.
    """

    def __init__(self, conn, buffer, imdsinfo, n_atoms):
        super(IMDv2Producer, self).__init__(conn, buffer, imdsinfo, n_atoms)
        self._data_elements = {"positions", "energies"}

    def _parse_and_insert_frame(self):
        with timeit() as parse_frame:
            self._expect_header(
                expected_type=IMDType.IMD_ENERGIES, expected_value=1
            )
            logger.debug("IMDProducer: Received energies packet")
            self._recv_data(self._byte_dict["energies"])
            self._expect_header(
                expected_type=IMDType.IMD_FCOORDS,
                expected_value=self.n_atoms,
            )
            logger.debug("IMDProducer: Received positions packet")
            self._recv_data(self._byte_dict["positions"])
            self._full_frames = self._buffer.insert(
                self._byte_dict, self._data_elements
            )

        self.parsed_frames += 1
        # Use the longest parse frame time to calculate
        # the timeout for the next frame
        self._parse_frame_time = max(
            parse_frame.elapsed, self._parse_frame_time
        )
        logger.debug(
            f"IMDProducer: Added frame #{self.parsed_frames - 1} to buffer in {parse_frame.elapsed} seconds"
        )

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
                    "{} due to closed connection".format(self.parsed_frames - 1)
                )
                self.running = False
                self._buffer.producer_finished = True
            # case 2: server has not closed the connection and we have no more data to read
        except socket.timeout:
            logger.debug(
                "IMDProducer: Assuming simulation is over at frame "
                "#{} due to read timeout,".format(self.parsed_frames - 1)
            )
            self.running = False
            self._buffer.producer_finished = True

        self._conn.settimeout(None)


class IMDv3Producer(AbstractIMDProducer):
    def _parse_and_insert_frame(self):
        with timeit() as parse_frame:

            data_elements = set()
            if (
                "energies" in self._byte_dict
                and (self.parsed_frames + 1 % self.imdsinfo.energies) == 0
            ):
                data_elements.add("energies")
                self._expect_header(
                    expected_type=IMDType.IMD_ENERGIES, expected_value=1
                )
                self._recv_data(self._byte_dict["energies"])

            if (
                "dimensions" in self._byte_dict
                and (self.parsed_frames + 1 % self.imdsinfo.dimensions) == 0
            ):
                data_elements.add("dimensions")
                self._expect_header(
                    expected_type=IMDType.IMD_BOX,
                )

            if (
                "positions" in self._byte_dict
                and (self.parsed_frames + 1 % self.imdsinfo.positions) == 0
            ):
                data_elements.add("positions")
                self._expect_header(
                    expected_type=IMDType.IMD_FCOORDS,
                    expected_value=self.n_atoms,
                )
                self._recv_data(self._byte_dict["positions"])

            if (
                "velocities" in self._byte_dict
                and (self.parsed_frames + 1 % self.imdsinfo.velocities) == 0
            ):
                data_elements.add("velocities")
                self._expect_header(
                    expected_type=IMDType.IMD_VELS,
                    expected_value=self.n_atoms,
                )

            if (
                "forces" in self._byte_dict
                and (self.parsed_frames + 1 % self.imdsinfo.forces) == 0
            ):
                data_elements.add("forces")
                self._expect_header(
                    expected_type=IMDType.IMD_FORCES,
                    expected_value=self.n_atoms,
                )

            self._full_frames = self._buffer.insert(
                self._byte_dict, data_elements
            )

        self.parsed_frames += 1
        # Use the longest parse frame time to calculate
        # the timeout for the next frame
        self._parse_frame_time = max(
            parse_frame.elapsed, self._parse_frame_time
        )
        logger.debug(
            f"IMDProducer: Added frame #{self.parsed_frames - 1} to buffer in {parse_frame.elapsed} seconds"
        )

    def _check_end_of_simulation(self):
        # Peek for an EOS Header packet
        header_bytes = self._conn.recv(8, socket.MSG_PEEK)
        header = parse_header_bytes(header_bytes)
        if header.type == IMDType.IMD_EOS:
            logger.debug(
                "IMDProducer: Received end of simulation packet, stopping producer"
            )
            self.running = False
            self._buffer.producer_finished = True


class CircularByteBuf:
    """
    Acts as interface between producer and consumer threads
    """

    # NOTE: Use 1 buffer for pos, vel, force rather than 3
    def __init__(self, buffer_size, n_atoms, timestep, imdsinfo):
        self._n_atoms = n_atoms
        self._buffer_size = buffer_size
        # Syncing reader and producer
        self._producer_finished = False

        self._ts = timestep
        self.imdsinfo = imdsinfo

        self._frame_size = 1
        # framesize is the number of bytes needed to hold 1 simulation frame in the buffer
        # one byte is used to hold flags for the type of data in the frame

        # Even though every simulation step might not contain all data,
        # we allocate space for all data in every frame. This is memory-innefficient,
        # but works for the majority of use cases, so we can optimize later if needed
        # this could be optimized by creating a timestep buffer, some kind of byte queue,
        # or a circular buffer with wrapped inserts and reads
        if self.imdsinfo.energies > 0:
            self._frame_size += 40
        if self.imdsinfo.dimensions > 0:
            self._frame_size += 36
        if self.imdsinfo.positions > 0:
            self._frame_size += self._n_atoms * 12
        if self.imdsinfo.velocities > 0:
            self._frame_size += self._n_atoms * 12
        if self.imdsinfo.forces > 0:
            self._frame_size += self._n_atoms * 12

        self._ts_buf = bytearray(self._frame_size)
        self._data_bytes = self._n_atoms * 12

        # offsets for each data element in the frame
        self._energy_offset = 1
        self._dim_offset = 41
        self._pos_offset = 77
        self._vel_offset = self._pos_offset + (self._data_bytes)
        self._force_offset = self._vel_offset + (self._data_bytes)

        if self._frame_size > self._buffer_size:
            raise MemoryError(
                f"Requested buffer size of {self._buffer_size} "
                + f"doesn't meet memory requirement of {self._frame_size}"
            )

        self._buf_frames = (self._buffer_size) // self._frame_size
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

    def insert(self, byte_dict, data_elements):
        """byte_dict is a dictionary of bytearrays
        where the keys are the names of the data elements
        and the values are the bytearrays of the data elements for that frame

        data_elements is a set of the names of the data elements that are present in the frame
        """
        with self._not_full:
            while not self._empty:
                self._not_full.wait()

            flags = 0
            self._memview[self._fill] = flags
            self._fill += 1
            if "energies" in data_elements:
                self._memview[
                    self._fill : self._fill + IMDENERGYPACKETLENGTH
                ] = byte_dict["energies"][:]
                self._fill += IMDENERGYPACKETLENGTH
                flags |= 1 << 4

            if "dimensions" in data_elements:
                self._memview[self._fill : self._fill + 36] = byte_dict[
                    "dimensions"
                ][:]
                self._fill += 36
                flags |= 1 << 3

            if "positions" in data_elements:
                self._memview[self._fill : self._fill + self._data_bytes] = (
                    byte_dict["positions"][:]
                )
                self._fill += self._data_bytes
                flags |= 1 << 2

            if "velocities" in data_elements:
                self._memview[self._fill : self._fill + self._data_bytes] = (
                    byte_dict["velocities"][:]
                )
                self._fill += self._data_bytes
                flags |= 1 << 1

            if "forces" in data_elements:
                self._memview[self._fill : self._fill + self._data_bytes] = (
                    byte_dict["forces"][:]
                )
                self._fill += self._data_bytes
                flags |= 1

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

        self._frame += 1

        with self._not_empty:
            while not self._full and not self.producer_finished:
                self._not_empty.wait()

            # Buffer is responsible for stopping iteration
            if self.producer_finished and not self._full:
                raise IOError from None

            # quickly unload memory to free mutex and later parse
            self._ts_buf[:] = self._memview[
                self._use : self._use + self._frame_size
            ]
            self._use = (self._use + self._frame_size) % len(self._buf)

            self._full -= 1
            self._empty += 1
            self._not_full.notify()

        flags = self._ts_buf[0]

        if flags >> 4 & 1:
            self._ts.data["step"] = np.frombuffer(
                self._ts_buf,
                dtype=f"{self.imdinfo.endianness}i4",
                offset=self._energy_offset,
            )[0]
            energydata = (
                "temperature",
                "total_energy",
                "potential_energy",
                "van_der_walls_energy",
                "coulomb_energy",
                "bonds_energy",
                "angles_energy",
                "dihedrals_energy",
                "improper_dihedrals_energy",
            )
            for i, name in enumerate(energydata):
                self._ts.data[name] = np.frombuffer(
                    self._ts_buf,
                    dtype=f"{self.imdinfo.endianness}f4",
                    offset=(self._energy_offset + 4 + (i * 4)),
                    count=1,
                )[0]

        if flags >> 3 & 1:
            dim = np.frombuffer(
                self._ts_buf, dtype=f"{self._end}f4", offset=self._dim_offset
            ).reshape((3, 3))
            self._ts.dimensions = core.triclinic_box(*dim)
        else:
            self._ts.dimensions = None

        if flags >> 2 & 1:
            self._ts.positions = np.frombuffer(
                self._ts_buf, dtype=f"{self._end}f4", offset=self._pos_offset
            ).reshape((self._n_atoms, 3))
        else:
            self._ts.has_positions = False

        if flags >> 1 & 1:
            self._ts.velocities = np.frombuffer(
                self._ts_buf, dtype=f"{self._end}f4", offset=self._vel_offset
            ).reshape((self._n_atoms, 3))
        else:
            self._ts.has_velocities = False

        if flags & 1:
            self._ts.forces = np.frombuffer(
                self._ts_buf, dtype=f"{self._end}f4", offset=self._force_offset
            ).reshape((self._n_atoms, 3))
        else:
            self._ts.has_forces = False

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
