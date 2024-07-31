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

import queue
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

        self.n_atoms = n_atoms
        logger.debug(f"IMDReader: n_atoms: {self.n_atoms}")

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
                self._buffer = CircularByteBuf(
                    self._buffer_size, imdsinfo, self._Timestep, self.n_atoms
                )
                self._producer = IMDProducer(
                    conn,
                    self._buffer,
                    imdsinfo,
                    self.n_atoms,
                )
            # Producer responsible for sending go packet
            self._producer.start()

        return self._read_frame(self._frame + 1)

    def _read_frame(self, frame):

        self._ts = self._buffer.consume_next_timestep()

        logger.debug(f"self._ts.positions: {self._ts.positions}")
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

        h_buf = bytearray(IMDHEADERSIZE)
        handshake_recieved = read_into_buf(conn, h_buf)
        if not handshake_recieved:
            raise ConnectionError("IMDReader: No handshake received.")

        header = IMDHeader(h_buf)

        if header.type != IMDHeaderType.IMD_HANDSHAKE:
            raise ValueError(
                f"Expected header type `IMD_HANDSHAKE`, got {header.type}"
            )

        if header.length not in IMDVERSIONS:
            # Try swapping endianness
            swapped = struct.unpack("<i", struct.pack(">i", header.length))[0]
            if swapped not in IMDVERSIONS:
                err_version = min(swapped, header.length)
                # NOTE: Add test for this
                raise ValueError(
                    f"Incompatible IMD version. Expected version in {IMDVERSIONS}, got {err_version}"
                )
            else:
                end = "<"
                ver = swapped
        else:
            ver = header.length

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
        elif ver == 3:
            sinfo = parse_imdv3_session_info(conn, end)

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
            self._buffer.notify_consumer_finished()
            # NOTE: is join necessary here?
            # self._producer.join()
        # NOTE: removeme after testing
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

    # NOTE: prevent auxiliary iteration methods from being called


class IMDProducer(threading.Thread):

    def __init__(self, conn, buffer, sinfo, n_atoms):
        super(IMDProducer, self).__init__()
        self._conn = conn
        self.sinfo = sinfo
        self._should_stop = False
        self._paused = False

        self._buffer = buffer

        self._frame = 0
        self._parse_frame_time = 0

        # The body of an x/v/f packet should contain
        # (4 bytes per float * 3 atoms * n_atoms) bytes
        self.n_atoms = n_atoms
        xvf_bytes = 12 * n_atoms

        # NOTE: does 0 mean every frame?
        self._header = bytearray(IMDHEADERSIZE)
        if self.sinfo.energies > 0:
            self._energies = bytearray(40)
        if self.sinfo.dimensions > 0:
            self._dimensions = bytearray(36)
        if self.sinfo.positions > 0:
            self._positions = bytearray(xvf_bytes)
        if self.sinfo.velocities > 0:
            self._velocities = bytearray(xvf_bytes)
        if self.sinfo.forces > 0:
            self._forces = bytearray(xvf_bytes)

    def _go(self):
        """
        Send a go packet to the client to start the simulation
        and begin receiving data.
        """
        # NOTE: removeme after testing
        print("sending go packet...")
        go = create_header_bytes(IMDHeaderType.IMD_GO, 0)
        self._conn.sendall(go)
        logger.debug("IMDProducer: Sent go packet to server")

    def _pause(self):
        """
        Block the simulation until the buffer has more space.
        """
        self._conn.settimeout(0)
        pause = create_header_bytes(IMDHeaderType.IMD_PAUSE, 0)
        try:
            self._conn.sendall(pause)
        except ConnectionResetError as e:
            # NOTE: test to ensure this is the correct error type to except
            # Simulation has already ended by the time we paused
            return False
        return True

    def _unpause(self):
        self._conn.settimeout(None)
        unpause = create_header_bytes(IMDHeaderType.IMD_PAUSE, 0)
        try:
            self._conn.sendall(unpause)
        except ConnectionResetError as e:
            # NOTE: test to ensure this is the correct error type to except
            # This shouldn't ever happen
            logger.debug("Unpause failed, connection reset")
            return False
        return True

    def run(self):
        self._go()

        while not self._should_stop:
            # NOTE: imdv3 may have a simulation end packet rather than relying on disconnect
            if sock_disconnected(self._conn):
                break
            # NOTE: add timeout

            logger.debug(f"IMDProducer: Attempting to get timestep")
            ts = self._buffer.get_timestep()
            # Reader is closed
            if ts is None:
                break

            logger.debug(f"IMDProducer: Got timstep")

            # This value is approximate
            n_empty_ts = self._buffer.get_empty_qsize()

            logger.debug(f"IMDProducer: Got empty qsize")

            # If buffer is more than 50% full, pause the simulation
            if not self._paused and n_empty_ts < self._buffer.capacity // 2:
                pause_success = self._pause()
                if not pause_success:
                    # Don't need to disconnect, pause failed due to simulation end
                    break
                self._paused = True

            # If buffer is less than 25% full, unpause the simulation
            if (
                self._paused
                and sock_empty(self._conn)
                and n_empty_ts >= self._buffer.capacity // 4
            ):
                unpause_success = self._unpause()
                if not unpause_success:
                    break
                self._paused = False

            logger.debug(f"IMDProducer: Attempting to read nrg and pos")
            # NOTE: This can be replaced with a simple parser if
            # the server doesn't send the final frame with all data
            e_header_success = self._expect_header(
                IMDHeaderType.IMD_ENERGIES, expected_value=1
            )
            if not e_header_success:
                break
            logger.debug(f"IMDProducer: Expected header")

            energies_successs = read_into_buf(self._conn, self._energies)
            if not energies_successs:
                break

            logger.debug(f"IMDProducer: Read nrg, reading pos")

            p_header_success = self._expect_header(
                IMDHeaderType.IMD_FCOORDS, expected_value=self.n_atoms
            )
            if not p_header_success:
                break
            logger.debug(f"IMDProducer: Expected header")
            positions_success = read_into_buf(self._conn, self._positions)

            if not positions_success:
                break

            logger.debug(f"IMDProducer: attempting to load ts")

            ts.frame = self._frame
            ts.positions = np.frombuffer(
                self._positions, dtype=f"{self.sinfo.endianness}f4"
            ).reshape((self.n_atoms, 3))

            logger.debug(f"IMDProducer: ts loaded- inserting it")

            self._buffer.insert(ts)

            logger.debug(f"IMDProducer: ts inserted")

        logger.debug("IMDProducer: break occuurred")

        # Tell reader not to expect more frames to be added
        self._buffer.notify_producer_finished()
        # MUST disconnect before stopping run loop
        # if simulation already ended, this method will do nothing
        self._disconnect()

        return

    def _expect_header(self, expected_type, expected_value=None):

        recv_success = read_into_buf(self._conn, self._header)
        logger.debug(f"IMDProducer: recv success: {recv_success}")
        if not recv_success:
            return False

        logger.debug(f"IMDProducer: header: {self._header}")
        header = IMDHeader(self._header)

        logger.debug(f"IMDProducer: header parsed")

        if header.type != expected_type:
            return False

        if expected_value is not None and header.length != expected_value:
            return False

        return True

    def _disconnect(self):
        if self._conn:
            try:
                disconnect = create_header_bytes(
                    IMDHeaderType.IMD_DISCONNECT, 0
                )
                self._conn.sendall(disconnect)
                self._conn.close()
                logger.debug("IMDProducer: Disconnected from server")
            # NOTE: verify this is the correct error type to except
            except ConnectionResetError:
                logger.debug(
                    f"IMDProducer: Server already terminated the connection"
                )


class CircularByteBuf:
    """
    Acts as interface between producer and consumer threads
    """

    # NOTE: Use 1 buffer for pos, vel, force rather than 3
    def __init__(self, buffer_size, imdsinfo, ts_class, n_atoms):
        self._buffer_size = buffer_size

        # Syncing reader and producer
        self._producer_finished = False
        self._consumer_finished = False

        self._prev_empty_ts = None

        self.imdsinfo = imdsinfo

        self._empty_q = queue.Queue()
        self._full_q = queue.Queue()

        self._full_ts_avail = threading.Condition(threading.Lock())
        self._empty_ts_avail = threading.Condition(threading.Lock())

        # NOTE: hardcoded for testing
        self._total_ts = 100
        for i in range(100):
            self._empty_q.put(ts_class(n_atoms, positions=True))

        # Timing for analysis
        self._t1 = None
        self._t2 = None
        self._start = True
        self._analyze_frame_time = None

        self._frame = 0

    def get_empty_qsize(self):
        return self._empty_q.qsize()

    def get_timestep(self):
        with self._empty_ts_avail:
            while self._empty_q.qsize() == 0 and not self._consumer_finished:
                self._empty_ts_avail.wait()

            if self._consumer_finished:
                return None

            ts = self._empty_q.get()
        return ts

    def insert(self, ts):
        self._full_q.put(ts)
        with self._full_ts_avail:
            self._full_ts_avail.notify()

    def consume_next_timestep(self):
        """Put empty_ts in the empty_q and get the next full timestep"""
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

        # Return the processed timestep
        if self._prev_empty_ts is not None:
            self._empty_q.put(self._prev_empty_ts)
            with self._empty_ts_avail:
                self._empty_ts_avail.notify()

        # Get the next timestep
        with self._full_ts_avail:
            while self._full_q.qsize() == 0 and not self._producer_finished:
                self._full_ts_avail.wait()

            # Buffer is responsible for stopping iteration
            if self._producer_finished and self._full_q.qsize() == 0:
                raise StopIteration from None

            ts = self._full_q.get()

        self._prev_empty_ts = ts

        return ts

    def notify_producer_finished(self):
        self._producer_finished = True
        with self._full_ts_avail:
            self._full_ts_avail.notify()

    def notify_consumer_finished(self):
        self._consumer_finished = True
        with self._empty_ts_avail:
            # noop if producer isn't waiting
            self._empty_ts_avail.notify()

    @property
    def analyze_frame_time(self):
        if self._analyze_frame_time is not None:
            return self._analyze_frame_time
        else:
            return None

    @property
    def capacity(self):
        return self._total_ts
