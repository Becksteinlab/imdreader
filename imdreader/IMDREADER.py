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
from .IMDClient import *
from .util import *
import socket
import threading
import numpy as np
import signal
import logging
import time
import select
import warnings

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
        **kwargs,
    ):
        """
        Parameters
        ----------
        filename : a string of the form "host:port" where host is the hostname
            or IP address of the listening GROMACS server and port
            is the port number.
        convert_units : bool (optional)
            convert units to MDAnalysis units [``True``]
        n_atoms : int (optional)
            number of atoms in the system. defaults to number of atoms
            in the topology. don't set this unless you know what you're doing.
        """
        self._producer = None
        super(IMDReader, self).__init__(filename, **kwargs)

        logger.debug("Reader initializing")

        if n_atoms is None:
            raise ValueError("IMDReader: n_atoms must be specified")
        self.n_atoms = n_atoms

        host, port = parse_host_port(filename)

        # This starts the simulation
        self._imdclient = IMDClient(host, port, n_atoms, **kwargs)

        self._imdsinfo = self._imdclient.get_imdsessioninfo()

        self.convert_units = convert_units
        # NOTE: changme after deciding how imdreader will handle units
        self.units = {
            "time": "ps",
            "length": "nm",
            "force": "kJ/(mol*nm)",
            "velocity": "nm/ps",
        }
        self.ts = self._Timestep(
            self.n_atoms,
            positions=(self._imdsinfo.positions > 0),
            velocities=(self._imdsinfo.velocities > 0),
            forces=(self._imdsinfo.forces > 0),
            **self._ts_kwargs,
        )

        self._frame = -1
        self._init_scope = True
        self._reopen_called = False
        self._read_next_timestep()

    def _read_next_timestep(self):
        # No rewinding- to both load the first frame on __init__
        # and access it during iteration, we need to store first ts in mem
        if not self._init_scope and self._frame == -1:
            self._frame += 1
            return self.ts

        return self._read_frame(self._frame + 1)

    def _read_frame(self, frame):

        try:
            imdf = self._imdclient.get_imdframe()
        except EOFError:
            # Not strictly necessary, but for clarity
            raise StopIteration

        self._load_imdframe_into_ts(imdf)

        if self.convert_units:
            self._convert_units()

        self._frame = frame

        if self._init_scope:
            self._init_scope = False

        logger.debug(f"IMDReader: Loaded frame {self._frame}")
        return self.ts

    def _load_imdframe_into_ts(self, imdf):
        self.ts.frame = self._frame
        # NOTE: need time.
        if imdf.energies is not None:
            self.ts.data.update(imdf.energies)
        if imdf.dimensions is not None:
            self.ts.dimensions = core.triclinic_box(*imdf.dimensions)
        if imdf.positions is not None:
            self.ts.positions = imdf.positions
        if imdf.velocities is not None:
            self.ts.velocities = imdf.velocities
        if imdf.forces is not None:
            self.ts.forces = imdf.forces

    def _convert_units(self):
        """converts time, position, velocity, and force values if they
        are not given in MDAnalysis standard units
        """

        self.ts.time = self.convert_time_from_native(self.ts.time)

        if self.ts.dimensions is not None:
            self.convert_pos_from_native(self.ts.dimensions[:3])

        if self.ts.has_positions:
            self.convert_pos_from_native(self.ts.positions)

        if self.ts.has_velocities:
            self.convert_velocities_from_native(self.ts.velocities)

        if self.ts.has_forces:
            self.convert_forces_from_native(self.ts.forces)

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
        # Don't stop client if only first ts was loaded
        # from __init__
        if self._init_scope:
            return
        self._imdclient.stop()
        # NOTE: removeme after testing
        logger.debug("IMDReader shut down gracefully.")

    # Incompatible methods
    def copy(self):
        raise NotImplementedError("IMDReader does not support copying")

    def _reopen(self):
        if self._reopen_called:
            raise RuntimeError("IMDReader: Cannot reopen IMD stream")
        self._frame = -1
        self._reopen_called = True

    def __getitem__(self, frame):
        """This method from ProtoReader must be overridden
        to prevent slicing that doesn't make sense in a stream.
        """
        raise RuntimeError("IMDReader: Trajectory can only be read in for loop")

    # NOTE: prevent auxiliary iteration methods from being called


# class IMDProducer(threading.Thread):

#     def __init__(
#         self,
#         conn,
#         buffer,
#         sinfo,
#         n_atoms,
#         pause_empty_proportion,
#         unpause_empty_proportion,
#     ):
#         super(IMDProducer, self).__init__()
#         self._conn = conn
#         self.sinfo = sinfo
#         self._paused = False

#         self._pause_empty_proportion = pause_empty_proportion
#         self._unpause_empty_proportion = unpause_empty_proportion

#         # Timeout for first frame should be longer
#         # than rest of frames
#         self._timeout = 5
#         self._conn.settimeout(self._timeout)

#         self._buffer = buffer

#         self._frame = 0
#         self._parse_frame_time = 0

#         # The body of an x/v/f packet should contain
#         # (4 bytes per float * 3 atoms * n_atoms) bytes
#         self.n_atoms = n_atoms
#         xvf_bytes = 12 * n_atoms

#         self._header = bytearray(IMDHEADERSIZE)
#         if self.sinfo.energies > 0:
#             self._energies = bytearray(40)
#         if self.sinfo.dimensions > 0:
#             self._dimensions = bytearray(36)
#         if self.sinfo.positions > 0:
#             self._positions = bytearray(xvf_bytes)
#         if self.sinfo.velocities > 0:
#             self._velocities = bytearray(xvf_bytes)
#         if self.sinfo.forces > 0:
#             self._forces = bytearray(xvf_bytes)

#     def _go(self):
#         """
#         Send a go packet to the client to start the simulation
#         and begin receiving data.
#         """
#         # NOTE: removeme after testing
#         print("sending go packet...")
#         go = create_header_bytes(IMDHeaderType.IMD_GO, 0)
#         self._conn.sendall(go)
#         logger.debug("IMDProducer: Sent go packet to server")

#     def _pause(self):
#         """
#         Block the simulation until the buffer has more space.
#         """
#         self._conn.settimeout(0)
#         logger.debug(
#             "IMDProducer: Pausing simulation because buffer is almost full"
#         )
#         pause = create_header_bytes(IMDHeaderType.IMD_PAUSE, 0)
#         try:
#             self._conn.sendall(pause)
#         except ConnectionResetError as e:
#             # Simulation has already ended by the time we paused
#             raise IndexError
#         # Edge case: pause occured in the time between server sends its last frame
#         # and closing socket
#         # Simulation is not actually paused but is over, but we still want to read remaining data
#         # from the socket

#     def _unpause(self):
#         self._conn.settimeout(self._timeout)
#         logger.debug("IMDProducer: Unpausing simulation, buffer has space")
#         unpause = create_header_bytes(IMDHeaderType.IMD_PAUSE, 0)
#         try:
#             self._conn.sendall(unpause)
#         except ConnectionResetError as e:
#             # Edge case: pause occured in the time between server sends its last frame
#             # and closing socket
#             # Simulation was never actually paused in this case and is now over
#             raise IndexError
#         # Edge case: pause & unpause occured in the time between server sends its last frame and closing socket
#         # in this case, the simulation isn't actually unpaused but over

#     def run(self):
#         self._go()

#         try:
#             while True:

#                 logger.debug(f"IMDProducer: Got timstep")

#                 # This value is approximate, doesn't acquire lock
#                 n_empty_ts = self._buffer.get_empty_qsize()

#                 logger.debug(f"IMDProducer: Got empty qsize {n_empty_ts}")

#                 logger.debug(
#                     f"IMDProducer: {n_empty_ts} // {self._buffer.n_ts} = {n_empty_ts // self._buffer.n_ts}"
#                 )
#                 logger.debug(f"IMDProducer: {self._pause_empty_proportion}")
#                 # If buffer is more than 50% full, pause the simulation
#                 if (
#                     not self._paused
#                     and n_empty_ts / self._buffer.n_ts
#                     <= self._pause_empty_proportion
#                 ):
#                     # if pause succeeds, simulation may still have ended
#                     self._pause()
#                     self._paused = True

#                 if self._paused:
#                     # If buffer is less than 25% full, unpause the simulation
#                     if (
#                         n_empty_ts / self._buffer.n_ts
#                         >= self._unpause_empty_proportion
#                     ):
#                         self._unpause()
#                         self._paused = False
#                     # If the buffer is still full but we've run out
#                     # of frames to read, wait until the buffer is less full
#                     elif not sock_contains_data(self._conn, 0):
#                         self._buffer.wait_for_space(
#                             self._unpause_empty_proportion
#                         )
#                         self._unpause()
#                         self._paused = False

#                 logger.debug(f"IMDProducer: Attempting to get timestep")
#                 ts = self._buffer.get_timestep()

#                 logger.debug(f"IMDProducer: Attempting to read nrg and pos")
#                 # NOTE: This can be replaced with a simple parser if
#                 # the server doesn't send the final frame with all data
#                 # as in xtc
#                 self._expect_header(
#                     IMDHeaderType.IMD_ENERGIES, expected_value=1
#                 )
#                 logger.debug(f"IMDProducer: Expected header")

#                 read_into_buf(self._conn, self._energies)

#                 self._load_energies(ts)
#                 logger.debug(f"IMDProducer: Read nrg, reading pos")

#                 self._expect_header(
#                     IMDHeaderType.IMD_FCOORDS, expected_value=self.n_atoms
#                 )

#                 logger.debug(f"IMDProducer: Expected header")
#                 read_into_buf(self._conn, self._positions)

#                 logger.debug(f"IMDProducer: attempting to load ts")

#                 ts.frame = self._frame
#                 ts.positions = np.frombuffer(
#                     self._positions, dtype=f"{self.sinfo.endianness}f"
#                 ).reshape((self.n_atoms, 3))

#                 logger.debug(f"IMDProducer: ts loaded- inserting it")

#                 self._buffer.insert(ts)

#                 logger.debug(f"IMDProducer: ts inserted")

#                 if self._frame == 0:
#                     self._conn.settimeout(1)

#                 self._frame += 1
#         except IndexError:
#             # Don't raise error if simulation ended in a way
#             # that we expected
#             pass
#         finally:

#             logger.debug("IMDProducer: simluation ended")

#             # Tell reader not to expect more frames to be added
#             self._buffer.notify_producer_finished()
#             # MUST disconnect before stopping run loop
#             # if simulation already ended, this method will do nothing
#             self._disconnect()

#             return

#     def _expect_header(self, expected_type, expected_value=None):

#         read_into_buf(self._conn, self._header)

#         logger.debug(f"IMDProducer: header: {self._header}")
#         header = IMDHeader(self._header)

#         logger.debug(f"IMDProducer: header parsed")

#         if header.type != expected_type:
#             raise RuntimeError

#         if expected_value is not None and header.length != expected_value:
#             raise RuntimeError

#     def _disconnect(self):
#         try:
#             disconnect = create_header_bytes(IMDHeaderType.IMD_DISCONNECT, 0)
#             self._conn.sendall(disconnect)
#             logger.debug("IMDProducer: Disconnected from server")
#         except (ConnectionResetError, BrokenPipeError):
#             logger.debug(
#                 f"IMDProducer: Attempted to disconnect but server already terminated the connection"
#             )
#         finally:
#             self._conn.close()

#     def _load_energies(self, ts):

#         energy_dict = IMDEnergyPacket(
#             self._energies, self.sinfo.endianness
#         ).data
#         logger.debug(f"IMDProducer: Loaded energies {energy_dict}")
#         ts.data.update(energy_dict)
#         logger.debug(f"IMDProducer: Updated ts with energies")


# class TimestepBuffer:
#     """
#     Acts as interface between producer and consumer threads
#     """

#     def __init__(self, buffer_size, imdsinfo, ts_class, n_atoms, ts_kwargs):

#         # Syncing reader and producer
#         self._producer_finished = False
#         self._consumer_finished = False

#         self._prev_empty_ts = None

#         self._empty_q = queue.Queue()
#         self._full_q = queue.Queue()
#         self._empty_ts_avail = threading.Condition(threading.Lock())
#         self._full_ts_avail = threading.Condition(threading.Lock())

#         # Allocate timesteps with all of xvf present in imdsinfo
#         # even if they aren't sent every frame. Can be optimized if needed
#         ts_memsize = approximate_timestep_memsize(
#             n_atoms,
#             (imdsinfo.energies > 0),
#             (imdsinfo.dimensions > 0),
#             (imdsinfo.positions > 0),
#             (imdsinfo.velocities > 0),
#             (imdsinfo.forces > 0),
#         )
#         self._total_ts = buffer_size // ts_memsize
#         logger.debug(
#             f"Timestepbuffer: Total timesteps allocated: {self._total_ts}"
#         )
#         for i in range(self._total_ts):
#             self._empty_q.put(
#                 ts_class(
#                     n_atoms,
#                     positions=(imdsinfo.positions > 0),
#                     velocities=(imdsinfo.velocities > 0),
#                     forces=(imdsinfo.forces > 0),
#                     **ts_kwargs,
#                 )
#             )

#         # Timing for analysis
#         self._t1 = None
#         self._t2 = None
#         self._start = True
#         self._analyze_frame_time = None

#         self._frame = 0

#     def get_empty_qsize(self):
#         return self._empty_q.qsize()

#     def get_timestep(self):
#         with self._empty_ts_avail:
#             while self._empty_q.qsize() == 0 and not self._consumer_finished:
#                 self._empty_ts_avail.wait()

#         if self._consumer_finished:
#             raise IndexError

#         ts = self._empty_q.get()

#         return ts

#     def wait_for_space(self, unpause_empty_proportion):
#         with self._empty_ts_avail:
#             while (
#                 self._empty_q.qsize() / self.n_ts < unpause_empty_proportion
#             ) and not self._consumer_finished:
#                 self._empty_ts_avail.wait()

#         if self._consumer_finished:
#             raise IndexError

#     def insert(self, ts):
#         self._full_q.put(ts)
#         with self._full_ts_avail:
#             self._full_ts_avail.notify()

#     def consume_next_timestep(self):
#         """Put empty_ts in the empty_q and get the next full timestep"""
#         # Start timer- one frame of analysis is starting (including removal
#         # from buffer)

#         self._t1 = self._t2
#         self._t2 = time.time()
#         if self._t1 is not None:
#             logger.debug(
#                 f"IMDReader: Frame #{self._frame} analyzed in {self._t2 - self._t1} seconds"
#             )
#             self._analyze_frame_time = self._t2 - self._t1

#         self._frame += 1

#         # Return the processed timestep
#         if self._prev_empty_ts is not None:
#             self._empty_q.put(self._prev_empty_ts)
#             with self._empty_ts_avail:
#                 self._empty_ts_avail.notify()

#         # Get the next timestep
#         with self._full_ts_avail:
#             while self._full_q.qsize() == 0 and not self._producer_finished:
#                 self._full_ts_avail.wait()

#         # Buffer is responsible for stopping iteration
#         if self._producer_finished and self._full_q.qsize() == 0:
#             raise StopIteration from None

#         ts = self._full_q.get()

#         self._prev_empty_ts = ts

#         return ts

#     def notify_producer_finished(self):
#         self._producer_finished = True
#         with self._full_ts_avail:
#             self._full_ts_avail.notify()

#     def notify_consumer_finished(self):
#         self._consumer_finished = True
#         with self._empty_ts_avail:
#             # noop if producer isn't waiting
#             self._empty_ts_avail.notify()

#     @property
#     def analyze_frame_time(self):
#         if self._analyze_frame_time is not None:
#             return self._analyze_frame_time
#         else:
#             return None

#     @property
#     def n_ts(self):
#         return self._total_ts


# def read_into_buf(sock, buf) -> bool:
#     """Receives len(buf) bytes into buf from the socket sock"""
#     view = memoryview(buf)
#     total_received = 0
#     while total_received < len(view):
#         try:
#             received = sock.recv_into(view[total_received:])
#             if received == 0:
#                 # Server called close()
#                 logger.debug(
#                     "IMDProducer: recv excepting due to server calling close()"
#                 )
#                 raise IndexError
#         except TimeoutError:
#             # Server is *likely* done sending frames
#             logger.debug("IMDProducer: recv excepting due to timeout")
#             raise IndexError
#         except BlockingIOError:
#             # Server is done sending frames
#             logger.debug("IMDProducer: recv excepting due to blocking")
#             raise IndexError
#         total_received += received
#     return True


# def sock_contains_data(sock, timeout) -> bool:
#     ready_to_read, ready_to_write, in_error = select.select(
#         [sock], [], [], timeout
#     )
#     return sock in ready_to_read
