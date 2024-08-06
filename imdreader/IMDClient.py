import socket
import threading
from .IMDProtocol import *
import logging
import queue
import select
import time
import numpy as np
from typing import Union, Dict

logger = logging.getLogger("imdreader.IMDREADER")


class IMDClient:
    def __init__(
        self,
        host,
        port,
        n_atoms,
        socket_bufsize=None,
        buffer_size=(10 * 1024**2),
        pause_empty_proportion=0.25,
        unpause_empty_proportion=0.5,
        **kwargs,
    ):

        conn = self._connect_to_server(host, port, socket_bufsize)
        self._imdsinfo = self._await_IMD_handshake(conn)
        self._buf = IMDFrameBuffer(
            buffer_size,
            self._imdsinfo,
            n_atoms,
            pause_empty_proportion,
            unpause_empty_proportion,
        )
        producer = IMDProducer(
            conn,
            self._buf,
            self._imdsinfo,
            n_atoms,
        )
        self._go(conn)
        producer.start()

    def get_imdframe(self):
        return self._buf.pop_full_imdframe()

    def get_imdsessioninfo(self):
        return self._imdsinfo

    def stop(self):
        self._buf.notify_consumer_finished()

    def _connect_to_server(self, host, port, socket_bufsize):
        """
        Establish connection with the server, failing out if this
        does not occur within 5 seconds.
        """
        conn = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        if socket_bufsize is not None:
            conn.setsockopt(
                socket.SOL_SOCKET, socket.SO_RCVBUF, self._socket_bufsize
            )
        conn.settimeout(60)
        try:
            conn.connect((host, port))
        except ConnectionRefusedError:
            raise ConnectionRefusedError(
                f"IMDReader: Connection to {host}:{port} refused"
            )
        return conn

    def _await_IMD_handshake(self, conn) -> IMDSessionInfo:
        """
        Wait for the server to send a handshake packet, then determine
        IMD session information.
        """
        end = ">"
        ver = None

        h_buf = bytearray(IMDHEADERSIZE)
        try:
            read_into_buf(conn, h_buf)
        except IndexError:
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

    def _go(self, conn):
        """
        Send a go packet to the client to start the simulation
        and begin receiving data.
        """
        go = create_header_bytes(IMDHeaderType.IMD_GO, 0)
        conn.sendall(go)
        logger.debug("IMDProducer: Sent go packet to server")


class IMDProducer(threading.Thread):

    def __init__(
        self,
        conn,
        buffer,
        sinfo,
        n_atoms,
    ):
        super(IMDProducer, self).__init__()
        self._conn = conn
        self._imdsinfo = sinfo
        self._paused = False

        # Timeout for first frame should be longer
        # than rest of frames
        self._timeout = 5
        self._conn.settimeout(self._timeout)

        self._buf = buffer

        self._frame = 0
        self._parse_frame_time = 0

        # The body of an x/v/f packet should contain
        # (4 bytes per float * 3 atoms * n_atoms) bytes
        self._n_atoms = n_atoms
        xvf_bytes = 12 * n_atoms

        self._header = bytearray(IMDHEADERSIZE)
        if self._imdsinfo.energies > 0:
            self._energies = bytearray(40)
        if self._imdsinfo.dimensions > 0:
            self._dimensions = bytearray(36)
        if self._imdsinfo.positions > 0:
            self._positions = bytearray(xvf_bytes)
        if self._imdsinfo.velocities > 0:
            self._velocities = bytearray(xvf_bytes)
        if self._imdsinfo.forces > 0:
            self._forces = bytearray(xvf_bytes)

    def _pause(self):
        """
        Block the simulation until the buffer has more space.
        """
        self._conn.settimeout(0)
        logger.debug(
            "IMDProducer: Pausing simulation because buffer is almost full"
        )
        pause = create_header_bytes(IMDHeaderType.IMD_PAUSE, 0)
        try:
            self._conn.sendall(pause)
        except ConnectionResetError as e:
            # Simulation has already ended by the time we paused
            raise IndexError
        # Edge case: pause occured in the time between server sends its last frame
        # and closing socket
        # Simulation is not actually paused but is over, but we still want to read remaining data
        # from the socket

    def _unpause(self):
        self._conn.settimeout(self._timeout)
        logger.debug("IMDProducer: Unpausing simulation, buffer has space")
        unpause = create_header_bytes(IMDHeaderType.IMD_PAUSE, 0)
        try:
            self._conn.sendall(unpause)
        except ConnectionResetError as e:
            # Edge case: pause occured in the time between server sends its last frame
            # and closing socket
            # Simulation was never actually paused in this case and is now over
            raise IndexError
        # Edge case: pause & unpause occured in the time between server sends its last frame and closing socket
        # in this case, the simulation isn't actually unpaused but over

    def run(self):
        try:
            while True:
                if not self._paused:
                    if self._buf.is_full():
                        self._pause()
                        self._paused = True

                if self._paused:
                    # wait for socket to empty before unpausing
                    if not sock_contains_data(self._conn, 0):
                        self._buf.wait_for_space()
                        self._unpause()
                        self._paused = False

                logger.debug(f"IMDProducer: Attempting to get timestep")
                imdf = self._buf.pop_empty_imdframe()

                logger.debug(f"IMDProducer: Attempting to read nrg and pos")
                # NOTE: This can be replaced with a simple parser if
                # the server doesn't send the final frame with all data
                # as in xtc
                self._expect_header(
                    IMDHeaderType.IMD_ENERGIES, expected_value=1
                )
                read_into_buf(self._conn, self._energies)

                logger.debug("read energy data")
                imdf.energies.update(
                    parse_energy_bytes(
                        self._energies, self._imdsinfo.endianness
                    )
                )
                logger.debug(f"IMDProducer: added energies to {imdf.energies}")

                logger.debug(f"IMDProducer: added energies to imdf")

                self._expect_header(
                    IMDHeaderType.IMD_FCOORDS, expected_value=self._n_atoms
                )

                logger.debug(f"IMDProducer: Expected header")
                read_into_buf(self._conn, self._positions)

                logger.debug(f"IMDProducer: attempting to load ts")

                imdf.positions = np.frombuffer(
                    self._positions, dtype=f"{self._imdsinfo.endianness}f"
                ).reshape((self._n_atoms, 3))

                logger.debug(f"IMDProducer: ts loaded- inserting it")

                self._buf.push_full_imdframe(imdf)

                logger.debug(f"IMDProducer: ts inserted")

                self._frame += 1
        except EOFError:
            # Don't raise error if simulation ended in a way
            # that we expected
            # i.e. consumer stopped or read_into_buf didn't find
            # full token of data
            pass
        finally:

            logger.debug("IMDProducer: simulation ended")

            # Tell reader not to expect more frames to be added
            self._buf.notify_producer_finished()
            # MUST disconnect before stopping run loop
            # if simulation already ended, this method will do nothing
            self._disconnect()

            return

    def _expect_header(self, expected_type, expected_value=None):

        read_into_buf(self._conn, self._header)

        logger.debug(f"IMDProducer: header: {self._header}")
        header = IMDHeader(self._header)

        logger.debug(f"IMDProducer: header parsed")

        if header.type != expected_type:
            raise RuntimeError

        if expected_value is not None and header.length != expected_value:
            raise RuntimeError

    def _disconnect(self):
        try:
            disconnect = create_header_bytes(IMDHeaderType.IMD_DISCONNECT, 0)
            self._conn.sendall(disconnect)
            logger.debug("IMDProducer: Disconnected from server")
        except (ConnectionResetError, BrokenPipeError):
            logger.debug(
                f"IMDProducer: Attempted to disconnect but server already terminated the connection"
            )
        finally:
            self._conn.close()


class IMDFrameBuffer:
    """
    Acts as interface between producer and consumer threads
    """

    def __init__(
        self,
        buffer_size,
        imdsinfo,
        n_atoms,
        pause_empty_proportion,
        unpause_empty_proportion,
    ):

        # Syncing reader and producer
        self._producer_finished = False
        self._consumer_finished = False

        self._prev_empty_imdf = None

        self._empty_q = queue.Queue()
        self._full_q = queue.Queue()
        self._empty_imdf_avail = threading.Condition(threading.Lock())
        self._full_imdf_avail = threading.Condition(threading.Lock())

        if pause_empty_proportion < 0 or pause_empty_proportion > 1:
            raise ValueError("pause_empty_proportion must be between 0 and 1")
        self._pause_empty_proportion = pause_empty_proportion
        if unpause_empty_proportion < 0 or unpause_empty_proportion > 1:
            raise ValueError("unpause_empty_proportion must be between 0 and 1")
        self._unpause_empty_proportion = unpause_empty_proportion

        if buffer_size <= 0:
            raise ValueError("Buffer size must be positive")
        # Allocate IMDFrames with all of xvf present in imdsinfo
        # even if they aren't sent every frame. Can be optimized if needed
        imdf_memsize = imdframe_memsize(n_atoms, imdsinfo)
        self._total_imdf = buffer_size // imdf_memsize
        logger.debug(
            f"IMDFRAMEBuffer: Total timesteps allocated: {self._total_imdf}"
        )
        if self._total_imdf == 0:
            raise ValueError(
                "Buffer size is too small to hold a single IMDFrame"
            )
        for i in range(self._total_imdf):
            self._empty_q.put(IMDFrame(n_atoms, imdsinfo))

        # Timing for analysis
        self._t1 = None
        self._t2 = None

        self._frame = 0

    def is_full(self):
        if (
            self._empty_q.qsize() / self._total_imdf
            <= self._pause_empty_proportion
        ):
            return True
        return False

    def wait_for_space(self):
        with self._empty_imdf_avail:
            while (
                self._empty_q.qsize() / self._total_imdf
                < self._unpause_empty_proportion
            ) and not self._consumer_finished:
                self._empty_imdf_avail.wait()

        if self._consumer_finished:
            raise EOFError

    def pop_empty_imdframe(self):
        with self._empty_imdf_avail:
            while self._empty_q.qsize() == 0 and not self._consumer_finished:
                self._empty_imdf_avail.wait()

        if self._consumer_finished:
            raise EOFError

        imdf = self._empty_q.get()

        return imdf

    def push_full_imdframe(self, imdf):
        self._full_q.put(imdf)
        with self._full_imdf_avail:
            self._full_imdf_avail.notify()

    def pop_full_imdframe(self):
        """Put empty_ts in the empty_q and get the next full timestep"""
        # Start timer- one frame of analysis is starting (including removal
        # from buffer)

        self._t1 = self._t2
        self._t2 = time.time()
        if self._t1 is not None:
            logger.debug(
                f"IMDReader: Frame #{self._frame} analyzed in {self._t2 - self._t1} seconds"
            )

        self._frame += 1

        # Return the processed IMDFrame
        if self._prev_empty_imdf is not None:
            self._empty_q.put(self._prev_empty_imdf)
            with self._empty_imdf_avail:
                self._empty_imdf_avail.notify()

        # Get the next IMDFrame
        logger.debug("IMDReader: Attempting to get next frame")
        with self._full_imdf_avail:
            while self._full_q.qsize() == 0 and not self._producer_finished:
                self._full_imdf_avail.wait()

        if self._producer_finished and self._full_q.qsize() == 0:
            logger.debug("IMDReader: Producer finished")
            raise EOFError

        imdf = self._full_q.get()

        self._prev_empty_imdf = imdf

        logger.debug(f"IMDReader: Got frame {self._frame}")

        return imdf

    def notify_producer_finished(self):
        self._producer_finished = True
        with self._full_imdf_avail:
            self._full_imdf_avail.notify()

    def notify_consumer_finished(self):
        self._consumer_finished = True
        with self._empty_imdf_avail:
            # noop if producer isn't waiting
            self._empty_imdf_avail.notify()


class IMDFrame:
    def __init__(self, n_atoms, imdsinfo):
        if imdsinfo.energies > 0:
            self.energies = {
                "step": 0,
                "temperature": 0.0,
                "total_energy": 0.0,
                "potential_energy": 0.0,
                "van_der_walls_energy": 0.0,
                "coulomb_energy": 0.0,
                "bonds_energy": 0.0,
                "angles_energy": 0.0,
                "dihedrals_energy": 0.0,
                "improper_dihedrals_energy": 0.0,
            }
        else:
            self.energies = None
        if imdsinfo.dimensions > 0:
            self.dimensions = np.empty((3, 3), dtype=np.float32)
        else:
            self.dimensions = None
        if imdsinfo.positions > 0:
            self.positions = np.empty((n_atoms, 3), dtype=np.float32)
        else:
            self.positions = None
        if imdsinfo.velocities > 0:
            self.velocities = np.empty((n_atoms, 3), dtype=np.float32)
        else:
            self.velocities = None
        if imdsinfo.forces > 0:
            self.forces = np.empty((n_atoms, 3), dtype=np.float32)
        else:
            self.forces = None


def imdframe_memsize(n_atoms, imdsinfo) -> int:
    """
    Calculate the memory size of an IMDFrame in bytes
    """
    memsize = 0
    if imdsinfo.energies > 0:
        memsize += 4 * 10
    if imdsinfo.dimensions > 0:
        memsize += 4 * 9
    if imdsinfo.positions > 0:
        memsize += 4 * 3 * n_atoms
    if imdsinfo.velocities > 0:
        memsize += 4 * 3 * n_atoms
    if imdsinfo.forces > 0:
        memsize += 4 * 3 * n_atoms

    return memsize


def read_into_buf(sock, buf) -> bool:
    """Receives len(buf) bytes into buf from the socket sock"""
    view = memoryview(buf)
    total_received = 0
    while total_received < len(view):
        try:
            received = sock.recv_into(view[total_received:])
            if received == 0:
                # Server called close()
                # Server is definitely done sending frames
                logger.debug(
                    "IMDProducer: recv excepting due to server calling close()"
                )
                raise EOFError
        except TimeoutError:
            # Server is *likely* done sending frames
            logger.debug("IMDProducer: recv excepting due to timeout")
            raise EOFError
        except BlockingIOError:
            # Occurs when timeout is 0 in place of a TimeoutError
            # Server is *likely* done sending frames
            logger.debug("IMDProducer: recv excepting due to blocking")
            raise EOFError
        total_received += received


def sock_contains_data(sock, timeout) -> bool:
    ready_to_read, ready_to_write, in_error = select.select(
        [sock], [], [], timeout
    )
    return sock in ready_to_read
