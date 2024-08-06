from imdreader.IMDProtocol import *
import MDAnalysis as mda
import numpy as np
import socket
import threading
import time
from imdreader.IMDProtocol import *
from imdreader.IMDREADER import read_into_buf, sock_contains_data
from MDAnalysisTests.datafiles import COORDINATES_TOPOLOGY, COORDINATES_TRR
import abc
import imdreader
import logging

logger = logging.getLogger(imdreader.IMDClient.__name__)


class Behavior(abc.ABC):
    """Abstract base class for behaviors for the DummyIMDServer to perform.
    Ensure that behaviors do not contain potentially infinite loops- they should be
    testing for a specific sequence of events"""

    @abc.abstractmethod
    def perform(self, *args, **kwargs):
        pass


class DefaultConnectionBehavior(Behavior):

    def perform(self, host, port, event_q):
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.bind((host, port))
        s.listen(60)
        conn, addr = s.accept()
        return (conn, addr)


class DefaultHandshakeV2Behavior(Behavior):
    def perform(self, conn, imdsessioninfo, event_q):
        header = struct.pack("!i", IMDHeaderType.IMD_HANDSHAKE.value)
        if imdsessioninfo.endianness == "<":
            header += struct.pack("<i", 2)
        else:
            header += struct.pack(">i", 2)
        conn.sendall(header)


class DefaultHandshakeV3Behavior(Behavior):
    def perform(self, conn, imdsessioninfo, event_q):
        pass


class DefaultAwaitGoBehavior(Behavior):
    def perform(self, conn, event_q):
        conn.settimeout(IMDAWAITGOTIME)
        head_buf = bytearray(IMDHEADERSIZE)
        read_into_buf(conn, head_buf)
        header = IMDHeader(head_buf)
        if header.type != IMDHeaderType.IMD_GO:
            raise ValueError("Expected IMD_GO packet, got something else")
        logger.debug("DummyIMDServer: Received IMD_GO")


class DefaultLoopV2Behavior(Behavior):
    """Default behavior doesn't allow pausing"""

    def perform(self, conn, traj, imdsessioninfo, event_q):
        conn.settimeout(1)
        headerbuf = bytearray(IMDHEADERSIZE)
        paused = False

        logger.debug("DummyIMDServer: Starting loop")

        for i in range(len(traj)):
            logger.debug(f"DummyIMDServer: generating frame {i}")

            energy_header = create_header_bytes(IMDHeaderType.IMD_ENERGIES, 1)

            energies = create_energy_bytes(
                i,
                i + 1,
                i + 2,
                i + 3,
                i + 4,
                i + 5,
                i + 6,
                i + 7,
                i + 8,
                i + 9,
                imdsessioninfo.endianness,
            )

            pos_header = create_header_bytes(
                IMDHeaderType.IMD_FCOORDS, traj.n_atoms
            )
            pos = np.ascontiguousarray(
                traj[i].positions, dtype=f"{imdsessioninfo.endianness}f"
            ).tobytes()

            conn.sendall(energy_header + energies)
            conn.sendall(pos_header + pos)
            logger.debug(f"Sent frame {i}")


class ExpectPauseLoopV2Behavior(DefaultLoopV2Behavior):
    """Waits for a pause & unpause in all frames after the first frame."""

    def perform(self, conn, traj, imdsessioninfo, event_q):
        conn.settimeout(10)
        headerbuf = bytearray(IMDHEADERSIZE)

        logger.debug("DummyIMDServer: Starting loop")

        for i in range(len(traj)):
            if i != 0:
                read_into_buf(conn, headerbuf)
                header = IMDHeader(headerbuf)
                if header.type != IMDHeaderType.IMD_PAUSE:
                    logger.debug(
                        f"DummyIMDServer: Expected IMD_PAUSE, got {header.type}"
                    )

                read_into_buf(conn, headerbuf)
                header = IMDHeader(headerbuf)
                if header.type != IMDHeaderType.IMD_PAUSE:
                    logger.debug(
                        f"DummyIMDServer: Expected IMD_PAUSE, got {header.type}"
                    )

            logger.debug(f"DummyIMDServer: generating frame {i}")

            energy_header = create_header_bytes(IMDHeaderType.IMD_ENERGIES, 1)

            energies = create_energy_bytes(
                i,
                i + 1,
                i + 2,
                i + 3,
                i + 4,
                i + 5,
                i + 6,
                i + 7,
                i + 8,
                i + 9,
                imdsessioninfo.endianness,
            )

            pos_header = create_header_bytes(
                IMDHeaderType.IMD_FCOORDS, traj.n_atoms
            )
            pos = np.ascontiguousarray(
                traj[i].positions, dtype=f"{imdsessioninfo.endianness}f"
            ).tobytes()

            conn.sendall(energy_header + energies)
            conn.sendall(pos_header + pos)
            logger.debug(f"Sent frame {i}")


class ExpectPauseUnpauseAfterLoopV2Behavior(DefaultLoopV2Behavior):
    def perform(self, conn, traj, imdsessioninfo, event_q):
        conn.settimeout(10)
        headerbuf = bytearray(IMDHEADERSIZE)

        logger.debug("DummyIMDServer: Starting loop")

        for i in range(len(traj)):
            logger.debug(f"DummyIMDServer: generating frame {i}")

            energy_header = create_header_bytes(IMDHeaderType.IMD_ENERGIES, 1)

            energies = create_energy_bytes(
                i,
                i + 1,
                i + 2,
                i + 3,
                i + 4,
                i + 5,
                i + 6,
                i + 7,
                i + 8,
                i + 9,
                imdsessioninfo.endianness,
            )

            pos_header = create_header_bytes(
                IMDHeaderType.IMD_FCOORDS, traj.n_atoms
            )
            pos = np.ascontiguousarray(
                traj[i].positions, dtype=f"{imdsessioninfo.endianness}f"
            ).tobytes()

            conn.sendall(energy_header + energies)
            conn.sendall(pos_header + pos)
            logger.debug(f"Sent frame {i}")

        read_into_buf(conn, headerbuf)
        header = IMDHeader(headerbuf)
        if header.type != IMDHeaderType.IMD_PAUSE:
            logger.debug(
                f"DummyIMDServer: Expected IMD_PAUSE, got {header.type}"
            )

        read_into_buf(conn, headerbuf)
        header = IMDHeader(headerbuf)
        if header.type != IMDHeaderType.IMD_PAUSE:
            logger.debug(
                f"DummyIMDServer: Expected IMD_PAUSE, got {header.type}"
            )


class DefaultDisconnectBehavior(Behavior):
    def perform(self, conn, event_q):
        # Gromacs uses the c equivalent of the SHUT_WR flag
        conn.shutdown(socket.SHUT_WR)
        conn.close()


def create_default_imdsinfo_v2():
    return IMDSessionInfo(
        version=2,
        endianness="<",
        imdterm=None,
        imdwait=None,
        imdpull=None,
        wrapped_coords=True,
        energies=1,
        dimensions=0,
        positions=1,
        velocities=0,
        forces=0,
    )


class DummyIMDServer(threading.Thread):
    """Performs the following steps in order:

    1. ConnectionBehavior.perform_connection()
    2. HandshakeBehavior.perform_handshake()
    3. AwaitGoBehavior.perform_await_go()
    4. LoopBehavior.perform_loop()
    5. DisconnectBehavior.perform_disconnect()

    Start the server by calling DummyIMDServer.start().
    """

    def __init__(
        self,
        traj,
        version,
    ):
        """
        If passing `traj` kwarg, ensure it is a copy of the trajectory used
        in the test to avoid moving the original trajectory "head" in the
        main thread.
        """
        super().__init__(daemon=True)

        logger.debug("DummyIMDServer: Initializing")

        self._traj = traj

        self._host = "localhost"
        self._port = 8888

        self.connection_behavior = DefaultConnectionBehavior()

        if version == 2:
            self._imdsessioninfo = create_default_imdsinfo_v2()
            self.handshake_behavior = DefaultHandshakeV2Behavior()
            self.loop_behavior = DefaultLoopV2Behavior()
        elif version == 3:
            # self.imdsessioninfo = create_default_imdsinfo_v3()
            # self.handshake_behavior = DefaultHandshakeV3Behavior()
            # self.loop_behavior = DefaultLoopBehaviorV3()
            pass

        self.await_go_behavior = DefaultAwaitGoBehavior()
        self.disconnect_behavior = DefaultDisconnectBehavior()

        self._event_q = []

    def run(self):
        conn = self.connection_behavior.perform(
            self.host, self.port, self._event_q
        )[0]
        self.handshake_behavior.perform(
            conn, self.imdsessioninfo, self._event_q
        )
        self.await_go_behavior.perform(conn, self._event_q)
        self.loop_behavior.perform(
            conn, self._traj, self.imdsessioninfo, self._event_q
        )
        self.disconnect_behavior.perform(conn, self._event_q)
        return

    @property
    def port(self):
        return self._port

    @port.setter
    def port(self, port):
        self._port = port

    @property
    def host(self):
        return self._host

    @host.setter
    def host(self, host):
        self._host = host

    @property
    def imdsessioninfo(self):
        return self._imdsessioninfo

    @imdsessioninfo.setter
    def imdsessioninfo(self, imdsessioninfo):
        self._imdsessioninfo = imdsessioninfo

    @property
    def connection_behavior(self):
        return self._connection_behavior

    @connection_behavior.setter
    def connection_behavior(self, connection_behavior):
        self._connection_behavior = connection_behavior

    @property
    def handshake_behavior(self):
        return self._handshake_behavior

    @handshake_behavior.setter
    def handshake_behavior(self, handshake_behavior):
        self._handshake_behavior = handshake_behavior

    @property
    def await_go_behavior(self):
        return self._await_go_behavior

    @await_go_behavior.setter
    def await_go_behavior(self, await_go_behavior):
        self._await_go_behavior = await_go_behavior

    @property
    def loop_behavior(self):
        return self._loop_behavior

    @loop_behavior.setter
    def loop_behavior(self, loop_behavior):
        self._loop_behavior = loop_behavior

    @property
    def disconnect_behavior(self):
        return self._disconnect_behavior

    @disconnect_behavior.setter
    def disconnect_behavior(self, disconnect_behavior):
        self._disconnect_behavior = disconnect_behavior

    @property
    def event_q(self):
        return self._event_q


def recvuntil(file_path, target_line, timeout):
    """
    Read from the file until the target line is found or the timeout occurs.

    Args:
        file_path (str): The path to the file to read from.
        target_line (str): The line to wait for.
        timeout (float): The timeout in seconds.

    Returns:
        str: The line containing the target line.

    Raises:
        TimeoutError: If the target line is not found within the timeout period.
    """
    end_time = time.time() + timeout
    buffer = ""

    while time.time() < end_time:
        time.sleep(0.1)  # Small delay to avoid busy-waiting
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines:
                buffer += line
                if target_line in line:
                    return line
    raise TimeoutError(
        f"Timeout after {timeout} seconds waiting for '{target_line}'"
    )


def get_free_port():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("", 0))
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        return s.getsockname()[1]
