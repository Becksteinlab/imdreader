from imdreader.IMDProtocol import *
import MDAnalysis as mda
import numpy as np
import socket
import threading
import time
from imdreader.IMDProtocol import *
from imdreader.IMDREADER import read_into_buf
from MDAnalysisTests.datafiles import COORDINATES_TOPOLOGY, COORDINATES_TRR


class Behavior(abc.ABC):
    @abc.abstractmethod
    def perform(self):
        pass


class DefaultConnectionBehavior(Behavior):
    def __init__(self, host, port):
        self.host = host
        self.port = port

    def perform(self):
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.bind((self.host, self.port))
        s.listen(60)
        conn, addr = s.accept()
        return (conn, addr)


class DefaultHandshakeV2Behavior(Behavior):
    def __init__(self, imdsessioninfo):
        self.version = imdsessioninfo.version
        self.endianness = imdsessioninfo.endianness

    def perform(self, conn):
        header = struct.pack("!i", IMDHeaderType.IMD_HANDSHAKE.value)
        if self.endianness == "<":
            header += struct.pack("<i", self.version)
        else:
            header += struct.pack(">i", self.version)
        conn.sendall(header)


class DefaultHandshakeV3Behavior(Behavior):
    def __init__(self, imdsessioninfo):
        pass

    def perform(self, conn):
        pass


class DefaultAwaitGoBehavior(Behavior):
    def __init__(self, timeout=1):
        self.timeout = timeout

    def perform(self, conn):
        conn.settimeout(self.timeout)
        head_buf = bytearray(IMDHEADERSIZE)
        read_into_buf(conn, head_buf)
        header = IMDHeader(head_buf)
        if header.type != IMDHeaderType.IMD_GO:
            raise ValueError("Expected IMD_GO packet, got something else")
        conn.settimeout(None)


class DefaultLoopBehaviorV2(Behavior):
    def __init__(self, trajectory, imdsessioninfo):
        self.traj = trajectory
        self.endianness = imdsessioninfo.endianness
        self.imdterm = imdsessioninfo.imdterm
        self.imdwait = imdsessioninfo.imdwait
        self.imdpull = imdsessioninfo.imdpull

    def perform(self, conn):
        conn.settimeout(1)
        headerbuf = bytearray(IMDHEADERSIZE)
        paused = False

        energies = np.zeros((len(self.traj), 10), dtype=np.float32)

        for i in range(len(self.traj)):
            while sock_contains_data(conn, 1) or paused:
                header_success = read_into_buf(conn, headerbuf)
                if header_success:
                    header = IMDHeader(headerbuf)
                    if header.type == IMDHeaderType.IMD_PAUSE:
                        paused = not paused

            energy_header = create_header_bytes(IMDHeaderType.IMD_ENERGIES, 1)
            energy = np.ascontiguousarray(
                energies[i], dtype=f"{self.endianness}f4"
            ).tobytes()
            pos_header = create_header_bytes(
                IMDHeaderType.IMD_FCOORDS, self.traj.n_atoms
            )
            pos = np.ascontiguousarray(
                self.traj[i].positions, dtype=f"{self.endianness}f4"
            ).tobytes()

            conn.sendall(energy_header + energy)
            conn.sendall(pos_header + pos)


class DefaultDisconnectBehavior(Behavior):
    def __init__(self, pre_shutdown_wait=0, pre_close_wait=0):
        self.pre_shutdown_wait = pre_shutdown_wait
        self.pre_close_wait = pre_close_wait

    def perform(self, conn):
        time.sleep(self.pre_shutdown_wait)
        # Gromacs uses the c equivalent of the SHUT_WR flag
        conn.shutdown(socket.SHUT_WR)
        time.sleep(self.pre_close_wait)
        conn.close()


class DummyIMDServer(threading.Thread):
    """Performs the following steps in order:

    1. ConnectionBehavior.perform_connection()
    2. HandshakeBehavior.perform_handshake()
    3. AwaitGoBehavior.perform_await_go()
    4. LoopBehavior.perform_loop()
    5. DisconnectBehavior.perform_disconnect()

    Each if these behaviors can be changed by calling DummyIMDServer.set_x_behavior(y) where x is the behavior name
    and y is the behavior object.

    Start the server by calling DummyIMDServer.start().
    """

    def __init__(
        self,
        host="localhost",
        port=8888,
        imdsessioninfo=None,
        traj=None,
    ):
        """
        If passing `traj` kwarg, ensure it is a copy of the trajectory used
        in the test to avoid modifying the original trajectory "head" in the
        main thread.
        """
        super().__init__(daemon=True)

        if traj is None:
            traj = mda.Universe(
                COORDINATES_TOPOLOGY, COORDINATES_TRR
            ).trajectory

        if imdsessioninfo is None:
            imdsessioninfo = IMDSessionInfo(
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

        self.connection_behavior = DefaultConnectionBehavior(host, port)

        if imdsessioninfo.version == 2:
            self.handshake_behavior = DefaultHandshakeV2Behavior(imdsessioninfo)
            self.loop_behavior = DefaultLoopBehaviorV2(traj, imdsessioninfo)
        elif imdsessioninfo.version == 3:
            self.connection_behavior = DefaultConnectionBehavior(
                "localhost", port
            )
            # self.handshake_behavior = DefaultHandshakeV3Behavior()

        self.await_go_behavior = DefaultAwaitGoBehavior()

        self.disconnect_behavior = DefaultDisconnectBehavior()

    def run(self):
        conn = self.connection_behavior.perform()[0]
        self.handshake_behavior.perform(conn)
        self.await_go_behavior.perform(conn)
        self.loop_behavior.perform(conn)
        self.disconnect_behavior.perform(conn)

    def set_connection_behavior(self, connection_behavior):
        self.connection_behavior = connection_behavior

    def set_handshake_behavior(self, handshake_behavior):
        self.handshake_behavior = handshake_behavior

    def set_await_go_behavior(self, await_go_behavior):
        self.await_go_behavior = await_go_behavior

    def set_loop_behavior(self, loop_behavior):
        self.loop_behavior = loop_behavior

    def set_disconnect_behavior(self, disconnect_behavior):
        self.disconnect_behavior = disconnect_behavior


def sock_contains_data(sock, timeout) -> bool:
    ready_to_read, ready_to_write, in_error = select.select(
        [sock], [], [], timeout
    )
    return sock in ready_to_read


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


def check_port_availability(port):
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        s.bind(("0.0.0.0", port))
    except socket.error as e:
        if e.errno == socket.errno.EADDRINUSE:
            print(f"Port {port} is already in use")
            return False
        else:
            raise
    finally:
        s.close()
    return True


def get_free_port():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("", 0))
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        return s.getsockname()[1]
