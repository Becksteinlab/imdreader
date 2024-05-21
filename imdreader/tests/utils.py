from typing import Tuple
from imdreader.IMDProtocol import *
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
import numpy as np
import socket
import threading
import time
import select


class DummyIMDServer(threading.Thread):
    def __init__(
        self,
        imdwait=True,
        imdpull=False,
        imdterm=True,
        port=8888,
        endianness="<",
        version=2,
        frames=None,
    ):
        self.port = port

        self.imdwait = imdwait
        self.imdpull = imdpull
        self.imdterm = imdterm
        self.endian = endianness
        self.version = version

        self.step = 0
        self.imdstep = -1

        self.energy = np.arange(9, dtype=np.float32)
        self.positions = np.array(
            [[0, 0, 0], [1, 1, 1], [2, 2, 2]], dtype=np.float32
        )

    def run(self):
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.socket.bind(("localhost", self.port))
        self._simulation_loop()

    def _send_handshake(self):
        type = struct.pack("!i", IMDType.IMD_HANDSHAKE.value)
        length = struct.pack(f"{self.endian}i", self.version)
        header = type + length
        self.conn.send(header)

    def disconnect(self):
        header = create_header_bytes(IMDType.IMD_DISCONNECT, 0)
        self.connection.close()
        self.socket.close()

    def _simulation_loop(self):
        if self.imdwait and self.imdstep == -1:
            self._handle_connect()

        # First, read all incoming packets from the client
        while self._peek_header():
            command = self._expect_header()
            if command.type == IMDType.IMD_KILL:
                pass
            if command.type == IMDType.IMD_DISCONNECT:
                pass
            if command.type == IMDType.IMD_MDCOMM:
                pass
            if command.type == IMDType.IMD_PAUSE:
                pass
            if command.type == IMDType.IMD_TRATE:
                pass

        # Then, send energy + position data

    def _await_go(self):
        self._expect_header(expected_type=IMDType.IMD_GO)

    def _expect_header(self, expected_type=None, expected_value=None):
        """
        Read a header packet from the socket.
        """
        header = parse_header_bytes(self._recv_n_bytes(IMDHEADERSIZE))
        if expected_type is not None and header.type != expected_type:
            raise ValueError(
                f"DummyIMDServer: Expected packet type {expected_type}, got {header.type}"
            )
        elif expected_value is not None and header.length != expected_value:
            raise ValueError(
                f"DummyIMDServer: Expected packet length {expected_value}, got {header.length}"
            )
        return header

    def _peek_header(self):
        """
        Peek at a header packet from the socket without consuming it
        """

        data = self.socket.recv(IMDHEADERSIZE, socket.MSG_PEEK)
        if data:
            header = parse_header_bytes(data)
            return header
        return None

    def _recv_n_bytes(self, num_bytes):
        """Receive an arbitrary number of bytes from the socket."""
        data = bytearray(num_bytes)
        view = memoryview(data)
        total_received = 0
        while total_received < num_bytes:
            chunk = self.conn.recv(num_bytes - total_received)
            if not chunk:
                raise ConnectionError("Socket connection was closed")
            view[total_received : total_received + len(chunk)] = chunk
            total_received += len(chunk)
        return data

    def _handle_connect(self):
        # GROMACS waits for a connection for an unlimited time here
        self.socket.listen(1)
        self.conn, self.address = self.socket.accept()

        self._send_handshake()
        # must recieve go within 1 second of sending handshake
        # if not, terminate connection
        self.conn.settimeout(1)
        self._await_go()
        self.conn.settimeout(None)


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
