from typing import Tuple
from imdreader.IMDProtocol import *
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
import numpy as np
import socket


class DummyIMDServer:
    def __init__(
        self,
        imdwait=True,
        imdpull=False,
        imdterm=True,
        port=8888,
        endianness="<",
        version=2,
    ):
        self.port = port
        self.imdstep = -1
        self.imdwait = imdwait
        self.imdpull = imdpull
        self.imdterm = imdterm
        self.endian = endianness
        self.version = version

        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.socket.bind(("localhost", self.port))
        # GROMACS waits for a connection for an unlimited time here
        self.socket.listen(1)
        self.conn, self.address = self.socket.accept()

        self._send_handshake()
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
            self._await_go()
        pass

    def _await_go(self):
        # must recieve go within 1 second of sending handshake
        # if not, terminate connection
        self.conn.settimeout(1)
        header = self._expect_header(expected_type=IMDType.IMD_GO)
        self.conn.settimeout(None)

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
