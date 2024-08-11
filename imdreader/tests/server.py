import socket
import struct
import numpy as np
import threading
import time
from imdreader.IMDProtocol import (
    IMDHeaderType,
    create_header_bytes,
    create_energy_bytes,
    IMDHEADERSIZE,
)
import logging
from imdreader.IMDClient import sock_contains_data, read_into_buf, IMDHeader

logger = logging.getLogger("imdreader.IMDREADER")


class TestIMDServer:

    def __init__(self, traj):
        self.traj = traj
        self.listen_socket = None
        self.conn = None
        self.accept_thread = None

    def listen_accept_handshake_send_ts(self, host, port, imdsinfo):
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.bind((host, port))
        logger.debug(f"TestIMDServer: Listening on {host}:{port}")
        s.listen(60)
        self.listen_socket = s
        self.accept_thread = threading.Thread(
            target=self._accept_and_handshake, args=(imdsinfo,)
        )
        self.accept_thread.start()

    def _accept_and_handshake(self, imdsinfo):
        logger.debug(f"TestIMDServer: Entering accept thread")
        waited = 0

        if sock_contains_data(self.listen_socket, 5):
            self.conn, _ = self.listen_socket.accept()
            if imdsinfo.version == 2:
                self._send_handshakeV2(imdsinfo.endianness)
            self._expect_go()
            self.send_frame(0, imdsinfo.endianness)
            return
        else:
            # IMDReader will fail out if it fails to connect
            return

    def join_accept_thread(self):
        logger.debug(f"TestIMDServer: Joining accept thread")
        self.accept_thread.join()

    def _send_handshakeV2(self, endianness="<"):
        header = struct.pack("!i", IMDHeaderType.IMD_HANDSHAKE.value)
        header += struct.pack(f"{endianness}i", 2)
        self.conn.sendall(header)

    def _expect_go(self):
        head_buf = bytearray(IMDHEADERSIZE)
        read_into_buf(self.conn, head_buf)
        header = IMDHeader(head_buf)
        if header.type != IMDHeaderType.IMD_GO:
            raise ValueError("Expected IMD_GO packet, got something else")

    def send_frame(self, i, endianness="<"):
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
            endianness,
        )
        pos_header = create_header_bytes(
            IMDHeaderType.IMD_FCOORDS, self.traj.n_atoms
        )

        pos = np.ascontiguousarray(
            self.traj[i].positions, dtype=f"{endianness}f"
        ).tobytes()

        self.conn.sendall(energy_header + energies)
        self.conn.sendall(pos_header + pos)

    def disconnect(self):
        self.conn.shutdown(socket.SHUT_RD)
        self.conn.close()

    def cleanup(self):
        try:
            self.conn.close()
        except:
            pass
        try:
            self.listen_socket.close()
        except:
            pass
