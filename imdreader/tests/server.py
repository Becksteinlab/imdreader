import socket
import struct
import numpy as np
import threading
import time
from imdreader.IMDProtocol import (
    IMDHeaderType,
    create_header_bytes,
    create_energy_bytes,
    create_imdv3_session_info_bytes,
    IMDHEADERSIZE,
)
import logging
from imdreader.IMDClient import sock_contains_data, read_into_buf, IMDHeader

logger = logging.getLogger("imdreader.IMDREADER")


class InThreadIMDServer:
    """Server operates in the same thread as the client (except for the method
    'listen_accept_handshake_send_ts', which spawns a new thread to get past the
    reader connecting to the server).

    Note that the trajectory units of the provided traj must be in MDAnalysis units
    """

    def __init__(self, traj):
        self.traj = traj
        self.listen_socket = None
        self.conn = None
        self.accept_thread = None

    def set_imdsessioninfo(self, imdsinfo):
        self.imdsinfo = imdsinfo

    def listen_accept_handshake_send_ts(self, host, port):
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.bind((host, port))
        logger.debug(f"InThreadIMDServer: Listening on {host}:{port}")
        s.listen(60)
        self.listen_socket = s
        self.accept_thread = threading.Thread(
            target=self._accept_handshake_send_ts,
        )
        self.accept_thread.start()

    def _accept_handshake_send_ts(self):
        logger.debug(f"InThreadIMDServer: Entering accept thread")
        waited = 0

        if sock_contains_data(self.listen_socket, 5):
            logger.debug(f"InThreadIMDServer: Accepting connection")
            self.conn, _ = self.listen_socket.accept()
            if self.imdsinfo.version == 2:
                self._send_handshakeV2()
            elif self.imdsinfo.version == 3:
                self._send_handshakeV3()
            self._expect_go()
            self.send_frame(0)
            return
        else:
            # IMDReader will fail out if it fails to connect
            return

    def _send_handshakeV2(self):
        header = struct.pack("!i", IMDHeaderType.IMD_HANDSHAKE.value)
        header += struct.pack(f"{self.imdsinfo.endianness}i", 2)
        self.conn.sendall(header)

    def _send_handshakeV3(self):
        logger.debug(f"InThreadIMDServer: Sending handshake V3")
        packet = struct.pack("!i", IMDHeaderType.IMD_HANDSHAKE.value)
        packet += struct.pack(f"{self.imdsinfo.endianness}i", 3)
        self.conn.sendall(packet)

        sinfo = struct.pack("!ii", IMDHeaderType.IMD_SESSIONINFO.value, 7)
        time = 1 if self.imdsinfo.time else 0
        energies = 1 if self.imdsinfo.energies else 0
        box = 1 if self.imdsinfo.box else 0
        positions = 1 if self.imdsinfo.positions else 0
        velocities = 1 if self.imdsinfo.velocities else 0
        forces = 1 if self.imdsinfo.forces else 0
        wrapped_coords = 0
        sinfo += struct.pack(
            f"{self.imdsinfo.endianness}BBBBBBB",
            time,
            energies,
            box,
            positions,
            velocities,
            forces,
            wrapped_coords,
        )
        logger.debug(f"InThreadIMDServer: Sending session info")
        self.conn.sendall(sinfo)

    def join_accept_thread(self):
        self.accept_thread.join()

    def _expect_go(self):
        logger.debug(f"InThreadIMDServer: Waiting for go")
        head_buf = bytearray(IMDHEADERSIZE)
        read_into_buf(self.conn, head_buf)
        header = IMDHeader(head_buf)
        if header.type != IMDHeaderType.IMD_GO:
            raise ValueError("Expected IMD_GO packet, got something else")

    def send_frames(self, start, end):
        for i in range(start, end):
            self.send_frame(i)

    def send_frame(self, i):
        endianness = self.imdsinfo.endianness

        if self.imdsinfo.time:
            time_header = create_header_bytes(IMDHeaderType.IMD_TIME, 1)
            time = struct.pack(
                f"{endianness}ff", self.traj[i].data["dt"], self.traj[i].time
            )

            self.conn.sendall(time_header + time)

        if self.imdsinfo.energies:
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
            self.conn.sendall(energy_header + energies)

        if self.imdsinfo.box:
            box_header = create_header_bytes(IMDHeaderType.IMD_BOX, 1)
            box = np.ascontiguousarray(
                self.traj[i].triclinic_dimensions, dtype=f"{endianness}f"
            ).tobytes()

            self.conn.sendall(box_header + box)

        if self.imdsinfo.positions:
            pos_header = create_header_bytes(
                IMDHeaderType.IMD_FCOORDS, self.traj.n_atoms
            )

            pos = np.ascontiguousarray(
                self.traj[i].positions, dtype=f"{endianness}f"
            ).tobytes()

            self.conn.sendall(pos_header + pos)

        if self.imdsinfo.velocities:
            vel_header = create_header_bytes(
                IMDHeaderType.IMD_VELOCITIES, self.traj.n_atoms
            )
            vel = np.ascontiguousarray(
                self.traj[i].velocities, dtype=f"{endianness}f"
            ).tobytes()

            self.conn.sendall(vel_header + vel)

        if self.imdsinfo.forces:
            force_header = create_header_bytes(
                IMDHeaderType.IMD_FORCES, self.traj.n_atoms
            )
            force = np.ascontiguousarray(
                self.traj[i].forces, dtype=f"{endianness}f"
            ).tobytes()

            self.conn.sendall(force_header + force)

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