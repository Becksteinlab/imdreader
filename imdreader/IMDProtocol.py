import struct
import logging
from enum import Enum, auto
from typing import Union
from dataclasses import dataclass

"""
IMD Packets have an 8 byte header and a variable length payload 
who's length is specified in the header.

struct IMDheader{
  int32 type;
  int32 length;
};

In a HANDSHAKE header packet, the length attribute is used to store the version 
of the IMD protocol
"""
IMDHEADERSIZE = 8
IMDENERGYPACKETLENGTH = 40
IMDBOXPACKETLENGTH = 36
IMDVERSIONS = {2, 3}
IMDAWAITGOTIME = 1


class IMDHeaderType(Enum):
    IMD_DISCONNECT = 0
    IMD_ENERGIES = 1
    IMD_FCOORDS = 2
    IMD_GO = 3
    IMD_HANDSHAKE = 4
    IMD_KILL = 5
    IMD_MDCOMM = 6
    IMD_PAUSE = 7
    IMD_TRATE = 8
    IMD_IOERROR = 9
    # New in IMD v3
    # IMD_BOX = 10
    # IMD_VELS = 11
    # IMD_FORCES = 12
    # IMD_EOS = 13


def parse_energy_bytes(data, endianness):
    keys = [
        "step",
        "temperature",
        "total_energy",
        "potential_energy",
        "van_der_walls_energy",
        "coulomb_energy",
        "bonds_energy",
        "angles_energy",
        "dihedrals_energy",
        "improper_dihedrals_energy",
    ]
    values = struct.unpack(f"{endianness}ifffffffff", data)
    return dict(zip(keys, values))


def create_energy_bytes(
    step,
    temperature,
    total_energy,
    potential_energy,
    van_der_walls_energy,
    coulomb_energy,
    bonds_energy,
    angles_energy,
    dihedrals_energy,
    improper_dihedrals_energy,
    endianness,
):
    return struct.pack(
        f"{endianness}ifffffffff",
        step,
        temperature,
        total_energy,
        potential_energy,
        van_der_walls_energy,
        coulomb_energy,
        bonds_energy,
        angles_energy,
        dihedrals_energy,
        improper_dihedrals_energy,
    )


class IMDHeader:
    """Convenience class to represent the header of an IMD packet"""

    def __init__(self, data):
        msg_type, length = struct.unpack("!ii", data)
        h_type = IMDHeaderType(msg_type)

        self.type = h_type
        self.length = length


@dataclass
class IMDSessionInfo:
    """Convenience class to represent the session information of an IMD connection

    '<' represents little endian and '>' represents big endian

    Data should be loaded into and out of buffers in the order of the fields in this class
    if present in the session for that step, i.e.
        1. energies,
        2. dimensions,
        etc.
    """

    version: int
    endianness: str
    imdterm: Union[bool, None]
    imdwait: Union[bool, None]
    imdpull: Union[bool, None]
    wrapped_coords: bool
    energies: int
    dimensions: int
    positions: int
    velocities: int
    forces: int


def parse_imdv3_session_info(sock, end) -> Union[IMDSessionInfo, None]:
    """Parses the session information packet of an IMD v3 connection"""
    pass


def create_header_bytes(msg_type: IMDHeaderType, length: int):
    # NOTE: add error checking for invalid packet msg_type here

    type = msg_type.value
    return struct.pack("!ii", type, length)


def parse_header_bytes(data):
    msg_type, length = struct.unpack("!ii", data)
    type = IMDHeaderType(msg_type)
    # NOTE: add error checking for invalid packet msg_type here
    return IMDHeader(type, length)
