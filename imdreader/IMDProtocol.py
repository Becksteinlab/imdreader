import socket
import struct
import logging
from enum import Enum, auto

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
IMDVERSION = 2


class IMDType(Enum):
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


class IMDHeader:
    """Convenience class to represent the header of an IMD packet"""

    def __init__(self, msg_type: IMDType, length: int):
        self.type = msg_type
        self.length = length


def create_header_bytes(msg_type: IMDType, length: int):
    # NOTE: add error checking for invalid packet msg_type here
    type = msg_type.value
    return struct.pack("!ii", type, length)


def parse_header_bytes(data):
    msg_type, length = struct.unpack("!ii", data)
    type = IMDType(msg_type)
    # NOTE: add error checking for invalid packet msg_type here
    return IMDHeader(type, length)
