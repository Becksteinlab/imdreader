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
    """
    Convenience class to represent the header of an IMD packet
    """
    def __init__(self, msg_type : IMDType, length : int):
        self.type = msg_type
        self.length = length

class IMDEnergies:
    def __init__(self, tstep, T, Etot, Epot, Evdw, Eelec, Ebond, Eangle, Edihe, Eimpr):
        self.tstep = tstep
        self.T = T
        self.Etot = Etot
        self.Epot = Epot
        self.Evdw = Evdw
        self.Eelec = Eelec
        self.Ebond = Ebond
        self.Eangle = Eangle
        self.Edihe = Edihe
        self.Eimpr = Eimpr

def create_header_bytes(msg_type : IMDType, length : int):
    # Network byte order is big-endian
    # NOTE: add error checking for invalid packet msg_type here
    type = msg_type.value
    return struct.pack('!ii', type, length)

def parse_header_bytes(data):
    msg_type, length = struct.unpack('!ii', data)
    type = IMDType(msg_type)
    # NOTE: add error checking for invalid packet msg_type here
    return IMDHeader(type, length)


class IMDConnection:
    
    
    def __init__(self, host, port):
        self.host = host
        self.port = port
        try:
            self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        except socket.error as e:
            logging.error("Failed to create socket: %s", e)
            return False
        return True

    def connect(self):
        try:
            self.sock.connect((self.host, self.port))
        except socket.error as e:
            logging.error("Connection failed: %s", e)
            return False
        return True
    
    def disconnect(self):
        header = create_header(IMDType.IMD_DISCONNECT, 0)
        self.sock.sendall(header)
        self.sock.close()
    
    def send_data(self, data):
        self.sock.sendall(data)
    
    def receive_data(self, size):
        return self.sock.recv(size)
    
    def handshake(self):
        header = create_header(IMDType.IMD_HANDSHAKE, IMDConnection.IMDVERSION)
        self.send_data(header)
        response = self.receive_data(IMDConnection.HEADERSIZE)
        _, version = parse_header(response)
        if version != IMDConnection.IMDVERSION:
            raise ValueError("Version mismatch")
    
    def send_energies(self, energies):
        data = struct.pack('!ifffffffffff', energies.tstep, energies.T, energies.Etot, energies.Epot,
                           energies.Evdw, energies.Eelec, energies.Ebond, energies.Eangle, energies.Edihe, energies.Eimpr)
        header = create_header(IMDType.IMD_ENERGIES, len(data))
        self.send_data(header + data)
    
    def receive_energies(self):
        header_data = self.receive_data(IMDConnection.HEADERSIZE)
        _, length = parse_header(header_data)
        data = self.receive_data(length)
        values = struct.unpack('!ifffffffffff', data)
        return IMDEnergies(*values)
