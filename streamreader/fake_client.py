import socket
import struct
from streamreader import *
from streamreader.IMDProtocol import *


client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
client_socket.connect(("localhost", 8000))

handshake = create_header_bytes(IMDType.IMD_HANDSHAKE, 2)
client_socket.sendall(handshake)

response = client_socket.recv(8)
go_signal = parse_header_bytes(response)
print(go_signal.type)

coordsheader = create_header_bytes(IMDType.IMD_FCOORDS, 12)
coordsbody = struct.pack('!fff', 1.1, 2.3, 3.4)
client_socket.sendall(coordsheader + coordsbody)
coordsbody = struct.pack('!fff', 1.2, 2.3, 3.4)
client_socket.sendall(coordsheader + coordsbody)

client_socket.close()