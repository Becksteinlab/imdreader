from MDAnalysis.coordinates.base import ReaderBase
import socket


class StreamReader(ReaderBase):
    """
    
    
    """
    def __init__(self, filename,
                 convert_units=True,
                 **kwargs):
        super(StreamReader, self).__init__(filename, **kwargs)
        # Start listening on provided port
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        # Set buffer size
        new_recv_buf_size = 512 * 1024  # Example: 512 KB
        new_send_buf_size = 512 * 1024  # Example: 512 KB
        s.setsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF, new_recv_buf_size)
        s.setsockopt(socket.SOL_SOCKET, socket.SO_SNDBUF, new_send_buf_size)