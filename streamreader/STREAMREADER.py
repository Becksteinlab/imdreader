from MDAnalysis.coordinates.base import ReaderBase
from .IMDProtocol import *
import socket
import threading
import numpy as np
from MDAnalysis.lib.util import store_init_arguments

class StreamReader(ReaderBase):
    """
    Reader for IMD protocol packets.

    Default buffer size is set to 8MB for testing
    Buffer_size kwarg is in bytes

    We are assuming the header will never be sent without the body as in the sample code.
    If this assumption is violated, the producer thread can cause a deadlock.
    """
    format = "STREAM"

    @store_init_arguments
    def __init__(self, filename,
                 convert_units=True,
                 forces=False,
                 num_atoms=None,
                 n_frames=None,
                 buffer_size=2**23,
                 **kwargs):
        super(StreamReader, self).__init__(filename, **kwargs)
        self._host, self._port = parse_host_port(filename)
        self._socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        # NOTE: Replace me with header packet which contains this information
        # OR get this information from the topology?
        if not n_frames:
            raise  ValueError("n_frames must be specified")
        self.n_frames = n_frames
        if not num_atoms:
            raise  ValueError("num_atoms must be specified")
        self.n_atoms = num_atoms
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)

        # NOTE: We need to know if we are getting forces or not to construct the correct buffers
        self._buffer = CircularNPBuf(buffer_size, self.n_atoms)
        # The body of a force or position packet should contain
        # (4 bytes per float * 3 atoms * n_atoms) bytes
        self._expected_data_bytes = 12 * self.n_atoms
        self._get_client_connection()
        self._await_IMD_handshake()
        self._start_producer_thread()
        self._frame = -1        
        
    def _get_client_connection(self):
        """
        Listen for a connection on the specified port.
        """
        self._socket.bind((self._host, self._port))
        self._socket.listen(1)
        print("waiting for connection...")
        self._conn, self._addr = self._socket.accept()

    def _await_IMD_handshake(self):
        """
        Wait for the client to send a handshake packet and
        check IMD Protocol version from sender.
        """
        print("waiting for handshake...")
        handshake = self._expect_header(IMDType.IMD_HANDSHAKE)
        if handshake.length != IMDVERSION:
            self._pause_simulation()
            raise ValueError(f"Incompatible IMD version. Expected {IMDVERSION}, got {handshake.length}")
        print("handshake received")

    def _pause_simulation(self):
        """
        Pause the simulation for buffer filling purposes
        or to exit gracefully.
        """
        stp = create_header_bytes(IMDType.IMD_PAUSE, 0)
        self._conn.sendall(stp)
        
    def _send_go_packet(self):
        """
        Send a go packet to the client to start the simulation
        and begin receiving data.
        """
        print("sending go packet...")
        header = create_header_bytes(IMDType.IMD_GO, 0)
        self._conn.sendall(header)

    def _start_producer_thread(self):
        """
        Starts the thread which reads packets from the socket and
        places them into a buffer for the parent thread to read.
        """
        self._producer_thread = threading.Thread(target=self._control_flow)
        self._producer_thread.daemon = True
        self._producer_thread.start()

    def _control_flow(self):
        """
        Producer thread method. Reads from the socket and
        sends a 'pause' signal if needed.
        """
        self._parsed_frames = 0

        while self._parsed_frames < self.n_frames:
            header = self._expect_header()
            if header.type == IMDType.IMD_ENERGIES and header.length == IMDENERGYPACKETLENGTH:
                self._recv_energies()
            elif header.type == IMDType.IMD_FCOORDS and header.length == self._expected_data_bytes:
                self._recv_fcoords()
                self._parsed_frames += 1
            else:
                raise ValueError("Unexpected packet type {}".format(header.type))

    def _expect_header(self, expected_type=None):
        """
        Read a header packet from the socket.
        """
        header = parse_header_bytes(self._recv_n_bytes(IMDHEADERSIZE))
        if expected_type is not None and header.type != expected_type:
            raise ValueError("Expected packet type {}, got {}".format(expected_type, header.type))
        return header

    def _recv_n_bytes(self, num_bytes):
        data = bytearray(num_bytes)
        view = memoryview(data)
        total_received = 0
        while total_received < num_bytes:
            chunk = self._conn.recv(num_bytes - total_received)
            if not chunk:
                raise ConnectionError("Socket connection was closed")
            view[total_received:total_received+len(chunk)] = chunk
            total_received += len(chunk)
        return data

    def _recv_fcoords(self):
        self._buffer.insert(self._recv_n_bytes(self._expected_data_bytes))

    def _recv_energies(self):
        self._recv_n_bytes(IMDENERGYPACKETLENGTH)

    def _read_next_timestep(self):
        if self._frame == -1:
            self._send_go_packet()
        return self._read_frame(self._frame + 1)

    def _read_frame(self, frame):
        self._frame = frame
        pos = self._buffer.get_frame()
        self.ts.positions = pos
        self.ts.frame = self._frame
        self.ts.dimensions = None
        return self.ts

    def _reopen(self):
        # NOTE: try setting self._frame to -1
        pass

    def rewind(self):
        pass

    @staticmethod
    def _format_hint(thing):
        try:
            parse_host_port(thing)
        except:
            return False
        return True
    
    def close(self):
        """Gracefully shut down the reader."""
        self._socket.close()
        print("StreamReader shut down gracefully.")
    
    # Incompatible methods
    def copy(self):
        raise NotImplementedError("StreamReader does not support copying")

    def __getitem__(self, frame):
        raise NotImplementedError("StreamReader does not support indexing")

def parse_host_port(filename):
        # Check if the format is correct
        parts = filename.split(':')
        if len(parts) == 2:
            host = parts[0]  # Hostname part
            try:
                port = int(parts[1])  # Convert the port part to an integer
                return (host, port)
            except ValueError:
                # Handle the case where the port is not a valid integer
                raise ValueError("Port must be an integer")
        else:
            # Handle the case where the format does not match "host:port"
            raise ValueError("Filename must be in the format 'host:port'")
        
class CircularNPBuf():
    """
    Thread-safe circular numpy buffer
    """
    def __init__(self, buffer_size, n_atoms):
        self._buf = np.empty(buffer_size // 4, dtype=np.float32)
        self._mutex = threading.Lock()
        self._not_empty = threading.Condition(self._mutex)
        self._not_full = threading.Condition(self._mutex)
        self._frame_size = (n_atoms * 3)
        # _empty refers to the number of empty
        # "frames" in the buffer
        # where a frame refers to space in a numpy array
        # for a full frame of positions or data
        self._empty = len(self._buf) // self._frame_size
        self._full = 0
        self._fill = 0
        self._use = 0

    def insert(self, data_bytes):
        with self._not_full:
            while not self._empty:
                self._not_full.wait()
            # Write that requires wrapping
            if self._fill + self._frame_size - 1 >= len(self._buf):
                n_first_write_elements = len(self._buf) - self._fill
                self._buf[self._fill : len(self._buf)] = np.frombuffer(data_bytes[:n_first_write_elements * 4], dtype='>f4')
                n_remaining_elements = self._frame_size - n_first_write_elements
                self._fill = 0
                self._buf[self._fill : self._fill + n_remaining_elements] = np.frombuffer(data_bytes[n_first_write_elements * 4:
                                                                                       (n_first_write_elements + 
                                                                                        n_remaining_elements) * 4],
                                                                                        dtype='>f4')
                self._fill += n_remaining_elements
            else:
                self._buf[self._fill : self._fill + self._frame_size] = np.frombuffer(data_bytes, dtype='>f4')
                self._fill += self._frame_size
            self._full += 1
            self._empty -= 1
            self._not_empty.notify()
    
    def get_frame(self):
        frame = np.empty(self._frame_size, dtype=np.float32)
        with self._not_empty:
            while not self._full:
                self._not_empty.wait()
            # Read that requires wrapping
            if self._use + self._frame_size - 1 >= len(self._buf):
                n_first_read_elements = len(self._buf) - self._use
                frame[0 : len(self._buf)] = self._buf[self._use : len(self._buf)]
                n_remaining_elements = self._frame_size - n_first_read_elements
                self._use = 0
                frame[n_first_read_elements : n_first_read_elements
                       + n_remaining_elements] = self._buf[self._use : self._use + n_remaining_elements]
                self._use += n_remaining_elements
            else:
                frame[:] = self._buf[self._use : self._use + self._frame_size]
                self._use += self._frame_size
            self._full -= 1
            self._empty += 1
            self._not_full.notify()
            return frame
