# cython: language_level=3
from libc.stdint cimport int32_t
from libc.stdlib cimport malloc, free
from libc.string cimport memcpy

# Declare the external C functions and structures
cdef extern from "imd.h":
    ctypedef enum IMDType:
        IMD_DISCONNECT, IMD_PAUSE, IMD_KILL, IMD_GO, IMD_HANDSHAKE, IMD_TRATE, IMD_MDCOMM, IMD_ENERGIES, IMD_FCOORDS, IMD_IOERROR

    ctypedef struct IMDheader:
        int32_t type
        int32_t length

    void fill_header(IMDheader *header, IMDType type, int32_t length)
    int imd_disconnect(void* s)
    int imd_pause(void* s)
    int imd_kill(void* s)
    int imd_go(void* s)
    int imd_handshake(void* s)
    int imd_trate(void* s, int32_t rate)
    int imd_send_mdcomm(void* s, int32_t n, const int32_t* indices, const float* forces)
    int imd_send_energies(void* s, const void* energies)  # Assuming IMDEnergies is defined elsewhere
    int imd_send_fcoords(void* s, int32_t n, const float* coords)
    int imd_recv_mdcomm(void* s, int32_t n, int32_t* indices, float* forces)
    int imd_recv_energies(void* s, void* energies)
    int imd_recv_fcoords(void* s, int32_t n, float* coords)
    int imd_recv_handshake(void* s)

cdef extern from "vmdsock.h":
    void* vmdsock_create()
    int vmdsock_bind(void* sock, int port)
    int vmdsock_listen(void* sock)
    void* vmdsock_accept(void* sock)
    int vmdsock_read(void* sock, char* buf, int nbytes)
    int vmdsock_write(void* sock, const char* buf, int nbytes)

# Provide Python wrapper functions
def create_socket():
    return vmdsock_create()

def bind_socket(sock, int port):
    return vmdsock_bind(sock, port)

def listen_socket(sock):
    return vmdsock_listen(sock)

def accept_socket(sock):
    return vmdsock_accept(sock)

def disconnect(sock):
    return imd_disconnect(sock)

def pause_simulation(sock):
    return imd_pause(sock)

def kill_simulation(sock):
    return imd_kill(sock)

def handshake(sock):
    return imd_handshake(sock)