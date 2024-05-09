# distutils: language = c

from libc.stdint cimport uint16_t, int32_t
from libc.sys.socket cimport sockaddr_in
from libc.sys.time cimport timeval
from libc.select cimport fd_set, FD_ZERO, FD_SET
from libc.stdlib cimport malloc, free

# Include necessary system headers for types used in the structs and functions
cdef extern from "sys/types.h":
    ctypedef int socklen_t

cdef extern from "vmdsock.h":
    # Define the vmdsocket struct if it's part of the exposed API
    ctypedef struct vmdsocket:
        sockaddr_in addr
        socklen_t addrlen
        int sd

    # Function declarations
    int vmdsock_init()
    void* vmdsock_create()
    int vmdsock_bind(void* s, int port)
    int vmdsock_listen(void* s)
    void* vmdsock_accept(void* s)
    int vmdsock_connect(void* s, const char* host, int port)
    int vmdsock_write(void* s, const void* buf, int len)
    int vmdsock_read(void* s, void* buf, int len)
    int vmdsock_selread(void* s, int timeout)
    int vmdsock_selwrite(void* s, int timeout)
    void vmdsock_shutdown(void* s)
    void vmdsock_destroy(void* s)