# imd.pxd
from libc.stdint cimport int32_t
from libc.stdfloat cimport float

# Assuming the definition of int32 as it depends on the architecture
ctypedef int32_t int32

# Enumeration for IMD message types
cdef enum IMDType:
    IMD_DISCONNECT   # close IMD connection, leaving sim running
    IMD_ENERGIES     # energy data block
    IMD_FCOORDS      # atom coordinates
    IMD_GO           # start the simulation
    IMD_HANDSHAKE    # endianism and version check message
    IMD_KILL         # kill the simulation job, shutdown IMD
    IMD_MDCOMM       # MDComm style force data
    IMD_PAUSE        # pause the running simulation
    IMD_TRATE        # set IMD update transmission rate
    IMD_IOERROR      # indicate an I/O error

# Structure for holding energy data
cdef struct IMDEnergies:
    int32 tstep      # integer timestep index
    float T          # Temperature in degrees Kelvin
    float Etot       # Total energy, in Kcal/mol
    float Epot       # Potential energy, in Kcal/mol
    float Evdw       # Van der Waals energy, in Kcal/mol
    float Eelec      # Electrostatic energy, in Kcal/mol
    float Ebond      # Bond energy, Kcal/mol
    float Eangle     # Angle energy, Kcal/mol
    float Edihe      # Dihedral energy, Kcal/mol
    float Eimpr      # Improper energy, Kcal/mol

# Declarations of external C functions
cdef extern from "imd.h":
    int imd_disconnect(void* s)
    int imd_pause(void* s)
    int imd_kill(void* s)
    int imd_handshake(void* s)
    int imd_trate(void* s, int32 rate)

    int imd_send_mdcomm(void* s, int32 n, const int32* indices, const float* forces)
    int imd_send_energies(void* s, const IMDEnergies* energies)
    int imd_send_fcoords(void* s, int32 n, const float* coords)

    int imd_recv_handshake(void* s)
    IMDType imd_recv_header(void* s, int32* length)
    int imd_recv_mdcomm(void* s, int32 n, int32* indices, float* forces)
    int imd_recv_energies(void* s, IMDEnergies* energies)
    int imd_recv_fcoords(void* s, int32 n, float* coords)