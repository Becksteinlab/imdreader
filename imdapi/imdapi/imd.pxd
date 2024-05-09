# distutils: language = c

cdef extern from "limits.h":
    pass

# Determine the appropriate integer type for int32 based on the size
cdef extern from *:
    """
    #if (INT_MAX == 2147483647)
    typedef int     int32;
    #else
    typedef short   int32;
    #endif
    """
    ctypedef int int32

cdef extern from "imd.h":
    # Enumeration for IMD control and data messages
    ctypedef enum IMDType:
        IMD_DISCONNECT,   # close IMD connection, leaving simulation running
        IMD_ENERGIES,     # energy data block
        IMD_FCOORDS,      # atom coordinates
        IMD_GO,           # start the simulation
        IMD_HANDSHAKE,    # endianism and version check message
        IMD_KILL,         # kill the simulation job, shutdown IMD
        IMD_MDCOMM,       # MDComm style force data
        IMD_PAUSE,        # pause the running simulation
        IMD_TRATE,        # set IMD update transmission rate
        IMD_IOERROR       # indicate an I/O error

    # Struct for energy data
    ctypedef struct IMDEnergies:
        int32 tstep       # integer timestep index
        float T           # Temperature in degrees Kelvin
        float Etot        # Total energy, in Kcal/mol
        float Epot        # Potential energy, in Kcal/mol
        float Evdw        # Van der Waals energy, in Kcal/mol
        float Eelec       # Electrostatic energy, in Kcal/mol
        float Ebond       # Bond energy, Kcal/mol
        float Eangle      # Angle energy, Kcal/mol
        float Edihe       # Dihedral energy, Kcal/mol
        float Eimpr       # Improper energy, Kcal/mol

    # Function declarations
    int imd_disconnect(void*)
    int imd_pause(void*)
    int imd_kill(void*)
    int imd_handshake(void*)
    int imd_trate(void*, int32)

    int imd_send_mdcomm(void*, int32, const int32*, const float*)
    int imd_send_energies(void*, const IMDEnergies*)
    int imd_send_fcoords(void*, int32, const float*)

    int imd_recv_handshake(void*)
    IMDType imd_recv_header(void*, int32*)
    int imd_recv_mdcomm(void*, int32, int32*, float*)
    int imd_recv_energies(void*, IMDEnergies*)
    int imd_recv_fcoords(void*, int32, float*)