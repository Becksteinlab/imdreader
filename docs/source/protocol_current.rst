
IMD Protocol as it is currently implemented in GROMACS
======================================================

Version: 2
Headersize: 8 bytes (two 32 bit integers)
Endianness: Little endian and big endian are both allowed for position and energy data. In headers, big endian is used except in 
            the handshake signal where the endianness of the server is preserved unswapped in the packet for the client to check.

Only one client is allowed per GROMACS server at a time

Header types
------------

GROMACS->Client Headers
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :widths: 10 30 90
   :header-rows: 1

   * - Type
     - Enum Integer Value
     - Description
   * - IMD_ENERGIES
     - 1
     - 
   * - IMD_FCOORDS
     - 2
     - 
   * - IMD_HANDSHAKE
     - 4
     -

Client->GROMACS Headers
^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :widths: 10 30 90
   :header-rows: 1

   * - Type
     - Enum Integer Value
     - Description
   * - IMD_DISCONNECT
     - 0
     - Tells the GROMACS server to disconnect its client socket and reset its IMD step to its default
   * - IMD_GO
     - 3
     - After receiving the handshake from GROMACS, the client has 1 second to send this signal to the server to begin
       receiving data. If the server does not receive this signal, it will disconnect the client. If the GROMACS server
       was started with ``-imdwait``, the simulation will not begin until this signal is received.
   * - IMD_KILL
     - 5
     - Tells GROMACS to stop the simulation (requires that the simulation was started with ``-imdterm``)
   * - IMD_MDCOMM
     - 6
     - 
   * - IMD_PAUSE
     - 7
     - After simulation has started, the client can send this signal to the GROMACS server to toggle the 
       simulation between paused and running states. This sends GROMACS into a blocking loop until
       the client sends another IMD_PAUSE signal.
   * - IMD_TRATE
     - 8
     - 

Unused headers
^^^^^^^^^^^^^^

.. list-table::
   :widths: 10 30 90
   :header-rows: 1

   * - Type
     - Enum Integer Value
     - Description
   * - IMD_IOERROR
     - 9
     - Used internally by GROMACS to signal an error in the I/O operations with the client. This is not sent by the client.
   * - Count
     - 10
     - Defined in GROMACS source code but unused

Protocol steps
--------------

1. Pre-connection
#################
GROMACS (Server) 1: Decide if simulation should wait for a client to connect to begin using option
GROMACS (Server) 2: Decide if IMD client should be able to terminate the simulation using option
GROMACS (Server) 3: Decide if pulling force data from the client is allowed using option
GROMACS (Server) 4: If waiting for client, block simulation run until GO signal received

2. Connection
#############
GROMACS (Server) 1: Create a TCP socket and bind it to a port
GROMACS (Server) 2: Listen for incoming connections, checking every 1 second
GROMACS (Server) 3: Accept the incoming connection, binding it to a new socket
GROMACS (Server) 4: Send the handshake signal to the client:

Header: 
    4 (IMD_HANDSHAKE)
    2 (Version, must be unswapped)

GROMACS (Server) 5: Wait up to 1 second to receive a GO signal from the client:

Header:
    3 (IMD_GO)
    <val> (Unused length attribute in header)

GROMACS (Server) 6: In every iteration of the md_do loop, first check for incoming packets from the client. These packets can be:

    Header:
        5 (IMD_KILL)
        <val> (Unused length attribute in header)

    Header:
        0 (IMD_DISCONNECT)
        <val> (Unused length attribute in header)

    Header:
        6 (IMD_MDCOMM)
        <val> (Number of forces that will be sent in the packet)

        Force packet:
            <count> of 32 byte force integers representing indices of atoms to apply force to
            <val> * 3 number of 32 byte force floats representing the force to apply to the atoms at 
            the corresponding indices

    Header:
        7 (IMD_PAUSE)
        <val> (Unused length attribute in header)

    Header:
        8 (IMD_TRATE)
        <val> (New transfer rate. Value of 0 means reset to default)

    Any other header sent will disconnect the client and print an error message but not stop the simulation.

Next, send energies and position data to the client IF this integration step lands in the rate step specified by the client (or the default, every step)

    Header:
        1 (IMD_ENERGIES)
        1 (Contstant value)
    Header:
        2 (IMD_FCOORDS)
        <val> (Number of atoms in the system)

        Position packet:




3. Disconnection
################
GROMACS (Server) 1: Shutdown the socket if it isn't already shutdown
GROMACS (Server) 2: Set the IMD frequency to default and begin listening for connections again

# If connecting a second time, handshake must be sent a second time (test this)
