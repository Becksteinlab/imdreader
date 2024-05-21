Proposed Protocol v3
====================

The suggested changes to the protocol are as follows:

1. Via an ``.mdp`` file setting, the user should be able to specify which of positions, forces, and velocities are sent.
   These do not need to be rate-adjustable like with xtc and trr output settings- they will always be output vai IMD for every frame if
   turned on.

2. The IMD_HANDSHAKE packet should have a 1-byte body which contains the configuration settings the simulation was setup with.
   This will allow the client to be choose appropriate packets to send and receive without redundant configuration.

   The modified handshake packet would look like this:

.. code-block:: none

    Header: 
        4 (int32) (IMD_HANDSHAKE)
        2 (int32) (Protocol version, must be unswapped so client can determine endianness)
    Body:
        <val> (bit) (imdpull: true or false)
        <val> (bit) (imdwait: true or false)
        <val> (bit) (imdterm: true or false)
        <val> (bit) (positions included: true or false)
        <val> (bit) (velocities included: true or false)
        <val> (bit) (forces included: true or false)
        <val> (bit) (dimensions included: true or false)
        <val> (bit) (unused bit)

3. The server should wait longer than 1 second (possibly up to 60s) for the go signal so that the client 
   has plenty of time to allocate memory buffers based on the endianness and data type information it received in the handshake packet.

4. In the simulation loop, the server will send the client data in this order (iff the configuration says to send it)
    
    i. Dimension data (IMD_DIM) in triclinic vectors

    ii. Position data (IMD_FCOORDS)
    
    iii. Velocity data (IMD_VELS) in the same manner as positions
    
    iv. Force data (IMD_FORCES) in the same manner as positions

5. The server will no longer send positions adjusted for visualiziation purposes, they will instead be sent in the same manner as .xtc or .trr
   positions are written to files. The client will be responsible for this adjustment since the dimension data is now available.

6. The server will send a new IMD_EOS (end of stream) packet after the last frame is sent unless the client initiates the disconnection with
   IMD_DISCONNECT.


Questions
^^^^^^^^^

- Should the server provide the client with the number of atoms?
