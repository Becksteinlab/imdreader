Proposed Protocol v3
====================

New packets
-----------

There will be 5 new packets:

1. (MDAnalysis -> Gromacs) IMD_V3SWITCH, sent immediately after IMD_GO to specify IMDv3 should be used

.. code-block:: none

   Header: 
      Value:
         10 (uint32) (IMD_V3SWITCH)
      Length:
         <val> (bit) Send IMD_TIME packets [bool]
         <val> (bit) Send IMD_ENERGIES packets [bool]
         <val> (bit) Send IMD_BOX packets [bool]
         <val> (bit) Send IMD_FCOORDS packets [bool]
         <val> (bit) Send IMD_VELOCITIES packets [bool]
         <val> (bit) Send IMD_FORCES packets [bool]
         <val> (bit) Wrap coords [bool, does nothing if FCOORDS not sent]
         0 (25 bits) Unused

2. (MDAnalysis -> Gromacs) IMD_RESUME, used to unpause a running simulation.

.. code-block:: none

   Header:
      Value: 
         11 (uint32) (IMD_RESUME)
      Length:
         <val> (unused length attribute)


3. (Gromacs -> MDAnalysis) IMD_TIME

.. code-block:: none

   Header:
      Value:
         12 (uint32) (IMD_TIME)
      Length:
         <val> (float32) (dt for the simulation)

   Body:
      <val> (float32) (Current time in picoseconds)


4. (Gromacs -> MDAnalysis) IMD_VELOCITIES

.. code-block:: none

   Header:
      Value:
         13 (uint32) (IMD_VELOCITIES)
      Length:
         <val> (uint32) (Number of atoms in the system)

   Body:
      <val> (float32) (n atoms * 3 values describing the velocities of each atom in the system in angstroms/picosecond)

5. (Gromacs -> MDAnalysis) IMD_FORCES

.. code-block:: none

   Header:
      Value:
         14 (uint32) (IMD_FORCES)
      Length:
         <val> (uint32) (Number of atoms in the system)
   Body:
      <val> (float32) (n atoms * 3 values describing the forces of each atom in the system in kilojoules/(mol*angstrom))



Units
-----

.. list-table::
   :widths: 10 20
   :header-rows: 1

   * - Measurement
     - Unit
   * - Length
     - Angstroms
   * - Velocity
     - Angstrom/picosecond
   * - Force
     - kilojoules/(mol*angstrom)

Packet order
------------

Data packets are always sent in this order, if present.

1. IMD_TIME
2. IMD_ENERGIES
3. IMD_BOX
4. IMD_FCOORDS
5. IMD_VELOCITIES
6. IMD_FORCES

For example, if the switch to send IMD_TIME was off in IMD_V3SWITCH, the resulting data packet order would be the same
except starting at 2.

Idempotency
-----------

If the IMD_V3SWITCH has been sent, making the simulation an IMDV3 simulation, IMD_PAUSE becomes an idempotent operation; 
sending it more than once has the same effect as sending it once. The only way to unpause a paused IMDV3 simulaton is to send
an IMD_RESUME packet, which is also idempotent.