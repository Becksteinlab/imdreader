import MDAnalysis as mda
from MDAnalysisTests.datafiles import COORDINATES_TOPOLOGY, COORDINATES_H5MD

u = mda.Universe(COORDINATES_TOPOLOGY, COORDINATES_H5MD)

print(u.atoms.positions)
u.trajectory[2]
print(u.atoms.positions)

for ts in u.trajectory:
    pass
