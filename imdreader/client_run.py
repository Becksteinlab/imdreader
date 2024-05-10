import imdreader
import MDAnalysis as mda
from MDAnalysisTests.datafiles import PSF
import os


u = mda.Universe(
    "example/imdexample/selected_atoms.pdb",
    "localhost:8888",
    n_frames=100,
    num_atoms=1789,
    format="IMD",
)
for ts in u.trajectory:
    print(ts.positions)
