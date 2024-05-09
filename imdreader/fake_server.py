import imdreader
import MDAnalysis as mda
from MDAnalysisTests.datafiles import PSF


connection = "localhost:8000"
u = mda.Universe(
    "imdapi/imdapi/testimd.pdb",
    connection,
    format="STREAM",
    n_frames=2,
    num_atoms=1,
)
for ts in u.trajectory:
    print(ts.positions)
