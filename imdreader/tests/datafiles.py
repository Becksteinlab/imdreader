"""
Location of data files
======================

Use as ::

    from imdreader.tests.datafiles import *

"""

__all__ = [
    "NPT_GRO",
    "MD_MDP",
    "MDOUT_MDP",
    "TOPOL_TOP",
    "TOPOL_TPR",
    "IMDGROUP_GRO",
    "OUT_TRR",
]

from importlib import resources
from pathlib import Path

_data_ref = resources.files("imdreader.data")

NPT_GRO = (_data_ref / "npt.gro").as_posix()
MD_MDP = (_data_ref / "md.mdp").as_posix()
MDOUT_MDP = (_data_ref / "md.mdp").as_posix()
TOPOL_TOP = (_data_ref / "topol.top").as_posix()
TOPOL_TPR = (_data_ref / "topol.tpr").as_posix()
IMDGROUP_GRO = (_data_ref / "imdgroup.gro").as_posix()
OUT_TRR = (_data_ref / "out.trr").as_posix()

del resources
