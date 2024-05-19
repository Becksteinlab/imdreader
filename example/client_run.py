import imdreader
from imdreader.tests.datafiles import *

import pytest
from numpy.testing import assert_allclose

import MDAnalysis as mda
from MDAnalysisTests.coordinates.base import assert_timestep_almost_equal

import sys
import threading
import logging
from pwn import *
import tempfile


def run_gmx():
    command = [
        "gmx",
        "mdrun",
        "-s",
        TOPOL_TPR,
        "-imdwait",
        "-imdpull",
        "-imdterm",
    ]

    p = process(
        command,
    )
    return p


def test_buffer_mgmt(run_gmx):
    run_gmx.readuntil(
        "IMD: Will wait until I have a connection and IMD_GO orders."
    )
    u = mda.Universe(
        IMDGROUP_GRO, "localhost:8888", num_atoms=100, buffer_size=12400
    )
    for ts in u.trajectory:
        sleep(0.1)
        print(ts)
    assert 1 == 1


run_gmx = run_gmx()
test_buffer_mgmt(run_gmx)
