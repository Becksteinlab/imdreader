"""
Unit and regression test for the imdreader package.
"""

# Import package, test suite, and other packages as needed
import imdreader
from imdreader.tests.datafiles import *
from .utils import *

import pytest
from numpy.testing import assert_allclose

import MDAnalysis as mda
from MDAnalysisTests.coordinates.base import assert_timestep_almost_equal

import sys
import threading
import logging
from pwn import *

"""
mdout generated using 
gmx grompp -f md.mdp -c argon_start.pdb -p argon.top

gmx mdrun -s topol.tpr

"""

"""
"-imdwait", "-imdpull", "-imdterm"
"""
"""
gmx output w imd:
IMD: Enabled. This simulation will accept incoming IMD connections.
IMD: Pausing simulation while no IMD connection present (-imdwait).
IMD: Allow termination of the simulation from IMD client (-imdterm).
IMD: Pulling from IMD remote is enabled (-imdpull).
IMD: Setting port for connection requests to 8888.
IMD: Setting up incoming socket.
IMD: Listening for IMD connection on port 8888.
IMD: Will wait until I have a connection and IMD_GO orders.
"""


@pytest.fixture()
def run_gmx(tmpdir):
    command = [
        "gmx",
        "mdrun",
        "-s",
        TOPOL_TPR,
        "-imdwait",
        "-imdpull",
        "-imdterm",
    ]
    with tmpdir.as_cwd():
        p = process(
            command,
        )
        try:
            # Yield the process to the test function; control resumes here after test ends
            yield p
        finally:
            # Ensure that the process is killed when the test ends, regardless of the result
            p.terminate()  # Sends a terminate signal (SIGTERM)
            try:
                p.wait(
                    timeout=10
                )  # Give some time for the process to terminate gracefully
            except TimeoutError:
                p.kill()  # Force kill if it does not terminate within the timeout
                p.wait()  # Wait again to ensure it's cleaned up


# NOTE: This test can't pass until dimensions are implemented in IMD 2.0
# def test_comp_imd_xtc(run_gmx):
#    # stdout, stderr = run_gmx.communicate()
#    run_gmx.readuntil(
#        "IMD: Will wait until I have a connection and IMD_GO orders."
#    )
#    u2 = mda.Universe(IMDGROUP_GRO, OUT_TRR)
#    u = mda.Universe(IMDGROUP_GRO, "localhost:8888", num_atoms=36688)
#    i = 0
#    streampos = np.empty((11, 36688, 3), dtype=np.float32)
#    for ts in u.trajectory:
#        u.atoms.wrap(
#            box=u2.trajectory[i].dimensions,
#            inplace=True,
#        )
#        streampos[i] = ts.positions[:]
#        i += 1
#
#    expected = streampos[0]
#    actual = u2.trajectory[0].positions
#    rtol = 1e-4
#    atol = 0
#    # Calculate the differences and the mask for elements that are not close
#    differences = np.abs(actual - expected)
#    not_close = (
#        np.isclose(actual, expected, rtol=rtol, atol=atol, equal_nan=True)
#        == False
#    )
#    print(u2.trajectory[0].dimensions)
#    if np.any(not_close):
#        actual_diffs = actual[not_close]
#        expected_diffs = expected[not_close]
#        diff_values = differences[not_close]
#        print("Differences found:")
#        print("Indices with differences:", np.nonzero(not_close))
#        print("Actual values:", actual_diffs)
#        print("Expected values:", expected_diffs)
#        print("Differences:", diff_values)
#        assert False, f"Arrays differ by more than atol={atol} and rtol={rtol}"
#    # print(run_gmx.recvlines(5))
#    # print(stdout)
#
#    # print(u)
#    # print(u.trajectory)
#    assert 1 == 1


def test_traj_len(run_gmx):
    run_gmx.readuntil(
        "IMD: Will wait until I have a connection and IMD_GO orders."
    )
    u2 = mda.Universe(IMDGROUP_GRO, OUT_TRR)
    u = mda.Universe(
        IMDGROUP_GRO,
        "localhost:8888",
        num_atoms=100,
    )
    for ts in u.trajectory:
        # Wait to simulate analysis
        sleep(0.1)

    assert len(u2.trajectory) == len(u.trajectory)


# NOTE: this test passes because the
# buffer threshold is set for 10% since the lennard jones simulation is so quick
# make sure to change this value in run() before actual use
def test_pause(run_gmx):

    # NOTE: assert this in output: Un-pause command received.
    # Provide a buffer small enough to force pausing the simulation
    run_gmx.readuntil(
        "IMD: Will wait until I have a connection and IMD_GO orders.",
        timeout=10,
    )
    u = mda.Universe(
        IMDGROUP_GRO,
        "localhost:8888",
        num_atoms=100,
        # 1240 bytes per frame
        buffer_size=62000,
    )
    for ts in u.trajectory:
        sleep(0.1)

    assert len(u.trajectory) == 101


# NOTE: This will fail before an exception queing system is setup
# def test_no_connection():
#
#    # assert this in output: Un-pause command received.
#    # Provide a buffer small enough to force pausing the simulation
#
#    u = mda.Universe(
#        IMDGROUP_GRO,
#        "localhost:8888",
#        num_atoms=100,
#        buffer_size=62000,
#    )
#    for ts in u.trajectory:
#        sleep(0.1)
#
#    assert len(u.trajectory) == 100

"""
import socket
import struct

conn = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

def connect():
    conn.connect(("localhost",8888))
    conn.recv(8)
    go = struct.pack("!ii", 3, 0)
    conn.sendall(go)

def pause():
    pause = struct.pack("!ii", 7, 0)
    conn.sendall(pause)

"""
