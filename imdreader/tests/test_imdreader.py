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
from pwn import *

"""
mdout generated using 
gmx grompp -f md.mdp -c npt.gro -maxwarn 3

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
        "-o",
        "out.trr",
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


@pytest.fixture(scope="function")
def server_thread():
    server = DummyIMDServer()
    thread = threading.Thread(target=server.run)
    thread.start()
    yield
    server.disconnect()  # Ensure proper shutdown
    thread.join()


def test_comp_imd_xtc(run_gmx):
    # stdout, stderr = run_gmx.communicate()
    run_gmx.readuntil(
        "IMD: Will wait until I have a connection and IMD_GO orders."
    )

    u = mda.Universe(
        IMDGROUP_GRO, "localhost:8888", n_frames=11, num_atoms=36688
    )
    i = 0
    streampos = np.empty((11, 36688, 3), dtype=np.float32)
    for ts in u.trajectory:
        streampos[i] = ts.positions[:]
        i += 1
    u2 = mda.Universe(IMDGROUP_GRO, "out.trr")

    expected = streampos[0]
    actual = u2.trajectory[0].positions
    rtol = 1e-7
    atol = 0
    # Calculate the differences and the mask for elements that are not close
    differences = np.abs(actual - expected)
    not_close = (
        np.isclose(actual, expected, rtol=rtol, atol=atol, equal_nan=True)
        == False
    )
    print(u2.trajectory[0].triclinic_dimensions[0])
    print(u2.trajectory[0].triclinic_dimensions[1])
    print(u2.trajectory[0].triclinic_dimensions[2])
    if np.any(not_close):
        actual_diffs = actual[not_close]
        expected_diffs = expected[not_close]
        diff_values = differences[not_close]
        print("Differences found:")
        print("Indices with differences:", np.nonzero(not_close))
        print("Actual values:", actual_diffs)
        print("Expected values:", expected_diffs)
        print("Differences:", diff_values)
        assert False, f"Arrays differ by more than atol={atol} and rtol={rtol}"
    # print(run_gmx.recvlines(5))
    # print(stdout)

    # print(u)
    # print(u.trajectory)
    assert 1 == 1


def test_buffer_mgmt(run_gmx):
    run_gmx.readuntil(
        "IMD: Will wait until I have a connection and IMD_GO orders."
    )
    u = mda.Universe(
        IMDGROUP_GRO,
        "localhost:8888",
        n_frames=11,
        num_atoms=36688,
        # 4 frames, each 440296 bytes
        buffer_size=3522368,
    )
    for ts in u.trajectory:
        sleep(0.1)
    assert 1 == 1
