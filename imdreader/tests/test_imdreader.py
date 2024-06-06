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
import subprocess

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


@pytest.fixture(autouse=True)
def log_config():
    logger = logging.getLogger("imdreader.IMDREADER")
    file_handler = logging.FileHandler("test.log")
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.setLevel(logging.DEBUG)
    yield
    logger.removeHandler(file_handler)


@pytest.fixture(scope="function")
def run_gmx(tmpdir):
    port = get_free_port()
    command = [
        "gmx",
        "mdrun",
        "-s",
        TOPOL_TPR,
        "-imdwait",
        "-imdpull",
        "-imdterm",
        "-imdport",
        str(port),
    ]
    with tmpdir.as_cwd():
        with open("gmx_output.log", "w") as f:
            p = subprocess.Popen(
                command,
                stdin=subprocess.PIPE,
                stdout=f,
                stderr=f,
                text=True,
                bufsize=1,
            )
            try:
                yield port
            finally:
                # Terminate the process
                p.terminate()
                try:
                    p.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    logger.error(
                        "Process did not terminate in time, killing it."
                    )
                    p.kill()
                    p.wait()

                # Ensure all file descriptors are closed
                f.close()


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
    port = run_gmx
    recvuntil(
        "gmx_output.log",
        "IMD: Will wait until I have a connection and IMD_GO orders.",
        60,
    )
    u2 = mda.Universe(IMDGROUP_GRO, OUT_TRR)
    u = mda.Universe(
        IMDGROUP_GRO,
        f"localhost:{port}",
    )
    for ts in u.trajectory:
        pass

    assert len(u2.trajectory) == len(u.trajectory)


def test_pause(run_gmx, caplog):
    port = run_gmx
    recvuntil(
        "gmx_output.log",
        "IMD: Will wait until I have a connection and IMD_GO orders.",
        60,
    )
    u = mda.Universe(
        IMDGROUP_GRO,
        f"localhost:{port}",
        # 1240 bytes per frame
        buffer_size=62000,
    )
    for ts in u.trajectory:
        time.sleep(0.05)

    assert len(u.trajectory) == 101
    assert (
        "IMDProducer: Pausing simulation because buffer is almost full"
        in caplog.text
    )
    assert "IMDProducer: Unpausing simulation, buffer has space" in caplog.text
    assert "data likely lost in frame" not in caplog.text


def test_no_connection(caplog):
    u = mda.Universe(
        IMDGROUP_GRO,
        "localhost:8888",
        buffer_size=62000,
    )
    for ts in u.trajectory:
        with pytest.raises(ConnectionError):
            pass
    # NOTE: assert this in output: No connection received. Pausing simulation.
    assert "IMDProducer: Connection to localhost:8888 refused" in caplog.text


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
