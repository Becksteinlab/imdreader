from MDAnalysisTests.datafiles import COORDINATES_TOPOLOGY, COORDINATES_TRR
import MDAnalysis as mda
from .utils import DummyIMDServer
from MDAnalysisTests.coordinates.base import (
    MultiframeReaderTest,
    BaseReference,
    BaseWriterTest,
    assert_timestep_almost_equal,
)
import numpy as np
import logging
import pytest


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


def test_traj_unchanged():
    u = mda.Universe(COORDINATES_TOPOLOGY, COORDINATES_TRR)
    u_imd = mda.Universe(COORDINATES_TOPOLOGY, "localhost:8888")
    server = DummyIMDServer()
    server.start()

    i = 0
    for ts in u_imd.trajectory:
        print(ts.dimensions)
        np.testing.assert_allclose(ts.positions, u.trajectory[i].positions)
        i += 1

    assert i == len(u.trajectory)
