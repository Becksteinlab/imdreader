from MDAnalysisTests.datafiles import (
    COORDINATES_TOPOLOGY,
    COORDINATES_TRR,
    COORDINATES_H5MD,
)
import MDAnalysis as mda
import imdreader
from imdreader.IMDClient import imdframe_memsize
from .utils import (
    IMDServerEventType,
    DummyIMDServer,
    get_free_port,
    ExpectPauseLoopV2Behavior,
    create_default_imdsinfo_v2,
)
from .server import TestIMDServer
from MDAnalysisTests.coordinates.base import (
    MultiframeReaderTest,
    BaseReference,
    BaseWriterTest,
    assert_timestep_almost_equal,
)
from MDAnalysisTests.coordinates.test_xdr import TRRReference
import numpy as np
import logging
import pytest
import time


# NOTE: removeme after initial testing
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


logger = logging.getLogger("imdreader.IMDREADER")

IMDENERGYKEYS = [
    "step",
    "temperature",
    "total_energy",
    "potential_energy",
    "van_der_walls_energy",
    "coulomb_energy",
    "bonds_energy",
    "angles_energy",
    "dihedrals_energy",
    "improper_dihedrals_energy",
]


class TestIMDReaderV2:

    @pytest.fixture
    def port(self):
        return get_free_port()

    @pytest.fixture
    def traj(self):
        return mda.coordinates.H5MD.H5MDReader(
            COORDINATES_H5MD, convert_units=False
        )

    @pytest.fixture
    def ref(self):
        return mda.coordinates.H5MD.H5MDReader(
            COORDINATES_H5MD, convert_units=False
        )

    @pytest.fixture
    def server(self, traj):
        server = DummyIMDServer(traj, 2)
        return server

    @pytest.fixture(params=[">", "<"])
    def setup_test_endianness_traj_unchanged(self, request, server, port):
        server.port = port
        server.imdsessioninfo.endianness = request.param
        server.start()
        server.wait_for_event(IMDServerEventType.LISTENING)
        return server, port

    def test_endianness_traj_unchanged(
        self, setup_test_endianness_traj_unchanged, ref
    ):
        _, port = setup_test_endianness_traj_unchanged

        reader = imdreader.IMDREADER.IMDReader(
            f"localhost:{port}",
            convert_units=False,
            n_atoms=ref.trajectory.n_atoms,
        )

        i = 0
        # Can't call assert in loop- this prevents reader's __exit__ from being called
        # if assert fails. Instead copy timesteps and then assert them
        timesteps = []

        for ts in reader:
            logger.debug(
                f"test_imdreader: positions for frame {i}: {ts.positions}"
            )
            timesteps.append(ts.copy())
            i += 1

        assert i == len(ref)

        for j in range(len(ref)):
            np.testing.assert_allclose(timesteps[j].positions, ref[j].positions)
            offset = 0
            for energy_key in IMDENERGYKEYS:
                assert timesteps[j].data[energy_key] == j + offset
                offset += 1

    @pytest.fixture
    def setup_test_pause_traj_unchanged(self, server, port):
        server.port = port
        server.loop_behavior = ExpectPauseLoopV2Behavior()
        server.start()
        server.wait_for_event(IMDServerEventType.LISTENING)
        return server, port

    def test_pause_traj_unchanged(self, setup_test_pause_traj_unchanged, ref):
        server, port = setup_test_pause_traj_unchanged

        # Give the buffer only 1 IMDFrame of memory
        # We expect the producer thread to have to
        # pause every frame (except the first)
        reader = imdreader.IMDREADER.IMDReader(
            f"localhost:{port}",
            convert_units=False,
            n_atoms=ref.trajectory.n_atoms,
            buffer_size=imdframe_memsize(
                ref.trajectory.n_atoms, server.imdsessioninfo
            ),
        )

        i = 0
        timesteps = []

        for ts in reader:
            time.sleep(1)
            timesteps.append(ts.copy())
            i += 1

        assert i == len(ref)

        for j in range(len(ref)):
            np.testing.assert_allclose(timesteps[j].positions, ref[j].positions)
            offset = 0
            for energy_key in IMDENERGYKEYS:
                assert timesteps[j].data[energy_key] == j + offset
                offset += 1

    def test_no_connection(self):
        with pytest.raises(ConnectionRefusedError):
            imdreader.IMDREADER.IMDReader("localhost:12345", n_atoms=1)


class TestIMDReaderWithBlockingServerV2:

    @pytest.fixture
    def port(self):
        return get_free_port()

    @pytest.fixture
    def traj(self):
        return mda.coordinates.H5MD.H5MDReader(
            COORDINATES_H5MD, convert_units=False
        )

    @pytest.fixture
    def ref(self):
        return mda.coordinates.H5MD.H5MDReader(
            COORDINATES_H5MD, convert_units=False
        )

    @pytest.fixture
    def server(self, traj):
        server = TestIMDServer(traj)
        yield server
        server.cleanup()

    @pytest.mark.parametrize("endianness", [">", "<"])
    def test_change_endianness_traj_unchanged(self, ref, server, endianness):
        imdsinfo = create_default_imdsinfo_v2()
        imdsinfo.endianness = endianness

        host = "localhost"
        port = get_free_port()

        # This also sends first frame to prevent blocking
        # in reader's init
        server.listen_accept_handshake_send_ts(host, port, imdsinfo)
        reader = imdreader.IMDREADER.IMDReader(
            f"localhost:{port}",
            n_atoms=ref.trajectory.n_atoms,
            convert_units=False,
        )

        i = 0
        timesteps = []
        for ts in reader:
            if i != 4:
                server.send_frame(i + 1, endianness=endianness)
            if i == 4:
                server.disconnect()
            timesteps.append(ts.copy())
            i += 1

        assert i == len(ref)
