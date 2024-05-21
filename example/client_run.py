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

OUTPUT_FILE = "gmx_output.log"

import os

os.chdir("/home/law/workspace/imdreader/example/imdexample/testrun")


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
    with open(OUTPUT_FILE, "w") as f:
        p = subprocess.Popen(
            command,
            stdin=subprocess.PIPE,
            stdout=f,
            stderr=subprocess.STDOUT,  # Redirect stderr to stdout
            text=True,
            bufsize=1,
        )
    return p


def stream_subprocess_output(process):
    while True:
        output = process.stdout.readline()
        if output == "" and process.poll() is not None:
            break
        if output:
            print(output.strip())
    rc = process.poll()
    return rc


run_gmx_process = run_gmx()

# Wait for the process to complete
run_gmx_process.wait()
