"""Tests for the sunsensor"""
import pytest
import os.path

import numpy as np

import hermes_core
import hermes_core.calibration as calib
from hermes_core.util.util import create_science_filename, parse_science_filename

SUNSENSOR_VECT_TEST_FILE = os.path.join(
    hermes_core.__path__[0], "tests/data/hermes_SS_l0_2022339-000000_v0.bin"
)
# check that the test file exists
os.path.exists(SUNSENSOR_VECT_TEST_FILE)

SUNSENSOR_NUM_VECTORS = 20


def test_sunsensor_filename():
    """ "Test that the test file has the correct format."""
    file_metadata = parse_science_filename(SUNSENSOR_VECT_TEST_FILE)
    assert isinstance(file_metadata, dict)


def test_parse_sensor_vector_packets():
    """Test that parsing function returns the expected data."""
    data = calib.parse_sunsensor_vector_packets(SUNSENSOR_VECT_TEST_FILE)

    expected_num_packets = 10
    # should have 10 packets
    assert len(data["SECONDS"]) == expected_num_packets

    # test the contents of the file
    assert (data["SECONDS"] == np.arange(expected_num_packets)).all()
    assert (data["SUBSECONDS"] == np.zeros(expected_num_packets)).all()
    assert (data["NUMVECTORSCONTAINED"] == 20 * np.ones(expected_num_packets)).all()
    assert (data["VECTORPOLLINTERVAL"] == np.zeros(expected_num_packets)).all()
    assert (data["TEMPERATURE"] == 25 * np.ones(expected_num_packets)).all()
    assert (data["PROCESSINGTIME"] == 1000 * np.ones(expected_num_packets)).all()

    vector_data = np.zeros(
        (3, SUNSENSOR_NUM_VECTORS, expected_num_packets), dtype="uint16"
    )

    for i in range(expected_num_packets):
        vector_data[0, :, i] = [
            data[f"X_INT_{j:02}"][i] for j in range(SUNSENSOR_NUM_VECTORS)
        ]
        vector_data[1, :, i] = [
            data[f"Y_INT_{j:02}"][i] for j in range(SUNSENSOR_NUM_VECTORS)
        ]
        vector_data[2, :, i] = [
            data[f"Z_INT_{j:02}"][i] for j in range(SUNSENSOR_NUM_VECTORS)
        ]

    for i in range(expected_num_packets):
        assert (vector_data[0, :, 0] == np.arange(SUNSENSOR_NUM_VECTORS)).all()
        assert (vector_data[1, :, 0] == np.arange(SUNSENSOR_NUM_VECTORS) + 1).all()
        assert (vector_data[2, :, 0] == np.arange(SUNSENSOR_NUM_VECTORS) + 2).all()
