"""Tests for Loading and Saving data from data containers"""

from collections import OrderedDict
from pathlib import Path
import pytest
import json
import numpy as np
from numpy.random import random
import tempfile
from astropy.timeseries import TimeSeries
from astropy.time import Time
from astropy.units import Quantity
from spacepy.pycdf import CDFError, CDF
from hermes_core.timedata import HermesData


def get_test_hermes_data():
    ts = TimeSeries()
    ts.meta.update(
        {
            "Descriptor": "EEA>Electron Electrostatic Analyzer",
            "Data_level": "l1>Level 1",
            "Data_version": "v0.0.1",
            "MODS": [
                "v0.0.0 - Original version.",
                "v1.0.0 - Include trajectory vectors and optics state.",
                "v1.1.0 - Update metadata: counts -> flux.",
                "v1.2.0 - Added flux error.",
                "v1.3.0 - Trajectory vector errors are now deltas.",
            ],
        }
    )

    # Create an astropy.Time object
    time = np.arange(10)
    time_col = Time(time, format="unix")
    ts["time"] = time_col

    # Add Measurement
    quant = Quantity(value=random(size=(10)), unit="m", dtype=np.uint16)
    ts["measurement"] = quant
    ts["measurement"].meta = OrderedDict(
        {
            "VAR_TYPE": "data",
            "CATDESC": "Test Data",
        }
    )
    hermes_data = HermesData(data=ts)
    return hermes_data


def test_cdf_io():
    """Test CDF IO Handler on Default Data"""
    # Get Test Datas
    td = get_test_hermes_data()

    with tempfile.TemporaryDirectory() as tmpdirname:
        # Convert HermesData the to a CDF File
        test_file_output_path = td.save(output_path=tmpdirname)

        # Load the CDF to a HermesData Object
        td_loaded = HermesData.load(test_file_output_path)

        assert td.shape == td_loaded.shape

        with pytest.raises(CDFError):
            td_loaded.save(output_path=tmpdirname)


def test_cdf_bad_file_path():
    """Test Loading CDF from a non-existant file"""
    with tempfile.TemporaryDirectory() as tmpdirname:
        # Try loading from non-existant_path
        with pytest.raises(FileNotFoundError):
            _ = HermesData.load(tmpdirname + "non_existant_file.cdf")


def test_cdf_nrv_support_data():
    """
    Test Loading Non-Record-Varying data with CDF IO Handler
    """
    # Get Test Datas
    td = get_test_hermes_data()

    with tempfile.TemporaryDirectory() as tmpdirname:
        # Convert HermesData the to a CDF File
        test_file_output_path = td.save(output_path=tmpdirname)

        # Load the JSON file as JSON
        with CDF(test_file_output_path, readonly=False) as cdf_file:
            # Add Non-Record-Varying Variable
            cdf_file["Test_NRV_Var"] = [1, 2, 3]

            # Add Support Data Variable
            cdf_file["Test_Support_Var"] = np.arange(10)
            cdf_file["Test_Support_Var"].meta["UNITS"] = "counts"

        # Make sure we can load the modified JSON
        td_loaded = HermesData.load(test_file_output_path)

        assert "Test_NRV_Var" in td_loaded.support
        assert "Test_Support_Var" in td_loaded.columns
