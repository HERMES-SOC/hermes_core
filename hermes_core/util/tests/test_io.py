"""Tests for Loading and Saving data from data containers"""

from collections import OrderedDict
from pathlib import Path
import datetime
import pytest
import numpy as np
from numpy.random import random
import tempfile
from astropy.timeseries import TimeSeries
from astropy.time import Time
from astropy.units import Quantity
from spacepy.pycdf import CDFError
from hermes_core.timedata import TimeData


def get_test_timedata():
    ts = TimeSeries()
    ts.meta.update(
        {
            "Descriptor": "EEA>Electron Electrostatic Analyzer",
            "Data_level": "l1>Level 1",
            "Data_version": "v0.0.1",
            "Start_time": datetime.datetime.now(),
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
            "VAR_TYPE": "metadata",
            "CATDESC": "Test Metadata",
        }
    )
    timedata = TimeData(data=ts)
    return timedata


def test_cdf_io():
    # Get Test Datas
    td = get_test_timedata()

    with tempfile.TemporaryDirectory() as tmpdirname:
        # Convert TimeData the to a CDF File
        test_file_output_path = td.save(output_path=tmpdirname, file_extension=".cdf")

        # Load the CDF to a TimeData Object
        td_loaded = TimeData.load(test_file_output_path)

        assert td.shape == td_loaded.shape

        with pytest.raises(CDFError):
            td_loaded.save(output_path=tmpdirname, file_extension=".cdf")


def test_json_io():
    # Get Test Datas
    td = get_test_timedata()

    with tempfile.TemporaryDirectory() as tmpdirname:
        # Convert TimeData the to a JSON File
        test_file_output_path = td.save(output_path=tmpdirname, file_extension=".json")

        # Load the JSON to a TimeData Object
        td_loaded = TimeData.load(test_file_output_path)

        assert td.shape == td_loaded.shape

        with pytest.raises(FileExistsError):
            td_loaded.save(output_path=tmpdirname, file_extension=".json")


def test_csv_io():
    # Get Test Datas
    td = get_test_timedata()

    with tempfile.TemporaryDirectory() as tmpdirname:
        # Convert TimeData the to a JSON File
        test_file_output_path = td.save(output_path=tmpdirname, file_extension=".csv")
        # Load the JSON to a TimeData Object
        td_loaded = TimeData.load(test_file_output_path)

        assert td.shape == td_loaded.shape

        print(test_file_output_path)
        with pytest.raises(FileExistsError):
            td_loaded.save(output_path=tmpdirname, file_extension=".csv")
