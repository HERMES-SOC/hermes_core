"""Tests for Loading and Saving data from data containers"""

from collections import OrderedDict
from pathlib import Path
import datetime
import pytest
import json
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
    timedata = TimeData(data=ts)
    return timedata


def test_cdf_io():
    """Test CDF IO Handler on Default Data"""
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
    """Test JSON IO Handler on Default Data"""
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


def test_json_bad_file_path():
    """Test Loading JSON from a non-existant file"""
    with tempfile.TemporaryDirectory() as tmpdirname:
        # Try loading from non-existant_path
        with pytest.raises(FileNotFoundError):
            _ = TimeData.load(tmpdirname + "non_existant_file.json")


def test_json_list_data():
    """
    Test Loading Data that is a List of Data objects.
    NOTE: This is how SPDF outputs it's JSON data.
    """
    # Get Test Datas
    td = get_test_timedata()

    with tempfile.TemporaryDirectory() as tmpdirname:
        # Convert TimeData the to a JSON File
        test_file_output_path = td.save(output_path=tmpdirname, file_extension=".json")

        json_data = {}
        # Load the JSON file as JSON
        with open(test_file_output_path) as json_file:
            json_data = json.load(json_file)

        # Induce JSON List of Data
        mod_json_data = [json_data]

        # Save the Modified JSON
        with open(test_file_output_path, "w") as json_file:
            json.dump(mod_json_data, json_file)

        # Make sure we can load the modified JSON
        td_loaded = TimeData.load(test_file_output_path)

        assert len(td) == len(td_loaded)
        assert td.shape == td_loaded.shape
        assert len(td.meta) == len(td_loaded.meta)

        for attr in td.meta:
            assert attr in td_loaded.meta

        for var in td.columns:
            assert var in td_loaded.columns
            assert len(td[var]) == len(td_loaded[var])
            assert len(td[var].meta) == len(td_loaded[var].meta)
            assert td[var].meta["VAR_TYPE"] == td_loaded[var].meta["VAR_TYPE"]


def test_json_nrv_data():
    """
    Test Loading Non-Record-Varying data with JSON IO Handler
    """
    # Get Test Datas
    td = get_test_timedata()

    with tempfile.TemporaryDirectory() as tmpdirname:
        # Convert TimeData the to a JSON File
        test_file_output_path = td.save(output_path=tmpdirname, file_extension=".json")

        json_data = {}
        # Load the JSON file as JSON
        with open(test_file_output_path) as json_file:
            json_data = json.load(json_file)

        # Add Non-Record-Varying Variable
        json_data["Test_NRV_Var"] = {"DAT": [1, 2, 3]}

        # Save the Modified JSON
        with open(test_file_output_path, "w") as json_file:
            json.dump(json_data, json_file)

        # Make sure we can load the modified JSON
        td_loaded = TimeData.load(test_file_output_path)

        assert "Test_NRV_Var" not in td_loaded


def test_json_multi_dim_data():
    """
    Test Loading Multi-Dimensional data with JSON IO Handler
    """
    # Get Test Datas
    td = get_test_timedata()

    with tempfile.TemporaryDirectory() as tmpdirname:
        # Convert TimeData the to a JSON File
        test_file_output_path = td.save(output_path=tmpdirname, file_extension=".json")

        json_data = {}
        # Load the JSON file as JSON
        with open(test_file_output_path) as json_file:
            json_data = json.load(json_file)

        # Add Non-Record-Varying Variable
        json_data["Test_multi_dim_Var"] = {"DAT": random(size=(10, 2)).tolist()}

        # Save the Modified JSON
        with open(test_file_output_path, "w") as json_file:
            json.dump(json_data, json_file)

        # Make sure we can load the modified JSON
        td_loaded = TimeData.load(test_file_output_path)

        assert "Test_multi_dim_Var" not in td_loaded


def test_csv_io():
    """Test CSV IO Handler on Default Data"""
    # Get Test Datas
    td = get_test_timedata()

    with tempfile.TemporaryDirectory() as tmpdirname:
        # Convert TimeData the to a JSON File
        test_file_output_path = td.save(output_path=tmpdirname, file_extension=".csv")
        # Load the JSON to a TimeData Object
        td_loaded = TimeData.load(test_file_output_path)

        assert td.shape == td_loaded.shape

        with pytest.raises(FileExistsError):
            td_loaded.save(output_path=tmpdirname, file_extension=".csv")


def test_csv_bad_file_path():
    """Test Loading CSV from a non-existant file"""
    with tempfile.TemporaryDirectory() as tmpdirname:
        # Try loading from non-existant_path
        with pytest.raises(FileNotFoundError):
            _ = TimeData.load(tmpdirname + "non_existant_file.csv")
