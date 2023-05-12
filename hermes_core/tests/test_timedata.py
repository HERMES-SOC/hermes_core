"""Tests for CDF Files to and from data containers"""

from typing import OrderedDict
from pathlib import Path
import datetime
import pytest
import numpy as np
from numpy.random import random

from astropy.timeseries import TimeSeries
from astropy.table import Column
from astropy.time import Time
from astropy.units import Quantity

import hermes_core
from hermes_core.timedata import TimeData
from hermes_core.util.validation import validate


def get_bad_timeseries():
    ts = TimeSeries()

    # Create an astropy.Time object
    time = np.arange(10)
    time_col = Time(time, format="unix").to_datetime()
    col = Column(data=time_col, name="time", meta={})
    ts.add_column(col)

    # Add Measurement
    col = Column(data=random(size=(10)), name="measurement", meta={"CATDESC": "Test Measurement"})
    ts.add_column(col)
    return ts


def get_test_timeseries():
    ts = TimeSeries()

    # Create an astropy.Time object
    time = np.arange(10)
    time_col = Time(time, format="unix")
    ts["time"] = time_col

    # Add Measurement
    quant = Quantity(value=random(size=(10)), unit="m")
    ts["measurement"] = quant
    ts["measurement"].meta = OrderedDict(
        {
            "VAR_TYPE": "metadata",
            "CATDESC": "Test Metadata",
            "DISPLAY_TYPE": "time_series",
            "LABLAXIS": "Label Axis",
        }
    )
    return ts


def test_time_data_empty_ts():
    with pytest.raises(ValueError):
        _ = TimeData(data=TimeSeries())


def test_time_data_bad_ts():
    ts = get_bad_timeseries()
    with pytest.raises(TypeError):
        _ = TimeData(data=ts)


def get_test_global_meta():
    global_attrs_template = TimeData.global_attribute_template()


def test_time_data_default():
    ts = get_test_timeseries()

    with pytest.raises(ValueError) as e:
        # We expect this to throw an error that the Instrument is not recognized.
        # The Instrument is one of the attributes required for generating the filename
        # Initialize a CDF File Wrapper
        test_data = TimeData(ts)

    # Test Deleting the Writer
    del ts


def test_time_data_valid_attrs():
    # fmt: off
    input_attrs = {
        "Descriptor": "EEA>Electron Electrostatic Analyzer",
        "Data_level": "l1>Level 1",
        "Data_version": "v0.0.1",
        "Start_time": datetime.datetime.now()
    }
    # fmt: on

    ts = get_test_timeseries()
    # Initialize a CDF File Wrapper
    test_data = TimeData(ts, meta=input_attrs)

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_data.save(output_path=test_cache)

    test_file_cache_path = Path(test_file_output_path)
    # Test the File Exists
    assert test_file_cache_path.exists()

    # Remove the Test File from Cache
    test_file_cache_path.unlink()

    # Test the File was Deleted
    assert not test_file_cache_path.exists()


def test_time_data_single_measurement():
    # fmt: off
    input_attrs = {
        "Descriptor": "EEA>Electron Electrostatic Analyzer",
        "Data_level": "l1>Level 1",
        "Data_version": "v0.0.1",
    }
    # fmt: on

    ts = get_test_timeseries()
    # Initialize a CDF File Wrapper
    test_data = TimeData(ts, meta=input_attrs)

    # Add Measurement
    test_data["test_var1"] = Quantity(value=random(size=(10)), unit="km")
    test_data["test_var1"].meta.update({"test_attr1": "test_value1"})

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_data.save(output_path=test_cache)

    test_file_cache_path = Path(test_file_output_path)
    # Test the File Exists
    assert test_file_cache_path.exists()

    # Remove the Test File from Cache
    test_file_cache_path.unlink()

    # Test the File was Deleted
    assert not test_file_cache_path.exists()


def test_time_data_generate_valid_cdf():
    # fmt: off
    input_attrs = {
        "DOI": "https://doi.org/<PREFIX>/<SUFFIX>",
        "Data_level": "L1>Level 1",  # NOT AN ISTP ATTR
        "Data_version": "0.0.1",
        "Descriptor": "EEA>Electron Electrostatic Analyzer",
        "Data_product_descriptor": "odpd",
        "HTTP_LINK": [
            "https://spdf.gsfc.nasa.gov/istp_guide/istp_guide.html",
            "https://spdf.gsfc.nasa.gov/istp_guide/gattributes.html",
            "https://spdf.gsfc.nasa.gov/istp_guide/vattributes.html"
        ],
        "Instrument_mode": "default",  # NOT AN ISTP ATTR
        "Instrument_type": "Electric Fields (space)",
        "LINK_TEXT": [
            "ISTP Guide",
            "Global Attrs",
            "Variable Attrs"
        ],
        "LINK_TITLE": [
            "ISTP Guide",
            "Global Attrs",
            "Variable Attrs"
        ],
        "MODS": [
            "v0.0.0 - Original version.",
            "v1.0.0 - Include trajectory vectors and optics state.",
            "v1.1.0 - Update metadata: counts -> flux.",
            "v1.2.0 - Added flux error.",
            "v1.3.0 - Trajectory vector errors are now deltas."
        ],
        "PI_affiliation": "HERMES",
        "PI_name": "HERMES SOC",
        "TEXT": "Valid Test Case",
    }
    # fmt: on

    ts = get_test_timeseries()
    # Initialize a CDF File Wrapper
    test_data = TimeData(ts, meta=input_attrs)

    # Add the Time column
    test_data["time"].meta.update(
        {
            "CATDESC": "TT2000 time tags",
            "VAR_TYPE": "support_data",
        }
    )

    # Add 'data' VAR_TYPE Attributes
    num_random_vars = 2
    for i in range(num_random_vars):
        # Add Measurement
        test_data.add_measurement(
            measure_name=f"test_var{i}",
            measure_data=Quantity(value=random(size=(10)), unit="km"),
            measure_meta={
                "VAR_TYPE": "data",
                "CATDESC": "Test Data",
                "DISPLAY_TYPE": "time_series",
                "LABLAXIS": "Label Axis",
            },
        )

    # Add 'support_data' VAR_TYPE Attributes
    num_random_vars = 2
    for i in range(num_random_vars):
        # Add Measurement
        test_data.add_measurement(
            measure_name=f"test_support{i}",
            measure_data=Quantity(value=random(size=(10)), unit="km"),
            measure_meta={
                "VAR_TYPE": "support_data",
                "CATDESC": "Test Support",
                "DISPLAY_TYPE": "time_series",
                "LABLAXIS": "Label Axis",
            },
        )

    # Add 'metadata' VAR_TYPE Attributes
    num_random_vars = 2
    for i in range(num_random_vars):
        # Add Measurement
        test_data.add_measurement(
            measure_name=f"test_metadata{i}",
            measure_data=Quantity(value=random(size=(10)), unit="km"),
            measure_meta={
                "VAR_TYPE": "metadata",
                "CATDESC": "Test Metadata",
                "DISPLAY_TYPE": "time_series",
                "LABLAXIS": "Label Axis",
            },
        )

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_data.save(output_path=test_cache)

    # Validate the generated CDF File
    result = validate(filepath=test_file_output_path)
    assert len(result) <= 2  # TODO Logical Source and File ID Do not Agree

    # Remove the File
    test_file_cache_path = Path(test_file_output_path)
    test_file_cache_path.unlink()


def test_time_data_from_cdf():
    # fmt: off
    input_attrs = {
        "DOI": "https://doi.org/<PREFIX>/<SUFFIX>",
        "Data_level": "L1>Level 1",  # NOT AN ISTP ATTR
        "Data_version": "0.0.1",
        "Descriptor": "EEA>Electron Electrostatic Analyzer",
        "Data_product_descriptor": "odpd",
        "HTTP_LINK": [
            "https://spdf.gsfc.nasa.gov/istp_guide/istp_guide.html",
            "https://spdf.gsfc.nasa.gov/istp_guide/gattributes.html",
            "https://spdf.gsfc.nasa.gov/istp_guide/vattributes.html"
        ],
        "Instrument_mode": "default",  # NOT AN ISTP ATTR
        "Instrument_type": "Electric Fields (space)",
        "LINK_TEXT": [
            "ISTP Guide",
            "Global Attrs",
            "Variable Attrs"
        ],
        "LINK_TITLE": [
            "ISTP Guide",
            "Global Attrs",
            "Variable Attrs"
        ],
        "MODS": [
            "v0.0.0 - Original version.",
            "v1.0.0 - Include trajectory vectors and optics state.",
            "v1.1.0 - Update metadata: counts -> flux.",
            "v1.2.0 - Added flux error.",
            "v1.3.0 - Trajectory vector errors are now deltas."
        ],
        "PI_affiliation": "HERMES",
        "PI_name": "HERMES SOC",
        "TEXT": "Valid Test Case",
    }
    # fmt: on

    ts = get_test_timeseries()
    # Initialize a CDF File Wrapper
    test_data = TimeData(ts, meta=input_attrs)

    # Add the Time column
    test_data["time"].meta.update(
        {
            "CATDESC": "TT2000 time tags",
            "VAR_TYPE": "support_data",
        }
    )

    # Add 'data' VAR_TYPE Attributes
    num_random_vars = 2
    for i in range(num_random_vars):
        # Add Measurement
        test_data.add_measurement(
            measure_name=f"test_var{i}",
            measure_data=Quantity(value=random(size=(10)), unit="km"),
            measure_meta={
                "VAR_TYPE": "data",
                "CATDESC": "Test Data",
                "DISPLAY_TYPE": "time_series",
                "LABLAXIS": "Label Axis",
            },
        )

    # Add 'support_data' VAR_TYPE Attributes
    num_random_vars = 2
    for i in range(num_random_vars):
        # Add Measurement
        test_data.add_measurement(
            measure_name=f"test_support{i}",
            measure_data=Quantity(value=random(size=(10)), unit="km"),
            measure_meta={
                "VAR_TYPE": "support_data",
                "CATDESC": "Test Support",
                "DISPLAY_TYPE": "time_series",
                "LABLAXIS": "Label Axis",
            },
        )

    # Add 'metadata' VAR_TYPE Attributes
    num_random_vars = 2
    for i in range(num_random_vars):
        # Add Measurement
        test_data.add_measurement(
            measure_name=f"test_metadata{i}",
            measure_data=Quantity(value=random(size=(10)), unit="km"),
            measure_meta={
                "VAR_TYPE": "metadata",
                "CATDESC": "Test Metadata",
                "DISPLAY_TYPE": "time_series",
                "LABLAXIS": "Label Axis",
            },
        )

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_data.save(output_path=test_cache)

    # Validate the generated CDF File
    result = validate(test_file_output_path)
    assert len(result) <= 2  # TODO Logical Source and File ID Do not Agree

    # Try to Load the CDF File in a new CDFWriter
    new_writer = TimeData.load(test_file_output_path)

    # Remove the Original File
    test_file_cache_path = Path(test_file_output_path)
    test_file_cache_path.unlink()

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path2 = new_writer.save(output_path=test_cache)
    assert test_file_output_path == test_file_output_path2

    # Validate the generated CDF File
    result2 = validate(test_file_output_path2)
    assert len(result2) <= 2  # TODO Logical Source and File ID Do not Agree
    assert len(result) == len(result2)

    # Remove the Second File
    test_file_cache_path2 = Path(test_file_output_path2)
    test_file_cache_path2.unlink()
    assert (not test_file_cache_path.exists()) and (not test_file_cache_path2.exists())
