"""Tests for cdf.py"""

from typing import OrderedDict
from pathlib import Path
import datetime
import pytest
import numpy as np
from numpy.random import random

import spacepy
from spacepy.pycdf import CDF, Var, gAttrList, zAttrList

from astropy.timeseries import TimeSeries
from astropy.table import Column
from astropy.time import Time
from astropy.units import Quantity

import hermes_core
from hermes_core.util.cdf import CDFWriter

import json


def get_bad_timeseries():
    ts = TimeSeries()

    # Create an astropy.Time object
    time = np.arange(10)
    time_col = Time(time, format="unix").to_datetime()
    col = Column(data=time_col, name="time", meta={})
    ts.add_column(col)

    # Add Variable
    col = Column(data=random(size=(10)), name="measurement", meta={"CATDESC": "Test Measurement"})
    ts.add_column(col)
    return ts


def get_test_timeseries():
    ts = TimeSeries()

    # Create an astropy.Time object
    time = np.arange(10)
    time_col = Time(time, format="unix")
    # col = Column(data=time_col, name="time", meta={})
    ts["time"] = time_col

    # Add Variable
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


def test_cdf_writer_empty_ts():
    with pytest.raises(ValueError):
        _ = CDFWriter(ts=TimeSeries())


def test_cdf_writer_bad_ts():
    ts = get_bad_timeseries()
    with pytest.raises(TypeError):
        _ = CDFWriter(ts=ts)


def test_cdf_writer_default_attrs():
    ts = get_test_timeseries()
    # Initialize a CDF File Wrapper
    test_writer = CDFWriter(ts)

    assert isinstance(test_writer.data, TimeSeries)
    assert isinstance(test_writer.meta, OrderedDict)
    assert len(test_writer.meta) > 0
    assert isinstance(test_writer.units, OrderedDict)
    assert len(test_writer.units) >= 0
    assert isinstance(test_writer.columns, OrderedDict)
    assert len(test_writer.columns) >= 2
    assert test_writer.time is not None
    assert test_writer.shape == (10, 2)
    assert str(test_writer) is not None
    assert test_writer.__repr__() is not None

    # Check "measurement"
    assert isinstance(test_writer["measurement"].meta, OrderedDict)
    assert "CATDESC" in test_writer["measurement"].meta

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    with pytest.raises(ValueError) as e:
        test_writer.write_cdf(output_path=test_cache)

    # Test Deleting the Writer
    del test_writer


def test_cdf_writer_valid_attrs():
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
    test_writer = CDFWriter(ts=ts)

    # Add Custom Data to the Wrapper
    test_writer.meta.update(input_attrs)

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_writer.write_cdf(output_path=test_cache)

    test_file_cache_path = Path(test_file_output_path)
    # Test the File Exists
    assert test_file_cache_path.exists()

    # Remove the Test File from Cache
    test_file_cache_path.unlink()

    # Test the File was Deleted
    assert not test_file_cache_path.exists()


def test_cdf_writer_invalid_derivation():
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
    test_writer = CDFWriter(ts=ts)

    # Add Custom Data to the Wrapper
    test_writer.meta.update(input_attrs)

    # Update the Global Attrs to be Invalid
    test_writer._global_attr_schema["New_Attr"] = {
        "derived": True,
        "required": True,
        "valid_check": None,
    }

    # Test Getting a Bad Variable
    with pytest.raises(KeyError) as e:
        test = test_writer["test"]

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    with pytest.raises(ValueError) as e:
        test_writer.write_cdf(output_path=test_cache)


def test_cdf_writer_overide_derived_attr():
    # fmt: off
    input_attrs = {
        "Descriptor": "EEA>Electron Electrostatic Analyzer",
        "Data_level": "l1>Level 1",
        "Data_version": "v0.0.1",
        "Source_name": None,
        "Data_product_descriptor": "odpd",
        "Data_type": "test>data_type",
        "Start_time": datetime.datetime.now()
    }
    # fmt: on

    ts = get_test_timeseries()
    # Initialize a CDF File Wrapper
    test_writer = CDFWriter(ts=ts)

    # Add Custom Data to the Wrapper
    test_writer.meta.update(input_attrs)

    # Derrive Attributes
    test_writer.derive_attributes()

    # Re-Overide Data
    input_attrs = {
        "Data_type": "test_override>data_type_override",
    }

    # Add Custom Data to the Wrapper
    test_writer.meta.update(input_attrs)

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_writer.write_cdf(output_path=test_cache)

    test_file_cache_path = Path(test_file_output_path)
    # Test the File Exists
    assert test_file_cache_path.exists()

    # Remove the Test File from Cache
    test_file_cache_path.unlink()

    # Test the File was Deleted
    assert not test_file_cache_path.exists()


def test_cdf_writer_single_variable():
    # fmt: off
    input_attrs = {
        "Descriptor": "EEA>Electron Electrostatic Analyzer",
        "Data_level": "l1>Level 1",
        "Data_version": "v0.0.1",
    }
    # fmt: on

    ts = get_test_timeseries()
    # Initialize a CDF File Wrapper
    test_writer = CDFWriter(ts=ts)

    # Add Custom Data to the Wrapper
    test_writer.meta.update(input_attrs)

    # Add Variable
    test_writer["test_var1"] = Quantity(value=random(size=(10)), unit="km")
    test_writer["test_var1"].meta.update({"test_attr1": "test_value1"})

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_writer.write_cdf(output_path=test_cache)

    test_file_cache_path = Path(test_file_output_path)
    # Test the File Exists
    assert test_file_cache_path.exists()

    # Remove the Test File from Cache
    test_file_cache_path.unlink()

    # Test the File was Deleted
    assert not test_file_cache_path.exists()


def test_cdf_writer_random_variable():
    # fmt: off
    input_attrs = {
        "Descriptor": "EEA>Electron Electrostatic Analyzer",
        "Data_level": "l1>Level 1",
        "Data_version": "v0.0.1",
    }
    # fmt: on

    ts = get_test_timeseries()
    # Initialize a CDF File Wrapper
    test_writer = CDFWriter(ts=ts)

    # Add Custom Data to the Wrapper
    test_writer.meta.update(input_attrs)

    # Add Variable Data
    N = 10  # Num Timesteps
    M = 2  # Num Columns

    num_random_vars = 10
    for i in range(num_random_vars):
        # Add Variable
        test_writer[f"test_var{i}"] = Quantity(value=random(size=(10)), unit="km")
        test_writer[f"test_var{i}"].meta.update({f"test_attr{i}": f"test_value{i}"})

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_writer.write_cdf(output_path=test_cache)

    test_file_cache_path = Path(test_file_output_path)
    # Test the File Exists
    assert test_file_cache_path.exists()

    # Remove the Test File from Cache
    test_file_cache_path.unlink()

    # Test the File was Deleted
    assert not test_file_cache_path.exists()


def test_cdf_writer_validate_missing_epoch_var_type():
    # fmt: off
    input_attrs = {
        "Descriptor": "EEA>Electron Electrostatic Analyzer",
        "Data_level": "l1>Level 1",
        "Data_version": "v0.0.1",
    }
    # fmt: on

    ts = get_test_timeseries()
    # Initialize a CDF File Wrapper
    test_writer = CDFWriter(ts=ts)

    # Add Custom Data to the Wrapper
    test_writer.meta.update(input_attrs)

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_writer.write_cdf(output_path=test_cache)

    # Validate the generated CDF File
    result = test_writer.validate_cdf(cdf_file_path=test_file_output_path, catch=True)

    assert "Variable: Epoch missing 'VAR_TYPE' attribute. Cannot Validate Variable." in result

    # Remove the File
    test_file_cache_path = Path(test_file_output_path)
    test_file_cache_path.unlink()


def test_cdf_writer_validate_present_epoch_var_type():
    # fmt: off
    input_attrs = {
        "Descriptor": "EEA>Electron Electrostatic Analyzer",
        "Data_level": "l1>Level 1",
        "Data_version": "v0.0.1",
    }
    # fmt: on

    ts = get_test_timeseries()
    # Initialize a CDF File Wrapper
    test_writer = CDFWriter(ts=ts)

    # Add Custom Data to the Wrapper
    test_writer.meta.update(input_attrs)

    # Add the Time column
    test_writer.time.meta.update({"VAR_TYPE": "time_series"})

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_writer.write_cdf(output_path=test_cache)

    # Validate the generated CDF File
    result = test_writer.validate_cdf(cdf_file_path=test_file_output_path, catch=True)

    assert result
    assert "Variable: Epoch missing 'VAR_TYPE' attribute. Cannot Validate Variable." not in result

    # Remove the File
    test_file_cache_path = Path(test_file_output_path)
    test_file_cache_path.unlink()


def test_cdf_writer_validate_multiple_var_type():
    # fmt: off
    input_attrs = {
        "Descriptor": "EEA>Electron Electrostatic Analyzer",
        "Data_level": "l1>Level 1",
        "Data_version": "v0.0.1",
    }
    # fmt: on

    ts = get_test_timeseries()
    # Initialize a CDF File Wrapper
    test_writer = CDFWriter(ts=ts)

    # Add Custom Data to the Wrapper
    test_writer.meta.update(input_attrs)

    # Add the Time column
    test_writer.time.meta.update({"VAR_TYPE": "time_series"})

    # Add Variable Data
    N = 10  # Num Timesteps

    # Add 'data' VAR_TYPE Attributes
    num_random_vars = 2
    for i in range(num_random_vars):
        # Add Variable
        test_writer[f"test_var{i}"] = Quantity(value=random(size=(10)), unit="km")
        test_writer[f"test_var{i}"].meta.update({f"VAR_TYPE": "data"})

    # Add 'support_data' VAR_TYPE Attributes
    num_random_vars = 2
    for i in range(num_random_vars):
        # Add Variable
        test_writer[f"test_support{i}"] = Quantity(value=random(size=(10)), unit="km")
        test_writer[f"test_support{i}"].meta.update({f"VAR_TYPE": "support_data"})

    # Add 'metadata' VAR_TYPE Attributes
    num_random_vars = 2
    for i in range(num_random_vars):
        # Add Variable
        test_writer[f"test_metadata{i}"] = Quantity(value=random(size=(10)), unit="km")
        test_writer[f"test_metadata{i}"].meta.update({f"VAR_TYPE": "metadata"})

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_writer.write_cdf(output_path=test_cache)

    # Validate the generated CDF File
    result = test_writer.validate_cdf(cdf_file_path=test_file_output_path, catch=True)

    assert result
    assert "Variable: Epoch missing 'VAR_TYPE' attribute. Cannot Validate Variable." not in result
    assert (
        "Variable: test_var0 missing 'VAR_TYPE' attribute. Cannot Validate Variable." not in result
    )

    # Remove the File
    test_file_cache_path = Path(test_file_output_path)
    test_file_cache_path.unlink()


def test_cdf_writer_generate_valid_cdf():
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
    test_writer = CDFWriter(ts=ts)

    # Add Custom Data to the Wrapper
    test_writer.meta.update(input_attrs)

    # Required Attributes
    default_attrs = test_writer._load_default_global_attr_schema()
    required_attrs = dict(filter(lambda item: item[-1]["required"], default_attrs.items()))

    # Add Variable Data
    N = 10  # Num Timesteps

    # Add the Time column
    test_writer["time"].meta.update(
        {
            "CATDESC": "TT2000 time tags",
            "VAR_TYPE": "support_data",
        }
    )

    # Add 'data' VAR_TYPE Attributes
    num_random_vars = 2
    for i in range(num_random_vars):
        # Add Variable
        test_writer[f"test_var{i}"] = Quantity(value=random(size=(10)), unit="km")
        test_writer[f"test_var{i}"].meta.update(
            {
                "VAR_TYPE": "data",
                "CATDESC": "Test Variable",
                "DISPLAY_TYPE": "time_series",
                "LABLAXIS": "Label Axis",
            }
        )

    # Add 'support_data' VAR_TYPE Attributes
    num_random_vars = 2
    for i in range(num_random_vars):
        # Add Variable
        test_writer[f"test_support{i}"] = Quantity(value=random(size=(10)), unit="km")
        test_writer[f"test_support{i}"].meta.update(
            {
                "VAR_TYPE": "support_data",
                "CATDESC": "Test Variable",
                "DISPLAY_TYPE": "time_series",
                "LABLAXIS": "Label Axis",
            },
        )

    # Add 'metadata' VAR_TYPE Attributes
    num_random_vars = 2
    for i in range(num_random_vars):
        # Add Variable
        test_writer[f"test_metadata{i}"] = Quantity(value=random(size=(10)), unit="km")
        test_writer[f"test_metadata{i}"].meta.update(
            {
                "VAR_TYPE": "metadata",
                "CATDESC": "Test Variable",
                "DISPLAY_TYPE": "time_series",
                "LABLAXIS": "Label Axis",
            }
        )

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_writer.write_cdf(output_path=test_cache)

    # Validate the generated CDF File
    result = test_writer.validate_cdf(cdf_file_path=test_file_output_path, catch=True)
    assert len(result) <= 2  # TODO Logical Source and File ID Do not Agree

    # Remove the File
    test_file_cache_path = Path(test_file_output_path)
    test_file_cache_path.unlink()


def test_cdf_writer_from_cdf():
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
    test_writer = CDFWriter(ts=ts)

    # Add Custom Data to the Wrapper
    test_writer.meta.update(input_attrs)

    # Add Variable Data
    N = 10  # Num Timesteps

    # Add the Time column
    test_writer["time"].meta.update(
        {
            "CATDESC": "TT2000 time tags",
            "VAR_TYPE": "support_data",
        }
    )

    # Add 'data' VAR_TYPE Attributes
    num_random_vars = 2
    for i in range(num_random_vars):
        # Add Variable
        test_writer[f"test_var{i}"] = Quantity(value=random(size=(10)), unit="km")
        test_writer[f"test_var{i}"].meta.update(
            {
                "VAR_TYPE": "data",
                "CATDESC": "Test Variable",
                "DISPLAY_TYPE": "time_series",
                "LABLAXIS": "Label Axis",
            }
        )

    # Add 'support_data' VAR_TYPE Attributes
    num_random_vars = 2
    for i in range(num_random_vars):
        # Add Variable
        test_writer[f"test_support{i}"] = Quantity(value=random(size=(10)), unit="km")
        test_writer[f"test_support{i}"].meta.update(
            {
                "VAR_TYPE": "support_data",
                "CATDESC": "Test Variable",
                "DISPLAY_TYPE": "time_series",
                "LABLAXIS": "Label Axis",
            },
        )

    # Add 'metadata' VAR_TYPE Attributes
    num_random_vars = 2
    for i in range(num_random_vars):
        # Add Variable
        test_writer[f"test_metadata{i}"] = Quantity(value=random(size=(10)), unit="km")
        test_writer[f"test_metadata{i}"].meta.update(
            {
                "VAR_TYPE": "metadata",
                "CATDESC": "Test Variable",
                "DISPLAY_TYPE": "time_series",
                "LABLAXIS": "Label Axis",
            }
        )

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_writer.write_cdf(output_path=test_cache)

    # Validate the generated CDF File
    result = test_writer.validate_cdf(cdf_file_path=test_file_output_path, catch=True)
    assert len(result) <= 2  # TODO Logical Source and File ID Do not Agree

    # Try to Load the CDF File in a new CDFWriter
    new_writer = CDFWriter.from_cdf(test_file_output_path)

    # Remove the Original File
    test_file_cache_path = Path(test_file_output_path)
    test_file_cache_path.unlink()

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path2 = new_writer.write_cdf(output_path=test_cache)
    assert test_file_output_path == test_file_output_path2

    # Validate the generated CDF File
    result2 = new_writer.validate_cdf(cdf_file_path=test_file_output_path2, catch=True)
    assert len(result2) <= 2  # TODO Logical Source and File ID Do not Agree
    assert len(result) == len(result2)

    # Remove the Second File
    test_file_cache_path2 = Path(test_file_output_path2)
    test_file_cache_path2.unlink()
    assert (not test_file_cache_path.exists()) and (not test_file_cache_path2.exists())
