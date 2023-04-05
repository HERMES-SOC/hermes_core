"""Tests for cdf.py"""

from pathlib import Path
import datetime
import pytest
import numpy as np
from numpy.random import random

import spacepy
from spacepy.pycdf import CDF, Var, gAttrList, zAttrList

from astropy.time import Time

import hermes_core
from hermes_core.util.cdf import CDFWriter

import json


def test_cdf_writer_default_attrs():
    # Initialize a CDF File Wrapper
    test_writer = CDFWriter()

    assert str(test_writer) is not None
    assert test_writer.__repr__() is not None

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    with pytest.raises(AttributeError) as e:
        test_writer.to_cdf(output_path=test_cache)

    # Test Deleting the Writer
    del test_writer


def test_cdf_writer_valid_attrs():
    # fmt: off
    input_attrs = [
        ("Descriptor", "EEA>Electron Electrostatic Analyzer"),
        ("Data_level", "l1>Level 1"),
        ("Data_version", "v0.0.1"),
        ("Start_time", datetime.datetime.now())
    ]
    # fmt: on

    # Initialize a CDF File Wrapper
    test_writer = CDFWriter()

    # Add Custom Data to the Wrapper
    test_writer.add_attributes_from_list(attributes=input_attrs)

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_writer.to_cdf(output_path=test_cache)

    # Save the CDF to a File
    test_writer.save_cdf()

    test_file_cache_path = Path(test_file_output_path)
    # Test the File Exists
    assert test_file_cache_path.exists()

    # Remove the Test File from Cache
    test_file_cache_path.unlink()

    # Test the File was Deleted
    assert not test_file_cache_path.exists()


def test_cdf_writer_invalid_derivation():
    # fmt: off
    input_attrs = [
        ("Descriptor", "EEA>Electron Electrostatic Analyzer"),
        ("Data_level", "l1>Level 1"),
        ("Data_version", "v0.0.1"),
        ("Start_time", datetime.datetime.now())
    ]
    # fmt: on

    # Initialize a CDF File Wrapper
    test_writer = CDFWriter()

    # Add Custom Data to the Wrapper
    test_writer.add_attributes_from_list(attributes=input_attrs)

    # Update the Global Attrs to be Invalid
    test_writer.global_attr_schema["New_Attr"] = {
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
        test_writer.to_cdf(output_path=test_cache)


def test_cdf_writer_overide_derived_attr():
    # fmt: off
    input_attrs = [
        ("Descriptor", "EEA>Electron Electrostatic Analyzer"),
        ("Data_level", "l1>Level 1"),
        ("Data_version", "v0.0.1"),
        ("Source_name", None),
        ("Data_product_descriptor", "odpd"),
        ("Data_type", "test>data_type"),
        ("Start_time", datetime.datetime.now())
    ]
    # fmt: on

    # Initialize a CDF File Wrapper
    test_writer = CDFWriter()

    # Add Custom Data to the Wrapper
    test_writer.add_attributes_from_list(attributes=input_attrs)

    # Derrive Attributes
    test_writer.derive_attributes()

    # Re-Overide Data
    input_attrs = {
        "Data_type": "test_override>data_type_override",
    }

    # Add Custom Data to the Wrapper
    test_writer.add_attributes_from_dict(attributes=input_attrs)

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_writer.to_cdf(output_path=test_cache)

    # Test the Value was not Derrived and used the Overriden Value
    assert str(test_writer.cdf.attrs["Data_type"]) == input_attrs["Data_type"]

    # Save the CDF to a File
    test_writer.save_cdf()

    test_file_cache_path = Path(test_file_output_path)
    # Test the File Exists
    assert test_file_cache_path.exists()

    # Remove the Test File from Cache
    test_file_cache_path.unlink()

    # Test the File was Deleted
    assert not test_file_cache_path.exists()


def test_cdf_writer_bad_variable():
    # Initialize a CDF File Wrapper
    test_writer = CDFWriter()

    with pytest.raises(ValueError) as e:
        # Add Variable
        test_writer.add_variable(
            var_name="test_var1",
            var_data=np.array(["test_data1"]),
            var_attrs={"test_attr1": "test_value1"},
        )


def test_cdf_writer_add_time():
    # Initialize a CDF File Wrapper
    test_writer = CDFWriter()

    # Create an astropy.Time object
    time = np.arange(50)
    time_col = Time(time, format="unix")

    # Add the Time column
    test_writer.add_time(time=time_col, time_attrs={})

    # Assert the Time Dimension in the TimeSeries Data matches the input
    assert len(test_writer) == 1
    assert "time" in test_writer
    assert (test_writer["time"].data.shape[0]) == len(time)


def test_cdf_writer_single_variable():
    # fmt: off
    input_attrs = [
        ("Descriptor", "EEA>Electron Electrostatic Analyzer"),
        ("Data_level", "l1>Level 1"),
        ("Data_version", "v0.0.1"),
    ]
    # fmt: on

    # Initialize a CDF File Wrapper
    test_writer = CDFWriter()

    # Add Custom Data to the Wrapper
    test_writer.add_attributes_from_list(attributes=input_attrs)

    # Create an astropy.Time object
    time = np.arange(50)
    time_col = Time(time, format="unix")

    # Add the Time column
    test_writer.add_time(time=time_col, time_attrs={})

    # Add Variable
    test_writer.add_variable(
        var_name="test_var1",
        var_data=np.array(["test_data1"]),
        var_attrs={"test_attr1": "test_value1"},
    )

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_writer.to_cdf(output_path=test_cache)

    # Assert Types of the Target CDF
    assert isinstance(test_writer.cdf, CDF)
    assert isinstance(test_writer.cdf["test_var1"], Var)
    assert isinstance(test_writer.cdf["test_var1"].attrs, zAttrList)
    assert isinstance(test_writer.cdf.attrs, gAttrList)

    # Test the Target CDF Containts the Variable Data
    assert test_writer.cdf["test_var1"][0] == "test_data1"

    # Test the Target CDF Contains the Variable Attributes
    assert test_writer.cdf["test_var1"].attrs["test_attr1"] == "test_value1"

    # Save the CDF to a File
    test_writer.save_cdf()

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

    # Initialize a CDF File Wrapper
    test_writer = CDFWriter()

    # Add Custom Data to the Wrapper
    test_writer.add_attributes_from_dict(attributes=input_attrs)

    # Add Variable Data
    N = 10  # Num Timesteps
    M = 2  # Num Columns

    # Create an astropy.Time object
    time = np.arange(N)
    time_col = Time(time, format="unix")

    # Add the Time column
    test_writer.add_time(time=time_col, time_attrs={})

    num_random_vars = 10
    for i in range(num_random_vars):
        data = random(size=(N, M))
        # Add Variable
        test_writer.add_variable(
            var_name=f"test_var{i}",
            var_data=data,
            var_attrs={f"test_attr{i}": f"test_value{i}"},
        )

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_writer.to_cdf(output_path=test_cache)

    # Assert Types of the Target CDF
    assert isinstance(test_writer.cdf, CDF)
    assert isinstance(test_writer.cdf["test_var1"], Var)
    assert isinstance(test_writer.cdf["test_var1"].attrs, zAttrList)
    assert isinstance(test_writer.cdf.attrs, gAttrList)

    # Test the Target CDF Containts the Variable Data
    assert test_writer.cdf["test_var1"].shape == (N, M)

    # Test the Target CDF Contains the Variable Attributes
    assert test_writer.cdf["test_var1"].attrs["test_attr1"] == "test_value1"

    # Save the CDF to a File
    test_writer.save_cdf()

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

    # Initialize a CDF File Wrapper
    test_writer = CDFWriter()

    # Add Custom Data to the Wrapper
    test_writer.add_attributes_from_dict(attributes=input_attrs)

    # Add Variable Data
    N = 10  # Num Timesteps

    # Create an astropy.Time object
    time = np.arange(N)
    time_col = Time(time, format="unix")

    # Add the Time column
    test_writer.add_time(time=time_col, time_attrs={})

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_writer.to_cdf(output_path=test_cache)

    # Validate the generated CDF File
    result = test_writer.validate_cdf(catch=True)

    assert "Variable: Epoch missing 'VAR_TYPE' attribute. Cannot Validate Variable." in result

    # Save the CDF to a File
    test_writer.save_cdf()
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

    # Initialize a CDF File Wrapper
    test_writer = CDFWriter()

    # Add Custom Data to the Wrapper
    test_writer.add_attributes_from_dict(attributes=input_attrs)

    # Add Variable Data
    N = 10  # Num Timesteps

    # Create an astropy.Time object
    time = np.arange(N)
    time_col = Time(time, format="unix")

    # Add the Time column
    test_writer.add_time(time=time_col, time_attrs={"VAR_TYPE": "time_series"})

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_writer.to_cdf(output_path=test_cache)

    # Validate the generated CDF File
    result = test_writer.validate_cdf(catch=True)

    assert result
    assert "Variable: Epoch missing 'VAR_TYPE' attribute. Cannot Validate Variable." not in result

    # Save the CDF to a File
    test_writer.save_cdf()
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

    # Initialize a CDF File Wrapper
    test_writer = CDFWriter()

    # Add Custom Data to the Wrapper
    test_writer.add_attributes_from_dict(attributes=input_attrs)

    # Add Variable Data
    N = 10  # Num Timesteps

    # Create an astropy.Time object
    time = np.arange(N)
    time_col = Time(time, format="unix")

    # Add the Time column
    test_writer.add_time(time=time_col, time_attrs={"VAR_TYPE": "time_series"})

    # Add 'data' VAR_TYPE Attributes
    num_random_vars = 2
    for i in range(num_random_vars):
        data = random(size=(N,))
        # Add Variable
        test_writer.add_variable(
            var_name=f"test_var{i}",
            var_data=data,
            var_attrs={f"VAR_TYPE": "data"},
        )

    # Add 'support_data' VAR_TYPE Attributes
    num_random_vars = 2
    for i in range(num_random_vars):
        data = random(size=(N,))
        # Add Variable
        test_writer.add_variable(
            var_name=f"test_support{i}",
            var_data=data,
            var_attrs={f"VAR_TYPE": "support_data"},
        )

    # Add 'metadata' VAR_TYPE Attributes
    num_random_vars = 2
    for i in range(num_random_vars):
        data = random(size=(N,))
        # Add Variable
        test_writer.add_variable(
            var_name=f"test_metadata{i}",
            var_data=data,
            var_attrs={f"VAR_TYPE": "metadata"},
        )

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_writer.to_cdf(output_path=test_cache)

    # Validate the generated CDF File
    result = test_writer.validate_cdf(catch=True)

    assert result
    assert "Variable: Epoch missing 'VAR_TYPE' attribute. Cannot Validate Variable." not in result

    # Save the CDF to a File
    test_writer.save_cdf()
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

    # Initialize a CDF File Wrapper
    test_writer = CDFWriter()

    # Required Attributes
    default_attrs = test_writer._load_default_global_attr_schema()
    required_attrs = dict(filter(lambda item: item[-1]["required"], default_attrs.items()))

    # Add Custom Data to the Wrapper
    test_writer.add_attributes_from_dict(attributes=input_attrs)

    # Add Variable Data
    N = 10  # Num Timesteps

    # Create an astropy.Time object
    time = np.arange(N)
    time_col = Time(time, format="unix")

    # Add the Time column
    test_writer["time"] = time_col
    test_writer["time"].meta = {
        "CATDESC": "TT2000 time tags",
        "FIELDNAM": "Epoch",
        # "FILLVAL": -9223372036854775808,
        "FILLVAL": spacepy.pycdf.lib.v_tt2000_to_datetime(-9223372036854775808),
        "VAR_TYPE": "time_series",
        "TIME_BASE": "J2000",
        "RESOLUTION": "1s",
        "TIME_SCALE": "Terrestrial Time (TT)",
        "REFERENCE_POSITION": "rotating Earth geoid",
    }

    # Add 'data' VAR_TYPE Attributes
    num_random_vars = 2
    for i in range(num_random_vars):
        data = random(size=(N,))
        # Add Variable
        test_writer.add_variable(
            var_name=f"test_var{i}",
            var_data=data,
            var_attrs={
                "VAR_TYPE": "data",
                "CATDESC": "Test Variable",
                "DEPEND_0": "Epoch",
                "DISPLAY_TYPE": "time_series",
                "FIELDNAM": f"test_var{i}",
                "FILLVAL": -1e31,
                "FORMAT": "F9.4",
                "LABLAXIS": "Label Axis",
                "SI_CONVERSION": "1.0e3>m",
                "UNITS": "km",
                "VALIDMIN": 0.0,
                "VALIDMAX": 1.0,
            },
        )

    # Add 'support_data' VAR_TYPE Attributes
    num_random_vars = 2
    for i in range(num_random_vars):
        data = random(size=(N,))
        # Add Variable
        test_writer.add_variable(
            var_name=f"test_support{i}",
            var_data=data,
            var_attrs={
                "VAR_TYPE": "support_data",
                "CATDESC": "Test Variable",
                "DEPEND_0": "Epoch",
                "DISPLAY_TYPE": "time_series",
                "FIELDNAM": f"test_support{i}",
                "FILLVAL": -1e31,
                "FORMAT": "F9.4",
                "LABLAXIS": "Label Axis",
                "SI_CONVERSION": "1.0e3>m",
                "UNITS": "km",
                "VALIDMIN": 0.0,
                "VALIDMAX": 1.0,
            },
        )

    # Add 'metadata' VAR_TYPE Attributes
    data = random(size=(N,))
    # Add Variable
    test_writer["test_metadata"] = data
    test_writer["test_metadata"].meta = {
        "VAR_TYPE": "metadata",
        "CATDESC": "Test Variable",
        "DEPEND_0": "Epoch",
        "DISPLAY_TYPE": "time_series",
        "FIELDNAM": f"test_metadata",
        "FILLVAL": -1e31,
        "FORMAT": "F9.4",
        "LABLAXIS": "Label Axis",
        "SI_CONVERSION": "1.0e3>m",
        "UNITS": "km",
        "VALIDMIN": 0.0,
        "VALIDMAX": 1.0,
    }

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_writer.to_cdf(output_path=test_cache)

    # Test number of Global Attrs in the generated CDF File (Result Data)
    assert len(test_writer.cdf.attrs) >= len(required_attrs.keys())

    # Validate the generated CDF File
    result = test_writer.validate_cdf(catch=True)
    assert len(result) <= 1  # TODO Logical Source and File ID Do not Agree

    # Save the CDF to a File
    test_writer.save_cdf()
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

    # Initialize a CDF File Wrapper
    test_writer = CDFWriter()

    # Add Custom Data to the Wrapper
    test_writer.add_attributes_from_dict(attributes=input_attrs)

    # Add Variable Data
    N = 10  # Num Timesteps

    # Create an astropy.Time object
    time = np.arange(N)
    time_col = Time(time, format="unix")

    # Add the Time column
    test_writer["time"] = time_col
    test_writer["time"].meta = {
        "CATDESC": "TT2000 time tags",
        "FIELDNAM": "Epoch",
        # "FILLVAL": -9223372036854775808,
        "FILLVAL": spacepy.pycdf.lib.v_tt2000_to_datetime(-9223372036854775808),
        "VAR_TYPE": "time_series",
        "TIME_BASE": "J2000",
        "RESOLUTION": "1s",
        "TIME_SCALE": "Terrestrial Time (TT)",
        "REFERENCE_POSITION": "rotating Earth geoid",
    }

    # Add 'data' VAR_TYPE Attributes
    num_random_vars = 2
    for i in range(num_random_vars):
        data = random(size=(N,))
        # Add Variable
        test_writer.add_variable(
            var_name=f"test_var{i}",
            var_data=data,
            var_attrs={
                "VAR_TYPE": "data",
                "CATDESC": "Test Variable",
                "DEPEND_0": "Epoch",
                "DISPLAY_TYPE": "time_series",
                "FIELDNAM": f"test_var{i}",
                "FILLVAL": -1e31,
                "FORMAT": "F9.4",
                "LABLAXIS": "Label Axis",
                "SI_CONVERSION": "1.0e3>m",
                "UNITS": "km",
                "VALIDMIN": 0.0,
                "VALIDMAX": 1.0,
            },
        )

    # Add 'support_data' VAR_TYPE Attributes
    num_random_vars = 2
    for i in range(num_random_vars):
        data = random(size=(N,))
        # Add Variable
        test_writer.add_variable(
            var_name=f"test_support{i}",
            var_data=data,
            var_attrs={
                "VAR_TYPE": "support_data",
                "CATDESC": "Test Variable",
                "DEPEND_0": "Epoch",
                "DISPLAY_TYPE": "time_series",
                "FIELDNAM": f"test_support{i}",
                "FILLVAL": -1e31,
                "FORMAT": "F9.4",
                "LABLAXIS": "Label Axis",
                "SI_CONVERSION": "1.0e3>m",
                "UNITS": "km",
                "VALIDMIN": 0.0,
                "VALIDMAX": 1.0,
            },
        )

    # Add 'metadata' VAR_TYPE Attributes
    data = random(size=(N,))
    # Add Variable
    test_writer["test_metadata"] = data
    test_writer["test_metadata"].meta = {
        "VAR_TYPE": "metadata",
        "CATDESC": "Test Variable",
        "DEPEND_0": "Epoch",
        "DISPLAY_TYPE": "time_series",
        "FIELDNAM": f"test_metadata",
        "FILLVAL": -1e31,
        "FORMAT": "F9.4",
        "LABLAXIS": "Label Axis",
        "SI_CONVERSION": "1.0e3>m",
        "UNITS": "km",
        "VALIDMIN": 0.0,
        "VALIDMAX": 1.0,
    }

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_writer.to_cdf(output_path=test_cache)

    # Validate the generated CDF File
    result = test_writer.validate_cdf(catch=True)
    assert len(result) <= 1  # TODO Logical Source and File ID Do not Agree

    # Save the CDF to a File
    test_writer.save_cdf()

    # Try to Load the CDF File in a new CDFWriter
    new_writer = CDFWriter.from_cdf(test_file_output_path)

    # Remove the Original File
    test_file_cache_path = Path(test_file_output_path)
    test_file_cache_path.unlink()

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path2 = new_writer.to_cdf(output_path=test_cache)
    assert test_file_output_path == test_file_output_path2

    # Validate the generated CDF File
    result2 = new_writer.validate_cdf(catch=True)
    assert len(result2) <= 1  # TODO Logical Source and File ID Do not Agree
    assert len(result) == len(result2)

    # Remove the Second File
    test_file_cache_path2 = Path(test_file_output_path2)
    test_file_cache_path2.unlink()
    assert (not test_file_cache_path.exists()) and (not test_file_cache_path2.exists())
