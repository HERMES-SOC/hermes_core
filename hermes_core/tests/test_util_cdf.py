"""Tests for cdf.py"""

from pathlib import Path
import pytest
import numpy as np
from numpy.random import random

from spacepy.pycdf import CDF, Var, gAttrList, zAttrList

from astropy.time import Time

import hermes_core
from hermes_core.util.cdf import CDFWriter

import json


def test_cdf_writer_default_attrs():
    # Initialize a CDF File Wrapper
    test_writer = CDFWriter()

    # Required Attributes
    default_attrs = test_writer._load_default_global_attr_schema()
    required_attrs = dict(filter(lambda item: item[-1]["required"], default_attrs.items()))

    # Test Number of Global Attrs in the Target Dict (Intermediate Data)
    assert len(test_writer.data.meta.keys()) == len(required_attrs.keys())

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    with pytest.raises(ValueError) as e:
        test_writer.to_cdf(output_path=test_cache)


def test_cdf_writer_valid_attrs():
    # fmt: off
    input_attrs = [
        ("Descriptor", "EEA>Electron Electrostatic Analyzer")
    ]
    # fmt: on

    # Initialize a CDF File Wrapper
    test_writer = CDFWriter()

    # Required Attributes
    default_attrs = test_writer._load_default_global_attr_schema()
    required_attrs = dict(filter(lambda item: item[-1]["required"], default_attrs.items()))

    # Add Custom Data to the Wrapper
    test_writer.add_attributes_from_list(attributes=input_attrs)

    # Test Number of Global Attrs in the Target Dict (Intermediate Data)
    assert len(test_writer.data.meta.keys()) == len(required_attrs.keys())

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_writer.to_cdf(output_path=test_cache)

    # # Test number of Global Attrs in the generated CDF File (Result Data)
    # assert len(test_writer.cdf.attrs) == len(required_attrs.keys())

    # Save the CDF to a File
    test_writer.save_cdf()

    test_file_cache_path = Path(test_file_output_path)
    # Test the File Exists
    assert test_file_cache_path.exists()

    # Remove the Test File from Cache
    test_file_cache_path.unlink()

    # Test the File was Deleted
    assert not test_file_cache_path.exists()


def test_cdf_writer_overide_derived_attr():
    # fmt: off
    input_attrs = [
        ("Descriptor", "EEA>Electron Electrostatic Analyzer"),
        ("Data_type", "test>data_type")
    ]
    # fmt: on

    # Initialize a CDF File Wrapper
    test_writer = CDFWriter()

    # Required Attributes
    default_attrs = test_writer._load_default_global_attr_schema()
    required_attrs = dict(filter(lambda item: item[-1]["required"], default_attrs.items()))

    # Add Custom Data to the Wrapper
    test_writer.add_attributes_from_list(attributes=input_attrs)

    # Test Number of Global Attrs in the Target Dict (Intermediate Data)
    assert len(test_writer.data.meta.keys()) == len(required_attrs.keys())

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_writer.to_cdf(output_path=test_cache)

    # Test number of Global Attrs in the generated CDF File (Result Data)
    # assert len(test_writer.cdf.attrs) == len(required_attrs.keys())

    # Test the Value was not Derrived and used the Overriden Value
    assert str(test_writer.cdf.attrs["Data_type"]) == input_attrs[-1][-1]

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
    assert (test_writer["time"].data.shape[0]) == len(time)


def test_cdf_writer_single_variable():
    # fmt: off
    input_attrs = [
        ("Descriptor", "EEA>Electron Electrostatic Analyzer"),
        ("Data_type", "test>data_type")
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
        "Data_type": "test>data_type",
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
        "Data_type": "test>data_type",
    }
    # fmt: on

    # Initialize a CDF File Wrapper
    test_writer = CDFWriter()
    # jsonify(test_writer.variable_attr_schema, "vattr_schema")

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
    # if len(result) > 0:
    #     jsonify(result, "validation_result")

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
        "Data_type": "test>data_type",
    }
    # fmt: on

    # Initialize a CDF File Wrapper
    test_writer = CDFWriter()
    # jsonify(test_writer.variable_attr_schema, "vattr_schema")

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
    # if len(result) > 0:
    #     jsonify(result, "validation_result")

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
        "Data_type": "test>data_type",
    }
    # fmt: on

    # Initialize a CDF File Wrapper
    test_writer = CDFWriter()
    # jsonify(test_writer.variable_attr_schema, "vattr_schema")

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
    if len(result) > 0:
        jsonify(result, "validation_result")

    assert result
    assert "Variable: Epoch missing 'VAR_TYPE' attribute. Cannot Validate Variable." not in result

    # Save the CDF to a File
    test_writer.save_cdf()
    # Remove the File
    test_file_cache_path = Path(test_file_output_path)
    test_file_cache_path.unlink()


def jsonify(obj, name):
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    result_output_path = str(Path(test_cache) / f"{name}.json")
    with open(result_output_path, "w") as f:
        json.dump(obj, f)
