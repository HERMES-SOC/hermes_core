"""Tests for cdf.py"""

from pathlib import Path
import pytest
import numpy as np
from numpy.random import random
from spacepy.pycdf import CDF, Var, gAttrList, zAttrList
import hermes_core
from hermes_core.util.cdf import CDFWriter


def test_cdf_writer_default_attrs():

    # Initialize a CDF File Wrapper
    test_writer = CDFWriter()

    # Required Attributes
    default_attrs = test_writer._load_default_global_attr_schema()
    required_attrs = dict(filter(lambda item: item[-1]["required"], default_attrs.items()))

    # Test Number of Global Attrs in the Target Dict (Intermediate Data)
    assert len(test_writer.global_attrs.keys()) == len(required_attrs.keys())

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
    assert len(test_writer.global_attrs.keys()) == len(required_attrs.keys())

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_writer.to_cdf(output_path=test_cache)

    # Test number of Global Attrs in the generated CDF File (Result Data)
    assert len(test_writer.target_cdf_file.attrs) == len(required_attrs.keys())

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
    assert len(test_writer.global_attrs.keys()) == len(required_attrs.keys())

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_writer.to_cdf(output_path=test_cache)

    # Test number of Global Attrs in the generated CDF File (Result Data)
    assert len(test_writer.target_cdf_file.attrs) == len(required_attrs.keys())

    # Test the Value was not Derrived and used the Overriden Value
    assert str(test_writer.target_cdf_file.attrs["Data_type"]) == input_attrs[-1][-1]

    # Save the CDF to a File
    test_writer.save_cdf()

    test_file_cache_path = Path(test_file_output_path)
    # Test the File Exists
    assert test_file_cache_path.exists()

    # Remove the Test File from Cache
    test_file_cache_path.unlink()

    # Test the File was Deleted
    assert not test_file_cache_path.exists()


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
    assert isinstance(test_writer.target_cdf_file, CDF)
    assert isinstance(test_writer.target_cdf_file["test_var1"], Var)
    assert isinstance(test_writer.target_cdf_file["test_var1"].attrs, zAttrList)
    assert isinstance(test_writer.target_cdf_file.attrs, gAttrList)

    # Test the Target CDF Containts the Variable Data
    assert test_writer.target_cdf_file["test_var1"][0] == "test_data1"

    # Test the Target CDF Contains the Variable Attributes
    assert test_writer.target_cdf_file["test_var1"].attrs["test_attr1"] == "test_value1"

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
    input_attrs = [
        ("Descriptor", "EEA>Electron Electrostatic Analyzer"),
        ("Data_type", "test>data_type"),
        ("Logical_file_id", "invalid_cdf.cdf")
    ]
    # fmt: on

    # Initialize a CDF File Wrapper
    test_writer = CDFWriter()

    # Add Custom Data to the Wrapper
    test_writer.add_attributes_from_list(attributes=input_attrs)

    num_random_vars = 10
    N = 200
    M = 20
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
    assert isinstance(test_writer.target_cdf_file, CDF)
    assert isinstance(test_writer.target_cdf_file["test_var1"], Var)
    assert isinstance(test_writer.target_cdf_file["test_var1"].attrs, zAttrList)
    assert isinstance(test_writer.target_cdf_file.attrs, gAttrList)

    # Test the Target CDF Containts the Variable Data
    assert test_writer.target_cdf_file["test_var1"].shape == (N, M)

    # Test the Target CDF Contains the Variable Attributes
    assert test_writer.target_cdf_file["test_var1"].attrs["test_attr1"] == "test_value1"

    # Save the CDF to a File
    test_writer.save_cdf()

    test_file_cache_path = Path(test_file_output_path)
    # Test the File Exists
    assert test_file_cache_path.exists()

    # Remove the Test File from Cache
    test_file_cache_path.unlink()

    # Test the File was Deleted
    assert not test_file_cache_path.exists()
