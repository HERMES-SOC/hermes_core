"""Tests for cdf.py"""

import pytest
from pathlib import Path
import hermes_core
from hermes_core.util.cdf import CDFWriter

DEFAULT_NUM_GLOBAL_ATTRS = 20


def test_cdf_writer_default_attrs():

    # Initialize a CDF File Wrapper
    test_writer = CDFWriter()

    # Test Number of Global Attrs in the Target Dict (Intermediate Data)
    assert len(test_writer.global_attrs.keys()) == DEFAULT_NUM_GLOBAL_ATTRS

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

    # Add Custom Data to the Wrapper
    test_writer.add_attributes_from_list(attributes=input_attrs)

    # Test Number of Global Attrs in the Target Dict (Intermediate Data)
    assert len(test_writer.global_attrs.keys()) == DEFAULT_NUM_GLOBAL_ATTRS

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_writer.to_cdf(output_path=test_cache)

    # Test number of Global Attrs in the generated CDF File (Result Data)
    assert len(test_writer.target_cdf_file.attrs) == DEFAULT_NUM_GLOBAL_ATTRS

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

    # Add Custom Data to the Wrapper
    test_writer.add_attributes_from_list(attributes=input_attrs)

    # Test Number of Global Attrs in the Target Dict (Intermediate Data)
    assert len(test_writer.global_attrs.keys()) == DEFAULT_NUM_GLOBAL_ATTRS

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_writer.to_cdf(output_path=test_cache)

    # Test number of Global Attrs in the generated CDF File (Result Data)
    assert len(test_writer.target_cdf_file.attrs) == DEFAULT_NUM_GLOBAL_ATTRS

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
