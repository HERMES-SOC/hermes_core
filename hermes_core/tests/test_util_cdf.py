"""Tests for cdf.py"""

import pytest
from pathlib import Path
import hermes_core
from hermes_core.util.cdf import CDFWriter

DEFAULT_NUM_GLOBAL_ATTRS = 20


def test_cdf_writer_default_attrs():

    # fmt: off
    input_data = {
        "gAttrList": {},
        "zAttrList": {
            "variable": {
                "value": None,
                "required": True,
                "valid_check": None,
                "derived": False
            },
        },
    }
    # fmt: on

    # Initialize a CDF File Wrapper
    test_writer = CDFWriter()

    # Add Custom Data to the Wrapper
    test_writer.add_data_from_dict(data=input_data)

    # Test Number of Global Attrs in the Target Dict (Intermediate Data)
    assert len(test_writer.target_dict["gAttrList"].keys()) == DEFAULT_NUM_GLOBAL_ATTRS

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    with pytest.raises(ValueError) as e:
        test_writer.to_cdf(output_path=test_cache)


def test_cdf_writer_valid_attrs():

    # fmt: off
    input_data = {
        "gAttrList": {
            "Descriptor": {
                "value": "EEA>Electron Electrostatic Analyzer"
            }
        },
        "zAttrList": {
            "test_attr_1": {
                "value": 'test_attr_1',
                "required": False,
                "valid_check": None,
                "derived": False
            },
            "test_attr_1": {
                "value": None,
                "required": True,
                "valid_check": None,
                "derived": False
            },
        },
    }
    # fmt: on

    # Initialize a CDF File Wrapper
    test_writer = CDFWriter()

    # Add Custom Data to the Wrapper
    test_writer.add_data_from_dict(data=input_data)

    # Test Number of Global Attrs in the Target Dict (Intermediate Data)
    assert len(test_writer.target_dict["gAttrList"].keys()) == DEFAULT_NUM_GLOBAL_ATTRS

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_writer.to_cdf(output_path=test_cache)

    # print("\nDEBUG THE CDF FILE OUTPUT")
    # print(f"Global Attrs: {len(test_writer.target_cdf_file.attrs)}")
    # print(test_writer.target_cdf_file.attrs)
    # print()
    # print(test_writer.target_cdf_file)

    # Test number of Global Attrs in the generated CDF File (Result Data)
    assert len(test_writer.target_cdf_file.attrs) == DEFAULT_NUM_GLOBAL_ATTRS

    # Test the Number of Variable Attrs in the Gnerated CDF File
    assert len(test_writer.target_cdf_file) == len(input_data["zAttrList"].keys())

    # Save the CDF to a File
    test_writer.save_cdf()

    test_file_cache_path = Path(test_file_output_path)
    # Test the File Exists
    assert test_file_cache_path.exists()

    # Remove the Test File from Cache
    test_file_cache_path.unlink()

    # Test the File was Deleted
    assert not test_file_cache_path.exists()


def test_cdf_writer_invalid_seed_data():

    # fmt: off
    seed_data = [1, 2, 3]
    # fmt: on

    with pytest.raises(AssertionError) as e:
        # Initialize a CDF File Wrapper
        test_writer = CDFWriter(seed_data=seed_data)


def test_cdf_writer_valid_seed_data():

    # fmt: off
    seed_data = {
        "gAttrList": {
            "Descriptor": {
                "value": "EEA>Electron Electrostatic Analyzer"
            }
        },
        "zAttrList": {
        },
    }
    # fmt: on

    # Initialize a CDF File Wrapper
    test_writer = CDFWriter(seed_data=seed_data)

    # Test Number of Global Attrs in the Target Dict (Intermediate Data)
    assert len(test_writer.target_dict["gAttrList"].keys()) == DEFAULT_NUM_GLOBAL_ATTRS


def test_cdf_writer_overide_derived_attr():

    # fmt: off
    seed_data = {
        "gAttrList": {
            "Descriptor": {
                "value": "EEA>Electron Electrostatic Analyzer"
            },
            "Data_type": {
                "value": 'test>data_type',
                "required": True,
                "valid_check": None,
                "derived": {
                    'method': 'get_data_type'
                }
            },
        },
        "zAttrList": {
        },
    }
    # fmt: on

    # Initialize a CDF File Wrapper
    test_writer = CDFWriter(seed_data=seed_data)

    # Test Number of Global Attrs in the Target Dict (Intermediate Data)
    assert len(test_writer.target_dict["gAttrList"].keys()) == DEFAULT_NUM_GLOBAL_ATTRS

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    test_file_output_path = test_writer.to_cdf(output_path=test_cache)

    # Test number of Global Attrs in the generated CDF File (Result Data)
    assert len(test_writer.target_cdf_file.attrs) == DEFAULT_NUM_GLOBAL_ATTRS

    # print(test_writer.target_cdf_file.attrs["Data_type"].strip())

    # Test the Value was not Derrived and used the Overriden Value
    assert (
        str(test_writer.target_cdf_file.attrs["Data_type"])
        == seed_data["gAttrList"]["Data_type"]["value"]
    )

    # Save the CDF to a File
    test_writer.save_cdf()

    test_file_cache_path = Path(test_file_output_path)
    # Test the File Exists
    assert test_file_cache_path.exists()

    # Remove the Test File from Cache
    test_file_cache_path.unlink()

    # Test the File was Deleted
    assert not test_file_cache_path.exists()


def test_cdf_writer_invalid_derive_method():

    # fmt: off
    seed_data = {
        "gAttrList": {
            "Descriptor": {
                "value": "EEA>Electron Electrostatic Analyzer"
            },
            "Data_type": {
                "value": 'test>data_type',
                "required": True,
                "valid_check": None,
                "derived": {
                    'method': 'invalid_method'
                }
            },
        },
        "zAttrList": {
        },
    }
    # fmt: on

    # Initialize a CDF File Wrapper
    test_writer = CDFWriter(seed_data=seed_data)

    # Test Number of Global Attrs in the Target Dict (Intermediate Data)
    assert len(test_writer.target_dict["gAttrList"].keys()) == DEFAULT_NUM_GLOBAL_ATTRS

    # Convert the Wrapper to a CDF File
    test_cache = Path(hermes_core.__file__).parent.parent / ".pytest_cache"
    with pytest.raises(ValueError) as e:
        test_writer.to_cdf(output_path=test_cache)

