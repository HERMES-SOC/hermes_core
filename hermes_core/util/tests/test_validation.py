from collections import OrderedDict
from pathlib import Path
import pytest
import tempfile
import numpy as np
from numpy.random import random
from astropy.time import Time
from astropy.timeseries import TimeSeries
import astropy.units as u
from spacepy.pycdf import CDF
import hermes_core
from hermes_core.timedata import HermesData
from hermes_core.util.schema import HermesDataSchema
from hermes_core.util.validation import validate

SAMPLE_CDF_FILE = "hermes_nms_default_l1_20160322_123031_v0.0.1.cdf"


def get_test_timeseries():
    """Get Test Data"""
    ts = TimeSeries()

    # Create an astropy.Time object
    time = np.arange(10)
    time_col = Time(time, format="unix")
    ts["time"] = time_col
    ts["time"].meta = OrderedDict({"CATDESC": "Epoch Time"})

    # Add Measurement
    quant = u.Quantity(value=random(size=(10)), unit="m", dtype=np.uint16)
    ts["measurement"] = quant
    ts["measurement"].meta = OrderedDict(
        {
            "VAR_TYPE": "data",
            "CATDESC": "Test Data",
        }
    )
    return ts


def test_non_cdf_file():
    """Function to Test a file using the CDFValidator that is not a CDF File"""
    invlid_path = str(Path(hermes_core.__file__).parent / "data" / "README.rst")
    with pytest.raises(ValueError):
        _ = validate(invlid_path)


def test_non_existant_file():
    """Function to Test a file using the CDFValidator that does not exist"""
    invlid_path = str(Path(hermes_core.__file__).parent / "data" / "test.cdf")
    result = validate(invlid_path)
    assert len(result) == 1
    assert "Could not open CDF File at path:" in result[0]


def test_missing_global_attrs():
    """Function to ensure missing global attributes are reported in validation"""

    # Create a Test HermesData
    ts = get_test_timeseries()
    template = HermesData.global_attribute_template("eea", "l2", "0.0.0")
    td = HermesData(data=ts, meta=template)

    # Convert to a CDF File and Validate
    with tempfile.TemporaryDirectory() as tmpdirname:
        out_file = td.save(tmpdirname)

        with CDF(out_file, readonly=False) as cdf:
            del cdf.meta["Descriptor"]

        # Validate
        result = validate(out_file)
        assert (
            "Required attribute (Descriptor) not present in global attributes."
            in result
        )
        assert "Required attribute (DOI) not present in global attributes." in result


def test_missing_var_type():
    """Function to ensure missing variable attributes are reported in validation"""

    # Create a Test HermesData
    ts = get_test_timeseries()
    template = HermesData.global_attribute_template("eea", "l2", "0.0.0")
    td = HermesData(data=ts, meta=template)

    # Convert to a CDF File and Validate
    with tempfile.TemporaryDirectory() as tmpdirname:
        out_file = td.save(tmpdirname)

        with CDF(out_file, readonly=False) as cdf:
            del cdf["measurement"].meta["VAR_TYPE"]

        # Validate
        result = validate(out_file)
        assert (
            "Variable: measurement missing 'VAR_TYPE' attribute. Cannot Validate Variable."
            in result
        )


def test_missing_variable_attrs():
    """Function to ensure missing variable attributes are reported in validation"""

    # Create a Test HermesData
    ts = get_test_timeseries()
    template = HermesData.global_attribute_template("eea", "l2", "0.0.0")
    td = HermesData(data=ts, meta=template)

    # Convert to a CDF File and Validate
    with tempfile.TemporaryDirectory() as tmpdirname:
        out_file = td.save(tmpdirname)

        with CDF(out_file, readonly=False) as cdf:
            del cdf["measurement"].meta["CATDESC"]
            del cdf["measurement"].meta["UNITS"]
            cdf["measurement"].meta["DISPLAY_TYPE"] = "bad_type"
            cdf["measurement"].meta["FORMAT"] = "bad_format"

        # Validate
        result = validate(out_file)
        assert "Variable: measurement missing 'CATDESC' attribute." in result
        assert (
            "Variable: measurement missing 'UNITS' attribute. Alternative: UNIT_PTR not found."
            in result
        )
        assert (
            "Variable: measurement Attribute 'DISPLAY_TYPE' not one of valid options.",
            "Was bad_type, expected one of time_series time_series>noerrorbars spectrogram stack_plot image",
        ) in result
        assert (
            "Variable: measurement Attribute 'FORMAT' value 'bad_format' does not match derrived format 'I5'"
            in result
        )
