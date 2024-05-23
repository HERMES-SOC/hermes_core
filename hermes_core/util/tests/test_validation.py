from collections import OrderedDict
from pathlib import Path
import pytest
import tempfile
import numpy as np
from numpy.random import random
import datetime
from astropy.time import Time
from astropy.timeseries import TimeSeries
import astropy.units as u
from spacepy.pycdf import CDF
import hermes_core
from hermes_core.timedata import HermesData
from hermes_core.util import const
from hermes_core.util.validation import validate, CDFValidator

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
    invlid_path = Path(hermes_core.__file__).parent / "data" / "README.rst"
    with pytest.raises(ValueError):
        _ = validate(invlid_path)


def test_non_existant_file():
    """Function to Test a file using the CDFValidator that does not exist"""
    invlid_path = Path(hermes_core.__file__).parent / "data" / "test.cdf"
    result = validate(invlid_path)
    assert len(result) == 1
    assert "Could not open CDF File at path:" in result[0]


def test_missing_global_attrs():
    """Function to ensure missing global attributes are reported in validation"""

    # Create a Test HermesData
    ts = get_test_timeseries()
    template = HermesData.global_attribute_template("eea", "l2", "0.0.0")
    td = HermesData(timeseries=ts, meta=template)

    # Convert to a CDF File and Validate
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmp_path = Path(tmpdirname)
        out_file = td.save(tmp_path)

        with CDF(str(out_file), readonly=False) as cdf:
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
    td = HermesData(timeseries=ts, meta=template)

    # Convert to a CDF File and Validate
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmp_path = Path(tmpdirname)
        out_file = td.save(tmp_path)

        with CDF(str(out_file), readonly=False) as cdf:
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
    td = HermesData(timeseries=ts, meta=template)

    # Convert to a CDF File and Validate
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmp_path = Path(tmpdirname)
        out_file = td.save(tmp_path)

        with CDF(str(out_file), readonly=False) as cdf:
            del cdf["measurement"].meta["CATDESC"]
            del cdf["measurement"].meta["UNITS"]
            cdf["measurement"].meta["DISPLAY_TYPE"] = "bad_type"

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


def test_valid_range_dimensioned():
    """Validmin/validmax with multiple elements"""
    with tempfile.TemporaryDirectory() as tmpdirname:
        # Create a Test CDF
        cdf = CDF(tmpdirname + "test.cdf", create=True)

        v = cdf.new("var1", data=[[1, 10], [2, 20], [3, 30]])
        v.attrs["VALIDMIN"] = [1, 20]
        v.attrs["VALIDMAX"] = [3, 30]
        v.attrs["FILLVAL"] = -100
        errs = CDFValidator()._validrange(v)
        assert 1 == len(errs)
        assert "Value 10 at index [0 1] under VALIDMIN [ 1 20]." == errs[0]
        v.attrs["VALIDMIN"] = [1, 10]
        errs = CDFValidator()._validrange(v)
        assert 0 == len(errs)
        v[0, 0] = -100
        errs = CDFValidator()._validrange(v)
        assert 0 == len(errs)


def test_valid_range_dimension_mismatch():
    """Validmin/validmax with something wrong in dimensionality"""
    with tempfile.TemporaryDirectory() as tmpdirname:
        # Create a Test CDF
        cdf = CDF(tmpdirname + "test.cdf", create=True)

        v = cdf.new("var1", data=[[1, 10], [2, 20], [3, 30]])
        v.attrs["VALIDMIN"] = [1, 10, 100]
        v.attrs["VALIDMAX"] = [3, 30, 127]
        errs = CDFValidator()._validrange(v)
        assert 2 == len(errs)
        assert (
            "VALIDMIN element count 3 does not match "
            "first data dimension size 2." == errs[0]
        )
        assert (
            "VALIDMAX element count 3 does not match "
            "first data dimension size 2." == errs[1]
        )


def test_valid_range_high_dimension():
    """Validmin/validmax with high-dimension variables"""
    with tempfile.TemporaryDirectory() as tmpdirname:
        # Create a Test CDF
        cdf = CDF(tmpdirname + "test.cdf", create=True)

        v = cdf.new(
            "var1",
            data=np.reshape(
                np.arange(27.0),
                (
                    3,
                    3,
                    3,
                ),
            ),
        )
        v.attrs["VALIDMIN"] = [1, 10, 100]
        v.attrs["VALIDMAX"] = [3, 30, 300]
        errs = CDFValidator()._validrange(v)
        assert 2 == len(errs)
        assert "Multi-element VALIDMIN only valid with 1D variable." == errs[0]
        assert "Multi-element VALIDMAX only valid with 1D variable." == errs[1]


def test_valid_range_wrong_type():
    """Validmin/validmax not matching variable type"""
    with tempfile.TemporaryDirectory() as tmpdirname:
        # Create a Test CDF
        cdf = CDF(tmpdirname + "test.cdf", create=True)

        v = cdf.new("var1", data=[1, 2, 3], type=const.CDF_INT4)
        v.attrs.new("VALIDMIN", data=1, type=const.CDF_INT2)
        v.attrs.new("VALIDMAX", data=3, type=const.CDF_INT2)
        errs = CDFValidator()._validrange(v)
        errs.sort()
        assert 0 == len(errs)


def test_valid_range_incompatible_type():
    """Validmin/validmax can't be compared to variable type"""
    with tempfile.TemporaryDirectory() as tmpdirname:
        # Create a Test CDF
        cdf = CDF(tmpdirname + "test.cdf", create=True)

        v = cdf.new("var1", data=[1, 2, 3], type=const.CDF_INT4)
        v.attrs.new("VALIDMIN", data="2")
        v.attrs.new("VALIDMAX", data="5")
        errs = CDFValidator()._validrange(v)
        errs.sort()
        assert 2 == len(errs)
        assert [
            "VALIDMAX type CDF_CHAR not comparable to variable type CDF_INT4.",
            "VALIDMIN type CDF_CHAR not comparable to variable type CDF_INT4.",
        ] == errs


def test_valid_range_nrv():
    """Validmin/validmax"""
    with tempfile.TemporaryDirectory() as tmpdirname:
        # Create a Test CDF
        cdf = CDF(tmpdirname + "test.cdf", create=True)

        v = cdf.new("var1", recVary=False, data=[1, 2, 3])
        v.attrs["VALIDMIN"] = 1
        v.attrs["VALIDMAX"] = 3
        assert 0 == len(CDFValidator()._validrange(v))
        v.attrs["VALIDMIN"] = 2
        errs = CDFValidator()._validrange(v)
        assert 1 == len(errs)
        assert "Value 1 at index 0 under VALIDMIN 2." == errs[0]
        v.attrs["VALIDMAX"] = 2
        errs = CDFValidator()._validrange(v)
        assert 2 == len(errs)
        assert "Value 3 at index 2 over VALIDMAX 2." == errs[1]


def test_valid_range_nrv_fillval():
    """Validmin/validmax with fillval set"""
    with tempfile.TemporaryDirectory() as tmpdirname:
        # Create a Test CDF
        cdf = CDF(tmpdirname + "test.cdf", create=True)

        v = cdf.new("var1", recVary=False, data=[1, 2, 3])
        v.attrs["VALIDMIN"] = 1
        v.attrs["VALIDMAX"] = 3
        v.attrs["FILLVAL"] = 99
        assert 0 == len(CDFValidator()._validrange(v))

        v.attrs["VALIDMIN"] = 2
        errs = CDFValidator()._validrange(v)
        assert 1 == len(errs)
        assert "Value 1 at index 0 under VALIDMIN 2." == errs[0]

        v.attrs["VALIDMAX"] = 2
        errs = CDFValidator()._validrange(v)
        assert 2 == len(errs)
        assert "Value 3 at index 2 over VALIDMAX 2." == errs[1]

        v.attrs["FILLVAL"] = 3
        errs = CDFValidator()._validrange(v)
        assert 1 == len(errs)
        assert "Value 1 at index 0 under VALIDMIN 2." == errs[0]

        v.attrs["FILLVAL"] = 1
        errs = CDFValidator()._validrange(v)
        assert 1 == len(errs)
        assert "Value 3 at index 2 over VALIDMAX 2." == errs[0]


def test_valid_range_fillval_float():
    """Validmin/validmax with fillval set, floating-point"""
    with tempfile.TemporaryDirectory() as tmpdirname:
        # Create a Test CDF
        cdf = CDF(tmpdirname + "test.cdf", create=True)

        v = cdf.new("var1", recVary=False, data=[1, 2, 3], type=const.CDF_DOUBLE)
        v.attrs["VALIDMIN"] = 0
        v.attrs["VALIDMAX"] = 10
        # This is a bit contrived to force a difference between attribute
        # and value that's only the precision of the float
        v.attrs.new("FILLVAL", -1e31, type=const.CDF_FLOAT)
        assert 0 == len(CDFValidator()._validrange(v))

        v[0] = -100
        errs = CDFValidator()._validrange(v)
        assert 1 == len(errs)
        assert "Value -100.0 at index 0 under VALIDMIN 0.0." == errs[0]

        v[0] = -1e31
        assert 0 == len(CDFValidator()._validrange(v))


def test_valid_range_fillval_float_wrong_type():
    """Validmin/validmax with fillval, floating-point, but fillval string"""
    with tempfile.TemporaryDirectory() as tmpdirname:
        # Create a Test CDF
        cdf = CDF(tmpdirname + "test.cdf", create=True)

        v = cdf.new("var1", recVary=False, data=[-1e31, 2, 3], type=const.CDF_DOUBLE)
        v.attrs["VALIDMIN"] = 0
        v.attrs["VALIDMAX"] = 10
        v.attrs.new("FILLVAL", b"badstuff", type=const.CDF_CHAR)
        expected = ["Value -1e+31 at index 0 under VALIDMIN 0.0."]
        errs = CDFValidator()._validrange(v)
        assert len(expected) == len(errs)
        for a, e in zip(sorted(errs), sorted(expected)):
            assert e == a


def test_valid_range_fillval_datetime():
    """Validmin/validmax with fillval set, Epoch var"""
    with tempfile.TemporaryDirectory() as tmpdirname:
        # Create a Test CDF
        cdf = CDF(tmpdirname + "test.cdf", create=True)

        v = cdf.new(
            "var1",
            data=[datetime.datetime(2010, 1, i) for i in range(1, 6)],
            type=const.CDF_EPOCH,
        )
        v.attrs["VALIDMIN"] = datetime.datetime(2010, 1, 1)
        v.attrs["VALIDMAX"] = datetime.datetime(2010, 1, 31)
        v.attrs["FILLVAL"] = datetime.datetime(9999, 12, 31, 23, 59, 59, 999000)
        assert 0 == len(CDFValidator()._validrange(v))

        v[-1] = datetime.datetime(2010, 2, 1)
        errs = CDFValidator()._validrange(v)
        assert 1 == len(errs)
        assert (
            "Value 2010-02-01 00:00:00 at index 4 over VALIDMAX "
            "2010-01-31 00:00:00." == errs[0]
        )

        v[-1] = datetime.datetime(9999, 12, 31, 23, 59, 59, 999000)
        assert 0 == len(CDFValidator()._validrange(v))


def test_valid_range_scalar():
    """Check validmin/max on a scalar"""
    with tempfile.TemporaryDirectory() as tmpdirname:
        # Create a Test CDF
        cdf = CDF(tmpdirname + "test.cdf", create=True)

        v = cdf.new("var1", recVary=False, data=1)
        v.attrs["VALIDMIN"] = 0
        v.attrs["VALIDMAX"] = 2
        v.attrs["FILLVAL"] = -100
        assert 0 == len(CDFValidator()._validrange(v))
        v.attrs["VALIDMIN"] = 2
        v.attrs["VALIDMAX"] = 3
        errs = CDFValidator()._validrange(v)
        assert 1 == len(errs)
        assert "Value 1 under VALIDMIN 2." == errs[0]


def test_valid_scale():
    """Check scale min and max."""
    with tempfile.TemporaryDirectory() as tmpdirname:
        # Create a Test CDF
        cdf = CDF(tmpdirname + "test.cdf", create=True)

        v = cdf.new("var1", recVary=False, data=[1, 2, 3])
        v.attrs["SCALEMIN"] = 1
        v.attrs["SCALEMAX"] = 3
        assert 0 == len(CDFValidator()._validscale(v))
        v.attrs["SCALEMIN"] = 5
        v.attrs["SCALEMAX"] = 3
        assert 1 == len(CDFValidator()._validscale(v))
        errs = CDFValidator()._validscale(v)
        assert "SCALEMIN > SCALEMAX." == errs[0]
        v.attrs["SCALEMIN"] = -200
        errs = CDFValidator()._validscale(v)
        assert 1 == len(errs)
        errs.sort()
        assert ["SCALEMIN (-200) outside valid data range (-128,127)."] == errs
        v.attrs["SCALEMIN"] = 200
        errs = CDFValidator()._validscale(v)
        assert 2 == len(errs)
        errs.sort()
        assert [
            "SCALEMIN (200) outside valid data range (-128,127).",
            "SCALEMIN > SCALEMAX.",
        ] == errs
        v.attrs["SCALEMAX"] = -200
        errs = CDFValidator()._validscale(v)
        assert 3 == len(errs)
        errs.sort()
        assert [
            "SCALEMAX (-200) outside valid data range (-128,127).",
            "SCALEMIN (200) outside valid data range (-128,127).",
            "SCALEMIN > SCALEMAX.",
        ] == errs
        v.attrs["SCALEMAX"] = 200
        errs = CDFValidator()._validscale(v)
        assert 2 == len(errs)
        errs.sort()
        assert [
            "SCALEMAX (200) outside valid data range (-128,127).",
            "SCALEMIN (200) outside valid data range (-128,127).",
        ] == errs


def test_valid_scale_dimensioned():
    """Validmin/validmax with multiple elements"""
    with tempfile.TemporaryDirectory() as tmpdirname:
        # Create a Test CDF
        cdf = CDF(tmpdirname + "test.cdf", create=True)

        v = cdf.new("var1", data=[[1, 10], [2, 20], [3, 30]])
        v.attrs["SCALEMIN"] = [2, 20]
        v.attrs["SCALEMAX"] = [300, 320]
        v.attrs["FILLVAL"] = -100
        errs = CDFValidator()._validscale(v)
        assert 1 == len(errs)
        errs.sort()
        assert [
            "SCALEMAX ([300 320]) outside valid data range (-128,127).",
        ] == errs
        v.attrs["SCALEMAX"] = [30, 32]
        errs = CDFValidator()._validscale(v)
        assert 0 == len(errs)


def test_valid_scale_dimension_mismatch():
    """Validmin/validmax with something wrong in dimensionality"""
    with tempfile.TemporaryDirectory() as tmpdirname:
        # Create a Test CDF
        cdf = CDF(tmpdirname + "test.cdf", create=True)

        v = cdf.new("var1", data=[[1, 10], [2, 20], [3, 30]])
        v.attrs["SCALEMIN"] = [1, 10, 100]
        v.attrs["SCALEMAX"] = [3, 30, 126]
        errs = CDFValidator()._validscale(v)
        assert 2 == len(errs)
        errs.sort()
        assert [
            "SCALEMAX element count 3 does not match " "first data dimension size 2.",
            "SCALEMIN element count 3 does not match " "first data dimension size 2.",
        ] == errs


def test_valid_scale_high_dimension():
    """scalemin/scalemax with high-dimension variables"""
    with tempfile.TemporaryDirectory() as tmpdirname:
        # Create a Test CDF
        cdf = CDF(tmpdirname + "test.cdf", create=True)

        v = cdf.new(
            "var1",
            data=np.reshape(
                np.arange(27.0),
                (
                    3,
                    3,
                    3,
                ),
            ),
        )
        v.attrs["SCALEMIN"] = [1, 10, 100]
        v.attrs["SCALEMAX"] = [3, 30, 300]
        errs = CDFValidator()._validscale(v)
        assert 2 == len(errs)
        assert "Multi-element SCALEMIN only valid with 1D variable." == errs[0]
        assert "Multi-element SCALEMAX only valid with 1D variable." == errs[1]
