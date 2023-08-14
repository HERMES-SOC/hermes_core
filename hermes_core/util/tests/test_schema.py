import pytest
import tempfile
import yaml
import datetime
from collections import OrderedDict
import numpy as np
from numpy.random import random
from astropy.time import Time
from astropy.timeseries import TimeSeries
from astropy.table import Table
import astropy.units as u
from spacepy.pycdf import CDF
from hermes_core.timedata import TimeData
from hermes_core.util.schema import HERMESDataSchema
from hermes_core.util import const


def test_hermes_data_schema():
    """Test Schema Template and Info Functions"""
    schema = HERMESDataSchema()

    # Global Attribute Schema
    assert schema.global_attribute_schema is not None
    assert isinstance(schema.global_attribute_schema, dict)

    # Variable Attribute Schema
    assert schema.variable_attribute_schema is not None
    assert isinstance(schema.variable_attribute_schema, dict)

    # Default Globla Attributes
    assert schema.default_global_attributes is not None
    assert isinstance(schema.default_global_attributes, dict)

    # Global Attribute Template
    assert schema.global_attribute_template() is not None
    assert isinstance(schema.global_attribute_template(), OrderedDict)

    # Measurement Attribute Template
    assert schema.measurement_attribute_template() is not None
    assert isinstance(schema.measurement_attribute_template(), OrderedDict)

    # Global Attribute Info
    assert schema.global_attribute_info() is not None
    assert isinstance(schema.global_attribute_info(), Table)
    assert isinstance(schema.global_attribute_info(attribute_name="Descriptor"), Table)
    with pytest.raises(KeyError):
        _ = schema.global_attribute_info(attribute_name="NotAnAttribute")

    # Measurement Attribute Info
    assert schema.measurement_attribute_info() is not None
    assert isinstance(schema.measurement_attribute_info(), Table)
    assert isinstance(
        schema.measurement_attribute_info(attribute_name="CATDESC"), Table
    )
    with pytest.raises(KeyError):
        _ = schema.measurement_attribute_info(attribute_name="NotAnAttribute")


def test_load_yaml_data():
    """Test Loading Yaml Data for Schema Files"""
    with tempfile.TemporaryDirectory() as tmpdirname:
        # This function writes invalid YAML content into a file
        invalid_yaml = """
        name: John Doe
        age 30
        """

        with open(tmpdirname + "test.yaml", "w") as file:
            file.write(invalid_yaml)

        # Load from an non-existant file
        yaml_data = HERMESDataSchema._load_yaml_data(tmpdirname + "test.yaml")
        assert yaml_data == {}


def test_global_attributes():
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
    template = TimeData.global_attribute_template("eea", "l2", "0.0.0")

    # Create TimeData
    td = TimeData(data=ts, meta=template)
    with pytest.raises(ValueError):
        _ = HERMESDataSchema()._derive_global_attribute(td, "bad_attribute")


def test_check_well_formed():
    """Test that the Data can be output to CDF"""

    # Badly Shaped Data
    data = [np.array([1, 2, 3]), np.array([[4, 5], [6, 7]])]
    with pytest.raises(ValueError):
        _ = HERMESDataSchema._check_well_formed(data)

    # Empty Data
    data = np.array([], dtype=object)
    d = HERMESDataSchema._check_well_formed(data)
    assert len(d) == 0

    # Badly Shaped Object
    data = np.array([np.array([1, 2, 3]), np.array([[4, 5], [6, 7]])], dtype=object)
    with pytest.raises(ValueError):
        _ = HERMESDataSchema._check_well_formed(data)


def test_types():
    """Function to test getting the CDF data types for different data types"""

    # String Type
    _, types, _ = HERMESDataSchema()._types("")
    assert types == [51, 52]

    # Time Types
    # data = np.array(
    #     ["2023-08-03T10:20:30.123456789", "2023-08-03T15:30:45.678901234"], dtype="datetime64[ns]"
    # )
    # _, types, _ = HERMESDataSchema()._types(data)
    # assert types == [const.CDF_TIME_TT2000, const.CDF_EPOCH16, const.CDF_EPOCH]


def test_type_guessing():
    """Guess CDF types based on input data"""
    samples = [
        [1, 2, 3, 4],
        [[1.2, 1.3, 1.4], [2.2, 2.3, 2.4]],
        ["hello", "there", "everybody"],
        datetime.datetime(2009, 1, 1),
        datetime.datetime(2009, 1, 1, 12, 15, 12, 1000),
        datetime.datetime(2009, 1, 1, 12, 15, 12, 1),
        [1.0],
        0.0,
        np.array([1, 2, 3], dtype=np.int32),
        np.array([1, 2, 4], dtype=np.float64),
        np.array([1, 2, 5], dtype=np.int64),
        2**62,
        -1.0,
        np.array([1, 2, 6], dtype="<u2"),
        np.array([1, 2, 7], dtype=">u2"),
        np.int64(-1 * 2**63),
        np.int32(-1 * 2**31),
        -1 * 2**31,
        np.array([5, 6, 7], dtype=np.uint8),
        [4611686018427387904],
        np.array([1], dtype=object),
        ["\U0001f600\U0001f600"],
    ]
    types = [
        (
            (4,),
            [
                const.CDF_BYTE,
                const.CDF_INT1,
                const.CDF_UINT1,
                const.CDF_INT2,
                const.CDF_UINT2,
                const.CDF_INT4,
                const.CDF_UINT4,
                const.CDF_INT8,
                const.CDF_FLOAT,
                const.CDF_REAL4,
                const.CDF_DOUBLE,
                const.CDF_REAL8,
            ],
            1,
        ),
        (
            (2, 3),
            [const.CDF_FLOAT, const.CDF_REAL4, const.CDF_DOUBLE, const.CDF_REAL8],
            1,
        ),
        ((3,), [const.CDF_CHAR, const.CDF_UCHAR], 9),
        ((), [const.CDF_TIME_TT2000, const.CDF_EPOCH, const.CDF_EPOCH16], 1),
        ((), [const.CDF_TIME_TT2000, const.CDF_EPOCH, const.CDF_EPOCH16], 1),
        ((), [const.CDF_TIME_TT2000, const.CDF_EPOCH16, const.CDF_EPOCH], 1),
        (
            (1,),
            [const.CDF_FLOAT, const.CDF_REAL4, const.CDF_DOUBLE, const.CDF_REAL8],
            1,
        ),
        ((), [const.CDF_FLOAT, const.CDF_REAL4, const.CDF_DOUBLE, const.CDF_REAL8], 1),
        ((3,), [const.CDF_INT4], 1),
        ((3,), [const.CDF_DOUBLE, const.CDF_REAL8], 1),
        ((3,), [const.CDF_INT8], 1),
        (
            (),
            [
                const.CDF_INT8,
                const.CDF_FLOAT,
                const.CDF_REAL4,
                const.CDF_DOUBLE,
                const.CDF_REAL8,
            ],
            1,
        ),
        ((), [const.CDF_FLOAT, const.CDF_REAL4, const.CDF_DOUBLE, const.CDF_REAL8], 1),
        ((3,), [const.CDF_UINT2], 1),
        ((3,), [const.CDF_UINT2], 1),
        ((), [const.CDF_INT8], 1),
        ((), [const.CDF_INT4], 1),
        (
            (),
            [
                const.CDF_INT4,
                const.CDF_INT8,
                const.CDF_FLOAT,
                const.CDF_REAL4,
                const.CDF_DOUBLE,
                const.CDF_REAL8,
            ],
            1,
        ),
        ((3,), [const.CDF_UINT1, const.CDF_UCHAR], 1),
        (
            (1,),
            [
                const.CDF_INT8,
                const.CDF_FLOAT,
                const.CDF_REAL4,
                const.CDF_DOUBLE,
                const.CDF_REAL8,
            ],
            1,
        ),
        (
            (1,),
            [
                const.CDF_BYTE,
                const.CDF_INT1,
                const.CDF_UINT1,
                const.CDF_INT2,
                const.CDF_UINT2,
                const.CDF_INT4,
                const.CDF_UINT4,
                const.CDF_INT8,
                const.CDF_FLOAT,
                const.CDF_REAL4,
                const.CDF_DOUBLE,
                const.CDF_REAL8,
            ],
            1,
        ),
        ((1,), [const.CDF_CHAR, const.CDF_UCHAR], 8),
    ]
    with pytest.raises(ValueError):
        HERMESDataSchema()._types([object()])
    for s, t in zip(samples, types):
        t = (t[0], [i.value for i in t[1]], t[2])
        assert t == HERMESDataSchema()._types(s)


def test_min_max_none():
    """Get min/max values for None types"""
    with pytest.raises(ValueError):
        HERMESDataSchema()._get_minmax(None)


def test_min_max_unknown():
    """Get min/max values for unknown types"""
    with pytest.raises(ValueError):
        HERMESDataSchema()._get_minmax("unknown_type")


def test_min_max_TT2000():
    """Get min/max values for TT2000 types"""
    minval, maxval = HERMESDataSchema()._get_minmax(const.CDF_TIME_TT2000)
    # Make sure the minimum isn't just plain invalid
    assert minval == datetime.datetime(1, 1, 1)
    assert maxval == datetime.datetime(2292, 4, 11, 11, 46, 7, 670776)


def test_min_max_Epoch16():
    """Get min/max values for Epoch16 types"""
    minval, maxval = HERMESDataSchema()._get_minmax(const.CDF_EPOCH16)
    # Make sure the minimum isn't just plain invalid
    assert minval == datetime.datetime(1, 1, 1)
    assert maxval == datetime.datetime(9999, 12, 31, 23, 59, 59)


def test_min_max_Epoch():
    """Get min/max values for EPOCH types"""
    minval, maxval = HERMESDataSchema()._get_minmax(const.CDF_EPOCH)
    # Make sure the minimum isn't just plain invalid
    assert minval == datetime.datetime(1, 1, 1)
    assert maxval == datetime.datetime(9999, 12, 31, 23, 59, 59)


def test_min_max_Float():
    """Get min/max values for a float"""
    minval, maxval = HERMESDataSchema()._get_minmax(const.CDF_FLOAT)
    np.allclose(-3.4028234663853e38, minval)
    np.allclose(3.4028234663853e38, maxval)


def test_min_max_Int():
    """Get min/max values for an integer"""
    minval, maxval = HERMESDataSchema()._get_minmax(const.CDF_INT1)
    assert -128 == minval
    assert 127 == maxval


def test_format():
    """Set the format"""
    # This is done by: type, expected, validmin, validmax (None okay)
    expected = (
        (const.CDF_EPOCH.value, "A24", None, None),
        (const.CDF_EPOCH16.value, "A36", None, None),
        (const.CDF_TIME_TT2000.value, "A29", None, None),
        (const.CDF_INT2.value, "I2", -2, 9),
        (const.CDF_UINT2.value, "I5", None, None),
        (const.CDF_INT2.value, "I6", None, None),
        (const.CDF_UINT2.value, "I2", 0, 10),
        (const.CDF_FLOAT.value, "F6.3", 1, 10),
        (const.CDF_FLOAT.value, "F6.2", 1, 100),
        (const.CDF_FLOAT.value, "F6.1", 1, 1000),
        (const.CDF_FLOAT.value, "G10.2E3", 1, -1),
    )

    with tempfile.TemporaryDirectory() as tmpdirname:
        cdf = CDF(tmpdirname + "test.cdf", create=True)

        for t, e, vmin, vmax in expected:
            v = cdf.new("var", type=t)
            if vmin is not None:
                v.attrs["VALIDMIN"] = vmin
            if vmin is not None:
                v.attrs["VALIDMAX"] = vmax
            format = HERMESDataSchema()._get_format(cdf["var"], t)
            assert e == format
            del cdf["var"]

        # Test Format Char
        v = cdf.new("var", data=["hi", "there"])
        format = HERMESDataSchema()._get_format(cdf["var"], const.CDF_CHAR.value)
        assert "A2" == format


def test_resolution():
    """Test Time Resolution"""
    ts = TimeSeries()
    # Create an astropy.Time object
    time = np.arange(1)
    time_col = Time(time, format="unix")
    ts["time"] = time_col
    ts["time"].meta = OrderedDict({"CATDESC": "Epoch Time"})

    # Add Measurement
    quant = u.Quantity(value=random(size=(1)), unit="m", dtype=np.uint16)
    ts["measurement"] = quant
    ts["measurement"].meta = OrderedDict(
        {
            "VAR_TYPE": "data",
            "CATDESC": "Test Data",
        }
    )
    template = TimeData.global_attribute_template("eea", "l2", "0.0.0")

    # Get Resolution
    with pytest.raises(ValueError):
        td = TimeData(data=ts, meta=template)


def test_reference_position():
    """Function to test time reference position"""
    assert (
        HERMESDataSchema()._get_reference_position(const.CDF_TIME_TT2000.value)
        == "rotating Earth geoid"
    )

    with pytest.raises(TypeError):
        HERMESDataSchema()._get_reference_position(const.CDF_EPOCH.value)


def test_time_base():
    """Function to test time base"""
    assert HERMESDataSchema()._get_time_base(const.CDF_TIME_TT2000.value) == "J2000"

    with pytest.raises(TypeError):
        HERMESDataSchema()._get_time_base(const.CDF_EPOCH.value)


def test_time_scale():
    """Function to test time scale"""
    assert (
        HERMESDataSchema()._get_time_scale(const.CDF_TIME_TT2000.value)
        == "Terrestrial Time (TT)"
    )

    with pytest.raises(TypeError):
        HERMESDataSchema()._get_time_scale(const.CDF_EPOCH.value)


def test_time_units():
    """Function to test time units"""
    assert HERMESDataSchema()._get_time_units(const.CDF_TIME_TT2000.value) == "ns"

    with pytest.raises(TypeError):
        HERMESDataSchema()._get_time_units(const.CDF_EPOCH.value)
