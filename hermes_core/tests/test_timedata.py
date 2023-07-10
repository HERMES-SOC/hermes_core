"""Tests for CDF Files to and from data containers"""

from collections import OrderedDict
from pathlib import Path
import datetime
import pytest
import numpy as np
from numpy.random import random
import tempfile
from astropy.timeseries import TimeSeries
from astropy.table import Column
from astropy.time import Time
from astropy.units import Quantity
import astropy.units as u
from spacepy.pycdf import CDFError
from matplotlib.axes import Axes
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
    col = Column(
        data=random(size=(10)), name="measurement", meta={"CATDESC": "Test Measurement"}
    )
    ts.add_column(col)
    return ts


def get_test_timeseries():
    ts = TimeSeries()

    # Create an astropy.Time object
    time = np.arange(10)
    time_col = Time(time, format="unix")
    ts["time"] = time_col

    # Add Measurement
    ts["measurement"] = Quantity(value=random(size=(10)), unit="m", dtype=np.uint16)
    ts["measurement"].meta = OrderedDict(
        {
            "VAR_TYPE": "metadata",
            "CATDESC": "Test Metadata",
        }
    )

    # Add 2D Variable
    ts["2d-var"] = Quantity(value=random(size=(10, 2)), unit="s", dtype=np.uint16)
    ts["2d-var"].meta = OrderedDict({"CATDESC": "Test 2-Dimensional Variable"})

    # Add 3D Variable
    ts["3d-var"] = Quantity(value=random(size=(10, 2, 2)), unit="s", dtype=np.uint16)
    ts["3d-var"].meta = OrderedDict({"CATDESC": "Test 3-Dimensional Variable"})

    return ts


def get_test_timedata():
    ts = TimeSeries(
        time_start="2016-03-22T12:30:31",
        time_delta=3 * u.s,
        data={"Bx": Quantity([1, 2, 3, 4], "gauss", dtype=np.uint16)},
    )
    input_attrs = TimeData.global_attribute_template("eea", "l1", "1.0.0")
    timedata = TimeData(data=ts, meta=input_attrs)
    timedata["Bx"].meta.update({"CATDESC": "Test"})
    return timedata


def get_test_global_meta():
    global_attrs_template = TimeData.global_attribute_template()


def test_non_timeseries():
    with pytest.raises(TypeError):
        _ = TimeData(data=[], meta={})


def test_timedata_empty_ts():
    with pytest.raises(ValueError):
        _ = TimeData(data=TimeSeries())


def test_timedata_bad_ts():
    ts = get_bad_timeseries()
    with pytest.raises(TypeError):
        _ = TimeData(data=ts)


def test_timedata_default():
    ts = get_test_timeseries()

    with pytest.raises(ValueError) as e:
        # We expect this to throw an error that the Instrument is not recognized.
        # The Instrument is one of the attributes required for generating the filename
        # Initialize a CDF File Wrapper
        test_data = TimeData(ts)

    # Test Deleting the Writer
    del ts


def test_timedata_valid_attrs():
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
    with tempfile.TemporaryDirectory() as tmpdirname:
        test_file_output_path = test_data.save(output_path=tmpdirname)
        test_file_cache_path = Path(test_file_output_path)
        # Test the File Exists
        assert test_file_cache_path.exists()


def test_global_attribute_template():
    # default
    assert isinstance(TimeData.global_attribute_template(), OrderedDict)

    # bad instrument
    with pytest.raises(ValueError):
        _ = TimeData.global_attribute_template(instr_name="test instrument")

    # bad Data Level
    with pytest.raises(ValueError):
        _ = TimeData.global_attribute_template(data_level="data level")

    # bad version
    with pytest.raises(ValueError):
        _ = TimeData.global_attribute_template(version="000")

    # good inputs
    template = TimeData.global_attribute_template(
        instr_name="eea", data_level="l2", version="1.3.6"
    )
    assert template["Descriptor"] == "EEA>Electron Electrostatic Analyzer"
    assert template["Data_level"] == "L2>Level 2"
    assert template["Data_version"] == "1.3.6"


def test_default_properties():
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

    # data
    assert isinstance(test_data.data, TimeSeries)

    # units
    assert isinstance(test_data.units, OrderedDict)

    # columns
    assert isinstance(test_data.columns, list)

    # time
    assert isinstance(test_data.time, Time)

    # time_range
    assert isinstance(test_data.time_range, tuple)

    # shape
    assert isinstance(test_data.shape, tuple)

    # __repr__
    assert isinstance(test_data.__repr__(), str)

    # __len__
    assert len(test_data) == test_data.shape[-1]

    # __contains__
    assert "time" in test_data
    with pytest.raises(KeyError):
        _ = test_data["test"]


def test_timedata_single_measurement():
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
    test_data["test_var1"].meta.update(
        {"test_attr1": "test_value1", "CATDESC": "Test data"}
    )

    # Convert the Wrapper to a CDF File
    with tempfile.TemporaryDirectory() as tmpdirname:
        test_file_output_path = test_data.save(output_path=tmpdirname)

        test_file_cache_path = Path(test_file_output_path)
        # Test the File Exists
        assert test_file_cache_path.exists()


def test_timedata_add_measurement():
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
    data_len = len(test_data)

    # Add non-Quantity
    with pytest.raises(TypeError):
        test_data.add_measurement(measure_name="test", data=[], meta={})

    # Good measurement
    q = Quantity(value=random(size=(10)), unit="s", dtype=np.uint16)
    q.meta = OrderedDict({"CATDESC": "Test Variable"})
    test_data.add_measurement(measure_name="test", data=q)
    assert len(test_data) == 1 + data_len

    # test remove_measurement
    test_data.remove_measurement("test")
    assert len(test_data) == data_len


def test_timedata_plot():
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
    q = Quantity(value=random(size=(10)), unit="m", dtype=np.uint16)
    q.meta = OrderedDict({"CATDESC": "Test Variable"})
    test_data.add_measurement(measure_name="test", data=q)

    # Plot All Columns
    ax = test_data.plot(subplots=True)
    assert isinstance(ax, np.ndarray)
    ax = test_data.plot(subplots=False)
    assert isinstance(ax, Axes)
    # Plot Single Column
    ax = test_data.plot(columns=["test"], subplots=True)
    assert isinstance(ax, Axes)
    ax = test_data.plot(columns=["test"], subplots=False)
    assert isinstance(ax, Axes)


def test_timedata_append():
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

    # Append Non-TimeSeries
    with pytest.raises(TypeError):
        test_data.append([])

    # Append Not-Enough Columns
    ts = TimeSeries()
    time = np.arange(start=10, stop=20)
    time_col = Time(time, format="unix").to_datetime()
    col = Column(data=time_col, name="time", meta={})
    ts.add_column(col)
    with pytest.raises(ValueError):
        test_data.append(ts)

    # Append Too-Many Columns
    ts = TimeSeries()
    time = np.arange(start=10, stop=20)
    time_col = Time(time, format="unix").to_datetime()
    col = Column(data=time_col, name="time", meta={})
    ts.add_column(col)
    ts["test1"] = Quantity(value=random(size=(10)), unit="m", dtype=np.uint16)
    ts["test2"] = Quantity(value=random(size=(10)), unit="m", dtype=np.uint16)
    with pytest.raises(ValueError):
        test_data.append(ts)

    # Append Good
    ts = TimeSeries()
    time = np.arange(start=10, stop=20)
    ts["time"] = Time(time, format="unix")
    ts["measurement"] = Quantity(value=random(size=(10)), unit="m", dtype=np.uint16)
    ts["2d-var"] = Quantity(value=random(size=(10, 2)), unit="s", dtype=np.uint16)
    ts["3d-var"] = Quantity(value=random(size=(10, 2, 2)), unit="s", dtype=np.uint16)
    test_data.append(ts)
    assert test_data.shape[0] == 20


def test_timedata_generate_valid_cdf():
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
            data=Quantity(value=random(size=(10)), unit="km"),
            meta={
                "CATDESC": "Test Data",
            },
        )

    # Add 'support_data' VAR_TYPE Attributes
    num_random_vars = 2
    for i in range(num_random_vars):
        # Add Measurement
        test_data.add_measurement(
            measure_name=f"test_support{i}",
            data=Quantity(value=random(size=(10)), unit="km"),
            meta={
                "VAR_TYPE": "support_data",
                "CATDESC": "Test Support",
            },
        )

    # Add 'metadata' VAR_TYPE Attributes
    num_random_vars = 2
    for i in range(num_random_vars):
        # Add Measurement
        test_data.add_measurement(
            measure_name=f"test_metadata{i}",
            data=Quantity(value=random(size=(10)), unit="km"),
            meta={
                "VAR_TYPE": "metadata",
                "CATDESC": "Test Metadata",
            },
        )

    # Convert the Wrapper to a CDF File
    with tempfile.TemporaryDirectory() as tmpdirname:
        # print out rights
        test_file_output_path = test_data.save(output_path=tmpdirname, overwrite=True)

        # Validate the generated CDF File
        result = validate(filepath=test_file_output_path)
        assert len(result) <= 1  # TODO Logical Source and File ID Do not Agree

        # Remove the File
        test_file_cache_path = Path(test_file_output_path)
        test_file_cache_path.unlink()


def test_timedata_from_cdf():
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
            data=Quantity(value=random(size=(10)), unit="km"),
            meta={
                "CATDESC": "Test Data",
            },
        )

    # Add 'support_data' VAR_TYPE Attributes
    num_random_vars = 2
    for i in range(num_random_vars):
        # Add Measurement
        test_data.add_measurement(
            measure_name=f"test_support{i}",
            data=Quantity(value=random(size=(10)), unit="km"),
            meta={
                "VAR_TYPE": "support_data",
                "CATDESC": "Test Support",
            },
        )

    # Add 'metadata' VAR_TYPE Attributes
    num_random_vars = 2
    for i in range(num_random_vars):
        # Add Measurement
        test_data.add_measurement(
            measure_name=f"test_metadata{i}",
            data=Quantity(value=random(size=(10)), unit="km"),
            meta={
                "VAR_TYPE": "metadata",
                "CATDESC": "Test Metadata",
            },
        )

    # Convert the Wrapper to a CDF File
    with tempfile.TemporaryDirectory() as tmpdirname:
        test_file_output_path = test_data.save(output_path=tmpdirname)

        # Validate the generated CDF File
        result = validate(test_file_output_path)
        assert len(result) <= 1  # TODO Logical Source and File ID Do not Agree

        # Try to Load the CDF File in a new CDFWriter
        new_writer = TimeData.load(test_file_output_path)

        # Remove the Original File
        test_file_cache_path = Path(test_file_output_path)
        test_file_cache_path.unlink()

        test_file_output_path2 = new_writer.save(output_path=tmpdirname)
        assert test_file_output_path == test_file_output_path2

        # Validate the generated CDF File
        result2 = validate(test_file_output_path2)
        assert len(result2) <= 1  # TODO Logical Source and File ID Do not Agree
        assert len(result) == len(result2)


@pytest.mark.parametrize(
    "bitlength",
    [
        np.uint8,
        np.uint16,
        np.uint32,
        np.uint64,
        np.int8,
        np.int16,
        np.int32,
        np.int64,
        np.float16,
        np.float32,
    ],
)
def test_bitlength_save_cdf(bitlength):
    """Check that it is possible to create a CDF file for all measurement bitlengths"""
    ts = TimeSeries(
        time_start="2016-03-22T12:30:31",
        time_delta=3 * u.s,
        data={"Bx": Quantity([1, 2, 3, 4], "gauss", dtype=bitlength)},
    )
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

    timedata = TimeData(data=ts, meta=input_attrs)
    timedata["Bx"].meta.update({"CATDESC": "Test"})
    with tempfile.TemporaryDirectory() as tmpdirname:
        test_file_output_path = timedata.save(output_path=tmpdirname)

        test_file_cache_path = Path(test_file_output_path)
        # Test the File Exists
        assert test_file_cache_path.exists()


def test_overwrite_save():
    """Test that when overwrite is set on save no error is generated when trying to create the same file twice"""
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
    td = get_test_timedata()
    td.meta.update(input_attrs)
    with tempfile.TemporaryDirectory() as tmpdirname:
        test_file_output_path = Path(td.save(output_path=tmpdirname))
        # Test the File Exists
        assert test_file_output_path.exists()
        # without overwrite set trying to create the file again should lead to an error
        with pytest.raises(CDFError):
            test_file_output_path = td.save(output_path=tmpdirname, overwrite=False)

        # with overwrite set there should be no error
        assert Path(td.save(output_path=tmpdirname, overwrite=True)).exists()
