"""Tests for CDF Files to and from data containers"""

from collections import OrderedDict
from pathlib import Path
import pytest
import numpy as np
from numpy.random import random
import tempfile
from astropy.timeseries import TimeSeries
from astropy.table import Column
from astropy.time import Time
from astropy.units import Quantity
import astropy.units as u
from astropy.nddata import NDData
from astropy.wcs import WCS
from ndcube import NDCube, NDCollection
from spacepy.pycdf import CDF, CDFError
from matplotlib.axes import Axes
from hermes_core.timedata import HermesData
from hermes_core.util.schema import HermesDataSchema
from hermes_core.util.validation import validate


def get_bad_timeseries():
    """
    TimeSeries returned does not contain Quantity Measurements
    """
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
    """
    Function to get test astropy.timeseries.TimeSeries to re-use in other tests
    """
    ts = TimeSeries()

    # Create an astropy.Time object
    time = np.arange(10)
    time_col = Time(time, format="unix")
    ts["time"] = time_col

    # Add Measurement
    quant = Quantity(value=random(size=(10)), unit="m", dtype=np.uint16)
    ts["measurement"] = quant
    ts["measurement"].meta = OrderedDict(
        {
            "VAR_TYPE": "data",
            "CATDESC": "Test Data",
        }
    )
    return ts


def get_test_hermes_data():
    """
    Function to get test hermes_core.timedata.HermesData objects to re-use in other tests
    """

    # Astropy TimeSeries
    ts = TimeSeries(
        time_start="2016-03-22T12:30:31",
        time_delta=3 * u.s,
        data={"Bx": Quantity([1, 2, 3, 4], "gauss", dtype=np.uint16)},
    )

    # Support Data / Non-Time Varying Data
    support = {"support_counts": NDData(data=[1])}

    # Spectra Data
    spectra = NDCollection(
        [
            (
                "test_spectra",
                NDCube(
                    data=random(size=(4, 10)),
                    wcs=WCS(naxis=2),
                    meta={"CATDESC": "Test Spectra Variable"},
                    unit="eV",
                ),
            )
        ]
    )

    # Global Metadata Attributes
    input_attrs = HermesData.global_attribute_template("eea", "l1", "1.0.0")

    # Create HermesData Object
    hermes_data = HermesData(
        timeseries=ts, support=support, spectra=spectra, meta=input_attrs
    )
    hermes_data.timeseries["Bx"].meta.update({"CATDESC": "Test"})
    return hermes_data


def test_non_timeseries():
    """
    Test Asserts the `timeseries` parameter must accept an astropy.timeseries.TimeSeries
    """
    with pytest.raises(TypeError):
        _ = HermesData(timeseries=[], meta={})


def test_hermes_data_empty_ts():
    """
    Test asserts the `timeseries` parameter must accept an astropy.timeseries.TimeSeries
    that includes `Time` and at least one other measurement.
    """
    with pytest.raises(ValueError):
        _ = HermesData(timeseries=TimeSeries())


def test_hermes_data_bad_ts():
    """
    Test asserts that all measurements in `timeseries` member must be of type
    `astropy.units.Quantity`
    """
    ts = get_bad_timeseries()
    with pytest.raises(TypeError):
        _ = HermesData(timeseries=ts)


def test_hermes_data_default():
    """
    Test asserts that a HermesDate object is created with the minimal set of required
    global metadata: Descriptor, Data_level, Data_version.
    """
    ts = get_test_timeseries()

    with pytest.raises(ValueError) as e:
        # We expect this to throw an error that the Instrument is not recognized.
        # The Instrument is one of the attributes required for generating the filename
        # Initialize a CDF File Wrapper
        test_data = HermesData(ts)

    # Test Deleting the Writer
    del ts


def test_multidimensional_timeseries():
    """
    Test asserts that HermesData cannot be created with multi-dimensional data in
    the astropy.timeseries.TimeSeries member
    """
    ts = get_test_timeseries()
    ts["var"] = Quantity(value=random(size=(10, 2)), unit="s", dtype=np.uint16)

    with pytest.raises(ValueError):
        _ = HermesData(ts)


def test_support_data():
    """
    Test asserts support / non-time-varying data is created properly
    """
    # fmt: off
    input_attrs = {
        "Descriptor": "EEA>Electron Electrostatic Analyzer",
        "Data_level": "l1>Level 1",
        "Data_version": "v0.0.1",
    }
    # fmt: on
    ts = get_test_timeseries()

    # Bad Support
    support = {"support_var": [1]}
    with pytest.raises(TypeError):
        _ = HermesData(ts, support=support, meta=input_attrs)

    # Support as NDData
    support = {}
    support["support_nddata"] = NDData(data=[1])

    # Support as Quantity
    support["support_quantity"] = Quantity(value=[1], unit="count")

    # Create HermesData
    test_data = HermesData(ts, support=support, meta=input_attrs)

    assert "support_nddata" in test_data.support
    assert test_data.support["support_nddata"].data[0] == 1
    assert "support_quantity" in test_data.support
    assert test_data.support["support_quantity"].data[0] == 1


def test_spectra_data():
    """
    Test asserts spectra / high-dimensional data is created properly
    through ndcube.NDCollection structures.
    """
    input_attrs = {
        "Descriptor": "EEA>Electron Electrostatic Analyzer",
        "Data_level": "l1>Level 1",
        "Data_version": "v0.0.1",
    }
    # fmt: on
    ts = get_test_timeseries()

    # Bad Spectra
    spectra = random(size=(10, 4))
    with pytest.raises(TypeError):
        _ = HermesData(ts, spectra=spectra, meta=input_attrs)

    # Good Spectra
    spectra = NDCollection(
        [
            (
                "test_spectra",
                NDCube(
                    data=random(size=(10, 10)),
                    wcs=WCS(naxis=2),
                    meta={"CATDESC": "Test Spectra Variable"},
                    unit="eV",
                ),
            )
        ]
    )

    # Create HermesData
    test_data = HermesData(ts, spectra=spectra, meta=input_attrs)

    assert "test_spectra" in test_data.spectra
    assert test_data.spectra["test_spectra"].data.shape == (10, 10)


def test_hermes_data_valid_attrs():
    """
    Test asserts that a minmally-defined HermesData object can be created
    and saved to a CDF file.
    """
    # fmt: off
    input_attrs = {
        "Descriptor": "EEA>Electron Electrostatic Analyzer",
        "Data_level": "l1>Level 1",
        "Data_version": "v0.0.1",
    }
    # fmt: on

    ts = get_test_timeseries()
    # Initialize a CDF File Wrapper
    test_data = HermesData(ts, meta=input_attrs)

    # Convert the Wrapper to a CDF File
    with tempfile.TemporaryDirectory() as tmpdirname:
        test_file_output_path = test_data.save(output_path=tmpdirname)
        test_file_cache_path = Path(test_file_output_path)
        # Test the File Exists
        assert test_file_cache_path.exists()


def test_global_attribute_template():
    """
    Test asserts that the HermesData.global_attribute_template()
    function can be used to create a minimal subset of required metadata
    """
    # default
    assert isinstance(HermesData.global_attribute_template(), OrderedDict)

    # bad instrument
    with pytest.raises(ValueError):
        _ = HermesData.global_attribute_template(instr_name="test instrument")

    # bad Data Level
    with pytest.raises(ValueError):
        _ = HermesData.global_attribute_template(data_level="data level")

    # bad version
    with pytest.raises(ValueError):
        _ = HermesData.global_attribute_template(version="000")

    # good inputs
    template = HermesData.global_attribute_template(
        instr_name="eea", data_level="ql", version="1.3.6"
    )
    assert template["Descriptor"] == "EEA>Electron Electrostatic Analyzer"
    assert template["Data_level"] == "QL>Quicklook"
    assert template["Data_version"] == "1.3.6"


def test_default_properties():
    """
    Test asserts the values of HermesData class attributes.
    """
    # Initialize a CDF File Wrapper
    test_data = get_test_hermes_data()

    # data
    assert isinstance(test_data.timeseries, TimeSeries)

    # support
    assert isinstance(test_data.support, dict)

    # data
    assert isinstance(test_data.data, dict)
    assert "timeseries" in test_data.data
    assert isinstance(test_data.data["timeseries"], TimeSeries)
    assert "support" in test_data.data
    assert isinstance(test_data.data["support"], dict)
    assert "spectra" in test_data.data
    assert isinstance(test_data.data["spectra"], NDCollection)

    # meta
    assert isinstance(test_data.meta, dict)

    # time
    assert isinstance(test_data.time, Time)

    # time_range
    assert isinstance(test_data.time_range, tuple)

    # __repr__
    assert isinstance(test_data.__repr__(), str)


def test_hermes_data_single_measurement():
    """
    Test assers that a HermesData object with a single added measurement
    can be created and saved to a CDF file.
    """
    # fmt: off
    input_attrs = {
        "Descriptor": "EEA>Electron Electrostatic Analyzer",
        "Data_level": "l1>Level 1",
        "Data_version": "v0.0.1",
    }
    # fmt: on

    ts = get_test_timeseries()
    # Initialize a CDF File Wrapper
    test_data = HermesData(ts, meta=input_attrs)

    # Add Measurement
    test_data.add_measurement("test_var1", Quantity(value=random(size=(10)), unit="km"))
    test_data.timeseries["test_var1"].meta.update(
        {"test_attr1": "test_value1", "CATDESC": "Test data"}
    )

    # Convert the Wrapper to a CDF File
    with tempfile.TemporaryDirectory() as tmpdirname:
        test_file_output_path = test_data.save(output_path=tmpdirname)

        test_file_cache_path = Path(test_file_output_path)
        # Test the File Exists
        assert test_file_cache_path.exists()


def test_hermes_data_add_measurement():
    """
    Asserts the HermesData.add_measurement() function adds data to the timeseries
    member as expected.
    """
    # fmt: off
    input_attrs = {
        "Descriptor": "EEA>Electron Electrostatic Analyzer",
        "Data_level": "l1>Level 1",
        "Data_version": "v0.0.1",
    }
    # fmt: on

    ts = get_test_timeseries()
    # Initialize a CDF File Wrapper
    test_data = HermesData(ts, meta=input_attrs)
    ts_len = len(test_data.timeseries.columns)

    # Add non-Quantity
    with pytest.raises(TypeError):
        test_data.add_measurement(measure_name="test", data=[], meta={})

    # Add multi-domensional data
    with pytest.raises(ValueError):
        q = Quantity(value=random(size=(10, 10, 10)), unit="s", dtype=np.uint16)
        test_data.add_measurement(measure_name="test", data=q, meta={})

    # Good measurement
    q = Quantity(value=random(size=(10)), unit="s", dtype=np.uint16)
    q.meta = OrderedDict({"CATDESC": "Test Variable"})
    test_data.add_measurement(measure_name="test", data=q)
    assert test_data.timeseries["test"].shape == (10,)

    # Add Dimensionless measurements (or Record-Varying) Data
    q = Quantity(value=random(size=(10)), unit=u.dimensionless_unscaled)
    test_data.add_measurement(
        measure_name="Test Dimensionless",
        data=q,
        meta={"CATDESC": "Test Dimensionless Data", "VAR_TYPE": "support_data"},
    )
    assert test_data.timeseries["Test Dimensionless"].shape == (10,)

    # Add Count-Based Record-Varying Data
    q = Quantity(value=random(size=(10)), unit=u.count)
    test_data.add_measurement(
        measure_name="Test Count",
        data=q,
        meta={"CATDESC": "Test Count Data", "VAR_TYPE": "support_data"},
    )
    assert test_data.timeseries["Test Count"].shape == (10,)

    # test remove_measurement
    test_data.remove("test")
    assert "test" not in test_data.timeseries.columns

    # Test non-existent variable
    with pytest.raises(ValueError):
        test_data.remove("bad_variable")


def test_hermes_data_add_support():
    """Function to Test Adding Support/ Non-Record-Varying Data"""
    # fmt: off
    input_attrs = {
        "Descriptor": "EEA>Electron Electrostatic Analyzer",
        "Data_level": "l1>Level 1",
        "Data_version": "v0.0.1",
    }
    # fmt: on

    ts = get_test_timeseries()
    # Initialize a CDF File Wrapper
    test_data = HermesData(ts, meta=input_attrs)

    # Add non-Quantity
    with pytest.raises(TypeError):
        test_data.add_support(name="test", data=[], meta={})

    # Add Test Metadata as NDData
    c = NDData(data=[1])
    test_data.add_support(
        name="Test Metadata",
        data=c,
        meta={"CATDESC": "Test Metadata Variable", "VAR_TYPE": "metadata"},
    )
    assert "Test Metadata" in test_data.support
    assert test_data.support["Test Metadata"].data[0] == 1

    # Add Test Metadata as Quantity
    test_data.add_support(
        name="Test Count",
        data=Quantity(value=[1], unit="count", dtype=np.uint16),
        meta={"CATDESC": "Test Metadata Count", "VAR_TYPE": "metadata"},
    )
    assert "Test Count" in test_data.support
    assert test_data.support["Test Count"].data[0] == 1

    # Test remove Support Data
    test_data.remove("Test Metadata")
    assert "Test Metadata" not in test_data.support


def test_hermes_data_add_spectra():
    """Function to Test Adding Spectra/ High-Dimensional Data"""
    # fmt: off
    input_attrs = {
        "Descriptor": "EEA>Electron Electrostatic Analyzer",
        "Data_level": "l1>Level 1",
        "Data_version": "v0.0.1",
    }
    # fmt: on

    ts = get_test_timeseries()
    # Initialize a CDF File Wrapper
    test_data = HermesData(ts, meta=input_attrs)

    # Add non-NDCube
    with pytest.raises(TypeError):
        test_data.add_spectra(name="test", data=[], meta={})

    # Add Test Data
    data = NDCube(
        data=random(size=(10, 10)),
        wcs=WCS(naxis=2),
        meta={"CATDESC": "Test Spectra Variable"},
        unit="eV",
    )
    test_data.add_spectra(
        name="Test Spectra",
        data=data,
        meta={"VAR_TYPE": "data"},
    )
    assert "Test Spectra" in test_data.spectra
    assert "CATDESC" in test_data.spectra["Test Spectra"].meta
    assert "VAR_TYPE" in test_data.spectra["Test Spectra"].meta
    assert test_data.spectra["Test Spectra"].data.shape == (10, 10)

    # Add a Second NDCube
    data = NDCube(
        data=random(size=(10, 10)),
        wcs=WCS(naxis=2),
        meta={"CATDESC": "Second Spectra Variable"},
        unit="eV",
    )
    test_data.add_spectra(
        name="Test2",
        data=data,
        meta={"VAR_TYPE": "data"},
    )
    assert "Test2" in test_data.spectra
    assert "CATDESC" in test_data.spectra["Test2"].meta
    assert "VAR_TYPE" in test_data.spectra["Test2"].meta
    assert test_data.spectra["Test2"].data.shape == (10, 10)

    # Test remove Support Data
    test_data.remove("Test Spectra")
    assert "Test Spectra" not in test_data.spectra


def test_hermes_data_plot():
    """
    Test asserts the HermesData.plot() function generates matplotlib
    images as expected.
    """
    # fmt: off
    input_attrs = {
        "Descriptor": "EEA>Electron Electrostatic Analyzer",
        "Data_level": "l1>Level 1",
        "Data_version": "v0.0.1",
    }
    # fmt: on

    ts = get_test_timeseries()
    # Initialize a CDF File Wrapper
    test_data = HermesData(ts, meta=input_attrs)
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


def test_hermes_data_append():
    """
    Test asserts the HermesData.append() function adds to the TimeSeries member
    as expected.
    """
    # fmt: off
    input_attrs = {
        "Descriptor": "EEA>Electron Electrostatic Analyzer",
        "Data_level": "l1>Level 1",
        "Data_version": "v0.0.1",
    }
    # fmt: on

    ts = get_test_timeseries()
    # Initialize a CDF File Wrapper
    test_data = HermesData(ts, meta=input_attrs)

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
    test_data.append(ts)
    assert len(test_data.timeseries) == 20


def test_hermes_data_generate_valid_cdf():
    """
    Test asserts the HermesData data container can create an ISTP compliant CDF based on
    the spacepy.pycdf.istp module.
    """
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
    support = {
        "nrv_var": NDData(
            data=[1, 2, 3],
            meta={"CATDESC": "Test Metadata Variable", "VAR_TYPE": "metadata"},
        )
    }
    # Initialize a CDF File Wrapper
    test_data = HermesData(ts, support=support, meta=input_attrs)

    # Add the Time column
    test_data.timeseries["time"].meta.update(
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
        assert len(result) <= 1  # Logical Source and File ID Do not Agree

        # Remove the File
        test_file_cache_path = Path(test_file_output_path)
        test_file_cache_path.unlink()


def test_hermes_data_from_cdf():
    """
    Test asserts that the HermesData class can be created by loading a CDF file.
    """
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
    support = {
        "nrv_var": NDData(
            data=[1, 2, 3],
            meta={"CATDESC": "Test Metadata Variable", "VAR_TYPE": "metadata"},
        )
    }
    # Initialize a CDF File Wrapper
    test_data = HermesData(ts, support=support, meta=input_attrs)

    # Add the Time column
    test_data.timeseries["time"].meta.update(
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
        assert len(result) <= 1  # Logical Source and File ID Do not Agree

        # Try to Load the CDF File in a new CDFWriter
        new_writer = HermesData.load(test_file_output_path)

        # Remove the Original File
        test_file_cache_path = Path(test_file_output_path)
        test_file_cache_path.unlink()

        test_file_output_path2 = new_writer.save(output_path=tmpdirname)
        assert test_file_output_path == test_file_output_path2

        # Validate the generated CDF File
        result2 = validate(test_file_output_path2)
        assert len(result2) <= 1  # Logical Source and File ID Do not Agree
        assert len(result) == len(result2)


def test_hermes_data_idempotency():
    """
    Test asserts that a HermesData object that is saved and loaded does not have any
    changes in it members, measurements, or metadata.
    """
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
    # Generate a base HermesData object
    test_data = get_test_hermes_data()
    test_data.meta.update(input_attrs)

    # Induce a Bad (Null) Global Attribute
    test_data.meta["Test Null Attr"] = ""

    # Induce an Non-Record-Varying Variable
    test_data.add_support(
        name="NRV_var", data=NDData(["Test NRV Data"]), meta={"CATDESC": "NRV Variable"}
    )

    # Induce a Variable with Bad UNITS
    test_data.add_support(
        name="Bad_units_var",
        data=NDData(
            [1, 2, 3, 4],
            meta={
                "UNITS": "Not A Unit",
                "CATDESC": "Test Variable with Incoherent UNITS",
            },
        ),
    )

    with tempfile.TemporaryDirectory() as tmpdirname:
        test_file_output_path = test_data.save(output_path=tmpdirname)

        # Try loading the *Invalid* CDF File
        loaded_data = HermesData.load(test_file_output_path)

        assert len(test_data.timeseries.columns) == len(loaded_data.timeseries.columns)
        assert len(test_data.support) == len(loaded_data.support)
        assert len(test_data.meta) == len(loaded_data.meta)

        for attr in test_data.meta:
            assert attr in loaded_data.meta

        for var in test_data.timeseries.columns:
            assert var in loaded_data.timeseries.columns
            assert len(test_data.timeseries[var]) == len(loaded_data.timeseries[var])
            assert len(test_data.timeseries[var].meta) == len(
                loaded_data.timeseries[var].meta
            )
            assert (
                test_data.timeseries[var].meta["VAR_TYPE"]
                == loaded_data.timeseries[var].meta["VAR_TYPE"]
            )

        for var in test_data.support:
            assert var in loaded_data.support
            assert (
                test_data.support[var].data.shape == loaded_data.support[var].data.shape
            )
            assert len(test_data.support[var].meta) == len(
                loaded_data.support[var].meta
            )
            assert (
                test_data.support[var].meta["VAR_TYPE"]
                == loaded_data.support[var].meta["VAR_TYPE"]
            )

        schema = HermesDataSchema()
        for var in test_data.spectra:
            assert var in loaded_data.spectra
            assert (
                test_data.spectra[var].data.shape == loaded_data.spectra[var].data.shape
            )
            assert len(test_data.spectra[var].meta) == len(
                loaded_data.spectra[var].meta
            )
            for _, prop, _ in schema.wcs_keyword_to_astropy_property:
                pval1 = getattr(test_data.spectra[var].wcs.wcs, prop)
                pval2 = getattr(loaded_data.spectra[var].wcs.wcs, prop)
                assert len(pval1) == len(pval2)
                for i in range(len(pval1)):
                    assert pval1[i] == pval2[i]


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

    hermes_data = HermesData(timeseries=ts, meta=input_attrs)
    hermes_data.timeseries["Bx"].meta.update({"CATDESC": "Test"})
    with tempfile.TemporaryDirectory() as tmpdirname:
        test_file_output_path = hermes_data.save(output_path=tmpdirname)

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
    td = get_test_hermes_data()
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


def test_without_cdf_lib():
    """Function to test HermesData Functions without the use of spacepy.pycdf libraries"""
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
    # Get Test TimeSeries
    ts = get_test_timeseries()

    # Disable CDF Libraries
    import spacepy.pycdf as pycdf

    pycdf.lib = None

    # Initialize a CDF File Wrapper
    test_data = HermesData(ts, meta=input_attrs)

    assert test_data.meta["CDF_Lib_version"] == "unknown version"

    # Reset/ Re-Enable CDF Libraries
    pycdf.lib = pycdf.Library(pycdf._libpath, pycdf._library)
