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
from hermes_core.timedata import HermesData

def get_test_timeseries(n=10):
    """
    Function to get test astropy.timeseries.TimeSeries to re-use in other tests
    """
    ts = TimeSeries()

    # Create an astropy.Time object
    time = np.arange(n)
    time_col = Time(time, format="unix")
    ts["time"] = time_col

    # Add Measurement
    quant = Quantity(value=random(size=(n)), unit="m", dtype=np.uint16)
    ts["measurement"] = quant
    ts["measurement"].meta = OrderedDict(
        {
            "VAR_TYPE": "data",
            "CATDESC": "Test Data",
        }
    )
    return ts


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
        tmp_path = Path(tmpdirname)
        test_file_output_path = test_data.save(output_path=tmp_path, overwrite=True)

        # Validate the generated CDF File
        result = HermesData.validate(file_path=test_file_output_path)
        assert len(result) <= 1  # Logical Source and File ID Do not Agree

        # Remove the File
        test_file_output_path.unlink()
