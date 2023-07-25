"""
This module provides an example for creating a CDF File using the `~hermes_core.timedata.TimeData`
class. This class is an abstraction of underlying data structures to make the handling of
measurement data easier when reading and writinf CDF data.
"""
from typing import OrderedDict
import numpy as np
import astropy.units as u
from astropy.timeseries import TimeSeries

# Import the `hermes_core` Package
from hermes_core.timedata import TimeData
from hermes_core.util.validation import validate


"""
# Creating a TimeData Container

A `TimeData` container must be created using:

- An `~astropy.timeseries.TimeSeries` containing at least two columns
    1. `time`
    2. At least one measurement column as an `~astropy.units.Quantity` object. 
- A global metadata dictionary that contains attribute information of the CDF to be created.

The following lines of code create these two pieces if required information.

"""

# Create a np.ndarray of example measurement data
bx1 = np.random.choice(a=[-1, 0, 1], size=1000).cumsum(0)
by1 = np.random.choice(a=[-1, 0, 1], size=1000).cumsum(0)
bz1 = np.random.choice(a=[-1, 0, 1], size=1000).cumsum(0)

# Create a TimeSeries with the example measurement and a Time column
ts = TimeSeries(
    time_start="2016-03-22T12:30:31",
    time_delta=3 * u.s,
    data={"Bx1": u.Quantity(value=bx1, unit="nanoTesla", dtype=np.int16)},
)

# You can also add new measurements to the TimeSeries directly
ts.add_column(col=u.Quantity(value=by1, unit="nanoTesla", dtype=np.int16), name="By1")
ts.add_column(col=u.Quantity(value=bz1, unit="nanoTesla", dtype=np.int16), name="Bz1")

# To make the creation of global metadata easier you can use the static
# `TimeData.global_attribute_template()` function.
global_attrs_template = TimeData.global_attribute_template()

"""
This gives an `OrderedDict` of empty attributes to be filled out

```
OrderedDict([('DOI', None),
             ('Data_level', None),
             ('Data_version', None),
             ('Descriptor', None),
             ('HTTP_LINK', None),
             ('Instrument_mode', None),
             ('Instrument_type', None),
             ('LINK_TEXT', None),
             ('LINK_TITLE', None),
             ('MODS', None),
             ('PI_affiliation', None),
             ('PI_name', None),
             ('TEXT', None)])
```

You can then provide values for each of these attributes 
"""

global_attrs_template["DOI"] = "https://doi.org/<PREFIX>/<SUFFIX>"
global_attrs_template["Data_level"] = "L1>Level 1"
global_attrs_template["Data_version"] = "0.0.1"
global_attrs_template[
    "Descriptor"
] = "nemisis>Noise Eliminating Magnetometer Instrument in a Small Integrated System"
global_attrs_template["Instrument_mode"] = "default"
global_attrs_template["Instrument_type"] = "Electric Fields (space)"
global_attrs_template["Data_product_descriptor"] = "odpd"

global_attrs_template["HTTP_LINK"] = [
    "https://spdf.gsfc.nasa.gov/istp_guide/istp_guide.html",
    "https://spdf.gsfc.nasa.gov/istp_guide/gattributes.html",
    "https://spdf.gsfc.nasa.gov/istp_guide/vattributes.html",
]
global_attrs_template["LINK_TEXT"] = ["ISTP Guide", "Global Attrs", "Variable Attrs"]
global_attrs_template["LINK_TITLE"] = ["ISTP Guide", "Global Attrs", "Variable Attrs"]

global_attrs_template["MODS"] = [
    "v0.0.1 - Original version.",
]
global_attrs_template["PI_affiliation"] = "HERMES"
global_attrs_template["PI_name"] = "HERMES SOC"
global_attrs_template["TEXT"] = "Sample CDF File"


"""
With bothe the `astropy.timeSeries.TimeSeries` aw well as the global `meta` defined, 
we can now create a `TimeData` object to manipulate the measuerments and metadata together. 
"""

timedata = TimeData(data=ts, meta=global_attrs_template)


"""
# Modifying a TimeData container

Once the TimeData container has been created you can add metadata and new measurements
"""

# To make the creation of variable metadata easier you can use the static
# `TimeData.measurement_attribute_template()` function.
template = TimeData.measurement_attribute_template()

"""
This gives an `OrderedDict` of empty attributes to be filled out

```
OrderedDict([('CATDESC', None)])
```

You can then provide values for each of these attributes 
"""

# Update the Metadata for each of the Measurements
timedata["Bx1"].meta.update(OrderedDict({"CATDESC": "Magnetic Field 1 X"}))
timedata["By1"].meta.update(OrderedDict({"CATDESC": "Magnetic Field 1 Y"}))
timedata["Bz1"].meta.update(OrderedDict({"CATDESC": "Magnetic Field 1 Z"}))

# You can also add new Measurements to the TimeData container
timedata.add_measurement(
    measure_name="Bx2",
    data=u.Quantity(
        value=np.random.choice(a=[-1, 0, 1], size=1000).cumsum(0), unit="nanoTesla", dtype=np.int16
    ),
    meta={
        "VAR_TYPE": "data",
        "CATDESC": "Magnetic Field 2 X",
    },
)

"""
# Saving the TimeData Container to a CDF File

The `TimeData` class delegetes to the appropriate IO module to save the CDF File. 
The `CDFHandler` uses the `spacepy.pycdf` module to convert all variable data and metadata
to a CDF format. 
"""

output_path = "./"
# Convert  to CDF File
cdf_file_path = timedata.save(output_path=output_path, overwrite=True)
