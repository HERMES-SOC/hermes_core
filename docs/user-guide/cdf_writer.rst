.. _cdf_writer:

**************************************
CDF File Writing and Interoperability
**************************************

Overview
========

The ``TimeData`` class provides a convenient and efficient way to work with science data files, including CDF 
files, by abstracting them as ``astropy.timeseries.TimeSeries`` tables with metadata attributes. This abstraction 
simplifies data management, enhances data discovery, and facilitates adherence to ISTP compliance standards 
when uploading CDF files to SPDF. By leveraging the ``TimeData`` class, scientists and researchers can 
focus more on their data analysis and scientific exploration.

CDF (Common Data Format) files are a binary file format commonly used in scientific research to store and 
exchange data. They provide a flexible structure for organizing and representing multidimensional datasets 
along with associated metadata. CDF files are widely used in space physics, astrophysics, and other 
scientific disciplines. Despite their versatility, CDF files have certain limitations. For example, working 
directly with CDF files can be complex and require knowledge of the CDF specifications. Additionally, 
managing metadata and accessing specific data elements can be cumbersome without proper abstractions.

To overcome the limitations of working with CDF files, the ``TimeData`` class facilitates the abstraction
of CDF files as ``astropy.timeseries.TimeSeries`` tables with metadata attributes. ``TimeSeries`` is a Python 
class for handling time series data that provides a convenient and familiar interface for working with 
tabular data. By representing a CDF file as an ``TimeSeries`` table, it becomes easier to manipulate, 
analyze, and visualize the data. Additionally, metadata attributes can be associated with the table, allowing 
for enhanced documentation and data discovery. It provides a simplified interface to write data and metadata 
to CDF files while automatically handling the complexities of the underlying CDF file format. The ``TimeData`` 
class encapsulates the necessary logic to convert the ``TimeSeries`` table into the appropriate CDF 
format, ensuring compatibility and data integrity.

ISTP (International Solar-Terrestrial Physics) compliance is a set of standards defined by the Space Physics 
Data Facility (SPDF) for submitting and sharing space physics data. ISTP compliance ensures that the data 
adheres to specific formatting requirements, quality control measures, and documentation standards. Uploading 
CDF files to SPDF typically requires conforming to the ISTP guidelines. The ``TimeData`` class plays a crucial 
role in ensuring ISTP compliance when uploading CDF files to SPDF. It handles the conversion of the data and 
metadata into the required format specified by ISTP guidelines, simplifying the process of generating 
compliant CDF files ready for upload to SPDF.


Creating a ``TimeData`` object from a ``TimeSeries``
====================================================

A ``TimeData`` must be created from a ``TimeSeries`` object containing at least two columns:
  1. `time` 
  2. at least one measurement
The `TimeData` container must also take in an input of global metadata for the CDF File to be created.

Import the required ``astropy`` dependencies and create an empty ``TimeSeries``::

    >>> from astropy.timeseries import TimeSeries
    >>> from astropy.time import Time
    >>> from astropy.units.quantity import Quantity
    >>> ts = TimeSeries()

Create a ``Time`` object that containts the time dimension and add it to the ``TimeSeries``::

    >>> import numpy as np
    >>> time_arr = np.arange(50)
    >>> time = Time(time_arr, format="unix")
    >>> ts["time"] = time

Create a ``Quantity`` object that contsins measurement data and add it to the ``TimeSeries``::

    >>> from numpy.random import random
    >>> data = random(size=(50))
    >>> quant = Quantity(value=data, unit="m")
    >>> ts["measurement"] = quant

Create a ``dict`` or ``OrderedDict`` containing global metadata required for generating the ``TimeData``.
To facilitate creating this dict see the notes for ``TimeData.global_attribute_template()`` below::

    >>> input_attrs = {
    ...     "DOI": "https://doi.org/<PREFIX>/<SUFFIX>",
    ...     "Data_level": "L1>Level 1",  # NOT AN ISTP ATTR
    ...     "Data_version": "0.0.1",
    ...     "Descriptor": "EEA>Electron Electrostatic Analyzer",
    ...     "Data_product_descriptor": "odpd",
    ...     "HTTP_LINK": [
    ...         "https://spdf.gsfc.nasa.gov/istp_guide/istp_guide.html",
    ...         "https://spdf.gsfc.nasa.gov/istp_guide/gattributes.html",
    ...         "https://spdf.gsfc.nasa.gov/istp_guide/vattributes.html"
    ...     ],
    ...     "Instrument_mode": "default",  # NOT AN ISTP ATTR
    ...     "Instrument_type": "Electric Fields (space)",
    ...     "LINK_TEXT": [
    ...         "ISTP Guide",
    ...         "Global Attrs",
    ...         "Variable Attrs"
    ...     ],
    ...     "LINK_TITLE": [
    ...         "ISTP Guide",
    ...         "Global Attrs",
    ...         "Variable Attrs"
    ...     ],
    ...     "MODS": [
    ...         "v0.0.0 - Original version.",
    ...         "v1.0.0 - Include trajectory vectors and optics state.",
    ...         "v1.1.0 - Update metadata: counts -> flux.",
    ...         "v1.2.0 - Added flux error.",
    ...         "v1.3.0 - Trajectory vector errors are now deltas."
    ...     ],
    ...     "PI_affiliation": "HERMES",
    ...     "PI_name": "HERMES SOC",
    ...     "TEXT": "Valid Test Case",
    ... }

You can now pass the ``TimeSeries`` and ``dict`` global metadata into a ``TimeData`` object,

    >>> from hermes_core.timedata import TimeData
    >>> timedata = TimeData(data=ts, meta=input_attrs)

The ``TimeData`` can the be updated, measurements added, metadata added, and written to a new CDF file. 


Creating a ``TimeData`` from an existing CDF File
===================================================

Given a current CDF File you can create a ``TimeData`` data container through passing a path to the CDF file::

    >>> from hermes_core.timedata import TimeData
    >>> timedata = TimeData.load("hermes_eea_default_ql_19700101_v0.0.1.cdf") # doctest: +SKIP  

The ``TimeData`` can the be updated, measurements added, metadata added, and written to a new CDF file. 


Adding data to a ``TimeData`` Container
=======================================

Data can be added to the ``TimeData`` using Astropy ``Quantity`` objects::

    >>> data = random(size=(50))
    >>> timedata["variable_name"] = Quantity(value=data, unit="m")

Variable metadata is derived for for the given measuerment automatically when adding a new ``Quantity``.

Measurement data and metadata can be aded together though the ``add_measurement()`` function::

    >>> timedata.add_measurement(
    ...     measure_name=f"test_metadata",
    ...     measure_data=Quantity(value=random(size=(50)), unit="km"),
    ...     measure_meta={
    ...         "VAR_TYPE": "metadata",
    ...         "CATDESC": "Test Metadata",
    ...         "DISPLAY_TYPE": "time_series",
    ...         "LABLAXIS": "Metadata Axis Label",
    ...     }
    ... )


Adding metadata attributes to a ``TimeData`` Container
======================================================

CDF file global metadata and variable metadata can be added through the ``TimeData`` data container. 

**Variable Metadata** can be updated for a ``TimeData`` variable using its ``.meta`` property 
which is an ``OrderdDict`` containing relevant information. The ``TimeData`` derives most
variable metadata required for ISTP compliance. However, there are a few pieces of metadaata
that must be supplied by users to generate ISTP-compliant CDF files:

* `CATDESC` : (Catalogue Description) This is a human readable desctiption of the data variable. 
* `DISPLAY_TYPE` : This tells the automated software, such as CDAWeb how the data should be 
    displayed.
* `LABLAXIS` : Used to label a plot axis or to provide a heading for data listing. 
* `VAR_TYPE` : Used in CDAWeb to indicate if the data should be used directly by users. 

For Example::

    >>> timedata["measurement"].meta.update({
    ...     "VAR_TYPE": "metadata",
    ...     "CATDESC": "Test Measurement",
    ...     "DISPLAY_TYPE": "time_series",
    ...     "LABLAXIS": "Measutement Label",
    ... })

A template of the required metadata can be obtained using the `TimeData.measurement_attribute_template()` function::

    >>> variable_attrs_template = TimeData.measurement_attribute_template()
    >>> variable_attrs_template
    OrderedDict([('CATDESC', None), ('DISPLAY_TYPE', None), ('LABLAXIS', None), ('VAR_TYPE', None)])

This can make the definition of variable metadata easier since instrument teams or users only need to supply
pieces of metadata that are in this template. Additional pieces of metadata can be added if desired.

**Global Metadata** can be updated for a ```TimeData``` object using the object's ``.meta`` parameter
which is an ``OrderdDict`` containing relevant information. The ``TimeData`` derives most global 
metadata required for ISTP compliance. However, there are a few pieces of metadata that must be 
supplied by users to successfuly generate ISTP-compliant CDF files:

* `Descriptor` : This attribute identifies the name of the instrument or sensor that collected
    the data. Both the a long name and a short name are given. For any data file, onle a
    single value is allowed.
* `Data_level` : This attribute identifies the level to which the data has been
    processed. Ex. "ql>Quick Look"
* `Data_version` : This attribute identifies the version (vX.Y.Z) of a particular CDF 
    data file.

For Example::

    >>> input_attrs = {
    ...     "Descriptor": "EEA>Electron Electrostatic Analyzer",
    ...     "Data_level": "l1>Level 1",
    ...     "Data_version": "v0.0.1",
    ... }
    >>> timedata.meta.update(input_attrs)

A template of the required metadata can be obtained using the `TimeData.global_attribute_template()` function::

    >>> global_attrs_template = TimeData.global_attribute_template()
    >>> global_attrs_template
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

This can make the definition of global metadata easier since instrument teams or users only need to supply
pieces of metadata that are in this template. Additional pieces of metadata can be added if desired.

Using ``TimeData`` to Write a CDF File
========================================

The ``TimeData`` class uses the ``spacepy.pycdf`` module to convert all variable data and metadata to 
a CDF format. Data cen be written to a CDF file using the ``save(...)`` method and by passing, 
as a parameter, a path to the folder where the CDF file should be saved. 

For example::

    >>> output_path = "./"
    >>> cdf_file_path = timedata.save(output_path) # doctest: +SKIP 

This returns the full path to the CDF file that was generated. From this you can validate and 
distribute your data as a CDF file.


Using ``TimeData`` to Validate a CDF File
===========================================

The ``TimeData`` uses the ``spacepy.pycdf.istp`` module for data validation, in addition to custom
tests for additional metadata. A CDF file can be validated using the ``validate(...)`` method
and by passing, as a parameter, a full path to the CDF file to be validated::

    >>> from hermes_core.util.validation import validate
    >>> validation_errors = validate(cdf_file_path) # doctest: +SKIP 

This returns a ``list[str]`` that contains any vlidation errors that were enountered when examining
the CDF file. If no validation errors were found the method will return an empty list. 