.. _cdf_writer:

*******************************
Opening and Writing HERMES Data
*******************************

Overview
========

The :py:class:`~hermes_core.timedata.TimeData` class provides a convenient and efficient way to work with HERMES science CDF data files.
The point of this class is to simplify data management, enhances data discovery, and facilitates adherence to CDF standards.

`CDF (Common Data Format) <https://cdf.gsfc.nasa.gov>`_ files are a binary file format commonly used by NASA scientific research to store and exchange data. They provide a flexible structure for organizing and representing multidimensional datasets along with associated metadata. CDF files are widely used in space physics. Because of their versatility, CDF files can be complex.
CDF standards exist to make it easier to work with these files.
`International Solar-Terrestrial Physics (ISTP) <https://spdf.gsfc.nasa.gov/istp_guide/vattributes.html#VAR_TYPE>`_ compliance is a set of standards defined by the Space Physics Data Facility (SPDF).
ISTP compliance ensures that the data adheres to specific formatting requirements, quality control measures, and documentation standards.
Uploading CDF files to the `NASA SPDF archive <https://spdf.gsfc.nasa.gov>`_ requires conforming to the ISTP guidelines.
In addition, HERMES maintains it's own standards in the CDF guide.

To make it easier to work with HERMES data, the :py:class:`~hermes_core.timedata.TimeData` class facilitates the abstraction of HERMES CDF files.
It allows users to read and write HERMES data and is compliant with `PyHC expectations <https://heliopython.org>`_.
The data is stored in a `~astropy.timeseries.TimeSeries` table while the metadata is stored in dictionaries.
`~astropy.timeseries.TimeSeries` is a Python class for handling scientific time series data that provides a convenient and familiar interface for working with tabular data.
By loading the contents of a CDF file into a `~astropy.timeseries.TimeSeries` table, it becomes easier to manipulate, analyze, and visualize the data.
Additionally, metadata attributes can be associated with the table, allowing for enhanced documentation and data discovery.
The :py:class:`~hermes_core.timedata.TimeData` class aims to provide a simplified interface to reading and writing HERMES data and metadata to CDF files while automatically handling the complexities of the underlying CDF file format.

Creating a ``TimeData`` object
==============================

A :py:class:`~hermes_core.timedata.TimeData` must be initialized by providing a `~astropy.timeseries.TimeSeries` object with at least one measurement.
There are many ways to initialize one but here is one example:

    >>> import astropy.units as u
    >>> from astropy.timeseries import TimeSeries
    >>> ts = TimeSeries(time_start='2016-03-22T12:30:31',
                        time_delta=3 * u.s,
                        data={'Bx': u.Quantity([1, 2, 3, 4], 'nanoTesla', dtype=np.uint16)})

Be mindful to set the right number of bits per measurement, in this case 16 bits.
If you do not, it will likely default to float64 and if you write a CDF file it will be larger than expected or needed.
The valid `~numpy.dtype` choices are uint8, uint16, uint32, uint64, int8, int16, int32, int64, float16, float32, float64, float164.
You can also create your time array directly

    >>> from astropy.time import Time
    >>> import astropy.units as u
    >>> from astropy.timeseries import TimeSeries
    >>> times = Time('2010-01-01 00:00:00', scale='utc') + TimeDelta(np.arange(100) * u.s)
    >>> ts = TimeSeries(time=times, data={'diff_e_flux': u.Quantity(np.arange(100) * 1e-3, '1/(cm**2 * s * eV * steradian)'), , dtype=np.float32)})

Note the use of `~astropy.time` and `astropy.units` which provide several advantages over using arrays of numbers and are required by :py:class:`~hermes_core.timedata.TimeData`.

Next create a `dict` or `~collections.OrderedDict` containing the required CDF global metadata.
The class function :py:func:`~hermes_core.timedata.TimeData.global_attribute_template`:: will provide you an empty version that you can fill in if needed.
Here is an example with filled in values.

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

You can now create the :py:class:`~hermes_core.timedata.TimeData` object,

    >>> from hermes_core.timedata import TimeData
    >>> timedata = TimeData(data=ts, meta=input_attrs)

The :py:class:`~hermes_core.timedata.TimeData` is mutable so you can edit it, add another measurement column or edit the metadata after the fact.
Your variable metadata can be found by querying the measurement column directly.

    >>> timedata['Bx'].meta

The class does its best to fill in metadata fields if it can and leaves others blank that it cannot.
Those should be filled in manually.
Be careful when editing metadata that was automatically generated as you might make the resulting CDF file non-compliant.

Putting it all together here is complete example

    >>> ts = TimeSeries(
        time_start="2016-03-22T12:30:31",
        time_delta=3 * u.s,
        data={"Bx": Quantity([1, 2, 3, 4], "gauss", dtype=bitlength)},
    >>> )
    >>> input_attrs = TimeData.global_attribute_template("eea", "l1", "1.0.0")
    >>> timedata = TimeData(data=ts, meta=input_attrs)
    >>> timedata['Bx'].meta.update({"CATDESC": "X component of the Magnetic field measured by HERMES"})

Creating a ``TimeData`` from an existing CDF File
=================================================

Given a current CDF File you can load it into a :py:class:`~hermes_core.timedata.TimeData` by providing a path to the CDF file::

    >>> from hermes_core.timedata import TimeData
    >>> timedata = TimeData.load("hermes_eea_default_ql_19700101_v0.0.1.cdf") # doctest: +SKIP

The :py:class:`~hermes_core.timedata.TimeData` can the be updated, measurements added, metadata added, and written to a new CDF file.

Adding data to a ``TimeData`` Container
=======================================

A new column of data can be added to an existing instance.
Remember that these new measurements must have the same time stamps as the existing ones and therefore the same number of measurements.
You can add the new column in one of two ways.
The more explicit approach is to use :py:func:`~hermes_core.timedata.TimeData.add_measurement` function::

    >>> timedata.add_measurement(
    ...     measure_name=f"By",
    ...     measure_data=u.Quantity(np.arange(len(timedata['Bx'])), 'Gauss', dtype=np.uint16),
    ...     measure_meta={
    ...         "CATDESC": "Test Metadata",
    ...     }
    ... )

Or you can just add the column directly.

    >>> timedata["By"] = u.Quantity(np.arange(len(timedata['Bx'])), 'Gauss', dtype=np.uint16)

Remember that you'll then have to fill in the meta data afterwards.

    >>> timedata['By'].meta.update(measure_meta)

Visualizing data in a ``TimeData`` Container
============================================
The :py:class:`~hermes_core.timedata.TimeData` provides a quick way to visualize its data through `~hermes_core.timedata.TimeData.plot`.


Writing a CDF File
==================

The :py:class:`~hermes_core.timedata.TimeData` class writes CDF files using the `~spacepy.pycdf` module.
This can be done using the :py:func:`~hermes_core.timedata.TimeData.save` method which only requires a path to the folder where the CDF file should be saved.
The filename is automatically derived consistent with HERMES filenaming requirements.
If no path is provided it writes the file to the current directory.
This function returns the full path to the CDF file that was generated.
From this you can validate and distribute your CDF file.

Validating a CDF File
=====================

The :py:class:`~hermes_core.timedata.TimeData` uses the `~spacepy.pycdf.istp` module for CDF validation, in addition to custom
tests for additional metadata. A CDF file can be validated using the :py:func:`~hermes_core.util.validation.validate` method
and by passing, as a parameter, the full path to the CDF file to be validated::

    >>> from hermes_core.util.validation import validate
    >>> validation_errors = validate(cdf_file_path) # doctest: +SKIP

This returns a `list[str]` that contains any validation errors that were encountered when examining the CDF file.
If no validation errors were found the method will return an empty list.

Adding metadata attributes
==========================

Additional CDF file global metadata and variable metadata can be easily added to a :py:class:`~hermes_core.timedata.TimeData` data container.
A list of ISTP-compliant attributes can be found `here <https://spdf.gsfc.nasa.gov/istp_guide/vattributes.html#VAR_TYPE>`_.

**Global Metadata** can be updated for a :py:class:`~hermes_core.timedata.TimeData` object using the object's :py:attr:`~hermes_core.timedata.TimeData.meta` parameter
which is an `~collections.OrderedDict` containing relevant information. The :py:class:`~hermes_core.timedata.TimeData` derives most global
metadata required for ISTP compliance. However, there are a few pieces of metadata that must be supplied by users to successfuly generate ISTP-compliant CDF files.
A template of the required metadata can be obtained using the :py:func:`~hermes_core.timedata.TimeData.global_attribute_template` function::

    >>> global_attrs_template = TimeData.global_attribute_template('eea', 'l1', '0.1.0')

This can make the definition of global metadata easier since instrument teams or users only need to supply pieces of metadata that are in this template.
Additional metadata items can be added if desired.

**Variable Metadata** can be updated for a :py:class:`~hermes_core.timedata.TimeData` variable using its :py:attr:`~hermes_core.timedata.TimeData.meta` property which is an `~collections.OrderedDict`.
The :py:class:`~hermes_core.timedata.TimeData` derives most variable metadata required for ISTP compliance.
However, there are a few pieces of metadata that must be supplied by users to generate ISTP-compliant CDF files:

* `CATDESC` : (Catalogue Description) This is a human readable description of the data variable.

A template of the required metadata can be obtained using the :py:func:`~hermes_core.timedata.TimeData.measurement_attribute_template` function::

    >>> variable_attrs_template = TimeData.measurement_attribute_template()
    >>> variable_attrs_template
    OrderedDict([('CATDESC', None),
                ('DISPLAY_TYPE', None),
                ('LABLAXIS', None),
                ('VAR_TYPE', None)])

If you use the :py:func:`~hermes_core.timedata.TimeData.add_measurement` function, it will automatically fill most of them in for you.
Additional pieces of metadata can be added if desired.

