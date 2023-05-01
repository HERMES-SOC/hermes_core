.. _cdf_writer:

**************************************
CDF File Writing and Interoperability
**************************************

Overview
========

The CDFWriter class provides a convenient and efficient way to work with CDF files by abstracting them as 
``astropy.TimeSeries`` tables with metadata attributes. This abstraction simplifies data management, enhances 
data discovery, and facilitates adherence to ISTP compliance standards when uploading CDF files to SPDF. By 
leveraging the CDFWriter class, scientists and researchers can focus more on their data analysis and 
scientific exploration.

CDF (Common Data Format) files are a binary file format commonly used in scientific research to store and exchange data. 
They provide a flexible structure for organizing and representing multidimensional datasets along with associated metadata.
CDF files are widely used in space physics, astrophysics, and other scientific disciplines. Despite their versatility, CDF files have certain limitations. For example, working directly with CDF files can be complex
and require knowledge of the CDF specifications. Additionally, managing metadata and accessing specific data elements
can be cumbersome without proper abstractions.

To overcome the limitations of working with CDF files, the CDFWriter class facilitates the abstraction
of CDF files as ``astropy.TimeSeries`` tables with metadata attributes. ``astropy.TimeSeries`` is a Python 
package for handling time series data that provides a convenient and familiar interface for working with 
tabular data. By representing a CDF file as an ``astropy.TimeSeries`` table, it becomes easier to manipulate, 
analyze, and visualize the data. Additionally, metadata attributes can be associated with the table, allowing 
for enhanced documentation and data discovery. It provides a simplified interface to write data and metadata 
to CDF files while automatically handling the complexities of the underlying CDF file format. The CDFWriter 
class encapsulates the necessary logic to convert the ``astropy.TimeSeries`` table into the appropriate CDF 
format, ensuring compatibility and data integrity.

ISTP (International Solar-Terrestrial Physics) compliance is a set of standards defined by the Space Physics 
Data Facility (SPDF) for submitting and sharing space physics data. ISTP compliance ensures that the data adheres 
to specific formatting requirements, quality control measures, and documentation standards. Uploading CDF files 
to SPDF typically requires conforming to the ISTP guidelines. The CDFWriter class plays a crucial role in ensuring 
ISTP compliance when uploading CDF files to SPDF. It handles the conversion of the data and metadata into the required 
format specified by ISTP guidelines, simplifying the process of generating compliant CDF files ready for upload to SPDF.


Creating a ``CDFWriter`` from a ``TimeSeries``
===============================================

A ``CDFWriter`` must be created from a ``TimeSeries`` object containing at least two columns:
  1. `time` 
  2. at least one measurement

Import the required ``astropy`` dependencies and create an empty ``TimeSeries``::

    >>> from `astropy.TimeSeries` import TimeSeries
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

You can now pass the ``TimeSeries`` into a ``CDFWriter``

    >>> frrom hermes_core.util.cdf import CDFWriter
    >>> writer = CDFWriter(ts=ts)

The ``CDFWriter`` can the be updated, measurements added, metadata added, and written to a new CDF file. 


Creating a ``CDFWriter`` from an existing CDF File
===================================================

Given a current CDF File you can create a ``CDFWriter`` data container through passing a path to the CDF file::

    >>> from hermes_core.util.cdf import CDFWriter
    >>> writer = CDFWriter.from_cdf("hermes_eea_default_ql_19700101_v0.0.1.cdf")

The ``CDFWriter`` can the be updated, measurements added, metadata added, and written to a new CDF file. 


Adding data to a ``CDFWriter``
===============================

Data can be added to the ``CDFWriter`` using Astropy ``Quantity`` objects::

    >>> writer["variable_name"] = Quantity(value=my_mumpy_data, unit="your_data_units")
    >>> wtiter

This will now display the ``CDFWriter`` object and a new column with name "variable_name".


Adding metadata attributes to a ``CDFWriter``
==============================================

CDF file global metadata and variable metadata can be added through the ``CDFWriter`` 

Variable Metadata can be updated for a ``CDFWriter`` variable using its ``.meta`` property 
which is an ``OrderdDict`` containing relevant information::

    >>> writer["variable_name"].meta.update({
      "VAR_TYPE": "data",
      "CATDESC": "Test Variable","
    })


Global Metadata cen be updated for a ```CDFWriter``` object using the object's ``.meta`` parameter
which is an ``OrderdDict`` containing relevant information::

    >>> writer.meta.update({
        "Descriptor": "EEA>Electron Electrostatic Analyzer",
        "Data_version": "0.0.1",
        "Data_level": "l1>Level 1",
    })


Using ``CDFWriter`` to Write a CDF File
========================================

The ``CDFWriter`` uses the ``spacepy.pycdf`` module to conver all variable data and metadata to 
a CDF format. Data cen be written to a CDF file using the ``write_cdf(...)`` method and by passing, 
as a parameter, a path to the folder where the CDF file should be saved. 

For example::

    >>> output_path = "./"
    >>> cdf_file_path = writer.write_cdf(output_path)

This returns the full path to the CDF file that was generated. From this you can validate and 
distribute your data as a CDF file.


Using ``CDFWriter`` to Validate a CDF File
===========================================

The ``CDFWriter`` uses the ``spacepy.pycdf.istp`` module for data validation, in addition to custom
tests for additional metadata. A CDF file can be validated using the ``validate_cdf(...)`` method
and by passong, as a parameter, a full path to the CDF file to be validated::

    >>> validation_errors = writer.validate_cdf(cdf_file_path)

This returns a ``list[str]`` that contains any vlidation errors that were enountered when examining
the CDF file. If no validation errors were found the method will return an empty list. 