.. _tutorial1:

Creating a CDF File
===================

This module provides an example for creating a CDF File using the `~hermes_core.timedata.HermesData`
class. This class is an abstraction of underlying data structures to make the handling of
measurement data easier when reading and writing CDF data.

    >>> from collections import OrderedDict
    >>> import numpy as np
    >>> import astropy.units as u
    >>> from astropy.timeseries import TimeSeries
    >>> from astropy.nddata import NDData
    >>> from astropy.wcs import WCS
    >>> from ndcube import NDCube, NDCollection
    >>> import tempfile
    >>> 
    >>> # Import the `hermes_core` Package
    >>> from hermes_core.timedata import HermesData
    >>> from hermes_core.util.validation import validate
    >>> 
    >>> # Create a np.ndarray of example measurement data
    >>> bx = np.random.choice(a=[-1, 0, 1], size=1000).cumsum(0)
    >>> by = np.random.choice(a=[-1, 0, 1], size=1000).cumsum(0)
    >>>
    >>> # Create a TimeSeries with the example measurement and a Time column
    >>> ts = TimeSeries(
    ...     time_start="2016-03-22T12:30:31",
    ...     time_delta=3 * u.s,
    ...     data={"Bx GSE": u.Quantity(value=bx, unit="nanoTesla", dtype=np.int16)},
    ... )
    >>> 
    >>> # You can also add new measurements to the TimeSeries directly
    >>> ts.add_column(col=u.Quantity(value=by, unit="nanoTesla", dtype=np.int16),
    ...     name="By GSE"
    ... )
    >>>
    >>> # Create support data or non-time-varying (time invariant) data
    >>> support_data = {
    ...     "data_mask": NDData(data=np.eye(100, 100, dtype=np.uint16))
    ... }
    >>> 
    >>> # Create high-dimensional data leveraging the API of NDCube
    >>> spectra = NDCollection(
    ...     [
    ...         (
    ...             "example_spectra",
    ...             NDCube(
    ...                 data=np.random.random(size=(4, 10)),
    ...                 wcs=WCS(naxis=2),
    ...                 meta={"CATDESC": "Example Spectra Variable"},
    ...                 unit="eV",
    ...             ),
    ...         )
    ...     ]
    ... )
    >>> 
    >>> # To make the creation of global metadata easier you can use the static
    >>> # `HermesData.global_attribute_template()` function.
    >>> global_attrs_template = HermesData.global_attribute_template()
    >>> 
    >>> global_attrs_template["DOI"] = "https://doi.org/<PREFIX>/<SUFFIX>"
    >>> global_attrs_template["Data_level"] = "L1>Level 2"
    >>> global_attrs_template["Data_version"] = "0.0.1"
    >>> global_attrs_template[
    ...     "Descriptor"
    ... ] = "nemisis>Noise Eliminating Magnetometer Instrument in a Small Integrated System"
    >>> global_attrs_template["Instrument_mode"] = "default"
    >>> global_attrs_template["Instrument_type"] = "Magnetic Fields (space)"
    >>> global_attrs_template["Data_product_descriptor"] = "odpd"
    >>> 
    >>> global_attrs_template["HTTP_LINK"] = [
    ...     "https://science.nasa.gov/missions/hermes",
    ...     "https://github.com/HERMES-SOC",
    ...     "https://github.com/HERMES-SOC/hermes_nemisis",
    ... ]
    >>> global_attrs_template["LINK_TEXT"] = ["HERMES homepage",
    ...     "HERMES SOC Github", "NEMISIS Analysis Tools"]
    >>> global_attrs_template["LINK_TITLE"] = ["HERMES homepage",
    ...     "HERMES SOC Github", "NEMISIS Analysis Tools"]
    >>> 
    >>> global_attrs_template["MODS"] = ["v0.0.1 - Original version."]
    >>> global_attrs_template["PI_affiliation"] = "NASA Goddard Space Flight Center"
    >>> global_attrs_template["PI_name"] = "Dr. Eftyhia Zesta"
    >>> global_attrs_template["TEXT"] = "Sample HERMES NEMISIS CDF File"
    >>> 
    >>> example_data = HermesData(
    ...     timeseries=ts, 
    ...     support=support_data, 
    ...     spectra=spectra, 
    ...     meta=global_attrs_template
    ... )
    >>> 
    >>> # To make the creation of variable metadata easier you can use the static
    >>> # `HermesData.measurement_attribute_template()` function.
    >>> template = HermesData.measurement_attribute_template()
    >>> 
    >>> # Update the Metadata for each of the Measurements
    >>> example_data.timeseries["Bx GSE"].meta.update(
    ...     OrderedDict({"CATDESC": "X component of magnetic Field GSE"}))
    >>> example_data.timeseries["By GSE"].meta.update(
    ...     OrderedDict({"CATDESC": "Y component of magnetic Field GSE"}))
    >>> 
    >>> # You can add new scalar time-variant measurements to the HermesData container
    >>> bz = np.random.choice(a=[-1, 0, 1], size=1000).cumsum(0)
    >>> example_data.add_measurement(
    ...     measure_name="Bz GSE",
    ...     data=u.Quantity(value=bz, unit="nanoTesla", dtype=np.int16),
    ...     meta={
    ...         "VAR_TYPE": "data",
    ...         "CATDESC": "Z component of magnetic Field GSE",
    ...     },
    ... )
    >>> 
    >>> # You can add new time-invariant data to the HermesData container
    >>> example_data.add_support(
    ...     name="calibration_const",
    ...     data=NDData(data=[1e-1]),
    ...     meta={
    ...         "CATDESC": "Calibration Factor", 
    ...         "VAR_TYPE": "metadata"
    ...     },
    ... )
    >>> 
    >>> # You can ass new spectral or high-dimensional data to the HermesData container
    >>> data = NDCube(
    ...     data=np.random.random(size=(1000, 10)),
    ...     wcs=WCS(naxis=2),
    ...     meta={"CATDESC": "Example Spectra Variable"},
    ...     unit="eV",
    ... )
    >>> example_data.add_spectra(
    ...     name="added_spectra",
    ...     data=data,
    ...     meta={"VAR_TYPE": "data"},
    ... )
    >>> 
    >>> # create the CDF File
    >>> DRYRUN=True
    >>> if DRYRUN:
    ...     with tempfile.TemporaryDirectory() as tmpdirname:
    ...         cdf_file_path = example_data.save(output_path=tmpdirname)
    ... else:
    ...     cdf_file_path = example_data.save(output_path="./", overwrite=True)

The file that this code generates is made available as a sample file in this
repository in :file:`hermes_core/data/sample/hermes_nms_default_l1_20160322T123031_v0.0.1.cdf`.