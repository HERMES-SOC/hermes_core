Creating a CDF File
===================

This module provides an example for creating a CDF File using the `~hermes_core.timedata.TimeData`
class. This class is an abstraction of underlying data structures to make the handling of
measurement data easier when reading and writing CDF data.

.. code-block:: python

    from collections import OrderedDict
    import numpy as np
    import astropy.units as u
    from astropy.timeseries import TimeSeries

    # Import the `hermes_core` Package
    from hermes_core.timedata import TimeData
    from hermes_core.util.validation import validate

    # Create a np.ndarray of example measurement data
    bx = np.random.choice(a=[-1, 0, 1], size=1000).cumsum(0)
    by = np.random.choice(a=[-1, 0, 1], size=1000).cumsum(0)

    # Create a TimeSeries with the example measurement and a Time column
    ts = TimeSeries(
        time_start="2016-03-22T12:30:31",
        time_delta=3 * u.s,
        data={"Bx GSE": u.Quantity(value=bx, unit="nanoTesla", dtype=np.int16)},
    )

    # You can also add new measurements to the TimeSeries directly
    ts.add_column(col=u.Quantity(value=by, unit="nanoTesla", dtype=np.int16),
        name="By GSE")

    # To make the creation of global metadata easier you can use the static
    # `TimeData.global_attribute_template()` function.
    global_attrs_template = TimeData.global_attribute_template()

    global_attrs_template["DOI"] = "https://doi.org/<PREFIX>/<SUFFIX>"
    global_attrs_template["Data_level"] = "L1>Level 2"
    global_attrs_template["Data_version"] = "0.0.1"
    global_attrs_template[
        "Descriptor"
    ] = "nemisis>Noise Eliminating Magnetometer Instrument in a Small Integrated System"
    global_attrs_template["Instrument_mode"] = "default"
    global_attrs_template["Instrument_type"] = "Magnetic Fields (space)"
    global_attrs_template["Data_product_descriptor"] = "odpd"

    global_attrs_template["HTTP_LINK"] = [
        "https://science.nasa.gov/missions/hermes",
        "https://github.com/HERMES-SOC",
        "https://github.com/HERMES-SOC/hermes_nemisis",
    ]
    global_attrs_template["LINK_TEXT"] = ["HERMES homepage",
        "HERMES SOC Github", "NEMISIS Analysis Tools"]
    global_attrs_template["LINK_TITLE"] = ["HERMES homepage",
        "HERMES SOC Github", "NEMISIS Analysis Tools"]

    global_attrs_template["MODS"] = ["v0.0.1 - Original version."]
    global_attrs_template["PI_affiliation"] = "NASA Goddard Space Flight Center"
    global_attrs_template["PI_name"] = "Dr. Eftyhia Zesta"
    global_attrs_template["TEXT"] = "Sample HERMES NEMISIS CDF File"

    timedata = TimeData(data=ts, meta=global_attrs_template)

    # To make the creation of variable metadata easier you can use the static
    # `TimeData.measurement_attribute_template()` function.
    template = TimeData.measurement_attribute_template()

    # Update the Metadata for each of the Measurements
    timedata["Bx GSE"].meta.update(
        OrderedDict({"CATDESC": "X component of magnetic Field GSE"}))
    timedata["By GSE"].meta.update(
        OrderedDict({"CATDESC": "Y component of magnetic Field GSE"}))

    # You can also add new Measurements to the TimeData container
    bz = np.random.choice(a=[-1, 0, 1], size=1000).cumsum(0)
    timedata.add_measurement(
        measure_name="Bz GSE",
        data=u.Quantity(value=bz, unit="nanoTesla", dtype=np.int16),
        meta={
            "VAR_TYPE": "data",
            "CATDESC": "Z component of magnetic Field GSE",
        },
    )

    # create the CDF File
    cdf_file_path = timedata.save(output_path="./", overwrite=True)

The file that this code generates is made available as a sample file in this
repository in :file:`hermes_core/data/sample/hermes_nms_default_l1_20160322_123031_v0.0.1.cdf`.