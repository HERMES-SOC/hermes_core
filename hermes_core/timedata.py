"""
Container class for Measurement Data.
"""

from typing import Optional, Union
import astropy
import ndcube
from swxsoc.swxdata import SWXData
from hermes_core.util.schema import HermesDataSchema

__all__ = ["HermesData"]


class HermesData(SWXData):
    """
    A generic object for loading, storing, and manipulating HERMES time series data.

    Parameters
    ----------
    timeseries :  `astropy.timeseries.TimeSeries`
        The time series of data. Columns must be `~astropy.units.Quantity` arrays.
    support : `Optional[dict[Union[astropy.units.Quantity, astropy.nddata.NDData]]]`
        Support data arrays which do not vary with time (i.e. Non-Record-Varying data).
    spectra : `Optional[ndcube.NDCollection]`
        One or more `ndcube.NDCube` objects containing spectral or higher-dimensional
        timeseries data.
    meta : `Optional[dict]`
        The metadata describing the time series in an ISTP-compliant format.

    Examples
    --------
    >>> import numpy as np
    >>> import astropy.units as u
    >>> from astropy.timeseries import TimeSeries
    >>> from ndcube import NDCube, NDCollection
    >>> from astropy.wcs import WCS
    >>> from astropy.nddata import NDData
    >>> from hermes_core.timedata import HermesData
    >>> # Create a TimeSeries structure
    >>> data = u.Quantity([1, 2, 3, 4], "gauss", dtype=np.uint16)
    >>> ts = TimeSeries(time_start="2016-03-22T12:30:31", time_delta=3 * u.s, data={"Bx": data})
    >>> # Create a Spectra structure
    >>> spectra = NDCollection(
    ...     [
    ...         (
    ...             "test_spectra",
    ...             NDCube(
    ...                 data=np.random.random(size=(4, 10)),
    ...                 wcs=WCS(naxis=2),
    ...                 meta={"CATDESC": "Test Spectra Variable"},
    ...                 unit="eV",
    ...             ),
    ...         )
    ...     ]
    ... )
    >>> # Create a Support Structure
    >>> support_data = {
    ...     "data_mask": NDData(data=np.eye(100, 100, dtype=np.uint16))
    ... }
    >>> # Create Global Metadata Attributes
    >>> input_attrs = HermesData.global_attribute_template("eea", "l1", "1.0.0")
    >>> # Create HermesData Object
    >>> hermes_data = HermesData(timeseries=ts, support=support_data, spectra=spectra, meta=input_attrs)

    Raises
    ------
    ValueError: If the number of columns is less than 2 or the required 'time' column is missing.
    TypeError: If any column, excluding 'time', is not an `astropy.units.Quantity` object with units.
    ValueError: If the elements of a `TimeSeries` column are multidimensional
    TypeError: If any `supoport` data elements are not type `astropy.nddata.NDData`
    TypeError: If `spectra` is not an `NDCollection` object.

    References
    ----------
    * `Astropy TimeSeries <https://docs.astropy.org/en/stable/timeseries/index.html/>`_
    * `Astropy Quantity and Units <https://docs.astropy.org/en/stable/units/index.html>`_
    * `Astropy Time <https://docs.astropy.org/en/stable/time/index.html>`_
    * `Astropy NDData <https://docs.astropy.org/en/stable/nddata/>`_
    * `Sunpy NDCube and NDCollection <https://docs.sunpy.org/projects/ndcube/en/stable/>`_
    * `Space Physics Guidelines for CDF (ISTP) <https://spdf.gsfc.nasa.gov/istp_guide/istp_guide.html>`_
    """

    def __init__(
        self,
        timeseries: astropy.timeseries.TimeSeries,
        support: Optional[
            dict[Union[astropy.units.Quantity, astropy.nddata.NDData]]
        ] = None,
        spectra: Optional[ndcube.NDCollection] = None,
        meta: Optional[dict] = None,
    ):
        # Call SwxSOC Init
        super().__init__(timeseries, support, spectra, meta)

        # Derive Metadata using custom Schema
        self.schema = HermesDataSchema()
        self._derive_metadata()
