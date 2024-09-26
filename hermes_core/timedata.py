"""
Container class for Measurement Data.
"""

from pathlib import Path
from collections import OrderedDict
from typing import Optional, Union
import astropy
import ndcube
from swxsoc.swxdata import SWXData
from swxsoc.util.util import VALID_DATA_LEVELS
from swxsoc.util.validation import CDFValidator
import hermes_core
from hermes_core.util.schema import HermesDataSchema

__all__ = ["HermesData"]


class HermesData(SWXData):
    """
    A generic object for loading, storing, and manipulating HERMES time series data.

    Parameters
    ----------
    timeseries :  `Union[astropy.timeseries.TimeSeries, Dict[str, astropy.timeseries.TimeSeries]]`
        The time-series data. This can be a single `astropy.timeseries.TimeSeries` object or a dictionary of `str` to `astropy.timeseries.TimeSeries` objects. If a dictionary, one key must be named 'epoch', the primary time axis. If non-index/time columns are included in any of the TimeSeries objects, they must be `~astropy.units.Quantity` arrays.
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
    ...     "data_mask": NDData(data=np.eye(100, 100, dtype=np.uint16), meta={"CATDESC": "Data Mask", "VAR_TYPE": "metadata"})
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
        timeseries: Union[
            astropy.timeseries.TimeSeries, dict[str, astropy.timeseries.TimeSeries]
        ],
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

    @staticmethod
    def global_attribute_template(
        instr_name: str = "", data_level: str = "", version: str = ""
    ) -> OrderedDict:
        """
        Function to generate a template of the required ISTP-compliant global attributes.

        Parameters
        ----------
        instr_name : `str`
            The instrument name. Must be "eea", "nemisis", "merit" or "spani".
        data_level : `str`
            The data level of the data. Must be "l0", "l1", "ql", "l2", "l3", "l4"
        version : `str`
            Must be of the form X.Y.Z.

        Returns
        -------
        template : `collections.OrderedDict`
            A template for required global attributes.
        """
        meta = HermesDataSchema().global_attribute_template()

        # Check the Optional Instrument Name
        if instr_name:
            if instr_name not in hermes_core.config["mission"]["inst_names"]:
                raise ValueError(
                    f"Instrument, {instr_name}, is not recognized. Must be one of {hermes_core.config['mission']['inst_names']}."
                )
            # Set the Property
            meta["Descriptor"] = (
                f"{instr_name.upper()}>{hermes_core.config['mission']['inst_to_fullname'][instr_name]}"
            )

        # Check the Optional Data Level
        if data_level:
            if data_level not in VALID_DATA_LEVELS:
                raise ValueError(
                    f"Level, {data_level}, is not recognized. Must be one of {VALID_DATA_LEVELS[1:]}."
                )
            # Set the Property
            if data_level != "ql":
                meta["Data_level"] = f"{data_level.upper()}>Level {data_level[1]}"
            else:
                meta["Data_level"] = f"{data_level.upper()}>Quicklook"

        # Check the Optional Data Version
        if version:
            # check that version is in the right format with three parts
            if len(version.split(".")) != 3:
                raise ValueError(
                    f"Version, {version}, is not formatted correctly. Should be X.Y.Z"
                )
            meta["Data_version"] = version
        return meta

    @staticmethod
    def measurement_attribute_template() -> OrderedDict:
        """
        Function to generate a template of the required measurement attributes.

        Returns
        -------
        template : `collections.OrderedDict`
            A template for required variable attributes that must be provided.
        """
        return HermesDataSchema().measurement_attribute_template()

    @staticmethod
    def validate(file_path: Path) -> list[str]:
        """
        Validate the CDF file.

        Parameters
        ----------
        file_path : `pathlib.Path`
            A fully specified file path of the CDF data file to validate.

        Returns
        -------
        errors : `list[str]`
            A list of validation errors returned. A valid file will result in an empty list being returned.
        """
        schema = HermesDataSchema()
        validator = CDFValidator(schema)
        return validator.validate(file_path)
