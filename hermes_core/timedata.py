"""
Container class for Measurement Data.
"""

from pathlib import Path
from collections import OrderedDict
from copy import deepcopy
from typing import Optional, Union
import numpy as np
import astropy
from astropy.time import Time
from astropy.timeseries import TimeSeries
from astropy.table import vstack
from astropy.nddata import NDData
from astropy import units as u
import ndcube
from ndcube import NDCube, NDCollection
import hermes_core
from hermes_core.util.schema import HermesDataSchema
from hermes_core.util.exceptions import warn_user
from hermes_core.util.util import VALID_DATA_LEVELS

__all__ = ["HermesData"]


class HermesData:
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
        # ================================================
        #               VALIDATE INPUTS
        # ================================================

        # Verify TimeSeries compliance
        if isinstance(timeseries, dict):
            for key, value in timeseries.items():
                self._validate_timeseries(value)
        else:
            self._validate_timeseries(timeseries)

        # Global Metadata Attributes are compiled from two places. You can pass in
        # global metadata throug the `meta` parameter or through the `TimeSeries.meta`
        # attribute.
        self._meta = {}
        if meta is not None and isinstance(meta, dict):
            self._meta.update(meta)
        if (
            isinstance(timeseries, TimeSeries)
            and timeseries.meta is not None
            and isinstance(timeseries.meta, dict)
        ):
            self._meta.update(timeseries.meta)

        # Check Global Metadata Requirements - Require Descriptor, Data_level, Data_Version
        if "Descriptor" not in self._meta or self._meta["Descriptor"] is None:
            raise ValueError(
                "'Descriptor' global meta attribute required for HERMES Instrument name"
            )
        if "Data_level" not in self._meta or self._meta["Data_level"] is None:
            raise ValueError(
                "'Data_level' global meta attribute required for HERMES data level"
            )
        if "Data_version" not in self._meta or self._meta["Data_version"] is None:
            raise ValueError(
                "'Data_version' global meta attribute is required for HERMES data version"
            )

        # Check NRV Data
        if support is not None:
            for key in support:
                if not (
                    isinstance(support[key], u.Quantity)
                    or isinstance(support[key], NDData)
                ):
                    raise TypeError(
                        f"Variable '{key}' must be an astropy.units.Quantity or astropy.nddata.NDData object"
                    )

        # Check Higher-Dimensional Spectra
        if spectra is not None:
            if not isinstance(spectra, NDCollection):
                raise TypeError(f"Spectra must be an ndcube.NDCollection object")

        # ================================================
        #         CREATE HERMES DATA STRUCTURES
        # ================================================

        if isinstance(timeseries, dict):
            self._timeseries = {}
            for key, value in timeseries.items():
                # Copy the TimeSeries
                self._timeseries[key] = TimeSeries(value, copy=True)
                # Add any Metadata from the original TimeSeries
                self._update_timeseries_measurement_meta(
                    timeseries=value, epoch_key=key
                )
        elif isinstance(timeseries, TimeSeries):
            self._timeseries = {
                hermes_core.DEFAULT_TIMESERIES_KEY: TimeSeries(timeseries, copy=True)
            }
            self._update_timeseries_measurement_meta(
                timeseries=timeseries, epoch_key=hermes_core.DEFAULT_TIMESERIES_KEY
            )

        # Copy the Non-Record Varying Data
        if support:
            self._support = deepcopy(support)
        else:
            self._support = {}

        # Add Support Metadata
        for key in self._support:
            self._support[key].meta = self.measurement_attribute_template()
            if hasattr(support[key], "meta"):
                self._support[key].meta.update(support[key].meta)

        # Copy the High-Dimensional Spectra
        if spectra:
            self._spectra = spectra
        else:
            self._spectra = NDCollection([])

        # ================================================
        #           DERIVE METADATA ATTRIBUTES
        # ================================================

        # Derive Metadata
        self.schema = HermesDataSchema()
        self._derive_metadata()

    @property
    def timeseries(self):
        """
        (`astropy.timeseries.TimeSeries` or `dict`) A `TimeSeries` representing one or more measurements as a function of time.
        If there are multiple `TimeSeries`, a dictionary is returned.
        """
        if len(self._timeseries) > 1:
            return {key: value for key, value in self._timeseries.items()}
        else:
            return self._timeseries[hermes_core.DEFAULT_TIMESERIES_KEY]

    @property
    def support(self):
        """
        (`dict[Union[astropy.units.Quantity, astropy.nddata.NDData]]`) A `dict` containing one or more non-time-varying support variables.
        """
        return self._support

    @property
    def spectra(self):
        """
        (`ndcube.NDCollection]`) A `NDCollection` object containing high-dimensional spectra data.
        """
        return self._spectra

    @property
    def data(self):
        """
        (`dict`) A `dict` containing each of `timeseries`, `spectra` and `support`.
        """
        return {
            "timeseries": self._timeseries,
            "spectra": self._spectra,
            "support": self._support,
        }

    @property
    def meta(self):
        """
        (`collections.OrderedDict`) Global metadata associated with the measurement data.
        """
        return self._meta

    @property
    def time(self):
        """
        (`astropy.time.Time`) The times of the measurements.
        """
        t = Time(self._timeseries[hermes_core.DEFAULT_TIMESERIES_KEY].time)
        # Set time format to enable plotting with astropy.visualisation.time_support()
        t.format = "iso"
        return t

    @property
    def time_range(self):
        """
        (`tuple`) The start and end times of the times.
        """
        return (self.time.min(), self.time.max())

    def __repr__(self):
        """
        Returns a representation of the `HermesData` class.
        """
        return self.__str__()

    def __str__(self):
        """
        Returns a string representation of the `HermesData` class.
        """
        str_repr = f"HermesData() Object:\n"
        # Global Attributes/Metedata
        str_repr += f"Global Attrs:\n"
        for attr_name, attr_value in self._meta.items():
            str_repr += f"\t{attr_name}: {attr_value}\n"
        # TimeSeries Data
        str_repr += f"TimeSeries Data:\n"
        for epoch_key, ts in self._timeseries.items():
            str_repr += f"\tTimeSeries: {epoch_key}\n"
            for var_name in ts.colnames:
                str_repr += f"\t\t{var_name}\n"
        # Support Data
        str_repr += f"Support Data:\n"
        for var_name in self._support.keys():
            str_repr += f"\t{var_name}\n"
        # Spectra Data
        str_repr += f"Spectra Data:\n"
        for var_name in self._spectra.keys():
            str_repr += f"\t{var_name}\n"
        return str_repr

    def __getitem__(self, var_name):
        """
        Get the data for a specific variable.

        Parameters
        ----------
        var_name : `str`
            The name of the variable to retrieve.

        Returns
        -------
        `astropy.timeseries.TimeSeries` or `astropy.nddata.NDData`
            The data for the variable.
        """
        for epoch_key, ts in self._timeseries.items():
            if var_name in ts.columns:
                return ts[var_name]
        if var_name in self.support:
            return self.support[var_name]
        if var_name in self.spectra:
            return self.spectra[var_name]
        else:
            raise KeyError(f"Variable {var_name} not found in HermesData object.")

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
        meta = HermesDataSchema.global_attribute_template()

        # Check the Optional Instrument Name
        if instr_name:
            if instr_name not in hermes_core.INST_NAMES:
                raise ValueError(
                    f"Instrument, {instr_name}, is not recognized. Must be one of {hermes_core.INST_NAMES}."
                )
            # Set the Property
            meta["Descriptor"] = (
                f"{instr_name.upper()}>{hermes_core.INST_TO_FULLNAME[instr_name]}"
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
        return HermesDataSchema.measurement_attribute_template()

    @staticmethod
    def get_timeseres_epoch_key(timeseries, var_data, var_meta: dict = None):
        """
        Function to determine the TimeSeries Epoch for a Record-Varying Variable.

        Parameters
        ----------
        timeseries : `dict[str, astropy.timeseries.TimeSeries]`
            A dictionary of `str` to `astropy.timeseries.TimeSeries` objects. Each `TimeSeries` object represents a different epoch.
        var_data : `astropy.units.Quantity`
            The variable data that we want to find the epoch for.
        var_meta : `dict`, optional
            The metadata associated with the variable data.
        """

        # Find the TimeSeries Epoch for this Record-Varying Variable
        if var_meta is not None and "DEPEND_0" in var_meta:
            epoch_key = var_meta["DEPEND_0"]
        else:
            # Check which epoch key to use
            potential_epoch_keys = []
            for key, ts in timeseries.items():
                if hasattr(var_data, "shape"):
                    if len(ts.time) == var_data.shape[0]:
                        potential_epoch_keys.append(key)
                elif hasattr(var_data, "data"):
                    if len(ts.time) == len(var_data.data):
                        potential_epoch_keys.append(key)
            if len(potential_epoch_keys) == 0:
                raise ValueError("No TimeSeries have the same length as the new data.")
            elif len(potential_epoch_keys) > 1:
                raise ValueError(
                    "Multiple TimeSeries have the same length as the new data."
                )
            epoch_key = potential_epoch_keys[0]
        return epoch_key

    def _validate_timeseries(self, timeseries: astropy.timeseries.TimeSeries):
        """
        Validate a timeseries.

        Parameters
        ----------
        timeseries : astropy.timeseries.TimeSeries
            The timeseries to validate.

        Raises
        ------
        TypeError
            If the timeseries is not a `astropy.timeseries.TimeSeries` object or a dictionary of `str` to `astropy.timeseries.TimeSeries` objects.
            If any column in the timeseries (other than 'time') is not an `astropy.units.Quantity` object.
        ValueError
            If the timeseries is empty.
            If any column in the timeseries is not a one-dimensional measurement.
        """
        if not isinstance(timeseries, astropy.timeseries.TimeSeries):
            raise TypeError(
                "timeseries must be a `astropy.timeseries.TimeSeries` object or a dictionary of `str` to `astropy.timeseries.TimeSeries` objects."
            )
        if (
            isinstance(timeseries, astropy.timeseries.TimeSeries)
            and len(timeseries) == 0
        ):
            raise ValueError(
                "timeseries cannot be empty, must include at least a 'time' column with valid times"
            )
        for colname in timeseries.columns:
            # Verify that all Measurements are `Quantity`
            if colname != "time" and not isinstance(timeseries[colname], u.Quantity):
                raise TypeError(
                    f"Column '{colname}' must be an astropy.units.Quantity object"
                )
            # Verify that the Column is only a single dimension
            if len(timeseries[colname].shape) > 1:  # If there is more than 1 Dimension
                raise ValueError(
                    f"Column '{colname}' must be a one-dimensional measurement. Split additional dimensions into unique measurements."
                )

    def _update_timeseries_measurement_meta(
        self, timeseries: TimeSeries, epoch_key: str
    ):
        """
        Update the metadata for a specific timeseries in the collection.

        This method updates the metadata for both the time attribute and the measurements
        in the timeseries. If the time attribute or a measurement has a `meta` attribute,
        its contents are added to the corresponding attribute in the stored timeseries.

        Parameters
        ----------
        timeseries : `astropy.timeseries.TimeSeries`
            The timeseries whose metadata is to be updated. This timeseries should already
            be part of the collection.
        epoch_key : str
            The key identifying the timeseries in the collection.
        """
        # Time Attributes
        self._timeseries[epoch_key]["time"].meta = OrderedDict()
        if hasattr(timeseries["time"], "meta"):
            self._timeseries[epoch_key]["time"].meta.update(timeseries["time"].meta)
        # Measurement Attributes
        for col in timeseries.columns:
            if col != "time":
                self._timeseries[epoch_key][
                    col
                ].meta = self.measurement_attribute_template()
                if hasattr(timeseries[col], "meta"):
                    self._timeseries[epoch_key][col].meta.update(timeseries[col].meta)

    def _derive_metadata(self):
        """
        Funtion to derive global and measurement metadata based on a HermesDataSchema
        """

        # Get Default Metadata
        for attr_name, attr_value in self.schema.default_global_attributes.items():
            self._update_global_attribute(attr_name, attr_value)

        # Global Attributes
        for attr_name, attr_value in self.schema.derive_global_attributes(self).items():
            self._update_global_attribute(attr_name, attr_value)

        for epoch_key, ts in self._timeseries.items():
            # Time Measurement Attributes
            for attr_name, attr_value in self.schema.derive_time_attributes(ts).items():
                self._update_timeseries_attribute(
                    epoch_key=epoch_key,
                    var_name="time",
                    attr_name=attr_name,
                    attr_value=attr_value,
                )

            # Other Measurement Attributes
            for col in [col for col in ts.columns if col != "time"]:
                for attr_name, attr_value in self.schema.derive_measurement_attributes(
                    self, col
                ).items():
                    self._update_timeseries_attribute(
                        epoch_key=epoch_key,
                        var_name=col,
                        attr_name=attr_name,
                        attr_value=attr_value,
                    )

        # Support/ Non-Record-Varying Data
        for col in self._support:
            for attr_name, attr_value in self.schema.derive_measurement_attributes(
                self, col
            ).items():
                self._update_support_attribute(
                    var_name=col, attr_name=attr_name, attr_value=attr_value
                )

        # Spectra/ High-Dimensional Data
        for col in self._spectra:
            for attr_name, attr_value in self.schema.derive_measurement_attributes(
                self, col
            ).items():
                self._update_spectra_attribute(
                    var_name=col, attr_name=attr_name, attr_value=attr_value
                )

    def _update_global_attribute(self, attr_name, attr_value):
        # If the attribute is set, check if we want to overwrite it
        if attr_name in self._meta and self._meta[attr_name] is not None:
            # We want to overwrite if:
            #   1) The actual value is not the derived value
            #   2) The schema marks this attribute to be overwriten
            if (
                self._meta[attr_name] != attr_value
                and self.schema.global_attribute_schema[attr_name]["overwrite"]
            ):
                warn_user(
                    f"Overriding Global Attribute {attr_name} : {self._meta[attr_name]} -> {attr_value}"
                )
                self._meta[attr_name] = attr_value
        # If the attribute is not set, set it
        else:
            self._meta[attr_name] = attr_value

    def _update_timeseries_attribute(self, epoch_key, var_name, attr_name, attr_value):
        if (
            attr_name in self._timeseries[epoch_key][var_name].meta
            and self._timeseries[epoch_key][var_name].meta[attr_name] is not None
            and attr_name in self.schema.variable_attribute_schema["attribute_key"]
        ):
            attr_schema = self.schema.variable_attribute_schema["attribute_key"][
                attr_name
            ]
            if (
                self._timeseries[epoch_key][var_name].meta[attr_name] != attr_value
                and attr_schema["overwrite"]
            ):
                warn_user(
                    f"Overriding TimeSeries {var_name} Attribute {attr_name} : {self._timeseries[epoch_key][var_name].meta[attr_name]} -> {attr_value}"
                )
                self._timeseries[epoch_key][var_name].meta[attr_name] = attr_value
        else:
            self._timeseries[epoch_key][var_name].meta[attr_name] = attr_value

    def _update_support_attribute(self, var_name, attr_name, attr_value):
        if (
            attr_name in self._support[var_name].meta
            and self._support[var_name].meta[attr_name] is not None
            and attr_name in self.schema.variable_attribute_schema["attribute_key"]
        ):
            attr_schema = self.schema.variable_attribute_schema["attribute_key"][
                attr_name
            ]
            if (
                self._support[var_name].meta[attr_name] != attr_value
                and attr_schema["overwrite"]
            ):
                warn_user(
                    f"Overriding Support {var_name} Attribute {attr_name} : {self._support[var_name].meta[attr_name]} -> {attr_value}"
                )
                self._support[var_name].meta[attr_name] = attr_value
        else:
            self._support[var_name].meta[attr_name] = attr_value

    def _update_spectra_attribute(self, var_name, attr_name, attr_value):
        if (
            attr_name in self._spectra[var_name].meta
            and self._spectra[var_name].meta[attr_name] is not None
            and attr_name in self.schema.variable_attribute_schema["attribute_key"]
        ):
            attr_schema = self.schema.variable_attribute_schema["attribute_key"][
                attr_name
            ]
            if (
                self._spectra[var_name].meta[attr_name] != attr_value
                and attr_schema["overwrite"]
            ):
                warn_user(
                    f"Overriding Spectra {var_name} Attribute {attr_name} : {self._spectra[var_name].meta[attr_name]} -> {attr_value}"
                )
                self._spectra[var_name].meta[attr_name] = attr_value
        else:
            self._spectra[var_name].meta[attr_name] = attr_value

    def add_measurement(self, measure_name: str, data: u.Quantity, meta: dict = None):
        """
        Add a new time-varying scalar measurement (column).

        Parameters
        ----------
        measure_name: `str`
            Name of the measurement to add.
        data: `astropy.units.Quantity`
            The data to add. Must have the same time stamps as the existing data.
        meta: `dict`, optional
            The metadata associated with the measurement.

        Raises
        ------
        TypeError: If var_data is not of type Quantity.
        ValueError: If data has more than one dimension
        """
        # Verify that all Measurements are `Quantity`
        if (not isinstance(data, u.Quantity)) or (not data.unit):
            raise TypeError(
                f"Measurement {measure_name} must be type `astropy.units.Quantity` and have `unit` assigned."
            )
        # Verify that the Column is only a single dimension
        if len(data.shape) > 1:  # If there is more than 1 Dimension
            raise ValueError(
                f"Column '{measure_name}' must be a one-dimensional measurement. Split additional dimensions into unique measurenents."
            )

        # Find the TimeSeries Epoch for this Record-Varying Variable
        epoch_key = HermesData.get_timeseres_epoch_key(self._timeseries, data, meta)

        # Add the new measurement
        self._timeseries[epoch_key][measure_name] = data
        # Add any Metadata from the original Quantity
        self._timeseries[epoch_key][
            measure_name
        ].meta = self.measurement_attribute_template()
        if hasattr(data, "meta"):
            self._timeseries[epoch_key][measure_name].meta.update(data.meta)
        if meta:
            self._timeseries[epoch_key][measure_name].meta.update(meta)

        # Derive Metadata Attributes for the Measurement
        self._derive_metadata()

    def add_timeseries(self, epoch_key: str, timeseries: TimeSeries):
        """
        Add a new TimeSeries object to the collection of epochs.

        Parameters
        ----------
        epoch_key: `str`
            The key to identify the new TimeSeries.
        timeseries: `astropy.timeseries.TimeSeries`
            The time-series data to add.
        """
        self._validate_timeseries(timeseries)

        # Check the epoch is not already used
        if epoch_key in self._timeseries:
            raise ValueError(f"Epoch key {epoch_key} is already in use.")
        else:
            # Add the TimeSeries
            self._timeseries[epoch_key] = TimeSeries(timeseries, copy=True)
            # Updata the Metadata
            self._update_timeseries_measurement_meta(
                timeseries=timeseries, epoch_key=epoch_key
            )

    def add_support(
        self,
        name: str,
        data: Union[astropy.units.Quantity, astropy.nddata.NDData],
        meta: Optional[dict] = None,
    ):
        """
        Add a new non-time-varying data array.

        Parameters
        ----------
        name: `str`
            Name of the data array to add.
        data: `Union[astropy.units.Quantity, astropy.nddata.NDData]`,
            The data to add.
        meta: `Optional[dict]`, optional
            The metadata associated for the data array.

        Raises
        ------
        TypeError: If var_data is not of type NDData.
        """
        # Verify that all Measurements are `NDData`
        if not (isinstance(data, u.Quantity) or isinstance(data, NDData)):
            raise TypeError(f"Measurement {name} must be type `astropy.nddata.NDData`.")

        self._support[name] = data
        # Add any Metadata from the original Quantity or NDData
        if hasattr(data, "meta"):
            self._support[name].meta.update(data.meta)
        else:
            self._support[name].meta = self.measurement_attribute_template()
        # Add any Metadata Passed not in the NDData
        if meta:
            self._support[name].meta.update(meta)

        # Derive Metadata Attributes for the Measurement
        self._derive_metadata()

    def add_spectra(self, name: str, data: NDCube, meta: dict = None):
        """
        Add a new time-varying vector measurement. This include higher-dimensional time-varying
        data.

        Parameters
        ----------
        name: `str`
            Name of the measurement to add.
        data: `ndcube.NDCube`
            The data to add. Must have the same time stamps as the existing data.
        meta: `dict`, optional
            The metadata associated with the measurement.

        Raises
        ------
        TypeError: If var_data is not of type NDCube.
        """
        # Verify that all Measurements are `NDCube`
        if not isinstance(data, NDCube):
            raise TypeError(f"Measurement {name} must be type `ndcube.NDCube`.")

        # Add the new measurement
        if len(self._spectra) == 0:
            aligned_axes = (0,)
            self._spectra = NDCollection([(name, data)], aligned_axes)
        else:
            # Check to see if we need to maintain the aligned axes
            if self._spectra.aligned_axes:
                first_aligned_axes = self._spectra.aligned_axes[
                    self._spectra._first_key
                ]
                aligned_axes = tuple(0 for _ in range(len(first_aligned_axes)))
                self._spectra.update([(name, data)], aligned_axes)
            else:
                self._spectra.update([(name, data)], self._spectra.aligned_axes)

        # Add any Metadata Passed not in the NDCube
        if meta:
            self._spectra[name].meta.update(meta)

        # Derive Metadata Attributes for the Measurement
        self._derive_metadata()

    def remove(self, measure_name: str):
        """
        Remove an existing measurement or support data array.

        Parameters
        ----------
        measure_name: `str`
            Name of the variable to remove.
        """
        found = False
        # Check TimeSeries
        for epoch_key, ts in self._timeseries.items():
            if measure_name in ts.columns:
                self._timeseries[epoch_key].remove_column(measure_name)
                found = True
        # Check Support
        if measure_name in self._support:
            self._support.pop(measure_name)
            found = True
        # Check Spectra
        elif measure_name in self._spectra:
            self._spectra.pop(measure_name)
            found = True
        # Otherwise Raise and Error
        if not found:
            raise ValueError(f"Data for Measurement {measure_name} not found.")

    def plot(self, axes=None, columns=None, subplots=True, **plot_args):
        """
        Plot the measurement data.

        Parameters
        ----------
        axes : `~matplotlib.axes.Axes`, optional
            If provided the image will be plotted on the given axes.
            Defaults to `None` and creates a new axis.
        columns : `list[str]`, optional
            If provided, only plot the specified measurements otherwise try to plot them all.
        subplots : `bool`
            If set, all columns are plotted in their own plot panel.
        **plot_args : `dict`, optional
            Additional plot keyword arguments that are handed to
            `~matplotlib.axes.Axes`.

        Returns
        -------
        `~matplotlib.axes.Axes`
            The plot axes.
        """
        from astropy.visualization import quantity_support, time_support
        from matplotlib.axes import Axes

        # Set up the plot axes based on the number of columns to plot
        axes, columns = self._setup_axes_columns(axes, columns, subplots=subplots)
        quantity_support()
        time_support()

        if subplots:
            i = 0
            if isinstance(axes, Axes):  # subplots is true but only one column given
                iter_axes = [axes]
            else:
                iter_axes = axes
            for this_ax, this_col in zip(iter_axes, columns):
                if i == 0:
                    this_ax.set_title(
                        f'{self.meta["Mission_group"]} {self.meta["Descriptor"]} {self.meta["Data_level"]}'
                    )
                    i += 1
                this_ax.plot(self.time, self.timeseries[this_col], **plot_args)
                this_ax.set_ylabel(self.timeseries[this_col].meta["LABLAXIS"])
        else:
            axes.set_title(
                f'{self.meta["Mission_group"]} {self.meta["Descriptor"]} {self.meta["Data_level"]}'
            )
            for this_col in columns:
                axes.plot(
                    self.time,
                    self.timeseries[this_col],
                    label=self.timeseries[this_col].meta["LABLAXIS"],
                    **plot_args,
                )
            axes.legend()
        # Setup the Time Axis
        self._setup_x_axis(axes)

        return axes

    def _setup_axes_columns(self, axes, columns, subplots=False):
        """
        Validate data for plotting, and get default axes/columns if not passed
        by the user.

        Code courtesy of sunpy.
        """
        import matplotlib.pyplot as plt

        # If no individual columns were input, try to plot all columns
        if columns is None:
            columns = list(self.timeseries.columns.copy())
            columns.remove("time")
        # Create Axes or Subplots for displaying the data
        if axes is None:
            if not subplots:
                axes = plt.gca()
            else:
                axes = plt.gcf().subplots(ncols=1, nrows=len(columns), sharex=True)

        return axes, columns

    @staticmethod
    def _setup_x_axis(ax):
        """
        Shared code to set x-axis properties.

        Code courtesy of sunpy.
        """
        import matplotlib.dates as mdates

        if isinstance(ax, np.ndarray):
            ax = ax[-1]

        locator = ax.xaxis.get_major_locator()
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))

    def append(self, timeseries: TimeSeries):
        """
        Add additional measurements to an existing column.

        Parameters
        ----------
        timeseries : `astropy.timeseries.TimeSeries`
            The data to be appended (rows) as a `TimeSeries` object.
        """
        # Verify TimeSeries compliance
        self._validate_timeseries(timeseries)

        # Check which epoch key to use
        selected_epoch_key = HermesData.get_timeseres_epoch_key(
            self._timeseries, timeseries.time
        )

        # Save Metadata since it is not carried over with vstack
        metadata_holder = {
            col: self._timeseries[selected_epoch_key][col].meta
            for col in self._timeseries[selected_epoch_key].columns
        }

        # Vertically Stack the TimeSeries
        self._timeseries[selected_epoch_key] = vstack(
            [self._timeseries[selected_epoch_key], timeseries]
        )

        # Add Metadata back to the Stacked TimeSeries
        for col in self._timeseries[selected_epoch_key].columns:
            self._timeseries[selected_epoch_key][col].meta = metadata_holder[col]

        # Re-Derive Metadata
        self._derive_metadata()

    def save(self, output_path: Path = None, overwrite: bool = False):
        """
        Save the data to a HERMES CDF file.

        Parameters
        ----------
        output_path : `pathlib.Path`, optional
            This can take two forms:

            1. A fully specified path to the directory where the file is to be saved. This will save the file with the default name of the Logical_file_id metadata attribute. This is the recommended approach to saving files.

            2. A fully specified path to the file to be saved, including the file name. This will ignore the Logical_file_id metadata attribute and save the file with the provided name. This is not the recommended approach, although it can be used for local development and testing.

            If not provided, saves to the current directory, using the file name of the Logical_file_id metadata attribute.
        overwrite : `bool`
            If set, overwrites existing file of the same name. This parameter should be used for development and testing purposes only. This parameter should not be set to `True` when saving data for archiving as a part of the production data processing pipeline.

        Returns
        -------
        path : `pathlib.Path`
            A path to the saved file.

        Raises
        ------
        FileNotFoundError: If the output_path points to a directory that does not exist.
        ValueError: If the output_path is a file and does not have a recognized extension.
        """
        from hermes_core.util.io import CDFHandler

        handler = CDFHandler()
        if not output_path:
            output_path = Path.cwd()
        return handler.save_data(data=self, file_path=output_path, overwrite=overwrite)

    @classmethod
    def load(cls, file_path: Path):
        """
        Load data from a file.

        Parameters
        ----------
        file_path : `pathlib.Path`
            A fully specified file path of the data file to load.

        Returns
        -------
        data : `HermesData`
            A `HermesData` object containing the loaded data.

        Raises
        ------
        ValueError: If the file type is not recognized as a file type that can be loaded.

        """
        from hermes_core.util.io import CDFHandler

        # Determine the file type
        file_extension = file_path.suffix

        # Create the appropriate handler object based on file type
        if file_extension == ".cdf":
            handler = CDFHandler()
        else:
            raise ValueError(f"Unsupported file type: {file_extension}")

        # Load data using the handler and return a HermesData object
        timeseries, support, spectra, meta = handler.load_data(file_path)
        return cls(timeseries=timeseries, support=support, spectra=spectra, meta=meta)
