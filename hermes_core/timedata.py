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
from hermes_core.util.io import CDFHandler
from hermes_core.util.schema import HermesDataSchema
from hermes_core.util.exceptions import warn_user
from hermes_core.util.util import VALID_DATA_LEVELS

__all__ = ["HermesData"]


class HermesData:
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
        # ================================================
        #               VALIDATE INPUTS
        # ================================================

        # Verify TimeSeries compliance
        if not isinstance(timeseries, TimeSeries):
            raise TypeError(
                "timeseries must be a `astropy.timeseries.TimeSeries` object."
            )

        if len(timeseries) == 0:
            raise ValueError(
                "timeseries cannot be empty, must include at least a 'time' column with valid times"
            )

        # Check individual Columns
        for colname in timeseries.columns:
            # Verify that all Measurements are `Quantity`
            if colname != "time" and not isinstance(timeseries[colname], u.Quantity):
                raise TypeError(
                    f"Column '{colname}' must be an astropy.units.Quantity object"
                )
            # Verify that the Column is only a single dimension
            if len(timeseries[colname].shape) > 1:  # If there is more than 1 Dimension
                raise ValueError(
                    f"Column '{colname}' must be a one-dimensional measurement. Split additional dimensions into unique measurenents."
                )

        # Global Metadata Attributes are compiled from two places. You can pass in
        # global metadata throug the `meta` parameter or through the `TimeSeries.meta`
        # attribute.
        _meta = {}
        if meta is not None and isinstance(meta, dict):
            _meta.update(meta)
        if timeseries.meta is not None and isinstance(timeseries.meta, dict):
            _meta.update(timeseries.meta)

        # Check Global Metadata Requirements - Require Descriptor, Data_level, Data_Version
        if "Descriptor" not in _meta or _meta["Descriptor"] is None:
            raise ValueError(
                "'Descriptor' global meta attribute required for HERMES Instrument name"
            )
        if "Data_level" not in _meta or _meta["Data_level"] is None:
            raise ValueError(
                "'Data_level' global meta attribute required for HERMES data level"
            )
        if "Data_version" not in _meta or _meta["Data_version"] is None:
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

        # Copy the TimeSeries
        self._timeseries = TimeSeries(timeseries, copy=True)

        # Add Input Metadata
        if meta is not None and isinstance(meta, dict):
            self._timeseries.meta.update(meta)

        # Add any Metadata from the original TimeSeries
        self._timeseries["time"].meta = OrderedDict()
        if hasattr(timeseries["time"], "meta"):
            self._timeseries["time"].meta.update(timeseries["time"].meta)

        # Add TimeSeries Measurement Metadata
        for col in self._timeseries.columns:
            if col != "time":
                self._timeseries[col].meta = self.measurement_attribute_template()
                if hasattr(timeseries[col], "meta"):
                    self._timeseries[col].meta.update(timeseries[col].meta)

        # Copy the Non-Record Varying Data
        if support:
            self._support = deepcopy(support)
        else:
            self._support = {}

        # Add Support Metadata
        for key in self._support:
            if hasattr(support[key], "meta"):
                self._support[key].meta.update(support[key].meta)
            else:
                self._support[key].meta = self.measurement_attribute_template()

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
        (`astropy.timeseries.TimeSeries`) A `TimeSeries` representing one or more measurements as a function of time.
        """
        return self._timeseries

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
        (`dict`) A `dict` containing each of `timeseries` and `support`.
        """
        return {
            "timeseries": self.timeseries,
            "spectra": self.spectra,
            "support": self.support,
        }

    @property
    def meta(self):
        """
        (`collections.OrderedDict`) Global metadata associated with the measurement data.
        """
        return self._timeseries.meta

    @property
    def time(self):
        """
        (`astropy.time.Time`) The times of the measurements.
        """
        t = Time(self._timeseries.time)
        # Set time format to enable plotting with astropy.visualisation.time_support()
        t.format = "iso"
        return t

    @property
    def time_range(self):
        """
        (`tuple`) The start and end times of the times.
        """
        return (self._timeseries.time.min(), self._timeseries.time.max())

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
        for attr_name, attr_value in self._timeseries.meta.items():
            str_repr += f"\t{attr_name}: {attr_value}\n"
        # TimeSeries Data
        str_repr += f"TimeSeries Data:\n"
        for var_name in self._timeseries.colnames:
            str_repr += f"\t{var_name}\n"
        # Support Data
        str_repr += f"Support Data:\n"
        for var_name in self._support.keys():
            str_repr += f"\t{var_name}\n"
        # Spectra Data
        str_repr += f"Spectra Data:\n"
        for var_name in self._spectra.keys():
            str_repr += f"\t{var_name}\n"
        return str_repr

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
            meta[
                "Descriptor"
            ] = f"{instr_name.upper()}>{hermes_core.INST_TO_FULLNAME[instr_name]}"

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

    def _derive_metadata(self):
        """
        Funtion to derive global and measurement metadata based on a HermesDataSchema
        """

        # Get Default Metadata
        for attr_name, attr_value in self.schema.default_global_attributes.items():
            self._update_global_attribute(attr_name, attr_value)

        # Global Attributes
        for attr_name, attr_value in self.schema.derive_global_attributes(
            self._timeseries
        ).items():
            self._update_global_attribute(attr_name, attr_value)

        # Time Measurement Attributes
        for attr_name, attr_value in self.schema.derive_time_attributes(
            self._timeseries
        ).items():
            self._update_timeseries_attribute(
                var_name="time", attr_name=attr_name, attr_value=attr_value
            )

        # Other Measurement Attributes
        for col in [col for col in self._timeseries.columns if col != "time"]:
            for attr_name, attr_value in self.schema.derive_measurement_attributes(
                self._timeseries, col
            ).items():
                self._update_timeseries_attribute(
                    var_name=col, attr_name=attr_name, attr_value=attr_value
                )

        # Support/ Non-Record-Varying Data
        for col in self._support:
            for attr_name, attr_value in self.schema.derive_measurement_attributes(
                self._support, col
            ).items():
                self._update_support_attribute(
                    var_name=col, attr_name=attr_name, attr_value=attr_value
                )

        # Spectra/ High-Dimensional Data
        for col in self._spectra:
            for attr_name, attr_value in self.schema.derive_measurement_attributes(
                self._spectra, col
            ).items():
                self._update_spectra_attribute(
                    var_name=col, attr_name=attr_name, attr_value=attr_value
                )

    def _update_global_attribute(self, attr_name, attr_value):
        # If the attribute is set, check if we want to overwrite it
        if (
            attr_name in self._timeseries.meta
            and self._timeseries.meta[attr_name] is not None
        ):
            # We want to overwrite if:
            #   1) The actual value is not the derived value
            #   2) The schema marks this attribute to be overwriten
            if (
                self._timeseries.meta[attr_name] != attr_value
                and self.schema.global_attribute_schema[attr_name]["overwrite"]
            ):
                warn_user(
                    f"Overriding Global Attribute {attr_name} : {self._timeseries.meta[attr_name]} -> {attr_value}"
                )
                self._timeseries.meta[attr_name] = attr_value
        # If the attribute is not set, set it
        else:
            self._timeseries.meta[attr_name] = attr_value

    def _update_timeseries_attribute(self, var_name, attr_name, attr_value):
        if (
            attr_name in self.timeseries[var_name].meta
            and self.timeseries[var_name].meta[attr_name] is not None
            and attr_name in self.schema.variable_attribute_schema["attribute_key"]
        ):
            attr_schema = self.schema.variable_attribute_schema["attribute_key"][
                attr_name
            ]
            if (
                self.timeseries[var_name].meta[attr_name] != attr_value
                and attr_schema["overwrite"]
            ):
                warn_user(
                    f"Overriding TimeSeries {var_name} Attribute {attr_name} : {self.timeseries[var_name].meta[attr_name]} -> {attr_value}"
                )
                self.timeseries[var_name].meta[attr_name] = attr_value
        else:
            self.timeseries[var_name].meta[attr_name] = attr_value

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

        self._timeseries[measure_name] = data
        # Add any Metadata from the original Quantity
        self._timeseries[measure_name].meta = self.measurement_attribute_template()
        if hasattr(data, "meta"):
            self._timeseries[measure_name].meta.update(data.meta)
        if meta:
            self._timeseries[measure_name].meta.update(meta)

        # Derive Metadata Attributes for the Measurement
        self._derive_metadata()

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
        if measure_name in self._timeseries.columns:
            self._timeseries.remove_column(measure_name)
        elif measure_name in self._support:
            self._support.pop(measure_name)
        elif measure_name in self._spectra:
            self._spectra.pop(measure_name)
        else:
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
        if not isinstance(timeseries, TimeSeries):
            raise TypeError("Data must be a TimeSeries object.")
        if len(timeseries.columns) < 2:
            raise ValueError("Data must have at least 2 columns")
        if len(self.timeseries.columns) != len(timeseries.columns):
            raise ValueError(
                (
                    f"Shape of curent TimeSeries ({len(self.timeseries.columns)}) does not match",
                    f"shape of data to add ({len(timeseries.columns)}).",
                )
            )

        # Check individual Columns
        for colname in self.timeseries.columns:
            if colname != "time" and not isinstance(
                self.timeseries[colname], u.Quantity
            ):
                raise TypeError(
                    f"Column '{colname}' must be an astropy.Quantity object"
                )

        # Save Metadata since it is not carried over with vstack
        metadata_holder = {
            col: self.timeseries[col].meta for col in self.timeseries.columns
        }

        # Vertically Stack the TimeSeries
        self._timeseries = vstack([self._timeseries, timeseries])

        # Add Metadata back to the Stacked TimeSeries
        for col in self.timeseries.columns:
            self.timeseries[col].meta = metadata_holder[col]

        # Re-Derive Metadata
        self._derive_metadata()

    def save(self, output_path: str = None, overwrite: bool = False):
        """
        Save the data to a HERMES CDF file.

        Parameters
        ----------
        output_path : `str`, optional
            A string path to the directory where file is to be saved.
            If not provided, saves to the current directory.
        overwrite : `bool`
            If set, overwrites existing file of the same name.
        Returns
        -------
        path : `str`
            A path to the saved file.
        """
        handler = CDFHandler()
        if not output_path:
            output_path = str(Path.cwd())
        if overwrite:
            cdf_file_path = Path(output_path) / (self.meta["Logical_file_id"] + ".cdf")
            if cdf_file_path.exists():
                cdf_file_path.unlink()
        return handler.save_data(data=self, file_path=output_path)

    @classmethod
    def load(cls, file_path: str):
        """
        Load data from a file.

        Parameters
        ----------
        file_path : `str`
            A fully specificed file path.

        Returns
        -------
        data : `HermesData`
            A `HermesData` object containing the loaded data.

        Raises
        ------
        ValueError: If the file type is not recognized as a file type that can be loaded.

        """
        # Determine the file type
        file_extension = Path(file_path).suffix

        # Create the appropriate handler object based on file type
        if file_extension == ".cdf":
            handler = CDFHandler()
        else:
            raise ValueError(f"Unsupported file type: {file_extension}")

        # Load data using the handler and return a HermesData object
        timeseries, support, spectra = handler.load_data(file_path)
        return cls(timeseries=timeseries, support=support, spectra=spectra)
