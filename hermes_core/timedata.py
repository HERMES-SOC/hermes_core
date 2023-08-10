"""
Container class for Measurement Data.
"""

from pathlib import Path
from collections import OrderedDict
import numpy as np
from astropy.time import Time
from astropy.timeseries import TimeSeries
from astropy.table import vstack, Column
from astropy import units as u
import hermes_core
from hermes_core.util.io import CDFHandler, NetCDFHandler
from hermes_core.util.schema import CDFSchema
from hermes_core.util.exceptions import warn_user
from hermes_core.util.util import VALID_DATA_LEVELS

__all__ = ["TimeData"]


class TimeData:
    """
    A generic object for loading, storing, and manipulating HERMES time series data.

    Parameters
    ----------
    data :  `astropy.timeseries.TimeSeries`
        The time series of data. Columns must be `~astropy.units.Quantity` arrays.
    meta : `dict`, optional
        The metadata describing the time series in an ISTP-compliant format.

    Examples
    --------
    >>> import astropy.units as u
    >>> from astropy.timeseries import TimeSeries
    >>> from hermes_core.timedata import TimeData
    >>> data = u.Quantity([1, 2, 3, 4], "gauss", dtype=np.uint16)
    >>> ts = TimeSeries(time_start="2016-03-22T12:30:31", time_delta=3 * u.s, data={"Bx": data})
    >>> input_attrs = TimeData.global_attribute_template("eea", "l1", "1.0.0")
    >>> timedata = TimeData(data=ts, meta=input_attrs)

    Raises
    ------
    ValueError: If the number of columns is less than 2 or the required 'time' column is missing.
    TypeError: If any column, excluding 'time', is not an astropy.Quantity object with units.
    ValueError: If the elements of a column are multidimensional

    References
    ----------
    * `Astropy TimeSeries <https://docs.astropy.org/en/stable/timeseries/index.html/>`_
    * `Astropy Quantity and Units <https://docs.astropy.org/en/stable/units/index.html>`_
    * `Astropy Time <https://docs.astropy.org/en/stable/time/index.html>`_
    * `Space Physics Guidelines for CDF (ISTP) <https://spdf.gsfc.nasa.gov/istp_guide/istp_guide.html>`_
    """

    def __init__(self, data, support_data=None, nrv_data=None, meta=None):
        # Verify TimeSeries compliance
        if not isinstance(data, TimeSeries):
            raise TypeError("Data must be a TimeSeries object.")
        if len(data.columns) < 2:
            raise ValueError("Data must have at least 2 columns")

        # Check individual Columns
        for colname in data.columns:
            # Verify that all Measurements are `Quantity`
            if colname != "time" and not isinstance(data[colname], u.Quantity):
                raise TypeError(
                    f"Column '{colname}' must be an astropy.Quantity object"
                )

        # Copy the TimeSeries
        self._data = TimeSeries(data, copy=True)

        # Add Input Metadata
        if meta is not None and isinstance(meta, dict):
            self._data.meta.update(meta)

        # Add any Metadata from the original TimeSeries
        self._data["time"].meta = OrderedDict()
        if hasattr(data["time"], "meta"):
            self._data["time"].meta.update(data["time"].meta)

        # Add Measurement Metadata
        for col in self._data.columns:
            if col != "time":
                self._data[col].meta = self.measurement_attribute_template()
                if hasattr(data[col], "meta"):
                    self._data[col].meta.update(data[col].meta)

        # Copy the Support Data
        if support_data:
            self.support_data = support_data
        else:
            self.support_data = {}

        # Copy the Non-Record Varying Data
        if nrv_data:
            self.nrv_data = nrv_data
        else:
            self.nrv_data = {}

        # Derive Metadata
        self.schema = CDFSchema()
        self._derive_metadata()

    @property
    def data(self):
        """
        (`astropy.timeseries.TimeSeries`) A `TimeSeries` representing one or more measurements as a function of time.
        """
        return self._data

    @property
    def meta(self):
        """
        (`collections.OrderedDict`) Global metadata associated with the measurement data.
        """
        return self._data.meta

    @property
    def units(self):
        """
        (`collections.OrderedDict`) The units of the measurement for each column in the `TimeSeries` table.
        """
        units = {}
        for name in self._data.columns:
            var_data = self._data[name]
            # Get the Unit
            if hasattr(var_data, "unit"):
                unit = var_data.unit
            else:
                unit = var_data.meta["UNITS"]
            units[name] = unit
        return OrderedDict(units)

    @property
    def columns(self):
        """
        (`list`) A list of all the names of the columns in data.
        """
        return self._data.colnames

    @property
    def time(self):
        """
        (`astropy.time.Time`) The times of the measurements.
        """
        t = Time(self._data.time)
        # Set time format to enable plotting with astropy.visualisation.time_support()
        t.format = "iso"
        return t

    @property
    def time_range(self):
        """
        (`tuple`) The start and end times of the times.
        """
        return (self._data.time.min(), self._data.time.max())

    @property
    def shape(self):
        """
        (`tuple`) The shape of the data, a tuple (nrows, ncols) including time
        """
        nrows = self._data.time.shape[0]
        ncols = len(self._data.columns)
        return (nrows, ncols)

    def __repr__(self):
        """
        Returns a representation of the `TimeData` class.
        """
        return self.__str__()

    def __str__(self):
        """
        Returns a string representation of the `TimeData` class.
        """
        str_repr = f"TimeData() Object:\n"
        # Global Attributes/Metedata
        str_repr += f"Global Attrs:\n"
        for attr_name, attr_value in self._data.meta.items():
            str_repr += f"\t{attr_name}: {attr_value}\n"
        # Measurement Data
        str_repr += f"Measurement Data:\n{self._data}\n"
        return str_repr

    def __len__(self):
        """
        Function to get the number of measurements.
        """
        return len(self._data.keys())

    def __getitem__(self, name):
        """
        Function to get a measurement.
        """
        if name in self._data.colnames:
            # Get the Data and Attrs for the named measurement
            var_data = self._data[name]
        elif name in self.support_data:
            var_data = self.support_data[name]
        elif name in self.nrv_data:
            var_data = self.nrv_data[name]
        else:
            raise KeyError(f"Can't find data measurement {name}")
        return var_data

    def __setitem__(self, name, data):
        """
        Function to add a new measurement.

        """
        # Set the Data for the named measurement
        self.add_measurement(measure_name=name, data=data, meta={})

    def __contains__(self, name):
        """
        Function to see whether a measurement is in the class.
        """
        return (
            name in self._data.columns
            or name in self.support_data
            or name in self.nrv_data
        )

    def __iter__(self):
        """
        Function to iterate over data measurements and attributes.
        """
        for name in self._data.columns:
            var_data = self._data[name]

            yield (name, var_data)

    @staticmethod
    def global_attribute_template(instr_name="", data_level="", version=""):
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
        meta = CDFSchema.global_attribute_template()

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
    def measurement_attribute_template():
        """
        Function to generate a template of the required measurement attributes.

        Returns
        -------
        template : `collections.OrderedDict`
            A template for required variable attributes that must be provided.
        """
        return CDFSchema.measurement_attribute_template()

    def _derive_metadata(self):
        """
        Funtion to derive global and measurement metadata based on a CDFSchema
        """

        # Get Default Metadata
        for attr_name, attr_value in self.schema.default_global_attributes.items():
            self._update_global_attribute(attr_name, attr_value)

        # Global Attributes
        for attr_name, attr_value in self.schema.derive_global_attributes(
            self._data
        ).items():
            self._update_global_attribute(attr_name, attr_value)

        # Time Measurement Attributes
        for attr_name, attr_value in self.schema.derive_time_attributes(
            self._data
        ).items():
            self._update_variable_attribute(
                var_name="time", attr_name=attr_name, attr_value=attr_value
            )

        # Other Measurement Attributes
        for col in [col for col in self._data.columns if col != "time"]:
            for attr_name, attr_value in self.schema.derive_measurement_attributes(
                self._data, col
            ).items():
                self._update_variable_attribute(
                    var_name=col, attr_name=attr_name, attr_value=attr_value
                )

        # Support Data
        for col in self.support_data:
            for attr_name, attr_value in self.schema.derive_measurement_attributes(
                self.support_data, col
            ).items():
                self._update_variable_attribute(
                    var_name=col, attr_name=attr_name, attr_value=attr_value
                )

        # NRV Data
        for col in self.nrv_data:
            for attr_name, attr_value in self.schema.derive_measurement_attributes(
                self.nrv_data, col
            ).items():
                self._update_variable_attribute(
                    var_name=col, attr_name=attr_name, attr_value=attr_value
                )

    def _update_global_attribute(self, attr_name, attr_value):
        # If the attribute is set, check if we want to overwrite it
        if attr_name in self._data.meta and self._data.meta[attr_name] is not None:
            # We want to overwrite if:
            #   1) The actual value is not the derived value
            #   2) The schema marks this attribute to be overwriten
            if (
                self._data.meta[attr_name] != attr_value
                and self.schema.global_attribute_schema[attr_name]["overwrite"]
            ):
                warn_user(
                    f"Overiding Global Attribute {attr_name} : {self._data.meta[attr_name]} -> {attr_value}"
                )
                self._data.meta[attr_name] = attr_value
        # If the attribute is not set, set it
        else:
            self._data.meta[attr_name] = attr_value

    def _update_variable_attribute(self, var_name, attr_name, attr_value):
        if (
            attr_name in self[var_name].meta
            and self[var_name].meta[attr_name] is not None
        ):
            attr_schema = self.schema.variable_attribute_schema["attribute_key"][
                attr_name
            ]
            if (
                self[var_name].meta[attr_name] != attr_value
                and attr_schema["overwrite"]
            ):
                warn_user(
                    f"Overiding {var_name} Attribute {attr_name} : {self[var_name].meta[attr_name]} -> {attr_value}"
                )
                self[var_name].meta[attr_name] = attr_value
        else:
            self[var_name].meta[attr_name] = attr_value

    def add_measurement(self, measure_name: str, data: u.Quantity, meta: dict = None):
        """
        Add a new measurement (column).

        Parameters
        ----------
        measure_name: `str`
            Name of the measurement to add.
        data: `astropy.units.Quantity`, `astropy.table.Column`
            The data to add. Must have the same time stamps as the existing data.
        meta: `dict`, optional
            The metadata associated with the measurement.

        Raises
        ------
        TypeError: If var_data is not of type Quantity.
        """
        # Verify that all Measurements are `Quantity`
        if isinstance(data, u.Quantity) and len(data) == len(self.time):
            self._data[measure_name] = data
            # Add any Metadata from the original Quantity
            self._data[measure_name].meta = self.measurement_attribute_template()
            if hasattr(data, "meta"):
                self._data[measure_name].meta.update(data.meta)
            if meta:
                self._data[measure_name].meta.update(meta)
        # Check Support Data
        elif (isinstance(data, Column)) and len(data) == len(self.time):
            self.support_data[measure_name] = data
            # Add any Metadata Passed not in the Column
            if meta:
                self.support_data[measure_name].meta.update(meta)
        # Consider it Metadata
        elif isinstance(data, Column):
            self.nrv_data[measure_name] = data
            # Add any Metadata Passed not in the Column
            if meta:
                self.nrv_data[measure_name].meta.update(meta)
        else:
            raise TypeError(
                f"Data for Measurement {measure_name} must be an astropy.table.Column or astropy.units.Quantity."
            )

        # Derive Metadata Attributes for the Measurement
        self._derive_metadata()

    def remove_measurement(self, measure_name: str):
        """
        Remove an existing measurement (column).

        Parameters
        ----------
        measure_name: `str`
            Name of the measurement to remove.
        """
        if measure_name in self._data.columns:
            self._data.remove_column(measure_name)
        elif measure_name in self.support_data:
            self.support_data.pop(measure_name)
        elif measure_name in self.nrv_data:
            self.nrv_data.pop(measure_name)
        else:
            raise ValueError(f"Data for Measurement {measure_name} not found.")

    def plot(self, axes=None, columns=None, subplots=True, **plot_args):
        """
        Plot the data.

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
                this_ax.plot(self.time, self.data[this_col], **plot_args)
                this_ax.set_ylabel(self.data[this_col].meta["LABLAXIS"])
        else:
            axes.set_title(
                f'{self.meta["Mission_group"]} {self.meta["Descriptor"]} {self.meta["Data_level"]}'
            )
            for this_col in columns:
                axes.plot(
                    self.time,
                    self.data[this_col],
                    label=self.data[this_col].meta["LABLAXIS"],
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
            columns = self.columns.copy()
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

    def append(self, data: TimeSeries):
        """
        Add additional measurements to an existing column.

        Parameters
        ----------
        data : `astropy.timeseries.TimeSeries`
            The data to be appended (rows) as a `TimeSeries` object.
        """
        # Verify TimeSeries compliance
        if not isinstance(data, TimeSeries):
            raise TypeError("Data must be a TimeSeries object.")
        if len(data.columns) < 2:
            raise ValueError("Data must have at least 2 columns")
        if len(self.data.columns) != len(data.columns):
            raise ValueError(
                (
                    f"Shape of curent TimeSeries ({self.shape}) does not match",
                    f"shape of data to add ({len(data.columns)}).",
                )
            )

        # Check individual Columns
        for colname in self.data.columns:
            if colname != "time" and not isinstance(self.data[colname], u.Quantity):
                raise TypeError(
                    f"Column '{colname}' must be an astropy.Quantity object"
                )

        # Save Metadata since it is not carried over with vstack
        metadata_holder = {col: self.data[col].meta for col in self.columns}

        # Vertically Stack the TimeSeries
        self._data = vstack([self._data, data])

        # Add Metadata back to the Stacked TimeSeries
        for col in self.columns:
            self.data[col].meta = metadata_holder[col]

        # Re-Derive Metadata
        self._derive_metadata()

    def save(self, output_path=None, overwrite=False):
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
    def load(cls, file_path):
        """
        Load data from a file.

        Parameters
        ----------
        file_path : `str`
            A fully specificed file path.

        Returns
        -------
        data : `TimeData`
            A `TimeData` object containing the loaded data.

        Raises
        ------
        ValueError: If the file type is not recognized as a file type that can be loaded.

        """
        # Determine the file type
        file_extension = Path(file_path).suffix

        # Create the appropriate handler object based on file type
        if file_extension == ".cdf":
            handler = CDFHandler()
        elif file_extension == ".nc":
            handler = NetCDFHandler()
        else:
            raise ValueError(f"Unsupported file type: {file_extension}")

        # Load data using the handler and return a TimeData object
        data, support_data, nrv_data = handler.load_data(file_path)
        return cls(data=data, support_data=support_data, nrv_data=nrv_data)
