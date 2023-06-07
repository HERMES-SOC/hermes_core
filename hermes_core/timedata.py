"""
Container class for Measurement Data.
"""

from pathlib import Path
from typing import OrderedDict
import numpy as np
from astropy.time import Time
from astropy.timeseries import TimeSeries
from astropy.table import vstack
from astropy import units as u
from hermes_core.util.io import CDFHandler, NetCDFHandler
from hermes_core.util.schema import CDFSchema
from hermes_core.util.exceptions import warn_user

__all__ = ["TimeData"]


class TimeData:
    """
    A generic object for loading, storing, and manipulating time series data for CDF files.

    Parameters
    ----------
    data : `~astropy.time.TimeSeries`
        A `astropy.time.TimeSeries` representing one or more measurements as a function of time.
    meta : `dict`, optional
        Meta data associated with the measurement data.
        Defaults to `None`.

    Attributes
    ----------
    data :
        A represntation of the undlerlying `astropy.timeseries.TimeSeries` object.
    meta : `dict`
        Metadata associated with the measurement data.
    units : `dict`
        A mapping from column names in `data` to the physical units of that measurement.
    columns : `list`
        A list of all the names of the columns in the data table.
    time : `astropy.time.Time`
        The times of the measurements.
    time_range : `tuple`
        The start and end times of the time axis.
    shape: `tuple`
        The shape of the data, a `tuple` `(nrows, ncols)`

    Examples
    --------
    >>> from hermes_core.timedata import TimeData
    >>>

    References
    ----------
    * `Astropy TimeSeries <https://docs.astropy.org/en/stable/timeseries/index.html/>`_

    """

    def __init__(self, data, meta=None, **kwargs):
        """
        Initialize a TimeData object.

        Parameters
        ----------
        data:  `astropy.time.TimeSeries`
            The time series data.
        meta: `dict`

        Raises
        ------
        ValueError: If the number of columns is less than 2 or the required 'time' column is missing.
        TypeError: If any column, excluding 'time', is not an astropy.Quantity object with units.
        """
        # Verify TimeSeries compliance
        if not isinstance(data, TimeSeries):
            raise TypeError("Data must be a TimeSeries object.")
        if len(data.columns) < 2:
            raise ValueError("Data must have at least 2 columns")
        if "time" not in data.columns:
            raise ValueError("Data must have a 'time' column")

        # Check individual Columns
        for colname in data.columns:
            if colname != "time" and not isinstance(data[colname], u.Quantity):
                raise TypeError(f"Column '{colname}' must be an astropy.Quantity object")

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
                self._data[col].meta = OrderedDict()
                if hasattr(data[col], "meta"):
                    self._data[col].meta.update(data[col].meta)

        # Derive Metadata
        self.schema = CDFSchema()
        self._derive_metadata()

    @property
    def data(self):
        """
        A `astropy.time.TimeSeries` representing one or more measurements as a function of time.
        """
        return self._data

    @property
    def meta(self):
        """
        Metadata associated with the measurement data.
        """
        return self._data.meta

    @property
    def units(self):
        """
        (OrderedDict) The units of the measurement for each column in the TimeSeries table.
        """
        units = {}
        for name in self._data.columns:
            var_data = self._data[name]

            # Get the Unit
            if hasattr(var_data, "unit"):
                unit = var_data.unit
            elif "UNITS" in var_data.meta and var_data.meta["UNITS"]:
                unit = var_data.meta["UNITS"]
            else:
                unit = None
            units[name] = unit
        return OrderedDict(units)

    @property
    def columns(self):
        """
        (List) A list of all the names of the columns in the data table.
        """
        return self._data.columns

    @property
    def time(self):
        """
        The times of the measurements.
        """
        t = Time(self._data.time)
        # Set time format to enable plotting with astropy.visualisation.time_support()
        t.format = "iso"
        return t

    @property
    def time_range(self):
        """
        The start and end times of the time axis.
        """
        return (self._data.time.min(), self._data.time.max())

    @property
    def shape(self):
        """
        The shape of the data, a tuple (nrows, ncols)
        """
        nrows = self._data.time.shape[0]
        ncols = len(self._data.columns)
        return (nrows, ncols)

    def __repr__(self):
        """
        Returns a representation of the CDFWriter class.
        """
        return self.__str__()

    def __str__(self):
        """
        Returns a string representation of the CDFWriter class.
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
        if name not in self._data.colnames:
            raise KeyError(f"Can't find data measurement {name}")
        # Get the Data and Attrs for the named measurement
        var_data = self._data[name]
        return var_data

    def __setitem__(self, name, data):
        """
        Function to add a new measurement.

        """
        # Set the Data for the named measurement
        self.add_measurement(measure_name=name, measure_data=data, measure_meta={})

    def __contains__(self, name):
        """
        Function to see whether a measurement is in the class.
        """
        return name in self._data.columns

    def __iter__(self):
        """
        Function to iterate over data measurements and attributes.
        """
        for name in self._data.columns:
            var_data = self._data[name]

            yield (name, var_data)

    @staticmethod
    def global_attribute_template():
        """
        Function to generate a template of the required global attributes
        that must be set for a valid CDF file to be generated from the TimeData container.

        Returns
        -------
        template : `OrderedDict`
            A template for required global attributes that must be provided.
        """
        return CDFSchema.global_attribute_template()

    @staticmethod
    def measurement_attribute_template():
        """
        Function to generate a template of the required measurement attributes
        that must be set for a valid CDF file to be generated from the TimeData container.

        Returns
        -------
        template : `OrderedDict`
            A template for required variable attributes that must be provided.
        """
        return CDFSchema.measurement_attribute_template()

    def _derive_metadata(self):
        """
        Funtion to derive global and measurement metadata based on a CDFSchema
        """

        # Get Default Metadata
        for attr_name, attr_value in self.schema.default_global_attributes.items():
            # If the attribute is set, check if we want to override it
            if attr_name in self._data.meta:
                # We want to override if:
                #   1) The actual value is not the derived value
                #   2) The schema marks this attribute to be overriden
                if (
                    self._data.meta[attr_name] != attr_value
                    and self.schema.global_attribute_schema[attr_name]["override"]
                ):
                    warn_user(
                        f"Overiding Global Attribute {attr_name} : {self._data.meta[attr_name]} -> {attr_value}"
                    )
                    self._data.meta[attr_name] = attr_value
            # If the attribute is not set, set it
            else:
                self._data.meta[attr_name] = attr_value

        # Global Attributes
        for attr_name, attr_value in self.schema.derive_global_attributes(self._data).items():
            if attr_name in self._data.meta:
                if (
                    self._data.meta[attr_name] != attr_value
                    and self.schema.global_attribute_schema[attr_name]["override"]
                ):
                    warn_user(
                        f"Overiding Global Attribute {attr_name} : {self._data.meta[attr_name]} -> {attr_value}"
                    )
                    self._data.meta[attr_name] = attr_value
            else:
                self._data.meta[attr_name] = attr_value

        # Time Measurement Attributes
        for attr_name, attr_value in self.schema.derive_time_attributes(self._data).items():
            if attr_name in self._data["time"].meta:
                attr_schema = self.schema.variable_attribute_schema["attribute_key"][attr_name]
                if self._data["time"].meta[attr_name] != attr_value and attr_schema["override"]:
                    warn_user(
                        f"Overiding Time Attribute {attr_name} : {self._data['time'].meta[attr_name]} -> {attr_value}"
                    )
                    self._data["time"].meta[attr_name] = attr_value
            else:
                self._data["time"].meta[attr_name] = attr_value

        # Other Measurement Attributes
        for col in [col for col in self._data.columns if col != "time"]:
            for attr_name, attr_value in self.schema.derive_measurement_attributes(
                self._data, col
            ).items():
                if attr_name in self._data[col].meta:
                    attr_schema = self.schema.variable_attribute_schema["attribute_key"][attr_name]
                    if self._data[col].meta[attr_name] != attr_value and attr_schema["override"]:
                        warn_user(
                            f"Overiding Measurement Attribute {attr_name} : {self._data[col].meta[attr_name]} -> {attr_value}"
                        )
                        self._data[col].meta[attr_name] = attr_value
                else:
                    self._data[col].meta[attr_name] = attr_value

    def add_measurement(self, measure_name: str, measure_data: u.Quantity, measure_meta: dict):
        """
        Function to add a new measurement.

        Parameters
        ----------
        measure_name: `str`
            Name of the measurement to add.

        measure_data: `Quantity`
            The data to add.

        measure_meta: `dict`
            The metadata associated with the measurement.

        Raises
        ------
        TypeError: If var_data is not of type Quantity.

        """
        # Verify that all Measurements are `Quantity`
        if (not isinstance(measure_data, u.Quantity)) or (not measure_data.unit):
            raise TypeError(
                f"Measurement {measure_name} must be type `astropy.units.Quantity` and have `unit` assigned."
            )

        self._data[measure_name] = measure_data
        # Add any Metadata from the original Quantity
        self._data[measure_name].meta = OrderedDict()
        if hasattr(measure_data, "meta"):
            self._data[measure_name].meta.update(measure_data.meta)
        self._data[measure_name].meta.update(measure_meta)

        # Derive Metadata Attributes for the Measurement
        self._derive_metadata()

    def remove_measurement(self, measure_name):
        """
        Remove a measuremt or measurement.
        """
        self._data.remove_column(measure_name)

    def plot(self, axes=None, columns=None, **plot_args):
        """
        Plot a plot of the data.

        Parameters
        ----------
        axes : `~matplotlib.axes.Axes`, optional
            If provided the image will be plotted on the given axes.
            Defaults to `None`, so the current axes will be used.
        columns : list[str], optional
            If provided, only plot the specified measurements.
        **plot_args : `dict`, optional
            Additional plot keyword arguments that are handed to
            :meth:`pandas.DataFrame.plot`.

        Returns
        -------
        `~matplotlib.axes.Axes`
            The plot axes.
        """
        axes, columns = self._setup_axes_columns(axes, columns)

        axes = self._data[columns].plot(ax=axes, **plot_args)

        units = set([self.units[col] for col in columns])
        if len(units) == 1:
            # If units of all columns being plotted are the same, add a unit
            # label to the y-axis.
            unit = u.Unit(list(units)[0])
            axes.set_ylabel(unit.to_string())

        self._setup_x_axis(axes)
        return axes

    def _setup_axes_columns(self, axes, columns, *, subplots=False):
        """
        Validate data for plotting, and get default axes/columns if not passed
        by the user.
        """
        # code from SunPy
        import matplotlib.pyplot as plt

        if columns is None:
            columns = self.columns
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
        """
        # code from SunPy
        import matplotlib.dates as mdates

        if isinstance(ax, np.ndarray):
            ax = ax[-1]

        locator = ax.xaxis.get_major_locator()
        formatter = ax.xaxis.get_major_formatter()
        if isinstance(formatter, mdates.AutoDateFormatter):
            # Override Matplotlib default date formatter (concise one is better)
            # but don't override any formatters pandas might have set
            ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))

    def append(self, data):
        """
        Function to add TimeSeries data to the end of the current data containers TimeSeries table.

        Parameters
        ----------
        data : `~astropy.time.TimeSeries`
            The data to be appended as a TimeSeries object.
        """
        # Verify TimeSeries compliance
        if not isinstance(data, TimeSeries):
            raise TypeError("Data must be a TimeSeries object.")
        if len(data.columns) < 2:
            raise ValueError("Data must have at least 2 columns")
        if "time" not in data.columns:
            raise ValueError("Data must have a 'time' column")
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
                raise TypeError(f"Column '{colname}' must be an astropy.Quantity object")

        # Save Metadata since it is not carried over with vstack
        metadata_holder = {col: self.data[col].meta for col in self.columns}

        # Vertically Stack the TimeSeries
        self._data = vstack([self._data, data])

        # Add Metadata back to the Stacked TimeSeries
        for col in self.columns:
            self.data[col].meta = metadata_holder[col]

        # Re-Derive Metadata
        self._derive_metadata()

    def save(self, output_path):
        """
        Save the data to a file using the specified handler.

        Parameters
        ----------
        output_path : `str`
            A string path to the directory where file is to be saved.

        Returns
        -------
        path : `str`
            A path to the saved file.

        Raises
        ------
        ValueError: If no handler is specified for saving data.

        """
        handler = CDFHandler()
        return handler.save_data(data=self, file_path=output_path)

    @classmethod
    def load(cls, file_path):
        """
        Load data from a file using the specified handler.

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
        ValueError: If the handler is not an instance of TimeDataIOHandler.

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
        data = handler.load_data(file_path)
        return cls(data)
