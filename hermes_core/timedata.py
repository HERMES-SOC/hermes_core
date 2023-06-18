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
import hermes_core
from hermes_core.util.io import CDFHandler, NetCDFHandler
from hermes_core.util.schema import CDFSchema
from hermes_core.util.exceptions import warn_user


__all__ = ["TimeData"]


class TimeData:
    """
    A generic object for loading, storing, and manipulating time series data for CDF files. This
    Data Container uses the `astropy.timeseries.TimeSeries` data structure to hold measurement data
    in addition to global and variable metadata for each masurement.

    Examples
    --------
    >>> from hermes_core.timedata import TimeData
    >>> time_series = ... # Define your `TimeSeries` here
    >>> global_meta = ... # Define your global metadata `dict` here
    >>> time_data = TimeData(data=time_series, meta=global_meta)

    References
    ----------
    * `Astropy TimeSeries <https://docs.astropy.org/en/stable/timeseries/index.html/>`_

    """

    def __init__(self, data, meta=None, **kwargs):
        """
        Initialize a TimeData object.

        Parameters
        ----------
        data:  `astropy.timeseries.TimeSeries`
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
            # Verify that all Measurements are `Quantity`
            if colname != "time" and not isinstance(data[colname], u.Quantity):
                raise TypeError(f"Column '{colname}' must be an astropy.Quantity object")
            # Verify that the Column is only a single dimension
            if len(data[colname].shape) > 1:  # If there is more than 1 Dimension
                raise ValueError(
                    f"Column '{colname}' must be a one-dimensional measurement. Split additional dimensions into unique measurenents."
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
                # TODO the following may belong in _derive_metadata
                self._data[col].meta = self.measurement_attribute_template()
                self._data[col].meta['LABLAXIS'] = f'{col} [{self._data[col].unit}]'
                self._data[col].meta['DISPLAY_TYPE'] = 'time_series'
                self._data[col].meta['VAR_TYPE'] = 'data'
                if hasattr(data[col], "meta"):
                    self._data[col].meta.update(data[col].meta)

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
            elif "UNITS" in var_data.meta and var_data.meta["UNITS"]:
                unit = var_data.meta["UNITS"]
            else:
                unit = None
            units[name] = unit
        return OrderedDict(units)

    @property
    def columns(self):
        """
        (`list`) A list of all the names of the columns in the data table.
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
        (`tuple`) The start and end times of the time axis.
        """
        return (self._data.time.min(), self._data.time.max())

    @property
    def shape(self):
        """
        (`tuple`) The shape of the data, a tuple (nrows, ncols)
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
    def global_attribute_template(instr_name, data_level, version):
        """
        Function to generate a template of the required global attributes
        that must be set for a valid CDF file to be generated from the `TimeData` container.

        Parameters
        ----------
        instr_name : str
            The instrument name. Must be "eea", "nemisis", "merit" or "spani".
        data_level : str
            The data level of the data. Must be "l1", "ql", "l3", "l4"
        version : str
            Must be of the form X.Y.Z.

        Returns
        -------
        template : `collections.OrderedDict`
            A template for required global attributes that must be provided.
        """
        meta = CDFSchema.global_attribute_template()

        meta['Descriptor'] = f"{instr_name.upper()}>{hermes_core.INST_TO_FULLNAME[instr_name]}"
        if data_level != 'ql':
            meta['Data_level'] = f"{data_level.upper()}>Level {data_level[1]}"
        else:
            meta['Data_level'] = f"{data_level.upper()}>Quicklook"
        meta['Data_version'] = version
        return meta

    @staticmethod
    def measurement_attribute_template():
        """
        Function to generate a template of the required measurement attributes
        that must be set for a valid CDF file to be generated from the `TimeData` container.

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
        Function to add a new measurement (column) to the `TimeData` data container.

        Parameters
        ----------
        measure_name: `str`
            Name of the measurement to add.

        measure_data: `astropy.units.Quantity`
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
        # Verify that the Column is only a single dimension
        if len(measure_data.shape) > 1:  # If there is more than 1 Dimension
            raise ValueError(
                f"Column '{measure_name}' must be a one-dimensional measurement. Split additional dimensions into unique measurenents."
            )

        self._data[measure_name] = measure_data
        # Add any Metadata from the original Quantity
        self._data[measure_name].meta = self.measurement_attribute_template()
        # TODO the following may belong in _derive_metadata
        self._data[measure_name].meta['LABLAXIS'] = f'{measure_name} [{measure_data.unit}]'
        self._data[measure_name].meta['DISPLAY_TYPE'] = 'time_series'
        self._data[measure_name].meta['VAR_TYPE'] = 'data'
        if hasattr(measure_data, "meta"):
            self._data[measure_name].meta.update(measure_data.meta)
        self._data[measure_name].meta.update(measure_meta)

        # Derive Metadata Attributes for the Measurement
        self._derive_metadata()

    def remove_measurement(self, measure_name: str):
        """
        Function to remove a measurement (column).

        Parameters
        ----------
        measure_name: `str`
            Name of the measurement to remove.
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
        columns : `list[str]`, optional
            If provided, only plot the specified measurements otherwise try to plot them all.
        **plot_args : `dict`, optional
            Additional plot keyword arguments that are handed to
            `~matplotlib.axes.Axes`.
        Returns
        -------
        `~matplotlib.axes.Axes`
            The plot axes.
        """
        # Set up the plot axes based on the number of columns to plot
        axes, columns = self._setup_axes_columns(axes, columns)
        from astropy.visualization import quantity_support
        with quantity_support():
            for col in columns:
                axes.plot(self.time, self._data[columns], **plot_args)

        #units = set([self.units[col] for col in columns])
        #if len(units) == 1:
            # If units of all columns being plotted are the same, add a unit
            # label to the y-axis.
        #    unit = u.Unit(list(units)[0])
        #    axes.set_ylabel(unit.to_string())

        # Setup the Time Axis
        self._setup_x_axis(axes)

        return axes

    def _setup_axes_columns(self, axes, columns, *, subplots=False):
        """
        Validate data for plotting, and get default axes/columns if not passed
        by the user.
        """
        import matplotlib.pyplot as plt

        # If no individual columns were input, try to plot all columns
        if columns is None:
            columns = self.columns
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
        """
        import matplotlib.dates as mdates

        if isinstance(ax, np.ndarray):
            ax = ax[-1]

        locator = ax.xaxis.get_major_locator()
        formatter = ax.xaxis.get_major_formatter()
        if isinstance(formatter, mdates.AutoDateFormatter):
            # Override Matplotlib default date formatter (concise one is better)
            # but don't override any formatters pandas might have set
            ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))

    def append(self, data: TimeSeries):
        """
        Function to add `TimeSeries` data to the end of the current data containers `TimeSeries` table.

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

    def save(self, output_path=None, overwrite=False):
        """
        Save the data to a CDF file in the directory specified by the `output_path`.

        Parameters
        ----------
        output_path : `str`, optional
            A string path to the directory where file is to be saved. If not provided, saves to the current directory.
        overwrite : `bool`
            If set, overwrite existing file of the same name.
        Returns
        -------
        path : `str`
            A path to the saved file.
        """
        handler = CDFHandler()
        if not output_path:
            output_path = str(Path.cwd())
        if overwrite:
            (Path(output_path) / (self.meta['Logical_file_id'] + '.cdf')).unlink()
        return handler.save_data(data=self, file_path=output_path)

    @classmethod
    def load(cls, file_path):
        """
        Load data from a file using one of the provided I/O handlers for recognized file types.

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
        data = handler.load_data(file_path)
        return cls(data)
