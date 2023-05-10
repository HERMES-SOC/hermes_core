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
from hermes_core.util.validation import CDFValidator, NetCDFValidator, FITSValidator
from hermes_core.util.io import TimeDataIOHandler, CDFHandler, NetCDFHandler, FITSHandler

__all__ = ["TimeData", "CDFMeta"]


def read(filepath):
    """
    A generic file reader.

    Parameters
    ----------
    filepath : `str`
        A fully specificed file path.

    Returns
    -------
        result : TimeData
    """
    # Determine the file type
    file_extension = Path(filepath).suffix

    # Create the appropriate handler object based on file type
    if file_extension == ".cdf":
        handler = CDFHandler()
    elif file_extension == ".nc":
        handler = NetCDFHandler()
    elif file_extension == ".fits":
        handler = FITSHandler()
    else:
        raise ValueError(f"Unsupported file type: {file_extension}")

    # Load data using the handler and return a TimeData object
    return TimeData.load(filepath, handler=handler)


def validate(filepath):
    """
    Validate a data file such as a CDF.

    Parameters
    ----------
    filepath : `str`
        A fully specificed file path.

    Returns
    -------
        result
    """
    # Determine the file type
    file_extension = Path(filepath).suffix

    # Create the appropriate validator object based on file type
    if file_extension == ".cdf":
        validator = CDFValidator()
    elif file_extension == ".nc":
        validator = NetCDFValidator()
    elif file_extension == ".fits":
        validator = FITSValidator()
    else:
        raise ValueError(f"Unsupported file type: {file_extension}")

    # Call the validate method of the validator object
    return validator.validate(filepath)


class TimeData:
    """
    A generic object for loading, storing, and manipulating time series data for cdf files.

    Parameters
    ----------
    data : `~astropy.time.TimeSeries`
        A `astropy.time.TimeSeries` representing one or more measurements as a function of time.
    meta : `dict`, optional
        Meta data associated with the measurement data.
        Defaults to `None`.
    units : `dict`, optional
        A mapping from column names in ``data`` to the physical units of that column.
        Defaults to `None`.
    handler (TimeDataIOHandler, optional): The handler for input/output operations. Defaults to None.


    Attributes
    ----------
    meta : `dict`
        Metadata associated with the measurement data.
    units : `dict`
        A mapping from column names in ``data`` to the physical units of that column.

    Examples
    --------
    >>> from hermes_core.timedata import TimeData
    >>>

    References
    ----------
    * `Astropy TimeSeries <https://docs.astropy.org/en/stable/timeseries/index.html/>`_

    """

    def __init__(self, data, meta=None, handler=None, **kwargs):
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
            raise ValueError("Science data must have at least 2 columns")
        if "time" not in data.columns:
            raise ValueError("Science data must have a 'time' column")

        # Check individual Columns
        for colname in data.columns:
            if colname != "time" and not isinstance(data[colname], u.Quantity):
                raise TypeError(f"Column '{colname}' must be an astropy.Quantity object")

        # Copy the TimeSeries
        self._data = TimeSeries(data, copy=True)

        # Add any Metadata from the original TimeSeries
        self._data.time.meta = OrderedDict()
        if hasattr(data["time"], "meta"):
            self._data.time.meta.update(data["time"].meta)
        for col in self._data.columns:
            if col != "time":
                self._data[col].meta = OrderedDict()
                if hasattr(data[col], "meta"):
                    self._data[col].meta.update(data[col].meta)

        # Add any Metadata from the original TimeSeries
        self.data.time.meta = OrderedDict()
        if hasattr(data["time"], "meta"):
            self.data.time.meta.update(data["time"].meta)
        self.data.time.meta.update(meta["time"])

        # Add Measurement Metadata
        for col in self.data.columns:
            if col != "time":
                self.data[col].meta = OrderedDict()
                if hasattr(data[col], "meta"):
                    self.data[col].meta.update(data[col].meta)
                self.data[col].meta.update(meta[col])

        self.handler = handler

        # Derive and add CDF metadata
        # the following is going to potentially overwrite metadata
        # we should just warn and just do it
        self._data.time.meta = CDFMeta.derive_time_attributes(self._data.time)
        for col in self._data.columns:
            self._data[col].meta.update(CDFMeta.derive_variable_attributes(self._data[col]))

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

        Returns
        -------
            dict: The metadata associated with the science data.
        """
        return self._data.meta

    @meta.setter
    def meta(self, value):
        """
        Set the metadata associated with the measurement data.

        Parameters:
            value (dict): The metadata to set.
        """
        self._data.meta = value

    @property
    def handler(self):
        """
        The handler for input/output operations.

        Returns:
            TimeDataIOHandler: The handler for input/output operations.

        """
        return self._handler

    @handler.setter
    def handler(self, value):
        """
        Set the handler for input/output operations.

        Parameters:
            value (TimeDataIOHandler): The handler to set.

        Raises:
            ValueError: If the handler is not an instance of TimeDataIOHandler.

        """
        if value is not None and not isinstance(value, TimeDataIOHandler):
            raise ValueError("Handler must be an instance of TimeDataIOHandler")
        self._handler = value

    @property
    def units(self):
        """
        (OrderedDict) The units of the columns of data.
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
        (List) A list of all the names of the columns in the data.
        """
        return self._data.columns

    @property
    def time(self):
        """
        The times of the measurements.
        """
        t = Time(self._data.index)
        # Set time format to enable plotting with astropy.visualisation.time_support()
        t.format = "iso"
        return t

    @property
    def shape(self):
        """
        The shape of the data, a tuple (nrows, ncols)
        """
        if "time" in self._data.columns:
            nrows = self._data.time.shape[0]
        else:
            nrows = 0
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
        # Variable Data
        str_repr += f"Variable Data:\n{self._data}\n"
        # Variable Attributes
        str_repr += f"Variable Attributes:\n"
        for col_name in self._data.columns:
            str_repr += f"\tVar: {col_name}\n"
            # for attr_name, attr_value in self._data[col_name].meta.items():

            #     str_repr += f"\t\t{attr_name}: {attr_value}\n"
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
        # Get the Data and Attrs for the named variable
        var_data = self._data[name]
        return var_data

    def __setitem__(self, name, data):
        """
        Function to add a new measurement.

        """
        # Set the Data for the named variable
        self.add_measurement(measure_name=name, measure_data=data, measure_meta={})

    def __contains__(self, name):
        """
        Function to see whether a variable is in the class.
        """
        return name in self._data.columns

    def __iter__(self):
        """
        Function to iterate over data variables and attributes.
        """
        for name in self._data.columns:
            var_data = self._data[name]

            yield (name, var_data)

    def add_measurement(self, measure_name: str, measure_data: u.Quantity, measure_meta: dict):
        """
        Function to add a new measurement or column.

        Parameters
        ----------
        measure_name: `str`
            Name of the measurement or column to add.

        measure_data: `Quantity` array
            The data to add.

        measure_meta: `dict`
            The metadata associated with the measurement.

        Raises
        ------
        TypeError: If var_data is not of type Quantity.

        """
        # Verify that all columns are `Quantity`
        if (not isinstance(measure_data, u.Quantity)) or (not measure_data.unit):
            raise TypeError(
                f"Column {measure_name} must be type `astropy.units.Quantity` and have `unit` assigned."
            )

        self.data[measure_name] = measure_data
        # Add any Metadata from the original Quantity
        self.data[measure_name].meta = OrderedDict()
        if hasattr(measure_data, "meta"):
            self.data[measure_name].meta.update(measure_data.meta)
        self.data[measure_name].meta.update(measure_meta)

    def remove_measurement(self, measure_name):
        """
        Remove a measuremt or column.
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

    @property
    def time_range(self):
        """
        The start and end times of the time axis.
        """
        return (self._data.index.min(), self._data.index.max())

    def to_cdf(self, filepath):
        """
        Write the data to a CDF.

        Parameters
        ----------
        filepath : `str`
            Fully specificied filepath of the output file.
        """
        pass

    def append(self, data):
        """
        Function to add TimeSeries data to the end of the current data containers TimeSeries table.

        Parameters:
            data (TimeSeries): The science data appended as a TimeSeries object.
        """
        # Verify TimeSeries compliance
        if not isinstance(data, TimeSeries):
            raise TypeError("Data must be a TimeSeries object.")
        if len(data.columns) < 2:
            raise ValueError("Science data must have at least 2 columns")
        if "time" not in data.columns:
            raise ValueError("Science data must have a 'time' column")
        if len(self.data.columns) != len(data.columns):
            raise ValueError(
                (
                    f"Shape of curent TimeSeries ({self.shape}) does not match",
                    "shape of data to add ({len(data.columns)}).",
                )
            )

        # Check individual Columns
        for colname in self.data.columns:
            if colname != "time" and not isinstance(self.data[colname], u.Quantity):
                raise TypeError(f"Column '{colname}' must be an astropy.Quantity object")

        # Save Metadata since it is not carried over with vstack
        metadata_holder = {col: self.data[col].meta for col in self.columns}

        # Vertically Stack the TimeSeries
        self.data = vstack([self.data, data])

        # Add Metadata back to the Stacked TimeSeries
        for col in self.columns:
            self.data[col].meta = metadata_holder[col]

    def save(self, output_path):
        """
        Save the science data to a file using the specified handler.

        Parameters:
            output_path (str): A string path to the directory where file is to be saved.

        Returns:
            str: A path to the saved file.

        Raises:
            ValueError: If no handler is specified for saving data.

        """
        if self.handler is None:
            raise ValueError("No handler specified for saving data")
        return self.handler.save_data(data=self, file_path=output_path)

    @classmethod
    def load(cls, file_path, handler):
        """
        Load science data from a file using the specified handler.

        Parameters:
            file_path (str): A fully specificed file path.
            handler (TimeDataIOHandler): The handler for input/output operations.

        Returns:
            TimeData: A TimeData object containing the loaded data.

        Raises:
            ValueError: If the handler is not an instance of TimeDataIOHandler.

        """
        if not isinstance(handler, TimeDataIOHandler):
            raise ValueError("Handler must be an instance of TimeDataIOHandler")
        data = handler.load_data(file_path)
        return cls(data, handler)


def CDFMeta():
    """An object to derive cdf meta data that is ISTP compliant."""

    def derive_global_attributes(self):
        """Function to derive global attributes"""
        # Loop through Global Attributes
        for attr_name, attr_schema in self._global_attr_schema.items():
            if attr_schema["derived"]:
                derived_value = self.derive_gattribute(attr_name=attr_name)
                # Only Derive Global Attributes if they have not been manually derived/overridden
                if (attr_name not in self.data.meta) or (not self.data.meta[attr_name]):
                    self.data.meta[attr_name] = derived_value
                else:
                    log.debug(
                        (
                            "Attribute: %s was marked for derivation (to be %s)"
                            "but was already overridden to %s"
                        ),
                        attr_name,
                        derived_value,
                        self.data.meta[attr_name],
                    )

    def derive_variable_attributes(column):
        """
        Given a TimeSeries column, derive CDF variable attributes
        """
        meta = OrderedDict()

        # Check the Attributes that can be derived
        meta["DEPEND_0"] = self._get_depend()
        meta["FIELDNAM"] = self._get_fieldnam(column)
        meta["FILLVAL"] = self._get_fillval(column)
        meta["FORMAT"] = self._get_format(column)
        meta["SI_CONVERSION"] = self._get_si_conversion(column)
        meta["UNITS"] = self._get_units(column)
        meta["VALIDMIN"] = self._get_validmin(column)
        meta["VALIDMAX"] = self._get_validmax(column)

        return meta

    def derive_time_attributes(time_column):
        """
        Given a TimeSeries time column, derive CDF time variable attributes

        Returns
        -------
        meta :  `dict`
        """
        meta = OrderedDict()

        # Check the Attributes that can be derived
        meta["REFERENCE_POSITION"] = self._get_reference_position()
        meta["RESOLUTION"] = self._get_resolution()
        meta["TIME_BASE"] = self._get_time_base()
        meta["TIME_SCALE"] = self._get_time_scale()
        meta["UNITS"] = self._get_time_units()

        return meta

    def derive_gattribute(self, attr_name):
        """
        Derive attributes based on other attribues in the CDF writer's internal dict state.
        """
        # SWITCH on the Derivation attr_name
        if attr_name == "Generation_date":
            return self._get_generation_date()
        elif attr_name == "Start_time":
            return self._get_start_time()
        elif attr_name == "Data_type":
            return self._get_data_type()
        elif attr_name == "Logical_file_id":
            return self._get_logical_file_id()
        elif attr_name == "Logical_source":
            return self._get_logical_source()
        elif attr_name == "Logical_source_description":
            return self._get_logical_source_description()
        else:
            raise ValueError(f"Derivation for Attribute ({attr_name}) Not Recognized")

    # =============================================================================================
    #                               VARIABLE ATTRIBUTE DERIVATIONS
    # =============================================================================================

    def _get_depend(self):
        return "Epoch"

    def _get_fieldnam(self, var_name):
        if var_name != "time":
            return deepcopy(var_name)
        else:
            return "Epoch"

    def _get_fillval(self, var_name):
        # Get the Variable Data
        var_data = self.data[var_name]
        if var_name == "time":
            # Guess the spacepy.pycdf.const CDF Data Type
            (guess_dims, guess_types, guess_elements) = _Hyperslice.types(var_data.to_datetime())
            # Get the FILLVAL for the gussed data type
            fillval = self._fillval_helper(cdf_type=guess_types[0])
            # guess_types[0] == spacepy.pycdf.const.CDF_TIME_TT2000.value:
            return spacepy.pycdf.lib.v_tt2000_to_datetime(fillval)
        else:
            # Guess the spacepy.pycdf.const CDF Data Type
            (guess_dims, guess_types, guess_elements) = _Hyperslice.types(var_data.value)
            # Get the FILLVAL for the gussed data type
            fillval = self._fillval_helper(cdf_type=guess_types[0])
            return fillval

    def _fillval_helper(self, cdf_type):
        # Fill value, indexed by the CDF type (numeric)
        fillvals = {}
        # Integers
        for i in (1, 2, 4, 8):
            fillvals[getattr(spacepy.pycdf.const, "CDF_INT{}".format(i)).value] = -(
                2 ** (8 * i - 1)
            )
            if i == 8:
                continue
            fillvals[getattr(spacepy.pycdf.const, "CDF_UINT{}".format(i)).value] = 2 ** (8 * i) - 1
        fillvals[spacepy.pycdf.const.CDF_EPOCH16.value] = (-1e31, -1e31)
        fillvals[spacepy.pycdf.const.CDF_REAL8.value] = -1e31
        fillvals[spacepy.pycdf.const.CDF_REAL4.value] = -1e31
        fillvals[spacepy.pycdf.const.CDF_CHAR.value] = " "
        fillvals[spacepy.pycdf.const.CDF_UCHAR.value] = " "
        # Equivalent pairs
        for cdf_t, equiv in (
            (spacepy.pycdf.const.CDF_TIME_TT2000, spacepy.pycdf.const.CDF_INT8),
            (spacepy.pycdf.const.CDF_EPOCH, spacepy.pycdf.const.CDF_REAL8),
            (spacepy.pycdf.const.CDF_BYTE, spacepy.pycdf.const.CDF_INT1),
            (spacepy.pycdf.const.CDF_FLOAT, spacepy.pycdf.const.CDF_REAL4),
            (spacepy.pycdf.const.CDF_DOUBLE, spacepy.pycdf.const.CDF_REAL8),
        ):
            fillvals[cdf_t.value] = fillvals[equiv.value]
        value = fillvals[cdf_type]
        return value

    def _get_format(self, var_name):
        # Get the Variable Data
        var_data = self.data[var_name]
        if var_name == "time":
            # Guess the spacepy.pycdf.const CDF Data Type
            (guess_dims, guess_types, guess_elements) = _Hyperslice.types(var_data.to_datetime())
            return self._format_helper(var_name, guess_types[0])
        else:
            # Guess the spacepy.pycdf.const CDF Data Type
            (guess_dims, guess_types, guess_elements) = _Hyperslice.types(var_data.value)
            return self._format_helper(var_name, guess_types[0])

    def _format_helper(self, var_name, cdftype):
        minn = "VALIDMIN"
        maxx = "VALIDMAX"
        # Get the Variable Data
        var_data = self.data[var_name]

        if cdftype in (
            spacepy.pycdf.const.CDF_INT1.value,
            spacepy.pycdf.const.CDF_INT2.value,
            spacepy.pycdf.const.CDF_INT4.value,
            spacepy.pycdf.const.CDF_INT8.value,
            spacepy.pycdf.const.CDF_UINT1.value,
            spacepy.pycdf.const.CDF_UINT2.value,
            spacepy.pycdf.const.CDF_UINT4.value,
            spacepy.pycdf.const.CDF_BYTE.value,
        ):
            if minn in var_data.meta:  # Just use validmin or scalemin
                minval = var_data.meta[minn]
            elif cdftype in (
                spacepy.pycdf.const.CDF_UINT1.value,
                spacepy.pycdf.const.CDF_UINT2.value,
                spacepy.pycdf.const.CDF_UINT4.value,
            ):  # unsigned, easy
                minval = 0
            elif cdftype == spacepy.pycdf.const.CDF_BYTE.value:
                minval = -(2**7)
            else:  # Signed, harder
                size = next(
                    (
                        i
                        for i in (1, 2, 4, 8)
                        if getattr(spacepy.pycdf.const, "CDF_INT{}".format(i)).value == cdftype
                    )
                )
                minval = -(2 ** (8 * size - 1))
            if maxx in var_data.meta:  # Just use max
                maxval = var_data.meta[maxx]
            elif cdftype == spacepy.pycdf.const.CDF_BYTE.value:
                maxval = 2**7 - 1
            else:
                size = next(
                    (
                        8 * i
                        for i in (1, 2, 4)
                        if getattr(spacepy.pycdf.const, "CDF_UINT{}".format(i)).value == cdftype
                    ),
                    None,
                )
                if size is None:
                    size = (
                        next(
                            (
                                8 * i
                                for i in (1, 2, 4, 8)
                                if getattr(spacepy.pycdf.const, "CDF_INT{}".format(i)).value
                                == cdftype
                            )
                        )
                        - 1
                    )
                maxval = 2**size - 1
            # Two tricks:
            # -Truncate and add 1 rather than ceil so get
            # powers of 10 (log10(10) = 1 but needs two digits)
            # -Make sure not taking log of zero
            if minval < 0:  # Need an extra space for the negative sign
                fmt = "I{}".format(int(math.log10(max(abs(maxval), abs(minval), 1))) + 2)
            else:
                fmt = "I{}".format(int(math.log10(maxval) if maxval != 0 else 1) + 1)
        elif cdftype == spacepy.pycdf.const.CDF_TIME_TT2000.value:
            fmt = "A{}".format(len("9999-12-31T23:59:59.999999999"))
        elif cdftype == spacepy.pycdf.const.CDF_EPOCH16.value:
            fmt = "A{}".format(len("31-Dec-9999 23:59:59.999.999.000.000"))
        elif cdftype == spacepy.pycdf.const.CDF_EPOCH.value:
            fmt = "A{}".format(len("31-Dec-9999 23:59:59.999"))
        elif cdftype in (
            spacepy.pycdf.const.CDF_REAL8.value,
            spacepy.pycdf.const.CDF_REAL4.value,
            spacepy.pycdf.const.CDF_FLOAT.value,
            spacepy.pycdf.const.CDF_DOUBLE.value,
        ):
            if "VALIDMIN" in var_data.meta and "VALIDMAX" in var_data.meta:
                range = var_data.meta["VALIDMAX"] - var_data.meta["VALIDMIN"]
            # If not, just use nothing.
            else:
                range = None
            # Find how many spaces we need for the 'integer' part of the number
            # (Use maxx-minn for this...effectively uses VALIDMIN/MAX for most
            # cases.)
            if range and (minn in var_data.meta and maxx in var_data.meta):
                if len(str(int(var_data.meta[maxx]))) >= len(str(int(var_data.meta[minn]))):
                    ln = str(int(var_data.meta[maxx]))
                else:
                    ln = str(int(var_data.meta[minn]))
            if range and ln and range < 0:  # Cover all our bases:
                range = None
            # Switch on Range
            if range and ln and range <= 11:  # If range <= 11, we want 2 decimal places:
                # Need extra for '.', and 3 decimal places (4 extra)
                fmt = "F{}.3".format(len([i for i in ln]) + 4)
            elif range and ln and 11 < range <= 101:
                # Need extra for '.' (1 extra)
                fmt = "F{}.2".format(len([i for i in ln]) + 3)
            elif range and ln and 101 < range <= 1000:
                # Need extra for '.' (1 extra)
                fmt = "F{}.1".format(len([i for i in ln]) + 2)
            else:
                # No range, must not be populated, copied from REAL4/8(s) above
                # OR we don't care because it's a 'big' number:
                fmt = "G10.2E3"
        elif cdftype in (spacepy.pycdf.const.CDF_CHAR.value, spacepy.pycdf.const.CDF_UCHAR.value):
            fmt = "A{}".format(len(var_data.value))
        else:
            raise ValueError(
                "Couldn't find FORMAT for {} of type {}".format(
                    var_name, spacepy.pycdf.lib.cdftypenames.get(cdftype, "UNKNOWN")
                )
            )
        return fmt

    def _get_reference_position(self):
        # Get the Variable Data
        var_data = self.time
        # Guess the spacepy.pycdf.const CDF Data Type
        (guess_dims, guess_types, guess_elements) = _Hyperslice.types(var_data.to_datetime())
        if guess_types[0] == spacepy.pycdf.const.CDF_TIME_TT2000.value:
            return "rotating Earth geoid"
        else:
            raise TypeError(f"Reference Position for Time type ({guess_types[0]}) not found.")

    def _get_resolution(self):
        # Get the Variable Data
        var_data = self.time
        times = len(var_data)
        if times < 2:
            raise ValueError(f"Can not derive Time Resolution, need 2 samples, found {times}.")
        # Calculate the Timedelta between two datetimes
        times = var_data.to_datetime()
        delta = times[1] - times[0]
        # Get the number of seconds between samples.
        delta_seconds = delta.total_seconds()
        return f"{delta_seconds}s"

    def _get_si_conversion(self, var_name):
        # Get the Variable Data
        var_data = self.data[var_name]
        if var_name == "time":
            time_unit_str = self._get_time_units()
            time_unit = u.s
            if time_unit_str == "ms":
                time_unit = u.ms
            if time_unit_str == "ns":
                time_unit = u.ns
            if time_unit_str == "ps":
                time_unit = u.ns
            conversion_rate = time_unit.to(u.s)
            si_conversion = f"{conversion_rate:e}>{u.s}"
        else:
            conversion_rate = var_data.unit.to(var_data.si.unit)
            si_conversion = f"{conversion_rate:e}>{var_data.si.unit}"
        return si_conversion

    def _get_time_base(self):
        # Get the Variable Data
        var_data = self.time
        # Guess the spacepy.pycdf.const CDF Data Type
        (guess_dims, guess_types, guess_elements) = _Hyperslice.types(var_data.to_datetime())
        if guess_types[0] == spacepy.pycdf.const.CDF_TIME_TT2000.value:
            return "J2000"
        else:
            raise TypeError(f"Time Base for Time type ({guess_types[0]}) not found.")

    def _get_time_scale(self):
        # Get the Variable Data
        var_data = self.time
        # Guess the spacepy.pycdf.const CDF Data Type
        (guess_dims, guess_types, guess_elements) = _Hyperslice.types(var_data.to_datetime())
        if guess_types[0] == spacepy.pycdf.const.CDF_TIME_TT2000.value:
            return "Terrestrial Time (TT)"
        else:
            raise TypeError(f"Time Scale for Time type ({guess_types[0]}) not found.")

    def _get_time_units(self):
        # Get the Variable Data
        var_data = self.time
        # Guess the spacepy.pycdf.const CDF Data Type
        (guess_dims, guess_types, guess_elements) = _Hyperslice.types(var_data.to_datetime())
        if guess_types[0] == spacepy.pycdf.const.CDF_EPOCH.value:
            return "ms"
        if guess_types[0] == spacepy.pycdf.const.CDF_TIME_TT2000.value:
            return "ns"
        if guess_types[0] == spacepy.pycdf.const.CDF_EPOCH16.value:
            return "ps"
        else:
            raise TypeError(f"Time Units for Time type ({guess_types[0]}) not found.")

    def _get_units(self, var_name):
        # Get the Variable Data
        var_data = self.data[var_name]
        unit = ""
        # Get the Unit from the TimeSeries Quantity if it exists
        if hasattr(var_data, "unit"):
            unit = var_data.unit.name
        return unit

    def _get_validmin(self, var_name):
        # Get the Variable Data
        var_data = self.data[var_name]
        if var_name == "time":
            # Guess the spacepy.pycdf.const CDF Data Type
            (guess_dims, guess_types, guess_elements) = _Hyperslice.types(var_data.to_datetime())
            # Get the Min Value
            minval, maxval = spacepy.pycdf.lib.get_minmax(guess_types[0])
            return minval + datetime.timedelta(seconds=1)
        else:
            # Guess the spacepy.pycdf.const CDF Data Type
            (guess_dims, guess_types, guess_elements) = _Hyperslice.types(var_data.value)
            # Get the Min Value
            minval, maxval = spacepy.pycdf.lib.get_minmax(guess_types[0])
            return minval

    def _get_validmax(self, var_name):
        # Get the Variable Data
        var_data = self.data[var_name]
        if var_name == "time":
            # Guess the spacepy.pycdf.const CDF Data Type
            (guess_dims, guess_types, guess_elements) = _Hyperslice.types(var_data.to_datetime())
            # Get the Max Value
            minval, maxval = spacepy.pycdf.lib.get_minmax(guess_types[0])
            return maxval - datetime.timedelta(seconds=1)
        else:
            # Guess the spacepy.pycdf.const CDF Data Type
            (guess_dims, guess_types, guess_elements) = _Hyperslice.types(var_data.value)
            # Get the Max Value
            minval, maxval = spacepy.pycdf.lib.get_minmax(guess_types[0])
            return maxval

    # =============================================================================================
    #                               GLOBAL ATTRIBUTE DERIVATIONS
    # =============================================================================================

    def _get_logical_file_id(self):
        """
        Function to get the `Logical_file_id` required global attribute.

        The attribute stores the name of the CDF File without the file
        extension (e.g. '.cdf'). This attribute is requires to avoid
        loss of the originial source in case of renaming.
        """
        attr_name = "Logical_file_id"
        if (attr_name not in self.data.meta) or (not self.data.meta[attr_name]):
            # Get Parts
            instrument_id = self._get_instrument_id()
            start_time = self._get_start_time()
            data_level = self._get_data_level()
            version = self._get_version()
            mode = self._get_instrument_mode()

            # Build Derivation
            science_filename = util.create_science_filename(
                instrument=instrument_id,
                time=start_time,
                level=data_level,
                version=version,
                mode=mode,
            )
            science_filename = science_filename.rstrip(util.FILENAME_EXTENSION)
        else:
            science_filename = self.data.meta[attr_name]
        return science_filename

    def _get_logical_source(self):
        """
        Function to get the `Logical_source` required global attribute.

        This attribute determines the file naming convention in the SKT Editor
        and is used by CDA Web.
        """
        attr_name = "Logical_source"
        if (attr_name not in self.data.meta) or (not self.data.meta[attr_name]):
            # Get Parts
            spacecraft_id = self._get_spacecraft_id()
            instrument_id = self._get_instrument_id()
            data_type = self._get_data_type()
            data_type_short_name, _ = data_type.split(">")

            # Build Derivation
            logical_source = f"{spacecraft_id}_{instrument_id}_{data_type_short_name}"
        else:
            logical_source = self.data.meta[attr_name]
        return logical_source

    def _get_logical_source_description(self):
        """
        Function to get the `Logical_source_description` required global attribute.

        This attribute writes out the full words associated with the encryped
        `Logical_source`  attribute.
        """
        attr_name = "Logical_source_description"
        if (attr_name not in self.data.meta) or (not self.data.meta[attr_name]):
            # Get Parts
            spacecraft_long_name = self._get_spacecraft_long_name()
            instrument_long_name = self._get_instrument_long_name()
            data_type = self._get_data_type()
            _, data_type_long_name = data_type.split(">")
            logical_source_description = (
                f"{spacecraft_long_name} {instrument_long_name} {data_type_long_name}"
            )
        else:
            logical_source_description = self.data.meta[attr_name]
        return logical_source_description

    def _get_data_type(self):
        """
        Function to get the `Data_type` required global attribute.

        This attribute is used by the CDF Writing software to create the filename.
        It is a combination of the following components:
            - mode
            - data_level
            - optional_data_product_descriptor
        """
        attr_name = "Data_type"
        if (attr_name not in self.data.meta) or (not self.data.meta[attr_name]):
            short_parts = []
            long_parts = []

            # Get `mode`
            mode_short_name = self._get_instrument_mode()
            mode_long_name = self._get_instrument_mode()
            if bool(mode_short_name and mode_long_name):
                short_parts.append(mode_short_name)
                long_parts.append(mode_long_name)

            # Get `data level`
            data_level_short_name = self._get_data_level()
            data_level_long_name = self._get_data_level_long_name()
            if bool(data_level_short_name and data_level_long_name):
                short_parts.append(data_level_short_name)
                long_parts.append(data_level_long_name)

            # Get `data product descriptor`
            odpd_short_name = self._get_data_product_descriptor()
            odpd_long_name = self._get_data_product_descriptor()
            if bool(odpd_short_name and odpd_long_name):
                short_parts.append(odpd_short_name)
                long_parts.append(odpd_long_name)

            # Build Derivation
            data_type = "_".join(short_parts) + ">" + " ".join(long_parts)
        else:
            data_type = self.data.meta[attr_name]
        return data_type

    def _get_spacecraft_id(self):
        """Function to get Spacecraft ID from Source_name Global Attribute"""
        attr_name = "Source_name"
        if (attr_name not in self.data.meta) or (not self.data.meta[attr_name]):
            # Get Module Default
            sc_id = hermes_core.MISSION_NAME
        else:
            sc_id = self.data.meta["Source_name"]
            # Formatting
            if ">" in sc_id:
                short_name, _ = sc_id.split(">")
                sc_id = short_name.lower()  # Makse sure its all lowercase
        return sc_id

    def _get_spacecraft_long_name(self):
        """Function to get Spacecraft ID from Source_name Global Attribute"""
        attr_name = "Source_name"
        if (attr_name not in self.data.meta) or (not self.data.meta[attr_name]):
            # Get Module Default
            sc_id = hermes_core.MISSION_NAME
        else:
            sc_id = self.data.meta["Source_name"]
            # Formatting
            if ">" in sc_id:
                _, long_name = sc_id.split(">")
                sc_id = long_name
        return sc_id

    def _get_instrument_id(self):
        """
        Function to get Instrument ID from Descriptor Global Attribute

        Instrument of investigation identifier shortened to three
        letter acronym.
        """
        attr_name = "Descriptor"
        if (attr_name not in self.data.meta) or (not self.data.meta[attr_name]):
            instr_id = None
        else:
            instr_id = self.data.meta["Descriptor"]
            # Formatting
            if ">" in instr_id:
                short_name, _ = instr_id.split(">")
                instr_id = short_name.lower()  # Makse sure its all lowercase
        return instr_id

    def _get_instrument_long_name(self):
        """
        Function to get Instrument ID from Descriptor Global Attribute

        Instrument of investigation identifier shortened to three
        letter acronym.
        """
        attr_name = "Descriptor"
        if (attr_name not in self.data.meta) or (not self.data.meta[attr_name]):
            instr_id = None
        else:
            instr_id = self.data.meta["Descriptor"]
            # Formatting
            if ">" in instr_id:
                _, long_name = instr_id.split(">")
                instr_id = long_name
        return instr_id

    def _get_data_level(self):
        """
        Function to get Data Level of CDF data

        The level to which the data product has been processed.
        """
        attr_name = "Data_level"
        if (attr_name not in self.data.meta) or (not self.data.meta[attr_name]):
            data_level = None
        else:
            data_level = self.data.meta["Data_level"]
            # Formatting
            if ">" in data_level:
                short_name, _ = data_level.split(">")
                data_level = short_name.lower()  # Makse sure its all lowercase
        return data_level

    def _get_data_level_long_name(self):
        """
        Function to get Data Level of CDF data

        The level to which the data product has been processed.
        """
        attr_name = "Data_level"
        if (attr_name not in self.data.meta) or (not self.data.meta[attr_name]):
            data_level = None
        else:
            data_level = self.data.meta["Data_level"]
            # Formatting
            if ">" in data_level:
                _, long_name = data_level.split(">")
                data_level = long_name
        return data_level

    def _get_data_product_descriptor(self):
        """
        Function to get the (Optional) Data Product Descriptor.

        This is an optional field that may not be needed for all products. Where it is used,
        identifier shouls be short (3-8 charachters) descriptors that are helpful to end users.
        If a descriptor contains multiple components, underscores are used top separate
        hose components.
        """
        attr_name = "Data_product_descriptor"
        if (attr_name not in self.data.meta) or (not self.data.meta[attr_name]):
            odpd = ""
        else:
            odpd = self.data.meta["Data_product_descriptor"]
        return odpd

    def _get_generation_date(self):
        """
        Function to get the date that the CDF was generated.
        """
        return datetime.datetime.now()

    def _get_start_time(self):
        """
        Function to get the start time of the data contained in the CDF
        given in format `YYYYMMDD_hhmmss`
        """
        gattr_name = "Start_time"
        vattr_name = "time"
        if (gattr_name in self.data.meta) and (self.data.meta[gattr_name]):
            start_time = self.data.meta[gattr_name]
        elif vattr_name not in self.data.columns:
            start_time = None
        else:
            # Get the Start Time from the TimeSeries
            start_time = self.data["time"].to_datetime()[0]
        return start_time

    def _get_version(self):
        """
        Function to get the 3-part version number of the data product.
        """
        attr_name = "Data_version"
        if (attr_name not in self.data.meta) or (not self.data.meta[attr_name]):
            version = None
        else:
            version_str = self.data.meta["Data_version"].lower()
            if "v" in version_str:
                _, version = version_str.split("v")
            else:
                version = version_str
        return version

    def _get_instrument_mode(self):
        """Function to get the mode attribute (TBS)"""
        attr_name = "Instrument_mode"
        if (attr_name not in self.data.meta) or (not self.data.meta[attr_name]):
            instr_mode = ""
        else:
            instr_mode = self.data.meta["Instrument_mode"]
        return instr_mode.lower()  # Makse sure its all lowercase
