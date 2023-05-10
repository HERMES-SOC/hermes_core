"""
Data Container class for Heliophysics Science Data.
"""

from pathlib import Path
from typing import OrderedDict
from astropy.timeseries import TimeSeries
from astropy.table import vstack
from astropy import units as u
from hermes_core.util.validation import CDFValidator, NetCDFValidator, FITSValidator
from hermes_core.util.io import TimeDataIOHandler, CDFHandler, NetCDFHandler, FITSHandler


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
    A class for storing and manipulating science data.

    Parameters:
        data (TimeSeries): The science data stored as a TimeSeries object.
        handler (TimeDataIOHandler, optional): The handler for input/output operations. Defaults to None.

    Attributes:
        meta (dict): Metadata associated with the science data.
        handler (TimeDataIOHandler): The handler for input/output operations.

    """

    def __init__(self, data, meta=None, handler=None, **kwargs):
        """
        Initialize a TimeData object.

        Parameters:
            data (TimeSeries): The science data stored as a TimeSeries object.
            meta: `dict`
            handler (TimeDataIOHandler, optional): The handler for input/output operations. Defaults to None.

        Raises:
            ValueError: If the number of columns is less than 2, the required 'time' column is missing,
                        or any column, excluding 'time', is not an astropy.Quantity object with units.

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
        self.data = TimeSeries(data, copy=True)

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

    @property
    def meta(self):
        """
        Metadata associated with the science data.

        Returns:
            dict: The metadata associated with the science data.

        """
        return self.data.meta

    @meta.setter
    def meta(self, value):
        """
        Set the metadata associated with the science data.

        Parameters:
            value (dict): The metadata to set.

        """
        self.data.meta = value

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
        (OrderedDict), returns the units of the columns of data
        """
        units = {}
        for name in self.data.columns:
            var_data = self.data[name]
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
        (OrderedDict), returns columns from data.columns
        """
        return self.data.columns

    @property
    def time(self):
        """
        returns the time array from data.time
        """
        if "time" in self.data.columns:
            return self.data.time
        else:
            return None

    @property
    def shape(self):
        """
        The shape of the data, a tuple (nrows, ncols)
        """
        if "time" in self.data.columns:
            nrows = self.data.time.shape[0]
        else:
            nrows = 0
        ncols = len(self.data.columns)
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
        for attr_name, attr_value in self.data.meta.items():
            str_repr += f"\t{attr_name}: {attr_value}\n"
        # Variable Data
        str_repr += f"Variable Data:\n{self.data}\n"
        # Variable Attributes
        str_repr += f"Variable Attributes:\n"
        for col_name in self.data.columns:
            str_repr += f"\tVar: {col_name}\n"
            # for attr_name, attr_value in self.data[col_name].meta.items():
            #     str_repr += f"\t\t{attr_name}: {attr_value}\n"
        return str_repr

    def __len__(self):
        """
        Function to get the number of variable data members in the CDFWriter class.
        """
        return len(self.data.keys())

    def __getitem__(self, name):
        """
        Function to get a data variable contained in the CDFWriter class.
        """
        if name not in self.data.colnames:
            raise KeyError(f"CDFWriter does not contain data variable {name}")
        # Get the Data and Attrs for the named variable
        var_data = self.data[name]
        return var_data

    def __setitem__(self, name, data):
        """
        Function to set a data variable conained in the CDFWriter class.
        """
        # Set the Data for the named variable
        self.add_measurement(measure_name=name, measure_data=data, measure_meta={})

    def __contains__(self, name):
        """
        Function to see whether a data variable is in the CDFWriter class.
        """
        return name in self.data.columns

    def __iter__(self):
        """
        Function to iterate over data variables and attributes in the CDFWriter class.
        """
        for name in self.data.columns:
            var_data = self.data[name]
            yield (name, var_data)

    def add_measurement(self, measure_name: str, measure_data: u.Quantity, measure_meta: dict):
        """
        Function to add variable data to the TimeData class. Variable data here is assumed
        to be array-like or matrix-like to be stored in the TimeData. Additionally, varaible
        attributes can be added though a native python `dict` of `(key, value)` pairs that
        are added to the TimeData variable.

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
        self.data.remove_column(measure_name)

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
