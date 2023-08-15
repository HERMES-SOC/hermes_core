from abc import ABC, abstractmethod
from pathlib import Path
from collections import OrderedDict
from astropy.timeseries import TimeSeries
from astropy.time import Time
from astropy.table import Column
import astropy.units as u
from hermes_core.util.exceptions import warn_user
from hermes_core.util.schema import HERMESDataSchema

__all__ = ["CDFHandler"]

# ================================================================================================
#                                   ABSTRACT HANDLER
# ================================================================================================


class TimeDataIOHandler(ABC):
    """
    Abstract base class for handling input/output operations of heliophysics data.
    """

    @abstractmethod
    def load_data(self, file_path):
        """
        Load data from a file.

        Parameters
        ----------
        file_path : `str`
            A fully specified file path.

        Returns
        -------
        data : `~astropy.time.TimeSeries`
            An instance of `TimeSeries` containing the loaded data.
        """
        pass

    @abstractmethod
    def save_data(self, data, file_path):
        """
        Save data to a file.

        Parameters
        ----------
        data : `hermes_core.timedata.TimeData`
            An instance of `TimeData` containing the data to be saved.
        file_path : `str`
            The fully specified file path to save into.
        """
        pass


# ================================================================================================
#                                   CDF HANDLER
# ================================================================================================


class CDFHandler(TimeDataIOHandler):
    """
    A concrete implementation of TimeDataIOHandler for handling heliophysics data in CDF format.

    This class provides methods to load and save heliophysics data from/to a CDF file.
    """

    def __init__(self):
        super().__init__()

        # CDF Schema
        self.schema = HERMESDataSchema()

    def load_data(self, file_path):
        """
        Load heliophysics data from a CDF file.

        Parameters
        ----------
        file_path : `str`
            The path to the CDF file.

        Returns
        -------
        data : `~astropy.time.TimeSeries`
            An instance of `TimeSeries` containing the loaded data.
        nrv_data : `dict`
            Non-record varying data contained in the file
        """
        from spacepy.pycdf import CDF

        if not Path(file_path).exists():
            raise FileNotFoundError(f"CDF Could not be loaded from path: {file_path}")

        # Create a new TimeSeries
        ts = TimeSeries()
        # Create a Data Structure for Support Data without Units
        support_data = {}
        # Create a Data Structure for Non-record Varying Data
        nrv_data = {}

        # Open CDF file with context manager
        with CDF(file_path) as input_file:
            # Add Global Attributes from the CDF file to TimeSeries
            input_global_attrs = {}
            for attr_name in input_file.attrs:
                if len(input_file.attrs[attr_name]) == 0:
                    # gAttr is not set
                    input_global_attrs[attr_name] = ""
                elif len(input_file.attrs[attr_name]) > 1:
                    # gAttr is a List
                    input_global_attrs[attr_name] = input_file.attrs[attr_name][:]
                else:
                    # gAttr is a single value
                    input_global_attrs[attr_name] = input_file.attrs[attr_name][0]
            ts.meta.update(input_global_attrs)

            # First Variable we need to add is time/Epoch
            if "Epoch" in input_file:
                time_data = Time(input_file["Epoch"][:].copy())
                time_attrs = {}
                for attr_name in input_file["Epoch"].attrs:
                    time_attrs[attr_name] = input_file["Epoch"].attrs[attr_name]
                # Create the Time object
                ts["time"] = time_data
                # Create the Metadata
                ts["time"].meta = OrderedDict()
                ts["time"].meta.update(time_attrs)

            # Add Variable Attributtes from the CDF file to TimeSeries
            for var_name in input_file:
                if var_name != "Epoch":  # Since we added this separately
                    # Extract the Variable's Metadata
                    var_attrs = {}
                    for attr_name in input_file[var_name].attrs:
                        var_attrs[attr_name] = input_file[var_name].attrs[attr_name]

                    if input_file[var_name].rv():
                        # Extract the Variable's Data
                        var_data = input_file[var_name][:].copy()
                        # See if it is `data` or `support_data`
                        if "UNITS" in var_attrs and len(var_data == len(ts["time"])):
                            # Load as Record-Varying `data`
                            try:
                                self._load_data_variable(
                                    ts, var_name, var_data, var_attrs
                                )
                            except ValueError:
                                warn_user(
                                    f"Cannot create Quantity for Variable {var_name} with UNITS {var_attrs['UNITS']}. Loading as Support Data."
                                )
                                self._load_support_data_variable(
                                    support_data, var_name, var_data, var_attrs
                                )
                        else:
                            # Load as `support_data`
                            self._load_support_data_variable(
                                support_data, var_name, var_data, var_attrs
                            )
                    else:
                        # Load NRV Data as `metadata`
                        self._load_metadata_variable(
                            nrv_data, input_file, var_name, var_attrs
                        )

        # Return the given TimeSeries, NRV Data
        return ts, support_data, nrv_data

    def _load_data_variable(self, ts, var_name, var_data, var_attrs):
        # Create the Quantity object
        var_data = u.Quantity(value=var_data, unit=var_attrs["UNITS"], copy=False)
        ts[var_name] = var_data
        # Create the Metadata
        ts[var_name].meta = OrderedDict()
        ts[var_name].meta.update(var_attrs)

    def _load_support_data_variable(self, support_data, var_name, var_data, var_attrs):
        # Create a Column entry for the variable
        support_data[var_name] = Column(data=var_data, meta=var_attrs)

    def _load_metadata_variable(self, nrv_data, input_file, var_name, var_attrs):
        # If the Var Data is Scalar
        if len(input_file[var_name].shape) == 0:
            # No Data
            var_data = input_file[var_name][...]
        # Otherwise its an array
        else:
            # Extract the Variable's Data
            var_data = input_file[var_name][:].copy()
        # Create a Column entry for the variable
        nrv_data[var_name] = Column(data=var_data, meta=var_attrs)

    def save_data(self, data, file_path):
        """
        Save heliophysics data to a CDF file.

        Parameters
        ----------
        data : `hermes_core.timedata.TimeData`
            An instance of `TimeData` containing the data to be saved.
        file_path : `str`
            The path to save the CDF file.

        Returns
        -------
        path : `str`
            A path to the saved file.
        """
        from spacepy.pycdf import CDF

        # Initialize a new CDF
        cdf_filename = f"{data.meta['Logical_file_id']}.cdf"
        output_cdf_filepath = str(Path(file_path) / cdf_filename)
        with CDF(output_cdf_filepath, masterpath="") as cdf_file:
            # Add Global Attriubtes to the CDF File
            self._convert_global_attributes_to_cdf(data, cdf_file)

            # Add zAttributes
            self._convert_variable_attributes_to_cdf(data, cdf_file)
        return output_cdf_filepath

    def _convert_global_attributes_to_cdf(self, data, cdf_file):
        # Loop though Global Attributes in target_dict
        for attr_name, attr_value in data.meta.items():
            # Make sure the Value is not None
            # We cannot add None Values to the CDF Global Attrs
            if attr_value is None:
                cdf_file.attrs[attr_name] = ""
            else:
                # Add the Attribute to the CDF File
                cdf_file.attrs[attr_name] = attr_value

    def _convert_variable_attributes_to_cdf(self, data, cdf_file):
        # Loop through Variable Attributes
        for var_name, var_data in data.__iter__():
            if var_name == "time":
                # Add 'time' in the TimeSeries as 'Epoch' within the CDF
                cdf_file["Epoch"] = var_data.to_datetime()
                # Add the Variable Attributes
                for var_attr_name, var_attr_val in var_data.meta.items():
                    cdf_file["Epoch"].attrs[var_attr_name] = var_attr_val
            else:
                # Add the Variable to the CDF File
                cdf_file[var_name] = var_data.value
                # Add the Variable Attributes
                for var_attr_name, var_attr_val in var_data.meta.items():
                    if var_attr_val is None:
                        raise ValueError(
                            f"Variable {var_name}: Cannot Add vAttr: {var_attr_name}. Value was {str(var_attr_val)}"
                        )
                    else:
                        # Add the Attribute to the CDF File
                        cdf_file[var_name].attrs[var_attr_name] = var_attr_val

        # Loop through Support Data
        for var_name, var_data in data.support_data.items():
            # Add the Variable to the CDF File
            cdf_file[var_name] = var_data.value
            # Add the Variable Attributes
            for var_attr_name, var_attr_val in var_data.meta.items():
                if var_attr_val is None:
                    raise ValueError(
                        f"Variable {var_name}: Cannot Add vAttr: {var_attr_name}. Value was {str(var_attr_val)}"
                    )
                else:
                    # Add the Attribute to the CDF File
                    cdf_file[var_name].attrs[var_attr_name] = var_attr_val

        # Loop through NRV Data
        for var_name, var_data in data.nrv_data.items():
            # Guess the data type to store
            # Documented in https://github.com/spacepy/spacepy/issues/707
            _, var_data_types, _ = self.schema._types(var_data)
            # Add the Variable to the CDF File
            cdf_file.new(
                name=var_name,
                data=var_data.value,
                type=var_data_types[0],
                recVary=False,
            )

            # Add the Variable Attributes
            for var_attr_name, var_attr_val in var_data.meta.items():
                if var_attr_val is None:
                    raise ValueError(
                        f"Variable {var_name}: Cannot Add vAttr: {var_attr_name}. Value was {str(var_attr_val)}"
                    )
                else:
                    # Add the Attribute to the CDF File
                    cdf_file[var_name].attrs[var_attr_name] = var_attr_val
