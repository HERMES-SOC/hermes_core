from abc import ABC, abstractmethod
from pathlib import Path
from collections import OrderedDict
from datetime import datetime
import numpy as np
from astropy.timeseries import TimeSeries
from astropy.time import Time
from hermes_core.util.exceptions import warn_user
from hermes_core.util.schema import CDFSchema

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
        self.schema = CDFSchema()

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
        """
        from spacepy.pycdf import CDF

        if not Path(file_path).exists():
            raise FileNotFoundError(f"CDF Could not be loaded from path: {file_path}")

        # Create a new TimeSeries
        ts = TimeSeries()

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

            # Get all the Keys for Measurement Variable Data
            # These are Keys where the underlying object is a `dict` that contains
            # additional data, and is not the `EPOCH` variable
            variable_keys = filter(lambda key: key != "Epoch", list(input_file.keys()))
            # Add Variable Attributtes from the CDF file to TimeSeries
            for var_name in variable_keys:
                # Make sure the Variable is record-varying
                if len(input_file[var_name]) != len(ts["time"]):
                    warn_user(
                        f"Measurement Variable {var_name} does not match the length of the EPOCH variable. Cannot add {var_name} to the TimeSeries"
                    )
                    # Skip to the next Variable
                    continue
                # Make sure the Variable is one-dimensional
                if isinstance(input_file[var_name][0], np.ndarray):
                    warn_user(
                        f"Measurement Variable {var_name} is Multi-Dimensional. Cannot add {var_name} to the TimeSeries"
                    )
                    # Skip to the next Variable
                    continue

                # If all checks pass, then add to the TimeSeries
                var_data = input_file[var_name][:].copy()
                var_attrs = {}
                for attr_name in input_file[var_name].attrs:
                    var_attrs[attr_name] = input_file[var_name].attrs[attr_name]
                # Create the Quantity object
                ts[var_name] = var_data
                ts[var_name].unit = var_attrs["UNITS"]
                # Create the Metadata
                ts[var_name].meta = OrderedDict()
                ts[var_name].meta.update(var_attrs)

        # Return the given TimeSeries
        return ts

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
                raise ValueError(f"Cannot Add gAttr: {attr_name}. Value was {str(attr_value)} ")
            else:
                # Add the Attribute to the CDF File
                cdf_file.attrs[attr_name] = attr_value

    def _convert_variable_attributes_to_cdf(self, data, cdf_file):
        # Loop through Variable Attributes in target_dict
        for var_name, var_data in data.__iter__():
            if var_name == "time":
                var_name = "Epoch"
                # Add 'time' in the TimeSeries as 'Epoch' within the CDF
                cdf_file["Epoch"] = var_data.to_datetime()
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


# ================================================================================================
#                                   NET CDF HANDLER
# ================================================================================================


class NetCDFHandler(TimeDataIOHandler):
    """
    A concrete implementation of TimeDataIOHandler for handling heliophysics data in NetCDF format.

    This class provides methods to load and save heliophysics data from/to a NetCDF file.
    """

    def load_data(self, file_path):
        """
        Load heliophysics data from a NetCDF file.

        Parameters
        ----------
        file_path : `str`
            The path to the NetCDF file.

        Returns
        -------
        data : `~astropy.time.TimeSeries`
            An instance of `TimeSeries` containing the loaded data.
        """
        pass

    def save_data(self, data, file_path):
        """
        Save heliophysics data to a NetCDF file.

        Parameters
        ----------
        data : `hermes_core.timedata.TimeData`
            An instance of `TimeData` containing the data to be saved.
        file_path : `str`
            The path to save the NetCDF file.

        Returns
        -------
        path : `str`
            A path to the saved file.
        """
        pass


# ================================================================================================
#                                   JSON DATA HANDLER
# ================================================================================================


class JSONDataHandler(TimeDataIOHandler):
    """
    A concrete implementation of TimeDataIOHandler for handling heliophysics data in JSON format.

    This class provides methods to load and save heliophysics data from/to a JSON file.
    """

    def load_data(self, file_path):
        """
        Load heliophysics data from a JSON file.

        Parameters
        ----------
        file_path : `str`
            The path to the JSON file.

        Returns
        -------
        data : `~astropy.time.TimeSeries`
            An instance of `TimeSeries` containing the loaded data.
        """
        import json

        if not Path(file_path).exists():
            raise FileNotFoundError(f"JSON Data Could not be loaded from path: {file_path}")

        # Create a new TimeSeries
        ts = TimeSeries()

        # Open CDF file with context manager
        with open(file_path) as input_file:
            data = json.load(input_file)

            # Data Downloaded from the SPDF sa JSON may be downloaded as a `list` object or as a `dict` object
            # If a `list` object then the data is in the first elemen of the `list` that should be a `dict`
            if isinstance(data, list):
                data = data[0]

            if "EPOCH_" in data:
                # `DAT` Appears to be the hardcoded Data Array
                time_data = Time(data["EPOCH_"]["DAT"])

                # Create the Time object
                ts["time"] = time_data
                # Create the Metadata
                ts["time"].meta = OrderedDict()

                # Get all the Metadata Attributs
                # This should be all the keys except for the "DAT" key which holds the measurement data
                metadata_attrs = filter(lambda key: key != "DAT", list(data["EPOCH_"].keys()))
                for attr_name in metadata_attrs:
                    ts["time"].meta[attr_name] = data["EPOCH_"][attr_name]

            # Get all the Keys for Measurement Variable Data
            # These are Keys where the underlying object is a `dict` that contains
            # additional data, and is not the `EPOCH` variable
            variable_keys = filter(
                lambda key: key != "EPOCH_" and isinstance(data[key], dict), list(data.keys())
            )
            # Loop Through the Variable Keys
            for key in variable_keys:
                # Get the Measurement Variable
                # This contains both the Metadata and the Data
                measurement_variable = data[key]
                # `DAT` Appears to be the hardcoded Data Array
                measurement_data = measurement_variable["DAT"]

                if len(measurement_data) != len(ts["time"]):
                    warn_user(
                        f"Measurement Variable {key} does not match the length of the EPOCH variable. Cannot Add {key} to the TimeSeries"
                    )
                    # Skip to the next Variable
                    continue
                if isinstance(measurement_data[0], list):
                    warn_user(
                        f"Measurement Variable {key} is Multi-Dimensional. Cannot Add {key} to the TimeSeries"
                    )
                    # Skip to the next Variable
                    continue

                # Add the Measurement to the TimeSeries as a Quantity
                ts[key] = measurement_data
                ts[key].unit = measurement_variable["UNITS"]
                # Create the Metadata
                ts[key].meta = OrderedDict()

                # Get all the Metadata Attributs
                # This should be all the keys except for the "DAT" key which holds the measurement data
                metadata_attrs = filter(lambda key: key != "DAT", list(measurement_variable.keys()))
                # Loop through the Metadata Attributes
                for attr_name in metadata_attrs:
                    # Add the Variable attribute to the Variable's Metadata in the TimeSeries
                    ts[key].meta[attr_name] = measurement_variable[attr_name]

            # Get all the Keys for Global Metadata
            # These are Keys where the underlying object is not a `dict`
            global_metadata_keys = filter(
                lambda key: not isinstance(data[key], dict), list(data.keys())
            )
            for key in global_metadata_keys:
                # Add to the TimeSeries meta
                ts.meta[key] = data[key]

        # Return the given TimeSeries
        return ts

    def save_data(self, data, file_path):
        """
        Save heliophysics data to a JSON file.

        Parameters
        ----------
        data : `hermes_core.timedata.TimeData`
            An instance of `TimeData` containing the data to be saved.
        file_path : `str`
            The path to save the JSON file.

        Returns
        -------
        path : `str`
            A path to the saved file.
        """
        import json

        # Initialize a new JSON File
        json_filename = f"{data.meta['Logical_file_id']}.json"
        output_json_filepath = str(Path(file_path) / json_filename)

        # Create a New JSON Object (dict) to Write to a JSON File
        output_json_data = {}

        # Add Global Attriubtes to the JSON File
        self._convert_global_attributes_to_json(data, output_json_data)

        # Add Measurement Data and Attributes to the JSON File
        self._convert_variable_attributes_to_json(data, output_json_data)

        with open(output_json_filepath, "w") as json_file:
            json.dump(output_json_data, json_file)

        return output_json_filepath

    def _convert_global_attributes_to_json(self, data, json_file: dict):
        # Loop though Global Attributes in target_dict
        for attr_name, attr_value in data.meta.items():
            # Make sure the Value is not None
            # We cannot add None Values to the CDF Global Attrs
            if attr_value is None:
                raise ValueError(f"Cannot Add gAttr: {attr_name}. Value was {str(attr_value)} ")
            else:
                # Add the Attribute to the JSON File
                if isinstance(attr_value, datetime):
                    json_file[attr_name] = attr_value.isoformat()
                else:
                    json_file[attr_name] = attr_value

    def _convert_variable_attributes_to_json(self, data, json_file: dict):
        # Loop through Variable Attributes in target_dict
        for var_name, var_data in data.__iter__():
            if var_name == "time":
                var_name = "EPOCH_"
                # Add 'time' in the TimeSeries as 'Epoch' within the JSON
                json_file["EPOCH_"] = {}
                # `DAT` Appears to be the hardcoded Data Array
                json_file["EPOCH_"]["DAT"] = var_data.to_value("isot").tolist()
            else:
                # Add the Variable to the JSON File
                json_file[var_name] = {}
                # `DAT` Appears to be the hardcoded Data Array
                json_file[var_name]["DAT"] = var_data.value.tolist()
            # Add the Variable Attributes
            for var_attr_name, var_attr_val in var_data.meta.items():
                if var_attr_val is None:
                    raise ValueError(
                        f"Variable {var_name}: Cannot Add vAttr: {var_attr_name}. Value was {str(var_attr_val)}"
                    )
                else:
                    # Add the Attribute to the JSON File
                    if isinstance(var_attr_val, datetime):
                        json_file[var_name][var_attr_name] = var_attr_val.isoformat()
                    else:
                        json_file[var_name][var_attr_name] = var_attr_val
