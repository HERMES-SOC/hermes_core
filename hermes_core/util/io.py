from abc import ABC, abstractmethod
from pathlib import Path
from typing import OrderedDict
from astropy.timeseries import TimeSeries
from astropy.time import Time
from astropy import units as u
from spacepy.pycdf import CDF, _Hyperslice
import hermes_core
from hermes_core import log
from hermes_core.util.exceptions import warn_user
from hermes_core.util.schema import CDFSchema

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
        data : TimeData
            An data container object containing the file data.
        """
        pass

    @abstractmethod
    def save_data(self, data, file_path):
        """
        Save data to a file.

        Parameters
        ----------
        data : TimeData
            An instance of TimeData containing the data to be saved.
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

        Parameters:
            file_path (str): The path to the CDF file.

        Returns:
            data (astropy.TimeSeries): An instance of astropy.TimeSeries containing the loaded data.
        """
        if not Path(file_path).exists():
            raise FileNotFoundError(f"CDF Could not be loaded from path: {file_path}")

        # Create a new TimeSeries
        ts = TimeSeries()

        # Open CDF file with context manager
        with CDF(file_path) as input_file:
            # Add Global Attributes from the CDF file to TimeSeries
            input_global_attrs = {}
            for attr_name in input_file.attrs:
                if len(input_file.attrs[attr_name]) > 1:
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

        Parameters:
            data (TimeData): An instance of TimeData containing the data to be saved.
            file_path (str): The path to save the CDF file.

        Returns:
            str: A path to the saved file.
        """

        # Initialize a new CDF
        cdf_filename = f"{data.meta['Logical_file_id']}.cdf"
        output_cdf_filepath = str(Path(file_path) / cdf_filename)
        with CDF(output_cdf_filepath, masterpath="") as cdf_file:
            # Add Global Attriubtes to the CDF File
            self._convert_global_attributes_to_cdf(data, cdf_file)

            # Add zAttributes
            self._convert_variable_attributes_to_cdf(data, cdf_file)
        return output_cdf_filepath

    def _convert_global_attributes_to_cdf(self, data, cdf_file: CDF):
        # Loop though Global Attributes in target_dict
        for attr_name, attr_value in data.meta.items():
            # Make sure the Value is not None
            # We cannot add None Values to the CDF Global Attrs
            if not attr_value:
                warn_user(f"Cannot Add gAttr: {attr_name}. Value was {str(attr_value)} ")
            else:
                # Add the Attribute to the CDF File
                cdf_file.attrs[attr_name] = attr_value

    def _convert_variable_attributes_to_cdf(self, data, cdf_file: CDF):
        # Loop through Variable Attributes in target_dict
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

        Parameters:
            file_path (str): The path to the NetCDF file.

        Returns:
            data (astropy.TimeSeries): An instance of astropy.TimeSeries containing the loaded data.
        """
        pass

    def save_data(self, data, file_path):
        """
        Save heliophysics data to a NetCDF file.

        Parameters:
            data (TimeData): An instance of TimeData containing the data to be saved.
            file_path (str): The path to save the NetCDF file.
        """
        pass


class FITSHandler(TimeDataIOHandler):
    """
    A concrete implementation of TimeDataIOHandler for handling heliophysics data in FITS format.

    This class provides methods to load and save heliophysics data from/to a FITS file.
    """

    def load_data(self, file_path):
        """
        Load heliophysics data from a FITS file.

        Parameters:
            file_path (str): The path to the FITS file.

        Returns:
            data (astropy.TimeSeries): An instance of astropy.astropy.TimeSeries containing the loaded data.
        """
        pass

    def save_data(self, data, file_path):
        """
        Save heliophysics data to a FITS file.

        Parameters:
            data (TimeData): An instance of TimeData containing the data to be saved.
            file_path (str): The path to save the FITS file.
        """
        pass
