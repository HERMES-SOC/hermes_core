from abc import ABC, abstractmethod
from pathlib import Path
from collections import OrderedDict
from astropy.timeseries import TimeSeries
from astropy.time import Time
from astropy.nddata import NDData
import astropy.units as u
from hermes_core.util.exceptions import warn_user
from hermes_core.util.schema import HermesDataSchema

__all__ = ["CDFHandler"]

# ================================================================================================
#                                   ABSTRACT HANDLER
# ================================================================================================


class HermesDataIOHandler(ABC):
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
        data : `hermes_core.timedata.HermesData`
            An instance of `HermesData` containing the data to be saved.
        file_path : `str`
            The fully specified file path to save into.
        """
        pass


# ================================================================================================
#                                   CDF HANDLER
# ================================================================================================


class CDFHandler(HermesDataIOHandler):
    """
    A concrete implementation of HermesDataIOHandler for handling heliophysics data in CDF format.

    This class provides methods to load and save heliophysics data from/to a CDF file.
    """

    def __init__(self):
        super().__init__()

        # CDF Schema
        self.schema = HermesDataSchema()

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
        support : `dict`
            Non-record-varying data contained in the file
        """
        from spacepy.pycdf import CDF

        if not Path(file_path).exists():
            raise FileNotFoundError(f"CDF Could not be loaded from path: {file_path}")

        # Create a new TimeSeries
        ts = TimeSeries()
        # Create a Data Structure for Non-record Varying Data
        support = {}

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
                # Make sure the Variable is not multi-dimensional
                if len(input_file[var_name].shape) > 1:
                    warn_user(
                        f"Measurement Variable {var_name} is multi-dimensional. Cannot add {var_name} to the TimeSeries"
                    )
                    # Skip to the next Variable
                    continue

                # Extract the Variable's Metadata
                var_attrs = {}
                for attr_name in input_file[var_name].attrs:
                    var_attrs[attr_name] = input_file[var_name].attrs[attr_name]

                # Extract the Variable's Data
                var_data = input_file[var_name][...]
                if input_file[var_name].rv():
                    # See if it is record-varying data with Units
                    if "UNITS" in var_attrs and len(var_data) == len(ts["time"]):
                        # Load as Record-Varying `data`
                        try:
                            self._load_data_variable(ts, var_name, var_data, var_attrs)
                        except ValueError:
                            warn_user(
                                f"Cannot create Quantity for Variable {var_name} with UNITS {var_attrs['UNITS']}. Creating Quantity with UNITS 'dimensionless_unscaled'."
                            )
                            # Swap Units
                            var_attrs["UNITS_DESC"] = var_attrs["UNITS"]
                            var_attrs["UNITS"] = u.dimensionless_unscaled.to_string()
                            self._load_data_variable(ts, var_name, var_data, var_attrs)
                    else:
                        # Load as `support`
                        self._load_support_variable(
                            support, var_name, var_data, var_attrs
                        )
                else:
                    # Load Non-Record-Varying Data as `support`
                    self._load_support_variable(support, var_name, var_data, var_attrs)

        # Return the given TimeSeries, NRV Data
        return ts, support

    def _load_data_variable(self, ts, var_name, var_data, var_attrs):
        # Create the Quantity object
        var_data = u.Quantity(value=var_data, unit=var_attrs["UNITS"], copy=False)
        ts[var_name] = var_data
        # Create the Metadata
        ts[var_name].meta = OrderedDict()
        ts[var_name].meta.update(var_attrs)

    def _load_support_variable(self, support, var_name, var_data, var_attrs):
        # Create a NDData entry for the variable
        support[var_name] = NDData(data=var_data, meta=var_attrs)

    def save_data(self, data, file_path):
        """
        Save heliophysics data to a CDF file.

        Parameters
        ----------
        data : `hermes_core.timedata.HermesData`
            An instance of `HermesData` containing the data to be saved.
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
        for var_name in data.timeseries.colnames:
            var_data = data.timeseries[var_name]
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

        # Loop through Non-Record-Varying Data
        for var_name, var_data in data.support.items():
            # Guess the data type to store
            # Documented in https://github.com/spacepy/spacepy/issues/707
            _, var_data_types, _ = self.schema._types(var_data.data)
            # Add the Variable to the CDF File
            cdf_file.new(
                name=var_name,
                data=var_data.data,
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
