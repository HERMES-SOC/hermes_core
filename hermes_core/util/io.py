from abc import ABC, abstractmethod
from pathlib import Path
from collections import OrderedDict
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
                    # Extract the Variable's Data and Metadata
                    var_data = input_file[var_name][:].copy()
                    var_attrs = {}
                    for attr_name in input_file[var_name].attrs:
                        var_attrs[attr_name] = input_file[var_name].attrs[attr_name]

                    # Check if the Variable is a Record-Variing Variable
                    if input_file[var_name].rv():
                        # Test the the Variable has `UNITS`
                        # If the variable has `UNITS` then it is likely a Measurement
                        # and should be stored with Measurements
                        if "UNITS" in var_attrs:
                            # Add to the main TimeSeries Table
                            ts[var_name] = var_data
                            ts[var_name].unit = var_attrs["UNITS"]
                            # Create the Metadata
                            ts[var_name].meta = OrderedDict()
                            ts[var_name].meta.update(var_attrs)
                        # If the variable does not have `UNITS` then it is most likely Support
                        # data and it should be stored with the Support Data
                        else:
                            # Create an empty dict entry for the variable
                            support_data[var_name] = {}
                            # Add the Variable Data
                            support_data[var_name]["data"] = var_data
                            # Add the Variable Metadata
                            support_data[var_name]["meta"] = var_attrs
                    # Add to the non-record varying data dictionary
                    else:
                        # Create an empty dict entry for the variable
                        nrv_data[var_name] = {}
                        # Add the Variable Data
                        nrv_data[var_name]["data"] = var_data
                        # Add the Variable Metadata
                        nrv_data[var_name]["meta"] = var_attrs

        # Return the given TimeSeries, NRV Data
        return ts, support_data, nrv_data

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
                raise ValueError(
                    f"Cannot Add gAttr: {attr_name}. Value was {str(attr_value)} "
                )
            else:
                # Add the Attribute to the CDF File
                cdf_file.attrs[attr_name] = attr_value

    def _convert_variable_attributes_to_cdf(self, data, cdf_file):
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
                    if var_attr_val is None:
                        raise ValueError(
                            f"Variable {var_name}: Cannot Add vAttr: {var_attr_name}. Value was {str(var_attr_val)}"
                        )
                    else:
                        # Add the Attribute to the CDF File
                        cdf_file[var_name].attrs[var_attr_name] = var_attr_val


# ================================================================================================
#                                   EPHEMERIS HANDLER
# ================================================================================================


class EphemDataHandler(TimeDataIOHandler):
    """
    A concrete implementation of TimeDataIOHandler for handling ephemeris data in `.txt` file format.
    """

    def load_data(self, file_path):
        """
        Load ephemeris data from a EphemData file.

        Parameters:
            file_path (str): The path to the EphemData file.

        Returns:
            data (astropy.TimeSeries): An instance of astropy.TimeSeries containing the loaded data.
        """
        if not Path(file_path).exists():
            raise FileNotFoundError(f"CDF Could not be loaded from path: {file_path}")

        file_name = Path(file_path).stem
        mission_name, _, _, satellite, _, phase = file_name.split("_")

        # Helpers
        generation_date = None
        title = None
        header_parsed = False
        separator_parsed = False

        # Intermittent Format
        columns = []
        data = {}
        units = {}

        # Open EphemData with a context manager
        with open(file_path) as input_file:
            for i, line in enumerate(input_file):
                if len(line) == 0:
                    continue
                # The Generation Date is the First line we're expecting from the EphemData
                # If we haven't parsed out the Generation Date from the File
                if len(line.strip()) > 0 and not generation_date:
                    # Parse the Generation Data from the EphemData file.
                    generation_date = datetime.strptime(
                        line.strip(), "%d %b %Y %H:%M:%S"
                    )
                # The EphemData Name is the Seconf line We're expecting from the EphemData
                elif len(line.strip()) > 0 and not title:
                    # Parse the Title from the EphemData file.
                    title = line.strip()
                elif len(line.strip()) > 0 and not header_parsed:
                    # Parse the Header
                    parts = line.strip().split()
                    i = 0
                    while i < len(parts):
                        column_name = parts[i]
                        if i + 1 < len(parts) and "(" in parts[i + 1]:
                            column_unit = parts[i + 1].strip("()")
                            if "/" in column_unit:
                                column_unit = " / ".join(column_unit.split("/"))
                            column_unit = column_unit.replace("sec", "s")
                            i += 1
                        else:
                            column_unit = 1

                        # Update Intermittent Structures
                        columns.append(column_name)
                        data[column_name] = []
                        units[column_name] = column_unit
                        i += 1
                    header_parsed = True
                elif len(line.strip()) > 0 and not separator_parsed:
                    line = line.strip().replace(" ", "")
                    if line == len(line) * line[0]:
                        separator_parsed = True
                # Now we're in the actual measurements part of the file
                elif len(line.strip()) > 0:
                    parts = line.strip().split()

                    # Time
                    time_str = " ".join(parts[:4])
                    time = datetime.strptime(time_str, "%d %b %Y %H:%M:%S.%f")
                    data[columns[0]].append(time)

                    for i, part in enumerate(parts[4:]):
                        # Float
                        value = float(part.strip())
                        data[columns[i + 1]].append(value)

        # Create a new TimeSeries
        ts = TimeSeries()

        # Add Time to the TimeSeries
        ts["time"] = data["Time"]

        # Add Columns to the TimeSeries
        for col in columns[1:]:
            # Create the Quantity object
            ts[col] = data[col]
            ts[col].unit = units[col]
            ts[col].meta = OrderedDict({"VAR_TYPE": "data"})

        # Update the Metadata of the Data Container
        ts.meta.update(
            {
                "Data_level": "L1>Level 1",
                "Data_version": "0.0.1",
                "Descriptor": "EEA>Electron Electrostatic Analyzer",
                "Generation_date": generation_date,
                "Mission_group": mission_name,
                "Phase": phase,
                "Satellite": satellite,
            }
        )

        # Return the given TimeSeries
        return ts

    def save_data(self, data, file_path):
        """
        Save ephemeris data to a EphemData file.

        Parameters:
            data (TimeData): An instance of TimeData containing the data to be saved.
            file_path (str): The path to save the EphemData file.
        """
        pass


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
