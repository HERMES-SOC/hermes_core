from abc import ABC, abstractmethod
from pathlib import Path
from typing import OrderedDict
from copy import deepcopy
import math
import datetime
from astropy.timeseries import TimeSeries
from astropy.time import Time
from astropy import units as u
import spacepy
from spacepy.pycdf import CDF, _Hyperslice
import hermes_core
from hermes_core import log
from hermes_core.util import util
from hermes_core.util.exceptions import warn_user
from hermes_core.util.schema import CDFSchema

# ================================================================================================
#                                   ABSTRACT HANDLER
# ================================================================================================


class ScienceDataIOHandler(ABC):
    """
    Abstract base class for handling input/output operations of heliophysics data.
    """

    @abstractmethod
    def load_data(self, file_path):
        """
        Load heliophysics data from a file.

        Parameters:
            file_path (str): The path to the data file.

        Returns:
            data (ScienceData): An instance of ScienceData containing the loaded data.
        """
        pass

    @abstractmethod
    def save_data(self, data, file_path):
        """
        Save heliophysics data to a file.

        Parameters:
            data (ScienceData): An instance of ScienceData containing the data to be saved.
            file_path (str): The path to save the data file.
        """
        pass


# ================================================================================================
#                                   CDF HANDLER
# ================================================================================================


class CDFHandler(ScienceDataIOHandler):
    """
    A concrete implementation of ScienceDataIOHandler for handling heliophysics data in CDF format.

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
            data (ScienceData): An instance of ScienceData containing the data to be saved.
            file_path (str): The path to save the CDF file.

        Returns:
            str: A path to the saved file.
        """

        # Update Default Attributes
        self._update_default_attributes(data)

        # Derive any Global Attributes
        self.derive_attributes(data)

        # Initialize a new CDF
        cdf_filename = f"{data.meta['Logical_file_id']}.cdf"
        output_cdf_filepath = str(Path(file_path) / cdf_filename)
        with CDF(output_cdf_filepath, masterpath="") as cdf_file:
            # Add Global Attriubtes to the CDF File
            self._convert_global_attributes_to_cdf(data, cdf_file)

            # Add zAttributes
            self._convert_variable_attributes_to_cdf(data, cdf_file)
        return output_cdf_filepath

    def _update_default_attributes(self, data):
        # Loop through the default Global Attributes
        for attr_name, attr_value in self.schema.default_global_attrs.items():
            if attr_name not in data.meta:
                data.meta[attr_name] = attr_value

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

    def derive_attributes(self, data):
        """Function to derive global attributes"""
        # Loop through Global Attributes
        for attr_name, attr_schema in self.schema.global_attribute_schema.items():
            if attr_schema["derived"]:
                derived_value = self._derive_global_attribute(data, attr_name=attr_name)
                # Only Derive Global Attributes if they have not been manually derived/overridden
                if (attr_name not in data.meta) or (not data.meta[attr_name]):
                    data.meta[attr_name] = derived_value
                else:
                    log.debug(
                        (
                            "Attribute: %s was marked for derivation (to be %s)"
                            "but was already overridden to %s"
                        ),
                        attr_name,
                        derived_value,
                        data.meta[attr_name],
                    )
        # Loop through Variable Attributes
        for name in data.columns:
            self._derive_variable_attributes(data, var_name=name)

    def _derive_variable_attributes(self, data, var_name):
        var_data = data[var_name]

        # Derive Time-Specific Attributes
        if var_name == "time":
            self._derive_time_attributes(data)

        # Check the Attributes that can be derived
        if "DEPEND_0" not in var_data.meta:
            var_data.meta["DEPEND_0"] = self._get_depend(data)
        if "FIELDNAM" not in var_data.meta:
            var_data.meta["FIELDNAM"] = self._get_fieldnam(data, var_name)
        if "FILLVAL" not in var_data.meta:
            var_data.meta["FILLVAL"] = self._get_fillval(data, var_name)
        if "FORMAT" not in var_data.meta:
            var_data.meta["FORMAT"] = self._get_format(data, var_name)
        if "SI_CONVERSION" not in var_data.meta:
            var_data.meta["SI_CONVERSION"] = self._get_si_conversion(data, var_name)
        if "UNITS" not in var_data.meta:
            var_data.meta["UNITS"] = self._get_units(data, var_name)
        if "VALIDMIN" not in var_data.meta:
            var_data.meta["VALIDMIN"] = self._get_validmin(data, var_name)
        if "VALIDMAX" not in var_data.meta:
            var_data.meta["VALIDMAX"] = self._get_validmax(data, var_name)

    def _derive_time_attributes(self, data):
        var_data = data.time

        # Check the Attributes that can be derived
        if "REFERENCE_POSITION" not in var_data.meta:
            var_data.meta["REFERENCE_POSITION"] = self._get_reference_position(data)
        if "RESOLUTION" not in var_data.meta:
            var_data.meta["RESOLUTION"] = self._get_resolution(data)
        if "TIME_BASE" not in var_data.meta:
            var_data.meta["TIME_BASE"] = self._get_time_base(data)
        if "TIME_SCALE" not in var_data.meta:
            var_data.meta["TIME_SCALE"] = self._get_time_scale(data)
        if "UNITS" not in var_data.meta:
            var_data.meta["UNITS"] = self._get_time_units(data)

    def _derive_global_attribute(self, data, attr_name):
        """
        Function to Derive Global Metadata Attributes
        """
        # SWITCH on the Derivation attr_name
        if attr_name == "Generation_date":
            return self._get_generation_date(data)
        elif attr_name == "Start_time":
            return self._get_start_time(data)
        elif attr_name == "Data_type":
            return self._get_data_type(data)
        elif attr_name == "Logical_file_id":
            return self._get_logical_file_id(data)
        elif attr_name == "Logical_source":
            return self._get_logical_source(data)
        elif attr_name == "Logical_source_description":
            return self._get_logical_source_description(data)
        else:
            raise ValueError(f"Derivation for Attribute ({attr_name}) Not Recognized")

    def _get_depend(self, data):
        return "Epoch"

    def _get_fieldnam(self, data, var_name):
        if var_name != "time":
            return deepcopy(var_name)
        else:
            return "Epoch"

    def _get_fillval(self, data, var_name):
        # Get the Variable Data
        var_data = data[var_name]
        if var_name == "time":
            # Guess the spacepy.pycdf.const CDF Data Type
            (guess_dims, guess_types, guess_elements) = _Hyperslice.types(var_data.to_datetime())
            # Get the FILLVAL for the gussed data type
            fillval = self._fillval_helper(data, cdf_type=guess_types[0])
            # guess_types[0] == spacepy.pycdf.const.CDF_TIME_TT2000.value:
            return spacepy.pycdf.lib.v_tt2000_to_datetime(fillval)
        else:
            # Guess the spacepy.pycdf.const CDF Data Type
            (guess_dims, guess_types, guess_elements) = _Hyperslice.types(var_data.value)
            # Get the FILLVAL for the gussed data type
            fillval = self._fillval_helper(data, cdf_type=guess_types[0])
            return fillval

    def _fillval_helper(self, data, cdf_type):
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

    def _get_format(self, data, var_name):
        # Get the Variable Data
        var_data = data[var_name]
        if var_name == "time":
            # Guess the spacepy.pycdf.const CDF Data Type
            (guess_dims, guess_types, guess_elements) = _Hyperslice.types(var_data.to_datetime())
            return self._format_helper(data, var_name, guess_types[0])
        else:
            # Guess the spacepy.pycdf.const CDF Data Type
            (guess_dims, guess_types, guess_elements) = _Hyperslice.types(var_data.value)
            return self._format_helper(data, var_name, guess_types[0])

    def _format_helper(self, data, var_name, cdftype):
        minn = "VALIDMIN"
        maxx = "VALIDMAX"
        # Get the Variable Data
        var_data = data[var_name]

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

    def _get_reference_position(self, data):
        # Get the Variable Data
        var_data = data.time
        # Guess the spacepy.pycdf.const CDF Data Type
        (guess_dims, guess_types, guess_elements) = _Hyperslice.types(var_data.to_datetime())
        if guess_types[0] == spacepy.pycdf.const.CDF_TIME_TT2000.value:
            return "rotating Earth geoid"
        else:
            raise TypeError(f"Reference Position for Time type ({guess_types[0]}) not found.")

    def _get_resolution(self, data):
        # Get the Variable Data
        var_data = data.time
        times = len(var_data)
        if times < 2:
            raise ValueError(f"Can not derive Time Resolution, need 2 samples, found {times}.")
        # Calculate the Timedelta between two datetimes
        times = var_data.to_datetime()
        delta = times[1] - times[0]
        # Get the number of second between samples.
        delta_seconds = delta.total_seconds()
        return f"{delta_seconds}s"

    def _get_si_conversion(self, data, var_name):
        # Get the Variable Data
        var_data = data[var_name]
        if var_name == "time":
            time_unit_str = self._get_time_units(data)
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

    def _get_time_base(self, data):
        # Get the Variable Data
        var_data = data.time
        # Guess the spacepy.pycdf.const CDF Data Type
        (guess_dims, guess_types, guess_elements) = _Hyperslice.types(var_data.to_datetime())
        if guess_types[0] == spacepy.pycdf.const.CDF_TIME_TT2000.value:
            return "J2000"
        else:
            raise TypeError(f"Time Base for Time type ({guess_types[0]}) not found.")

    def _get_time_scale(self, data):
        # Get the Variable Data
        var_data = data.time
        # Guess the spacepy.pycdf.const CDF Data Type
        (guess_dims, guess_types, guess_elements) = _Hyperslice.types(var_data.to_datetime())
        if guess_types[0] == spacepy.pycdf.const.CDF_TIME_TT2000.value:
            return "Terrestrial Time (TT)"
        else:
            raise TypeError(f"Time Scale for Time type ({guess_types[0]}) not found.")

    def _get_time_units(self, data):
        # Get the Variable Data
        var_data = data.time
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

    def _get_units(self, data, var_name):
        # Get the Variable Data
        var_data = data[var_name]
        unit = ""
        # Get the Unit from the TimeSeries Quantity if it exists
        if hasattr(var_data, "unit"):
            unit = var_data.unit.name
        return unit

    def _get_validmin(self, data, var_name):
        # Get the Variable Data
        var_data = data[var_name]
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

    def _get_validmax(self, data, var_name):
        # Get the Variable Data
        var_data = data[var_name]
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

    def _get_logical_file_id(self, data):
        """
        Function to get the `Logical_file_id` required global attribute.

        The attribute stores the name of the CDF File without the file
        extension (e.g. '.cdf'). This attribute is requires to avoid
        loss of the originial source in case of renaming.
        """
        attr_name = "Logical_file_id"
        if (attr_name not in data.meta) or (not data.meta[attr_name]):
            # Get Parts
            instrument_id = self._get_instrument_id(data)
            start_time = self._get_start_time(data)
            data_level = self._get_data_level(data)
            version = self._get_version(data)
            mode = self._get_instrument_mode(data)

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
            science_filename = data.meta[attr_name]
        return science_filename

    def _get_logical_source(self, data):
        """
        Function to get the `Logical_source` required global attribute.

        This attribute determines the file naming convention in the SKT Editor
        and is used by CDA Web.
        """
        attr_name = "Logical_source"
        if (attr_name not in data.meta) or (not data.meta[attr_name]):
            # Get Parts
            spacecraft_id = self._get_spacecraft_id(data)
            instrument_id = self._get_instrument_id(data)
            data_type = self._get_data_type(data)
            data_type_short_name, _ = data_type.split(">")

            # Build Derivation
            logical_source = f"{spacecraft_id}_{instrument_id}_{data_type_short_name}"
        else:
            logical_source = data.meta[attr_name]
        return logical_source

    def _get_logical_source_description(self, data):
        """
        Function to get the `Logical_source_description` required global attribute.

        This attribute writes out the full words associated with the encryped
        `Logical_source`  attribute.
        """
        attr_name = "Logical_source_description"
        if (attr_name not in data.meta) or (not data.meta[attr_name]):
            # Get Parts
            spacecraft_long_name = self._get_spacecraft_long_name(data)
            instrument_long_name = self._get_instrument_long_name(data)
            data_type = self._get_data_type(data)
            _, data_type_long_name = data_type.split(">")
            logical_source_description = (
                f"{spacecraft_long_name} {instrument_long_name} {data_type_long_name}"
            )
        else:
            logical_source_description = data.meta[attr_name]
        return logical_source_description

    def _get_data_type(self, data):
        """
        Function to get the `Data_type` required global attribute.

        This attribute is used by the CDF Writing software to create the filename.
        It is a combination of the following components:
            - mode
            - data_level
            - optional_data_product_descriptor
        """
        attr_name = "Data_type"
        if (attr_name not in data.meta) or (not data.meta[attr_name]):
            short_parts = []
            long_parts = []

            # Get `mode`
            mode_short_name = self._get_instrument_mode(data)
            mode_long_name = self._get_instrument_mode(data)
            if bool(mode_short_name and mode_long_name):
                short_parts.append(mode_short_name)
                long_parts.append(mode_long_name)

            # Get `data level`
            data_level_short_name = self._get_data_level(data)
            data_level_long_name = self._get_data_level_long_name(data)
            if bool(data_level_short_name and data_level_long_name):
                short_parts.append(data_level_short_name)
                long_parts.append(data_level_long_name)

            # Get `data product descriptor`
            odpd_short_name = self._get_data_product_descriptor(data)
            odpd_long_name = self._get_data_product_descriptor(data)
            if bool(odpd_short_name and odpd_long_name):
                short_parts.append(odpd_short_name)
                long_parts.append(odpd_long_name)

            # Build Derivation
            data_type = "_".join(short_parts) + ">" + " ".join(long_parts)
        else:
            data_type = data.meta[attr_name]
        return data_type

    def _get_spacecraft_id(self, data):
        """Function to get Spacecraft ID from Source_name Global Attribute"""
        attr_name = "Source_name"
        if (attr_name not in data.meta) or (not data.meta[attr_name]):
            # Get Module Default
            sc_id = hermes_core.MISSION_NAME
        else:
            sc_id = data.meta["Source_name"]
            # Formatting
            if ">" in sc_id:
                short_name, _ = sc_id.split(">")
                sc_id = short_name.lower()  # Makse sure its all lowercase
        return sc_id

    def _get_spacecraft_long_name(self, data):
        """Function to get Spacecraft ID from Source_name Global Attribute"""
        attr_name = "Source_name"
        if (attr_name not in data.meta) or (not data.meta[attr_name]):
            # Get Module Default
            sc_id = hermes_core.MISSION_NAME
        else:
            sc_id = data.meta["Source_name"]
            # Formatting
            if ">" in sc_id:
                _, long_name = sc_id.split(">")
                sc_id = long_name
        return sc_id

    def _get_instrument_id(self, data):
        """
        Function to get Instrument ID from Descriptor Global Attribute

        Instrument of investigation identifier shortened to three
        letter acronym.
        """
        attr_name = "Descriptor"
        if (attr_name not in data.meta) or (not data.meta[attr_name]):
            instr_id = None
        else:
            instr_id = data.meta["Descriptor"]
            # Formatting
            if ">" in instr_id:
                short_name, _ = instr_id.split(">")
                instr_id = short_name.lower()  # Makse sure its all lowercase
        return instr_id

    def _get_instrument_long_name(self, data):
        """
        Function to get Instrument ID from Descriptor Global Attribute

        Instrument of investigation identifier shortened to three
        letter acronym.
        """
        attr_name = "Descriptor"
        if (attr_name not in data.meta) or (not data.meta[attr_name]):
            instr_id = None
        else:
            instr_id = data.meta["Descriptor"]
            # Formatting
            if ">" in instr_id:
                _, long_name = instr_id.split(">")
                instr_id = long_name
        return instr_id

    def _get_data_level(self, data):
        """
        Function to get Data Level of CDF data

        The level to which the data product has been processed.
        """
        attr_name = "Data_level"
        if (attr_name not in data.meta) or (not data.meta[attr_name]):
            data_level = None
        else:
            data_level = data.meta["Data_level"]
            # Formatting
            if ">" in data_level:
                short_name, _ = data_level.split(">")
                data_level = short_name.lower()  # Makse sure its all lowercase
        return data_level

    def _get_data_level_long_name(self, data):
        """
        Function to get Data Level of CDF data

        The level to which the data product has been processed.
        """
        attr_name = "Data_level"
        if (attr_name not in data.meta) or (not data.meta[attr_name]):
            data_level = None
        else:
            data_level = data.meta["Data_level"]
            # Formatting
            if ">" in data_level:
                _, long_name = data_level.split(">")
                data_level = long_name
        return data_level

    def _get_data_product_descriptor(self, data):
        """
        Function to get the (Optional) Data Product Descriptor.

        This is an optional field that may not be needed for all products. Where it is used,
        identifier shouls be short (3-8 charachters) descriptors that are helpful to end users.
        If a descriptor contains multiple components, underscores are used top separate
        hose components.
        """
        attr_name = "Data_product_descriptor"
        if (attr_name not in data.meta) or (not data.meta[attr_name]):
            odpd = ""
        else:
            odpd = data.meta["Data_product_descriptor"]
        return odpd

    def _get_generation_date(self, data):
        """
        Function to get the date that the CDF was generated.
        """
        return datetime.datetime.now()

    def _get_start_time(self, data):
        """
        Function to get the start time of the data contained in the CDF
        given in format `YYYYMMDD_hhmmss`
        """
        gattr_name = "Start_time"
        vattr_name = "time"
        if (gattr_name in data.meta) and (data.meta[gattr_name]):
            start_time = data.meta[gattr_name]
        elif vattr_name not in data.columns:
            start_time = None
        else:
            # Get the Start Time from the TimeSeries
            start_time = data["time"].to_datetime()[0]
        return start_time

    def _get_version(self, data):
        """
        Function to get the 3-part version number of the data product.
        """
        attr_name = "Data_version"
        if (attr_name not in data.meta) or (not data.meta[attr_name]):
            version = None
        else:
            version_str = data.meta["Data_version"].lower()
            if "v" in version_str:
                _, version = version_str.split("v")
            else:
                version = version_str
        return version

    def _get_instrument_mode(self, data):
        """Function to get the mode attribute (TBS)"""
        attr_name = "Instrument_mode"
        if (attr_name not in data.meta) or (not data.meta[attr_name]):
            instr_mode = ""
        else:
            instr_mode = data.meta["Instrument_mode"]
        return instr_mode.lower()  # Makse sure its all lowercase


# ================================================================================================
#                                   NET CDF HANDLER
# ================================================================================================


class NetCDFHandler(ScienceDataIOHandler):
    """
    A concrete implementation of ScienceDataIOHandler for handling heliophysics data in NetCDF format.

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
            data (ScienceData): An instance of ScienceData containing the data to be saved.
            file_path (str): The path to save the NetCDF file.
        """
        pass


class FITSHandler(ScienceDataIOHandler):
    """
    A concrete implementation of ScienceDataIOHandler for handling heliophysics data in FITS format.

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
            data (ScienceData): An instance of ScienceData containing the data to be saved.
            file_path (str): The path to save the FITS file.
        """
        pass
