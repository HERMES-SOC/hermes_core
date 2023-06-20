from abc import ABC, abstractmethod
from pathlib import Path
from collections import OrderedDict
from copy import deepcopy
import math
import datetime
import yaml
from astropy import units as u
import spacepy
from spacepy.pycdf import _Hyperslice
import hermes_core
from hermes_core import log
from hermes_core.util import util

__all__ = ["CDFSchema"]

DEFAULT_GLOBAL_CDF_ATTRS_SCHEMA_FILE = "hermes_default_global_cdf_attrs_schema.yaml"
DEFAULT_GLOBAL_CDF_ATTRS_FILE = "hermes_default_global_cdf_attrs.yaml"
DEFAULT_VARIABLE_CDF_ATTRS_SCHEMA_FILE = "hermes_default_variable_cdf_attrs_schema.yaml"


class FileTypeSchema(ABC):
    """Abstract class representing the schema of a file type."""

    @property
    @abstractmethod
    def global_attribute_schema(self):
        """Schema for global attributes of the file."""
        pass

    @property
    @abstractmethod
    def variable_attribute_schema(self):
        """Schema for variable attributes of the file."""
        pass


class CDFSchema(FileTypeSchema):
    """Schema for CDF files."""

    def __init__(self):
        super().__init__()

        # Data Validation, Complaiance, Derived Attributes
        self._global_attr_schema = CDFSchema._load_default_global_attr_schema()

        # Data Validation and Compliance for Variable Data
        self._variable_attr_schema = CDFSchema._load_default_variable_attr_schema()

        # Load Default Global Attributes
        self.default_global_attributes = CDFSchema._load_default_attributes()

    @property
    def global_attribute_schema(self):
        """Schema for variable attributes of the file."""
        return self._global_attr_schema

    @property
    def variable_attribute_schema(self):
        """Schema for variable attributes of the file."""
        return self._variable_attr_schema

    @staticmethod
    def _load_default_global_attr_schema() -> dict:
        # The Default Schema file is contained in the `hermes_core/data` directory
        default_schema_path = str(
            Path(hermes_core.__file__).parent / "data" / DEFAULT_GLOBAL_CDF_ATTRS_SCHEMA_FILE
        )
        # Load the Schema
        return CDFSchema._load_yaml_data(yaml_file_path=default_schema_path)

    @staticmethod
    def _load_default_variable_attr_schema() -> dict:
        # The Default Schema file is contained in the `hermes_core/data` directory
        default_schema_path = str(
            Path(hermes_core.__file__).parent / "data" / DEFAULT_VARIABLE_CDF_ATTRS_SCHEMA_FILE
        )
        # Load the Schema
        return CDFSchema._load_yaml_data(yaml_file_path=default_schema_path)

    @staticmethod
    def _load_default_attributes():
        # The Default Attributes file is contained in the `hermes_core/data` directory
        default_attributes_path = str(
            Path(hermes_core.__file__).parent / "data" / DEFAULT_GLOBAL_CDF_ATTRS_FILE
        )
        return CDFSchema._load_yaml_data(yaml_file_path=default_attributes_path)

    @staticmethod
    def _load_yaml_data(yaml_file_path):
        """
        Function to load data from a Yaml file.

        Parameters
        ----------
        yaml_file_path: `str`
            Path to schem file to be used for CDF formatting.

        """
        assert isinstance(yaml_file_path, str)
        assert Path(yaml_file_path).exists()
        # Load the Yaml file to Dict
        yaml_data = {}
        with open(yaml_file_path, "r") as f:
            try:
                yaml_data = yaml.safe_load(f)
            except yaml.YAMLError as exc:
                log.critical(exc)
        return yaml_data

    @staticmethod
    def global_attribute_template():
        """
        Function to generate a template of required global attributes
        that must be set for a valid CDF.

        Returns
        -------
        template : `OrderedDict`
            A template for required global attributes that must be provided.
        """
        template = OrderedDict()
        global_attribute_schema = CDFSchema._load_default_global_attr_schema()
        default_global_attributes = CDFSchema._load_default_attributes()
        for attr_name, attr_schema in global_attribute_schema.items():
            if (
                attr_schema["required"]
                and not attr_schema["derived"]
                and attr_name not in default_global_attributes
            ):
                template[attr_name] = None
        return template

    @staticmethod
    def measurement_attribute_template():
        """
        Function to generate a template of required measurement attributes
        that must be set for a valid CDF measurement variable.

        Returns
        -------
        template: `OrderedDict`
            A template for required variable attributes that must be provided.
        """
        template = OrderedDict()
        measurement_attribute_schema = CDFSchema._load_default_variable_attr_schema()
        for attr_name, attr_schema in measurement_attribute_schema["attribute_key"].items():
            if attr_schema["required"] and not attr_schema["derived"]:
                template[attr_name] = None
        return template

    def derive_measurement_attributes(self, data, var_name):
        """
        Function to derive metadata for the given measurement.

        Parameters
        ----------
        data : `hermes_core.timedata.TimeData`
            An instance of `TimeData` to derive metadata from
        var_name : `str`
            The name of the measurement to derive metadata for

        Returns
        -------
        attributes: `OrderedDict`
            A dict containing `key: value` pairs of metadata attributes.
        """
        measurement_attributes = OrderedDict()

        # Check the Attributes that can be derived
        if not var_name == "time":
            measurement_attributes["DEPEND_0"] = self._get_depend(data)
        measurement_attributes["DISPLAY_TYPE"] = self._get_display_type(data, var_name)
        measurement_attributes["FIELDNAM"] = self._get_fieldnam(data, var_name)
        measurement_attributes["FILLVAL"] = self._get_fillval(data, var_name)
        measurement_attributes["FORMAT"] = self._get_format(data, var_name)
        measurement_attributes["LABLAXIS"] = self._get_lablaxis(data, var_name)
        measurement_attributes["SI_CONVERSION"] = self._get_si_conversion(data, var_name)
        measurement_attributes["UNITS"] = self._get_units(data, var_name)
        measurement_attributes["VALIDMIN"] = self._get_validmin(data, var_name)
        measurement_attributes["VALIDMAX"] = self._get_validmax(data, var_name)
        measurement_attributes["VAR_TYPE"] = self._get_var_type(data, var_name)
        return measurement_attributes

    def derive_time_attributes(self, data):
        """
        Function to derive metadata for the time measurement.

        Parameters
        ----------
        data : `hermes_core.timedata.TimeData`
            An instance of `TimeData` to derive metadata from.

        Returns
        -------
        attributes : `OrderedDict`
            A dict containing `key: value` pairs of time metadata attributes.
        """
        time_attributes = self.derive_measurement_attributes(data, "time")
        # Check the Attributes that can be derived
        time_attributes["REFERENCE_POSITION"] = self._get_reference_position(data)
        time_attributes["RESOLUTION"] = self._get_resolution(data)
        time_attributes["TIME_BASE"] = self._get_time_base(data)
        time_attributes["TIME_SCALE"] = self._get_time_scale(data)
        time_attributes["UNITS"] = self._get_time_units(data)
        return time_attributes

    def derive_global_attributes(self, data):
        """
        Function to derive global attributes for the given measurement data.

        Parameters
        ----------
        data : `hermes_core.timedata.TimeData`
            An instance of `TimeData` to derive metadata from.

        Returns
        -------
        attributes : `OrderedDict`
            A dict containing `key: value` pairs of global metadata attributes.
        """
        global_attributes = OrderedDict()
        # Loop through Global Attributes
        for attr_name, attr_schema in self.global_attribute_schema.items():
            if attr_schema["derived"]:
                derived_value = self._derive_global_attribute(data, attr_name=attr_name)
                global_attributes[attr_name] = derived_value
        return global_attributes

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

    # =============================================================================================
    #                             VARIABLE METADATA DERIVATIONS
    # =============================================================================================

    def _get_depend(self, data):
        return "Epoch"

    def _get_display_type(self, data, var_name):
        return "time_series"

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
            fmt = "A{}".format(len(var_data))
        else:
            raise ValueError(
                "Couldn't find FORMAT for {} of type {}".format(
                    var_name, spacepy.pycdf.lib.cdftypenames.get(cdftype, "UNKNOWN")
                )
            )
        return fmt

    def _get_lablaxis(self, data, var_name):
        # return f"{var_name} [{data[var_name].unit}]"
        return f"{var_name} [{self._get_units(data, var_name)}]"

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
            unit = var_data.unit.to_string()
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

    def _get_var_type(self, data, var_name):
        return "data"

    # =============================================================================================
    #                             GLOBAL METADATA DERIVATIONS
    # =============================================================================================

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
        day = datetime.date.today()
        return datetime.datetime(day.year, day.month, day.day)

    def _get_start_time(self, data):
        """
        Function to get the start time of the data contained in the CDF
        given in format `YYYYMMDD_hhmmss`
        """
        # Get the Start Time from the TimeSeries
        return data["time"].to_datetime()[0]

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


class NetCDFSchema(FileTypeSchema):
    """Schema for NetCDF files."""

    @property
    def global_attribute_schema(self):
        """Schema for global attributes of the file."""
        return {
            "attribute1": {"type": "string", "required": True},
            "attribute2": {"type": "int", "required": False},
            # Define more global attributes and their schemas
        }

    @property
    def variable_attribute_schema(self):
        """Schema for variable attributes of the file."""
        return {
            "variable1": {
                "attribute1": {"type": "float", "required": True},
                "attribute2": {"type": "string", "required": False},
                # Define more attributes for variable1 and their schemas
            },
            "variable2": {
                "attribute1": {"type": "int", "required": True},
                "attribute2": {"type": "string", "required": False},
                # Define more attributes for variable2 and their schemas
            },
            # Define more variables and their attribute schemas
        }
