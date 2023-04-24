"""
This module provides CDF file functionality.

"""
from typing import Union, OrderedDict
import datetime
from pathlib import Path
import yaml
from copy import deepcopy
import numpy as np
import spacepy
from spacepy.pycdf import CDF, _Hyperslice
from spacepy.pycdf.istp import FileChecks

from astropy.timeseries import TimeSeries
from astropy.time import Time
from astropy.table import Column
from astropy.units import Quantity

import hermes_core
from hermes_core import log
from hermes_core.util import util
from hermes_core.util.exceptions import warn_user

DEFAULT_GLOBAL_CDF_ATTRS_SCHEMA_FILE = "hermes_default_global_cdf_attrs_schema.yaml"
DEFAULT_GLOBAL_CDF_ATTRS_FILE = "hermes_default_global_cdf_attrs.yaml"
DEFAULT_VARIABLE_CDF_ATTRS_SCHEMA_FILE = "hermes_default_variable_cdf_attrs_schema.yaml"


class CDFWriter:
    """
    Python Class to provide an intermediate interafce to generate CDF files
    vlidated against ISTP standards. One limitation of the `spacepy.pycdf` library is
    that you can not create CDF-like data structures in memory without generating and saving
    a local `.cdf` File. This class provides an interface to manipulate CDF-like data
    structures in memory before the generation and saving of a local `.cdf` file.

    The class uses an `atropy.timeseries.TimeSeries` table to provide an intermediate
    interfce for modifying CDF-like data structure. The `TimeSeries` object can contain
    global metadata, single or multidimensional time series data, and variable metadata
    about each column/variable member in the `TimeSeries`. This intermediate structure is
    converted to a `spacepy.pycdf.CDF` data structure through the `to_cdf()` function. This
    function executes attribute derivations for required global attributes, if they have
    not already been overwritten and the individual sub-attributes are present. Global
    metadata is extracted from the `OrderedDict()` `.meta` property of the `TimeSeries`
    and added as `(key, value)` pairs in the `CDF` data structure. Variable data and
    attributes are converted from the `Column` `.meta` and `.data` members respectively to
    the `CDF` variable data and attributes respectively.



    """

    def __init__(self, ts: TimeSeries):
        # Verify that the input TimeSeries has all needed Columns
        if len(ts.columns) == 0:
            raise ValueError(
                (
                    "Input TimeSeries cannot be empty. ",
                    "Must contain at minimum 'time' and one measurement column.",
                )
            )
        # Verify that the `time` column is present
        if "time" not in ts.columns:
            raise ValueError(
                "TimeSeries must contain at minimum 'time' and one measurement column."
            )
        # Verify that the `time` column is a `Time` type
        if not isinstance(ts.time, Time):
            raise TypeError("TimeSeries 'time' column must be type `astropy.time.Time`.")
        # Verify that all columns are `Quantity`
        for col in ts.columns:
            if col != "time":
                if (not isinstance(ts[col], Quantity)) or (not ts[col].unit):
                    raise TypeError(
                        f"Column {col} must be type `astropy.units.Quantity` and have `unit` assigned."
                    )

        # Intermediate Data Structure
        self.data = ts.copy()

        # Load Default Global Attributes
        self._load_default_attributes()

        # Add any Metadata from the original TimeSeries
        self.data.time.meta = OrderedDict()
        if hasattr(ts.time, "meta"):
            self.data.time.meta.update(ts.time.meta)
        for col in self.data.columns:
            if col != "time":
                self.data[col].meta = OrderedDict()
                if hasattr(ts[col], "meta"):
                    self.data[col].meta.update(ts[col].meta)

        # Data Validation, Complaiance, Derived Attributes
        self._global_attr_schema = self._load_default_global_attr_schema()

        # Data Validation and Compliance for Variable Data
        self._variable_attr_schema = self._load_default_variable_attr_schema()

    @property
    def meta(self):
        """
        (OrderedDict), returns the metadata contained in the TimeSeries
        """
        return self.data.meta

    @meta.setter
    def meta(self, value):
        self.data.meta = value

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
        str_repr = f"CDFWriter() Object:\n"
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
        self._add_variable(var_name=name, var_data=data, var_attrs={})

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

    def _load_default_global_attr_schema(self) -> dict:
        # The Default Schema file is contained in the `hermes_core/data` directory
        default_schema_path = str(
            Path(hermes_core.__file__).parent / "data" / DEFAULT_GLOBAL_CDF_ATTRS_SCHEMA_FILE
        )
        # Load the Schema
        return self._load_yaml_data(yaml_file_path=default_schema_path)

    def _load_default_variable_attr_schema(self) -> dict:
        # The Default Schema file is contained in the `hermes_core/data` directory
        default_schema_path = str(
            Path(hermes_core.__file__).parent / "data" / DEFAULT_VARIABLE_CDF_ATTRS_SCHEMA_FILE
        )
        # Load the Schema
        return self._load_yaml_data(yaml_file_path=default_schema_path)

    def _load_default_attributes(self):
        # The Default Attributes file is contained in the `hermes_core/data` directory
        default_attributes_path = str(
            Path(hermes_core.__file__).parent / "data" / DEFAULT_GLOBAL_CDF_ATTRS_FILE
        )
        self.add_attributes_from_yaml(attributes_path=default_attributes_path)

    def _load_yaml_data(self, yaml_file_path):
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

    def add_attributes_from_yaml(self, attributes_path: str):
        """
        Function to add attributes to the CDFWriter's internal Dict
        state through a yaml file containing (key, value) pairs of
        attribute names and values.

        Parameters
        ----------
        attributes_path: `str`
            Path to a yaml file containing (key, value) pairs to be
            added as attributes.

        """
        # Load the Yaml file to Dict
        attribute_data = self._load_yaml_data(attributes_path)

        # Add the Attributes through the Dict
        self.data.meta.update(attribute_data)

    def _add_variable(self, var_name: str, var_data: Quantity, var_attrs: dict):
        """
        Function to add variable data to the CDFWriter class. Variable data here is assumed
        to be array-like or matrix-like to be stored in the CDF. Additionally, varaible
        attributes can be added though a native python `dict` of `(key, value)` pairs that
        are added to the CDF variable.

        Parameters
        ----------
        var_name: `str`
            Name of the variable to add to the CDFWriter.

        var_data: `Quantity`
                the data added to the internal `TimeSeries` with the `var_name` as the
                column name and `var_attrs` as the column metadata. Data contained in the
                `Quantity.info` member is combined with the `var_attrs` to add additional
                attributes or update and override existing attributes.

        var_attrs: `dict`
            A collection of `(key,value)` pairs to add as attributes of the CDF varaible.

        """
        # Verify that all columns are `Quantity`
        if (not isinstance(var_data, Quantity)) or (not var_data.unit):
            raise TypeError(
                f"Column {var_name} must be type `astropy.units.Quantity` and have `unit` assigned."
            )

        self.data[var_name] = var_data
        # Add any Metadata from the original TimeSeries
        self.data[var_name].meta = OrderedDict()
        if hasattr(var_data, "meta"):
            self.data[var_name].meta.update(var_attrs)

        # Derive any Attributes that can be Derived
        self._derive_variable_attributes(var_name=var_name)

    # =============================================================================================
    #                         CDF FILE LOADING FUNCTIONS
    # =============================================================================================

    @staticmethod
    def from_cdf(pathname: str):
        """
        Function to create a CDFWriter object from an existing CDF File. This allows
        someone to add data to a skeleton CDF and save the resulting data struture to
        a new, validated CDF file.

        Parameters
        ----------
        pathname: `str`
            A string path of the file to open

        """
        if not Path(pathname).exists():
            raise FileNotFoundError(f"CDF Could not be loaded from path: {pathname}")

        # Create a new TimeSeries
        ts = TimeSeries()

        # Open CDF file with context manager
        with CDF(pathname) as input_file:
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

        # Initialize a new CDFWriter object
        writer = CDFWriter(ts)

        # Return the created CDFWriter objet
        return writer

    # =============================================================================================
    #                         CDF FILE GENERATION FUNCTIONS
    # =============================================================================================

    def write_cdf(self, output_path="./"):
        """
        Function to convert members of the CDFWriter's internal dict
        state to a CDF File using the pyCDF module from spacepy. This function
        saves and closes the file pointer to the CDF file after converting all
        data from the internal `TimeSeries` to the `spacepy.pycdf.CDF` format.

        Parameters
        ----------
        output_path: `str`
            A string path to the directory where the resulting CSF file is to be saved.
            By default the CDF file is saved to the current working directory.

        """
        # Derive any Global Attributes
        self.derive_attributes()

        # Initialize a new CDF
        cdf_filename = f"{self.data.meta['Logical_file_id']}.cdf"
        log.info("Generating CDF: %s", cdf_filename)
        output_cdf_filepath = str(Path(output_path) / cdf_filename)
        with CDF(output_cdf_filepath, masterpath="") as cdf_file:
            # Add Global Attriubtes to the CDF File
            self._convert_global_attributes_to_cdf(cdf_file)

            # Add zAttributes
            self._convert_variable_attributes_to_cdf(cdf_file)

        return output_cdf_filepath

    def _convert_global_attributes_to_cdf(self, cdf_file: CDF):
        # Loop though Global Attributes in target_dict
        for attr_name, attr_value in self.data.meta.items():
            # Make sure the Value is not None
            # We cannot add None Values to the CDF Global Attrs
            if not attr_value:
                warn_user(f"Cannot Add gAttr: {attr_name}. Value was {str(attr_value)} ")
            else:
                # Add the Attribute to the CDF File
                cdf_file.attrs[attr_name] = attr_value

    def _convert_variable_attributes_to_cdf(self, cdf_file: CDF):
        # Loop through Variable Attributes in target_dict
        for var_name, var_data in self.__iter__():
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

    # =============================================================================================
    #                         DATA VALIDATION / VERIFICATION FUNCTIONS
    # =============================================================================================

    def validate_cdf(self, cdf_file_path, catch=True):
        """
        Function to validate a CDF file.
        The function executes the `pycdf.istp.FileChecks` tests to determine any
        of the metadata properties contained in the file are incorrect or if there
        are missing metadata properties in the CDF file.

        Parameters
        ----------
        cdf_file_path: `str`
            A string path to the CDF File to be validated.

        catch: `bool`
            A bool value for whether to catch errors in checking the validation of the
            generated CDF file. If `True` any excetions in validation checking will
            throw an error. If `False` any exceptions will not be raised, but the test
            failure will be recorded and returned.

        Returns
        -------
        result : `list`
            A `list` containing a descritption of each vlidation failure from the
            `FileChecks` tests. Returns an empty list if compliant.

        """
        # Initialize Validation Errrors
        validation_errors = []

        try:
            # Open CDF file with context manager
            with CDF(cdf_file_path) as cdf_file:
                # Verify that all `required` global attributes in the schema are present
                global_attr_validation_errors = self.validate_global_attr_schema(cdf_file=cdf_file)
                validation_errors.extend(global_attr_validation_errors)

                # Verify that all `required` variable attributes in the schema are present
                variable_attr_validation_errors = self.validate_variable_attr_schema(
                    cdf_file=cdf_file
                )
                validation_errors.extend(variable_attr_validation_errors)

                # Validate the CDF Using ISTP Module
                istp_validation_errors = FileChecks.all(f=cdf_file, catch=catch)
                validation_errors.extend(istp_validation_errors)

        except IOError:
            validation_errors.append(f"Could not open CDF File at path: {cdf_file_path}")

        return validation_errors

    def validate_global_attr_schema(self, cdf_file: CDF):
        """
        Function to ensure all required global attributes in the schema are present
        in the generated CDF File.
        """
        global_attr_validation_errors = []
        # Loop for each attribute in the schema
        for attr_name, attr_schema in self._global_attr_schema.items():
            # If it is a required attribute and not present
            if attr_schema["required"] and (attr_name not in cdf_file.attrs):
                global_attr_validation_errors.append(
                    f"Required attribute ({attr_name}) not present in global attributes.",
                )
        return global_attr_validation_errors

    def validate_variable_attr_schema(self, cdf_file: CDF):
        """
        Function to ensure all required variable attributes in the schema are present
        in the generated CDF file.
        """
        variable_attr_validation_errors = []

        # Loop for each Variable in the CDF File
        for var_name in cdf_file:
            # Get the `Var()` Class for the Variable
            var_data = cdf_file[var_name]

            # Get the Variable Type to compare the required attributes
            var_type = ""
            if "VAR_TYPE" in var_data.attrs:
                var_type = var_data.attrs["VAR_TYPE"]
                variable_errors = self._validate_variable(cdf_file, var_name, var_type)
                variable_attr_validation_errors.extend(variable_errors)
            else:
                variable_attr_validation_errors.append(
                    f"Variable: {var_name} missing 'VAR_TYPE' attribute. Cannot Validate Variable."
                )

        return variable_attr_validation_errors

    def _validate_variable(self, cdf_file: CDF, var_name: str, var_type: str):
        """
        Function to Validate an individual Variable.
        """
        variable_errors = []
        # Get the Expected Attributes for the Variable Type
        var_type_attrs = self._variable_attr_schema[var_type]

        # Get the `Var()` Class for the Variable
        var_data = cdf_file[var_name]

        # Loop for each Variable Attribute in the schema
        for attr_name in var_type_attrs:
            attr_schema = self._variable_attr_schema["attribute_key"][attr_name]
            # If it is a required attribute and not present
            if attr_schema["required"] and attr_name not in var_data.attrs:
                variable_errors.append(f"Variable: {var_name} missing '{attr_name}' attribute.")
            else:
                # If the Var Data can be Validated
                if "valid_values" in attr_schema:
                    attr_valid_values = attr_schema["valid_values"]
                    attr_value = var_data.attrs[attr_name]
                    if attr_value not in attr_valid_values:
                        variable_errors.append(
                            (
                                f"Variable: {var_name} Attribute '{attr_name}' not one of valid options.",
                                f"Was {attr_value}, expected one of {attr_valid_values}",
                            )
                        )
        return variable_errors

    # =============================================================================================
    #                               ATTRIBUTE DERIVATION FUNCTIONS
    # =============================================================================================

    def derive_attributes(self):
        """Function to derive global attributes"""
        # Loop through Global Attributes
        for attr_name, attr_schema in self._global_attr_schema.items():
            if attr_schema["derived"]:
                derived_value = self._derive_global_attribute(attr_name=attr_name)
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
        # Loop through Variable Attributes
        for name in self.data.columns:
            self._derive_variable_attributes(var_name=name)

    def _derive_variable_attributes(self, var_name):
        var_data = self.data[var_name]

        # Derive Time-Specific Attributes
        if var_name == "time":
            self._derive_time_attributes()

        # Check the Attributes that can be derived
        if "DEPEND_0" not in var_data.meta:
            var_data.meta["DEPEND_0"] = self._get_depend()
        if "FIELDNAM" not in var_data.meta:
            var_data.meta["FIELDNAM"] = self._get_fieldnam(var_name)
        if "FILLVAL" not in var_data.meta:
            var_data.meta["FILLVAL"] = self._get_fillval(var_name)
        if "UNITS" not in var_data.meta:
            var_data.meta["UNITS"] = self._get_units(var_name)
        if "VALIDMIN" not in var_data.meta:
            var_data.meta["VALIDMIN"] = self._get_validmin(var_name)
        if "VALIDMAX" not in var_data.meta:
            var_data.meta["VALIDMAX"] = self._get_validmax(var_name)

    def _derive_time_attributes(self):
        var_data = self.data.time

        # Check the Attributes that can be derived
        if "REFERENCE_POSITION" not in var_data.meta:
            var_data.meta["REFERENCE_POSITION"] = self._get_reference_position()
        if "RESOLUTION" not in var_data.meta:
            var_data.meta["RESOLUTION"] = self._get_resolution()
        if "TIME_BASE" not in var_data.meta:
            var_data.meta["TIME_BASE"] = self._get_time_base()
        if "TIME_SCALE" not in var_data.meta:
            var_data.meta["TIME_SCALE"] = self._get_time_scale()
        if "UNITS" not in var_data.meta:
            var_data.meta["UNITS"] = self._get_time_units()

    def _derive_global_attribute(self, attr_name):
        """
        Function to Derive Attributes based on other attribues in the CDF
        writer's internal dict state.
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
