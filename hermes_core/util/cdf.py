"""
This module provides CDF file functionality.

"""
import datetime
from pathlib import Path
import yaml
import numpy as np
from spacepy.pycdf import CDF
from spacepy.pycdf.istp import FileChecks

from astropy.timeseries import TimeSeries
from astropy.time import Time
from astropy.table import Column

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

    def __init__(self):
        # Current Working CDF File
        self.cdf = None

        # Intermediate Data Structure
        self.data = TimeSeries()

        # Load Default Global Attributes
        self._load_default_attributes()

        # Data Validation, Complaiance, Derived Attributes
        self.global_attr_schema = self._load_default_global_attr_schema()

        # Data Validation and Compliance for Variable Data
        self.variable_attr_schema = self._load_default_variable_attr_schema()

    def __del__(self):
        """
        Deconstructor for CDFWriter class.
        """
        if self.cdf:
            self.cdf.close()
        self.cdf = None

    def __repr__(self):
        """
        Returns a representation of the CDFWriter class.
        """
        return self.__str__()

    def __str__(self):
        """
        Returns a string representation of the CDFWriter class.
        """
        str_repr = (
            f"CDFWriter: Converted: {self.cdf is not None}"
            f"\nGlobal Attrs:\n{self.data.meta}"
            f"\nVariable Data:\n{self.data}"
            # f"\nVariable Attributes: {self.variable_attrs}"
        )
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
        if name == "time":
            self.add_time(time=data, time_attrs={})
        else:
            self.add_variable(var_name=name, var_data=data, var_attrs={})

    def __contains__(self, name):
        """
        Function to see whether a data variable is in the CDFWriter class.
        """
        return name in self.data.columns

    def __iter__(self):
        """
        Function to iterate over data variables and attributes in the CDFWriter class.
        """
        for var_name in self.data.keys():
            var_data = self.data[var_name]
            yield (var_name, var_data)

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

    def add_attributes_from_list(self, attributes: list):
        """
        Function to add data to the CDFWriter's internal Dict
        state by passing a native list type Object.

        Parameters
        ----------
        attributes: `list`
            A native python `list` object containing `tuples` of `(key, value)`
            attributes.

        Examples
        -------
        ```
        example_data = [
            ("attribute_example_1", "attribute_value_1"),
            ("attribute_example_2", "attribute_value_2")
        ]
        example_writer = CDFWriter()
        # Add the Example Data
        example_writer.add_data_from_list(attributes=example_data)
        ```
        """
        # Loop through attributes to add
        for attr_key, attr_value in attributes:
            # Check to see if the attribute is already present
            if attr_key in self.data.meta.keys():
                # Log Debug that we're overriding
                log.debug(
                    "In `add_attributes_from_list()` Overriding value for attr %s from %s to %s",
                    attr_key,
                    self.data.meta[attr_key],
                    attr_value,
                )
            # Override or set the vale of the attribute in state
            self.data.meta[attr_key] = attr_value

    def add_attributes_from_dict(self, attributes: dict):
        """
        Function to add data to the CDFWriter's internal Dict
        state by passing a native list type Object.

        Parameters
        ----------
        attributes: `dict`
            A native python `dict` object containing `(key, value)`
            attributes.

        Examples
        -------
        ```
        example_data = {
            "attribute_example_1": "attribute_value_1",
            "attribute_example_2": "attribute_value_2"
        }
        example_writer = CDFWriter()
        # Add the Example Data
        example_writer.add_data_from_dict(attributes=example_data)
        ```
        """
        # Loop through attributes to add
        for attr_key, attr_value in attributes.items():
            # Check to see if the attribute is already present
            if attr_key in self.data.meta.keys():
                # Log Debug that we're overriding
                log.debug(
                    "In `add_attributes_from_list()` Overriding value for attr %s from %s to %s",
                    attr_key,
                    self.data.meta[attr_key],
                    attr_value,
                )
            # Override or set the vale of the attribute in state
            self.data.meta[attr_key] = attr_value

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
        self.add_attributes_from_dict(attributes=attribute_data)

    def add_time(self, time: Time, time_attrs: dict):
        """
        Function to add required `time` column to the CDFWriter class.

        Parameters
        ----------
        time: `astropy.time.Time`
            An `astropy.time.Time` object that contains the time dimension for the CDF
            object to be created.

        time_attrs: `dict`
            A collection of `(key,value)` pairs to add as attributes of the CDF varaible.
        """
        # Add the Time Coumn
        self.add_variable(var_name="time", var_data=time.to_datetime(), var_attrs=time_attrs)

    def add_variable(self, var_name: str, var_data: np.ndarray, var_attrs: dict):
        """
        Function to add variable data to the CDFWriter class. Variable data here is assumed
        to be array-like or matrix-like to be stored in the CDF. Additionally, varaible
        attributes can be added though a native python `dict` of `(key, value)` pairs that
        are added to the CDF variable.

        Parameters
        ----------
        var_name: `str`
            Name of the variable to add to the CDFWriter.

        var_data: `np.ndarray`
            Array-like or matrix-like data to add to the CDFWriter.

        var_attrs: `dict`
            A collection of `(key,value)` pairs to add as attributes of the CDF varaible.

        """
        # Add the Variable as a new Column
        var_column = Column(data=var_data, name=var_name, meta=var_attrs)
        self.data.add_column(col=var_column)

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

        # Initialize a new CDFWriter object
        writer = CDFWriter()

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
            writer.add_attributes_from_dict(input_global_attrs)

            # First Variable we need to add is time/Epoch
            if "Epoch" in input_file:
                time_data = Time(input_file["Epoch"][:].copy())
                time_attrs = {}
                for attr_name in input_file["Epoch"].attrs:
                    time_attrs[attr_name] = input_file["Epoch"].attrs[attr_name]
                writer.add_time(time_data, time_attrs)

            # Add Variable Attributtes from the CDF file to TimeSeries
            for var_name in input_file:
                if var_name != "Epoch":  # Since we added this separately
                    var_data = input_file[var_name][:].copy()
                    var_attrs = {}
                    for attr_name in input_file[var_name].attrs:
                        var_attrs[attr_name] = input_file[var_name].attrs[attr_name]
                    writer.add_variable(var_name, var_data, var_attrs)

        # Return the created CDFWriter objet
        return writer

    @staticmethod
    def from_timeseries(ts: TimeSeries):
        """
        Function to create a CDFWriter object from an existing CDF File. This allows
        someone to add data to a skeleton CDF and save the resulting data struture to
        a new, validated CDF file.

        Parameters
        ----------
        ts: `TimeSeries`
            A TimeSeries object to create a CDFWriter from

        """
        # Initialize a new CDFWriter object
        writer = CDFWriter()

        # Add Global Attributes
        writer.add_attributes_from_dict(ts.meta)

        # Add the Time Data
        writer.add_variable("time", ts["time"].data, ts["time"].meta)

        # Add Columns
        for col in ts.itercols():
            if col.name != "time":
                writer.add_variable(col.name, col.data, col.meta)

        # Return the new Writer
        return writer

    # =============================================================================================
    #                         CDF FILE GENERATION FUNCTIONS
    # =============================================================================================

    def to_cdf(self, output_path="./"):
        """
        Function to convert members of the CDFWriter's internal dict
        state to a CDF File using the pyCDF module from spacepy

        Parameters
        ----------
        output_path: `str`
            A string path to the directory where the resulting CSF file is to be saved.
            By default the CDF file is saved to the current working directory.

        """
        assert self.cdf is None
        assert Path(output_path).exists()

        # Derive any Global Attributes
        self.derive_attributes()

        # Initialize a new CDF
        cdf_filename = f"{self.data.meta['Logical_file_id']}.cdf"
        log.info("Generating CDF: %s", cdf_filename)
        output_cdf_filepath = str(Path(output_path) / cdf_filename)
        self.cdf = CDF(output_cdf_filepath, masterpath="")

        # Add Global Attriubtes to the CDF File
        self._convert_global_attributes_to_cdf()

        # Add zAttributes
        self._convert_variable_attributes_to_cdf()

        return output_cdf_filepath

    def derive_attributes(self):
        """Function to derive global attributes in `gAttrList` member"""
        # Loop through Global Attributes
        for attr_name, attr_schema in self.global_attr_schema.items():
            if attr_schema["derived"]:
                derived_value = self._derive_attribute(attr_name=attr_name)
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

    def _derive_attribute(self, attr_name):
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

    def _convert_global_attributes_to_cdf(self):
        # Loop though Global Attributes in target_dict
        for attr_name, attr_value in self.data.meta.items():
            # Make sure the Value is not None
            # We cannot add None Values to the CDF Global Attrs
            if not attr_value:
                warn_user(f"Cannot Add gAttr: {attr_name}. Value was {str(attr_value)} ")
                # Set to a intentionally invalid value
                # attr_value = ""
            else:
                # Add the Attribute to the CDF File
                self.cdf.attrs[attr_name] = attr_value

    def _convert_variable_attributes_to_cdf(self):
        # Loop through Variable Attributes in target_dict
        for var_name, var_data in self.__iter__():
            if var_name == "time":
                # Add 'time' in the TimeSeries as 'Epoch' within the CDF
                self.cdf["Epoch"] = var_data.data
                # Add the Variable Attributes
                for var_attr_name, var_attr_val in var_data.meta.items():
                    self.cdf["Epoch"].attrs[var_attr_name] = var_attr_val
            else:
                # Add the Variable to the CDF File
                self.cdf[var_name] = var_data.data
                # Add the Variable Attributes
                for var_attr_name, var_attr_val in var_data.meta.items():
                    self.cdf[var_name].attrs[var_attr_name] = var_attr_val

    def save_cdf(self):
        """Function to save and close CDF File"""
        assert self.cdf is not None
        # Save the CDF
        self.cdf.save()
        self.cdf.close()
        # Reset the Local Member
        self.cdf = None

    # =============================================================================================
    #                         DATA VALIDATION / VERIFICATION FUNCTIONS
    # =============================================================================================

    def validate_cdf(self, catch=False):
        """
        Function to validate the CDF file generated after calling `to_cdf()`.
        The function executes the `pycdf.istp.FileChecks` tests to determine any
        of the metadata properties contained in the file are incorrect or if there
        are missing metadata properties in the CDF file.

        Parameters
        ----------
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

        # Verify that all `required` global attributes in the schema are present
        global_attr_validation_errors = self.validate_global_attr_schema()
        validation_errors.extend(global_attr_validation_errors)

        # Verify that all `required` variable attributes in the schema are present
        variable_attr_validation_errors = self.validate_variable_attr_schema()
        validation_errors.extend(variable_attr_validation_errors)

        # Validate the CDF Using ISTP Module
        istp_validation_errors = FileChecks.all(f=self.cdf, catch=catch)
        validation_errors.extend(istp_validation_errors)

        return validation_errors

    def validate_global_attr_schema(self):
        """
        Function to ensure all required global attributes in the schema are present
        in the generated CDF File.
        """
        global_attr_validation_errors = []
        # Loop for each attribute in the schema
        for attr_name, attr_schema in self.global_attr_schema.items():
            # If it is a required attribute and not present
            if attr_schema["required"] and (attr_name not in self.cdf.attrs):
                global_attr_validation_errors.append(
                    f"Required attribute ({attr_name}) not present in global attributes.",
                )
        return global_attr_validation_errors

    def validate_variable_attr_schema(self):
        """
        Function to ensure all required variable attributes in the schema are present
        in the generated CDF file.
        """
        variable_attr_validation_errors = []

        # Loop for each Variable in the CDF File
        for var_name in self.cdf:
            # Get the `Var()` Class for the Variable
            var_data = self.cdf[var_name]

            # Get the Variable Type to compare the required attributes
            var_type = ""
            if "VAR_TYPE" in var_data.attrs:
                var_type = var_data.attrs["VAR_TYPE"]
                variable_errors = self._validate_variable(var_name, var_type)
                variable_attr_validation_errors.extend(variable_errors)
            else:
                variable_attr_validation_errors.append(
                    f"Variable: {var_name} missing 'VAR_TYPE' attribute. Cannot Validate Variable."
                )

        return variable_attr_validation_errors

    def _validate_variable(self, var_name, var_type):
        """
        Function to Validate an individual Variable.
        """
        variable_errors = []
        # Get the Expected Attributes for the Variable Type
        var_type_attrs = self.variable_attr_schema[var_type]

        # Get the `Var()` Class for the Variable
        var_data = self.cdf[var_name]

        # Loop for each Variable Attribute in the schema
        for attr_name in var_type_attrs:
            attr_schema = self.variable_attr_schema["attribute_key"][attr_name]
            # If it is a required attribute and not present
            if attr_schema["required"] and attr_name not in var_data.attrs:
                variable_errors.append(f"Variable: {var_name} missing '{attr_name}' attribute.")

        return variable_errors

    # =============================================================================================
    #                               ATTRIBUTE DERIVATION FUNCTIONS
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
            start_time = self.data["time"][0]
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
