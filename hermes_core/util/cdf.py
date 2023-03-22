"""
This module provides CDF file functionality.

"""
import datetime
from pathlib import Path
import yaml
from spacepy import pycdf

import hermes_core
from hermes_core import log
from hermes_core.util import util
from hermes_core.util.exceptions import warn_user

DEFAULT_SCHEMA_FILE = "hermes_default_schema.yaml"
DEFAULT_ATTRIBUTES_FILE = "hermes_default_attributes.yaml"


class CDFWriter:
    """
    Python Class to Convert abstract Data Structures to CDF file
    formats compliant with SpacePy's pyCDF Module. These CDF Files
    can then be indexed and shared through SPDF.

    The purpose of this class is to provide an intermediate interafce to generate CDF Files.
    One limitation of the pyCDF library is that you can not create CDF-like data structures
    in memory without generating and saving a local `.cdf` File. This class provides an
    interface to manipulate CDF-like data structures in memory before the generation and saving
    of a local `.cdf` file.

    This class uses an internal native python `dict` data structure to provide a CDF-like
    data structure that can be manipulated in memory. The internal `dict` used `self.global_attrs`
    is formatted in the following way.

    ```py
    self.global_attrs = {
        'global_attribute_name1': 'global_attribute_value1',
        'global_attribute_name2': 'global_attribute_value2',
    }

    ```
    """

    def __init__(self):

        # Current Working CDF File
        self.target_cdf_file = None

        # Data Validation, Complaiance, Derived Attributes
        self.schema = {}
        # Load Default schema Data
        self._load_default_schema_data()

        # Intermediate Data Structure for Global / Metadata Attributes
        self.global_attrs = {}
        # Load Default Global Attributes
        self._load_default_attributes()

        # Intermediate Data Structure for Variable Data
        self.variable_attrs = {}

    def _load_default_schema_data(self):
        # The Default Schema file is contained in the `hermes_core/data` directory
        default_schema_path = str(Path(hermes_core.__file__).parent / "data" / DEFAULT_SCHEMA_FILE)
        # Load the Schema
        self.load_schema_data(schema_data_path=default_schema_path)

    def _load_default_attributes(self):
        # The Default Attributes file is contained in the `hermes_core/data` directory
        default_attributes_path = str(
            Path(hermes_core.__file__).parent / "data" / DEFAULT_ATTRIBUTES_FILE
        )
        # Load the Schema
        self.add_attributes_from_yaml(attributes_path=default_attributes_path)

    def load_schema_data(self, schema_data_path):
        """
        Function to load schema for CDF formatting data from a Yaml file.

        Parameters
        ----------
        schema_data_file: `str`
            Path to schem file to be used for CDF formatting.

        """
        assert isinstance(schema_data_path, str)
        assert Path(schema_data_path).exists()
        # Load the Yaml file to Dict
        schema_data = {}
        with open(schema_data_path, "r") as f:
            try:
                schema_data = yaml.safe_load(f)
            except yaml.YAMLError as exc:
                log.critical(exc)
        # Set the Current Schema
        self.schema = schema_data

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
        example_writer.add_data_from_dict(attributes=example_data)
        ```
        """
        # Loop through attributes to add
        for attr_key, attr_value in attributes:
            # Check to see if the attribute is already present
            if attr_key in self.global_attrs.keys():
                # Log Debug that we're overriding
                log.debug(
                    "In `add_attributes_from_list()` Overriding value for attr %s from %s to %s",
                    attr_key,
                    self.global_attrs[attr_key],
                    attr_value,
                )
            # Override or set the vale of the attribute in state
            self.global_attrs[attr_key] = attr_value

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
        assert isinstance(attributes_path, str)
        assert Path(attributes_path).exists()
        # Load the Yaml file to Dict
        attribute_data = {}
        with open(attributes_path, "r") as f:
            try:
                attribute_data = yaml.safe_load(f)
            except yaml.YAMLError as exc:
                log.critical(exc)

        # Loop through attributes to add
        for attr_key, attr_value in attribute_data.items():
            # Check to see if the attribute is already present
            if attr_key in self.global_attrs.keys():
                # Log Debug that we're overriding
                log.debug(
                    "In `add_attributes_from_list()` Overriding value for attri %s from %s to %s",
                    attr_key,
                    self.global_attrs[attr_key],
                    attr_value,
                )
            # Override or set the vale of the attribute in state
            self.global_attrs[attr_key] = attr_value

    def to_cdf(self, output_path="./"):
        """
        Function to convert members of the CDFWriter's internal dict
        state to a CDF File using the pyCDF module from spacepy
        """
        assert self.target_cdf_file is None
        assert Path(output_path).exists()

        # Derive any Global Attributes
        self.derive_attributes()

        # Verify that all `required` attributes in the schema are present
        self.verify_global_attr_schema()

        # Initialize a new CDF
        cdf_filename = f"{self.global_attrs['Logical_file_id']}.cdf"
        log.info("Generating CDF: %s", cdf_filename)
        output_cdf_filepath = str(Path(output_path) / cdf_filename)
        self.target_cdf_file = pycdf.CDF(output_cdf_filepath, masterpath="")

        # Add Global Attriubtes to the CDF File
        self._convert_global_attributes_to_cdf()

        # Add zAttributes
        if "zAttrList" in self.global_attrs:
            self._convert_variable_attributes_to_cdf()

        return output_cdf_filepath

    def derive_attributes(self):
        """Function to derive global attributes in `gAttrList` member"""
        # Loop through Global Attributes
        for attr_name, attr_schema in self.schema.items():
            if attr_schema["derived"]:
                derived_value = self._derive_attribute(attr_name=attr_name, attr_schema=attr_schema)
                # Only Derive Global Attributes if they have not been manually derived/overridden
                if not self.global_attrs[attr_name]:
                    self.global_attrs[attr_name] = derived_value
                else:
                    log.debug(
                        (
                            "Attribute: %s was marked for derivation (to be %s)"
                            "but was already overridden to %s"
                        ),
                        attr_name,
                        derived_value,
                        self.global_attrs[attr_name],
                    )

    def _derive_attribute(self, attr_name, attr_schema):
        """
        Function to Derive Attributes based on other attribues in the CDF
        writer's internal dict state.
        """
        # Cache for easier access
        method = attr_schema["derived"]

        # SWITCH on the Derivation Method
        if method == "get_generation_date":
            return self._get_generation_date()
        elif method == "get_start_time":
            return self._get_start_time()
        elif method == "get_data_type":
            return self._get_data_type()
        elif method == "get_logical_file_id":
            return self._get_logical_file_id()
        elif method == "get_logical_source":
            return self._get_logical_source()
        elif method == "get_logical_source_description":
            return self._get_logical_source_description()
        else:
            raise ValueError(f"Derivation Method ({method}) Not Recognized")

    def verify_global_attr_schema(self):
        """
        Function to ensure all required attributes in the schema are present
        in the `self.global_attr` member.

        Raises `ValueError` if there are required attributes in the schema that
        are not in the global attributes.
        """
        # Loop for each attribute in the schema
        for attr_name, attr_schema in self.schema.items():
            # If it is a required attribute and not present
            if attr_schema["required"] and (attr_name not in self.global_attrs):
                raise ValueError(
                    (
                        f"Required attribute ({attr_name}) not present in global attributes.",
                        "Make sure to add all required attributes to the gloabl attributes",
                        "before genereating a CDF file.",
                    )
                )

    def _convert_global_attributes_to_cdf(self):
        # Loop though Global Attributes in target_dict
        for attr_name, attr_value in self.global_attrs.items():
            # Make sure the Value is not None
            # We cannot add None Values to the CDF Global Attrs
            if not attr_value:
                warn_user((f"Cannot Add gAttr: {attr_name}. Value was {str(attr_value)}. "))
                # Set to a intentionally invalid value
                attr_value = f"{attr_name}-Not-Provided"
            # Add the Attribute to the CDF File
            self.target_cdf_file.attrs[attr_name] = attr_value

    def _convert_variable_attributes_to_cdf(self):
        # Loop through Variable Attributes in target_dict
        for attr_name, attr_value in self.variable_attrs.items():
            # Make sure the Value is not None
            # We cannot add None Values to the CDF Vars
            if not attr_value:
                warn_user((f"Cannot Add zAttr: {attr_name}. Value was {str(attr_value)}. "))
                # Set to a intentionally invalid value
                attr_value = f"{attr_name}-Not-Provided"
            # Add the Attribute to the CDF File
            self.target_cdf_file[attr_name] = attr_value

    def save_cdf(self):
        """Function to save and close CDF File"""
        assert self.target_cdf_file is not None
        # Save the CDF
        self.target_cdf_file.save()
        self.target_cdf_file.close()

    # ==============================================================================================
    #                               ATTRIBUTE DERIVATION FUNCTIONS
    # ==============================================================================================

    def _get_logical_file_id(self):
        """
        Function to get the `Logical_file_id` required global attribute.

        The attribute stores the name of the CDF File without the file
        extension (e.g. '.cdf'). This attribute is requires to avoid
        loss of the originial source in case of renaming.
        """
        if not self.global_attrs["Logical_file_id"]:
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
            science_filename = self.global_attrs["Logical_file_id"]
        return science_filename

    def _get_logical_source(self):
        """
        Function to get the `Logical_source` required global attribute.

        This attribute determines the file naming convention in the SKT Editor
        and is used by CDA Web.
        """
        if not self.global_attrs["Logical_source"]:
            # Get Parts
            spacecraft_id = self._get_spacecraft_id()
            instrument_id = self._get_instrument_id()
            data_type = self._get_data_type()
            data_type_short_name, _ = data_type.split(">")
            # Build Derivation
            logical_source = f"{spacecraft_id}_{instrument_id}_{data_type_short_name}"
        else:
            logical_source = self.global_attrs["Logical_source"]
        return logical_source

    def _get_logical_source_description(self):
        """
        Function to get the `Logical_source_description` required global attribute.

        This attribute writes out the full words associated with the encryped
        `Logical_source`  attribute.
        """
        if not self.global_attrs["Logical_source_description"]:
            # Get Parts
            spacecraft_long_name = self._get_spacecraft_long_name()
            instrument_long_name = self._get_instrument_long_name()
            data_level_long_name = self._get_data_level_long_name()
            logical_source_description = (
                f"{data_level_long_name} {spacecraft_long_name} {instrument_long_name}"
            )
        else:
            logical_source_description = self.global_attrs["Logical_source_description"]
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
        if not self.global_attrs["Data_type"]:
            # Get Short Parts
            mode = self._get_instrument_mode()
            data_level = self._get_data_level()
            odpd = self._get_data_product_descriptor()

            # Get Long Parts
            mode_long_name = self._get_instrument_mode().upper()  # NOTE Seems to be Upper case?
            data_level_long_name = self._get_data_level_long_name()
            odpd_long_name = (
                self._get_data_product_descriptor().upper()
            )  # NOTE Seems to be Uppar case?

            # Build Derivation (With ODPD)
            if odpd:
                data_type = (
                    f"{mode}_{data_level}_{odpd}>"
                    f"{mode_long_name} {data_level_long_name} {odpd_long_name}"
                )
            else:
                data_type = f"{mode}_{data_level}>{mode_long_name} {data_level_long_name}"
        else:
            data_type = self.global_attrs["Data_type"]
        return data_type

    def _get_spacecraft_id(self):
        """Function to get Spacecraft ID from Source_name Global Attribute"""
        sc_id = self.global_attrs["Source_name"]
        assert sc_id is not None
        # Formatting
        if ">" in sc_id:
            short_name, _ = sc_id.split(">")
            sc_id = short_name.lower()  # Makse sure its all lowercase
        return sc_id

    def _get_spacecraft_long_name(self):
        """Function to get Spacecraft ID from Source_name Global Attribute"""
        sc_id = self.global_attrs["Source_name"]
        assert sc_id is not None
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
        instr_id = self.global_attrs["Descriptor"]
        assert instr_id is not None
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
        instr_id = self.global_attrs["Descriptor"]
        assert instr_id is not None
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
        data_level = self.global_attrs["Data_level"]
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
        data_level = self.global_attrs["Data_level"]
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
        if "Data_product_descriptor" in self.global_attrs:
            odpd = self.global_attrs["Data_product_descriptor"]
        else:
            odpd = ""
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
        # TODO derrive from epoch
        return datetime.datetime.now()

    def _get_version(self):
        """
        Function to get the 3-part version number of the data product.
        """
        version_str = self.global_attrs["Data_version"].lower()
        assert version_str is not None
        if "v" in version_str:
            _, version = version_str.split("v")
        else:
            version = version_str
        return version

    def _get_instrument_mode(self):
        """Function to get the mode attribute (TBS)"""
        instr_mode = self.global_attrs["Instrument_mode"]
        assert instr_mode is not None
        return instr_mode.lower()  # Makse sure its all lowercase
