"""
This module provides CDF file functionality.

"""
import os
import datetime
import logging
from pathlib import Path
import yaml
from spacepy import pycdf
from flatten_dict import flatten, unflatten

import hermes_core
from hermes_core.util import util


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

    This class usses an internal native python `dict` data structure to provide a CDF-like
    data structure that can be maniputed in memory. The internal `dict` used `self.target_dict`
    is formatted in the following way.

    ```py
    self.target_dict = {
        'gAttrList`: {
            'global_attribute_name': {
                'derived' False,
                'required': False,
                'value': 'global_attribute_value'
            }
        }
        'zAttrList`: {

        }
    }
    ```

    Where:
    - `gAttrList` is a python `dict` object used to represent global attributes embedded in a
        CDF file. The global attributes list is keyed by names of global attributes that shoule
        be added to the CDF. Each attribute contains another python `dict` object that is used
        to configure options and metadata about that attribute.

        NOTE: This is intended to behave like the `pyCDF` `gAttrList` data structure.
    - `zAttrList` is a python `dict` object used to represent variable attributes that contain all
        telemetry or science data within the CDF file. the variable attributes list is keyed by the
        names of variable attributes where each attribute contains another python `dict` object that
        configures options and metadata about the variable attribute, as well as contain the actual
        telemetry or science data associated with the given variable.

        NOTE: This intended to behave like the `pyCDF` `zAttrList` data structure.

    `CDFWriter` objects are instantiated using:
        - `seed_data`: a python `dict` object of the format above, containing global and variable
            attributes to be added to a CDF file.
        - `template_data_file`: the path to a template dictionary-like, CDF-like data structure
            that is loaded from a `template.yaml` file. The file represents the nested dctionary
            data structures that are contained in the `self.target_dict` member.
            NOTE: Defaults to `cdf_template.yml` found in the module's `data` directory.
    """

    def __init__(self, seed_data=None, template_data_file="cdf_template.yaml"):

        # Current Working CDF File
        self.target_cdf_file = None

        # Inter Data Structure before Creating CDF File
        self.target_dict = {}

        # Load Default Template Data
        template_file = Path(hermes_core.__file__).parent / "data" / template_data_file
        default_template_data = {}
        with open(template_file, "r") as f:
            try:
                default_template_data = yaml.safe_load(f)
            except yaml.YAMLError as exc:
                logging.critical(exc)
        # Add TEMPLATE to target_dict
        self.add_data_from_dict(default_template_data)

        # If Present, add seed_data to target_dict
        if seed_data:
            assert type(seed_data) is dict
            self.add_data_from_dict(seed_data)

    def add_data_from_dict(self, data):
        """
        Function to add data to the CDFWriter's internal Dict
        state by passing a native Dict type Object.

        Parameters
        ----------
        data: `dict`
            A native python dict object containing global attributes and variables.

        Examples
        -------
        ```
        example_data = {
            "gAttrList": {
                "gAttrExample": {"value": None, "required": True, "valid_check": None, "derived": False},
            },
            "zAttrList": {
                "varExample": {"value": None, "required": True, "valid_check": None, "derived": False},
            },
        }
        example_writer = CDFWriter()
        # Add the Example Data
        example_writer.add_data_from_dict(data=example_data)
        ```
        """
        # Crete new Template
        compiled = {}
        # Add the Current Target Dict First
        compiled.update(flatten(self.target_dict))
        # Add the New Template Data
        compiled.update(flatten(data))
        # Reset the Target Dict
        self.target_dict = unflatten(compiled)

    def derive_attributes(self):
        """Function to derive global attributes in `gAttrList` member"""
        # Loop through Global Attributes
        for attr_name, info in self.target_dict["gAttrList"].items():
            # Only Derive Global Attributes if they have not been manually derived/overridden
            if info["derived"]:
                derived_value = self._derive_attribute(attr_name=attr_name, info=info)
                if not info["value"]:
                    info["value"] = derived_value
                else:
                    logging.debug(
                        (
                            "Attribute: %s was marked for derivation (to be %s)"
                            "but was already overridden to %s"
                        ),
                        attr_name,
                        derived_value,
                        info["value"],
                    )

    def _derive_attribute(self, attr_name, info):
        """
        Function to Derive Attributes based on other attribues in the CDF
        writer's internal dict state.
        """
        # Cache for easier access
        method = info["derived"]["method"]

        # SWITCH on the Derivation Method
        if method == "get_generation_date":
            return datetime.datetime.now().strftime("%Y%m%d")
        elif method == "get_data_type":
            return self.get_data_type()
        elif method == "get_logical_file_id":
            return self.get_logical_file_id()
        elif method == "get_logical_source":
            return self.get_logical_source()
        elif method == "get_logical_source_description":
            return self.get_logical_source_description()
        else:
            raise ValueError(f"Derivation Method ({method}) Not Recognized")

    def to_cdf(self, output_path="./"):
        """
        Function to convert members of the CDFWriter's internal dict
        state to a CDF File using the pyCDF module from spacepy
        """
        assert self.target_cdf_file is None

        # Derive any Global Attributes
        self.derive_attributes()

        # Initialize a new CDF
        cdf_filename = f"{self.target_dict['gAttrList']['Logical_file_id']['value']}.cdf"
        logging.info(f"Generating CDF: {cdf_filename}")
        output_cdf_filepath = os.path.join(output_path, cdf_filename)
        self.target_cdf_file = pycdf.CDF(output_cdf_filepath, masterpath="")

        # Add Global Attriubtes to the CDF File
        if "gAttrList" in self.target_dict:
            self._convert_global_attributes_to_cdf()

        # Add zAttributes
        if "zAttrList" in self.target_dict:
            self._convert_variable_attributes_to_cdf()

        return output_cdf_filepath

    def _convert_global_attributes_to_cdf(self):
        # Loop though Global Attributes in target_dict
        for attr_name, info in self.target_dict["gAttrList"].items():
            # Make sure the Value is not None
            # We cannot add None Values to the CDF Global Attrs
            if not info["value"]:
                logging.warning(
                    (
                        f"Cannot Add gAttr: {attr_name:30}. Value was {str(info['value']):5}. "
                        f"Required: {str(info['required']):5}. Derived: {str(info['derived']):5}"
                    )
                )
                # Set to a intentionally invalid value
                info["value"] = f"{attr_name}-Not-Provided"
            # Add the Attribute to the CDF File
            self.target_cdf_file.attrs[attr_name] = info["value"]

    def _convert_variable_attributes_to_cdf(self):
        # Loop through Variable Attributes in target_dict
        for attr_name, info in self.target_dict["zAttrList"].items():
            # Make sure the Value is not None
            # We cannot add None Values to the CDF Vars
            if not info["value"]:
                logging.warning(
                    (
                        f"Cannot Add zAttr: {attr_name:30}. Value was {str(info['value']):5}. "
                        f"Required: {str(info['required']):5}. Derived: {str(info['derived']):5}"
                    )
                )
                # Set to a intentionally invalid value
                info["value"] = f"{attr_name}-Not-Provided"
            # Add the Attribute to the CDF File
            self.target_cdf_file[attr_name] = [info["value"]]

    def save_cdf(self):
        """Function to save and close CDF File"""
        assert self.target_cdf_file is not None
        # Save the CDF
        self.target_cdf_file.save()
        self.target_cdf_file.close()

    # ==============================================================================================
    #                               ATTRIBUTE DERIVATION FUNCTIONS
    # ==============================================================================================

    def get_logical_file_id(self):
        """
        Function to get the `Logical_file_id` required global attribute.

        The attribute stores the name of the CDF File without the file
        extension (e.g. '.cdf'). This attribute is requires to avoid
        loss of the originial source in case of renaming.
        """
        if not self.target_dict["gAttrList"]["Logical_file_id"]["value"]:
            # Get Parts
            instrumentId = self.get_instrument_id()
            startTime = self.get_start_time()
            dataLevel = self.get_data_level()
            version = self.get_version()
            mode = self.get_mode()

            # Build Derivation
            science_filename = util.create_science_filename(
                instrument=instrumentId, time=startTime, level=dataLevel, version=version, mode=mode
            )
            science_filename = science_filename.rstrip(util.FILENAME_EXTENSION)
        else:
            science_filename = self.target_dict["gAttrList"]["Logical_file_id"]["value"]
        return science_filename

    def get_logical_source(self):
        """
        Function to get the `Logical_source` required global attribute.

        This attribute determines the file naming convention in the SKT Editor
        and is used by CDA Web.
        """
        if not self.target_dict["gAttrList"]["Logical_source"]["value"]:
            # Get Parts
            spacecraft_id = self.get_spacecraft_id()
            instrument_id = self.get_instrument_id()
            data_type = self.get_data_type()
            data_type_short_name, _ = data_type.split(">")
            # Build Derivation
            logical_source = f"{spacecraft_id}_{instrument_id}_{data_type_short_name}"
        else:
            logical_source = self.target_dict["gAttrList"]["Logical_source"]["value"]
        return logical_source

    def get_logical_source_description(self):
        """
        Function to get the `Logical_source_description` required global attribute.

        This attribute writes out the full words associated with the encryped
        `Logical_source`  attribute.
        """
        if not self.target_dict["gAttrList"]["Logical_source_description"]["value"]:
            # Get Parts
            spacecraft_long_name = self.get_spacecraft_long_name()
            instrument_long_name = self.get_instrument_long_name()
            data_level_long_name = self.get_data_level_long_name()
            logical_source_description = (
                f"{data_level_long_name} {spacecraft_long_name} {instrument_long_name}"
            )
        else:
            logical_source_description = self.target_dict["gAttrList"][
                "Logical_source_description"
            ]["value"]
        return logical_source_description

    def get_data_type(self):
        """
        Function to get the `Data_type` required global attribute.

        This attribute is used by the CDF Writing software to create the filename.
        It is a combination of the following components:
            - mode
            - data_level
            - optional_data_product_descriptor
        """
        if not self.target_dict["gAttrList"]["Data_type"]["value"]:
            # Get Short Parts
            mode = self.get_mode()
            data_level = self.get_data_level()
            odpd = self.get_optional_data_product_descriptor()

            # Get Long Parts
            mode_long_name = self.get_mode().upper()  # NOTE Seems to be Upper case?
            data_level_long_name = self.get_data_level_long_name()
            odpd_long_name = (
                self.get_optional_data_product_descriptor().upper()
            )  # NOTE Seems to be Uppar case?

            # Build Derivation
            data_type = (
                f"{mode}_{data_level}_{odpd}>"
                f"{mode_long_name} {data_level_long_name} {odpd_long_name}"
            )
        else:
            data_type = self.target_dict["gAttrList"]["Data_type"]["value"]
        return data_type

    def get_spacecraft_id(self):
        """Function to get Spacecraft ID from Source_name Global Attribute"""
        sc_id = self.target_dict["gAttrList"]["Source_name"]["value"]
        assert sc_id is not None
        # Formatting
        if ">" in sc_id:
            short_name, _ = sc_id.split(">")
            sc_id = short_name.lower()  # Makse sure its all lowercase
        return sc_id

    def get_spacecraft_long_name(self):
        """Function to get Spacecraft ID from Source_name Global Attribute"""
        sc_id = self.target_dict["gAttrList"]["Source_name"]["value"]
        assert sc_id is not None
        # Formatting
        if ">" in sc_id:
            _, long_name = sc_id.split(">")
            sc_id = long_name
        return sc_id

    def get_instrument_id(self):
        """
        Function to get Instrument ID from Descriptor Global Attribute

        Instrument of investigation identifier shortened to three
        letter acronym.
        """
        instr_id = self.target_dict["gAttrList"]["Descriptor"]["value"]
        assert instr_id is not None
        # Formatting
        if ">" in instr_id:
            short_name, _ = instr_id.split(">")
            instr_id = short_name.lower()  # Makse sure its all lowercase
        return instr_id

    def get_instrument_long_name(self):
        """
        Function to get Instrument ID from Descriptor Global Attribute

        Instrument of investigation identifier shortened to three
        letter acronym.
        """
        instr_id = self.target_dict["gAttrList"]["Descriptor"]["value"]
        assert instr_id is not None
        # Formatting
        if ">" in instr_id:
            _, long_name = instr_id.split(">")
            instr_id = long_name
        return instr_id

    def get_data_level(self):
        """
        Function to get Data Level of CDF data

        The level to which the data product has been processed.
        """
        # TODO
        data_level = "l1>Level 1"
        # Formatting
        if ">" in data_level:
            short_name, _ = data_level.split(">")
            data_level = short_name.lower()  # Makse sure its all lowercase
        return data_level

    def get_data_level_long_name(self):
        """
        Function to get Data Level of CDF data

        The level to which the data product has been processed.
        """
        # TODO
        data_level = "l1>Level 1"
        # Formatting
        if ">" in data_level:
            _, long_name = data_level.split(">")
            data_level = long_name
        return data_level

    def get_optional_data_product_descriptor(self):
        """
        Function to get the (Optional) Data Product Descriptor.

        This is an optional field that may not be needed for all products. Where it is used,
        identifier shouls be short (3-8 charachters) descriptors that are helpful to end users.
        If a descriptor contains multiple components, underscores are used top separate
        hose components.
        """
        # TODO
        return "odpd"

    def get_start_time(self):
        """
        Function to get the start time of the data contained in the CDF
        given in format `YYYYMMDD_hhmmss`
        """
        # TODO
        return datetime.datetime.now()

    def get_version(self):
        """
        Function to get the 3-part version number of the data product.
        """
        version_str = self.target_dict["gAttrList"]["Data_version"]["value"].lower()
        assert version_str is not None
        if "v" in version_str:
            _, version = version_str.split("v")
        else:
            version = version_str
        return version

    def get_mode(self):
        """Function to get the mode attribute (TBS)"""
        # TODO
        return "mode"
