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
    """

    def __init__(self, template_data_file="cdf_template.yaml", seed_data={}):

        # Load Default Template Data
        template_file = Path(hermes_core.__file__).parent / "data" / template_data_file
        default_template_data = {}
        with open(template_file, "r") as f:
            try:
                default_template_data = yaml.safe_load(f)
            except yaml.YAMLError as exc:
                logging.critical(exc)

        # Current Working CDF File
        self.target_cdf_file = None

        # Temp Conversion to CDF
        self.target_dict = {}
        # Add TEMPLATES to Target_Dict
        self.add_data_from_dict(default_template_data)
        self.add_data_from_dict(seed_data)

    def add_data_from_dict(self, data={}):
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

    def _derive_attribute(self, attr_name, info):
        """
        Function to Derive Attributes based on other attribues in the CDF
        writer's internal dict state.
        """
        # Cache for easier access
        method = info["derived"]["method"]

        # SWITCH on the Derivation Method
        if method == "get_current_time":
            return datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
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
        for attr_name, info in self.target_dict["gAttrList"].items():
            if info["derived"]:
                info["value"] = self._derive_attribute(attr_name=attr_name, info=info)

        # Initialize a new CDF
        cdf_filename = f"{self.target_dict['gAttrList']['Logical_file_id']['value']}.cdf"
        logging.info(f"Generating CDF: {cdf_filename}")
        output_cdf_filepath = os.path.join(output_path, cdf_filename)
        self.target_cdf_file = pycdf.CDF(output_cdf_filepath, masterpath="")

        # Add Global Attriubtes to the CDF File
        if "gAttrList" in self.target_dict:
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
                    info["value"] = "TEST"
                # Add the Attribute to the CDF File
                self.target_cdf_file.attrs[attr_name] = info["value"]

        # Add zAttributes
        if "zAttrList" in self.target_dict:
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
                    info["value"] = "TEST"
                # Add the Attribute to the CDF File
                self.target_cdf_file[attr_name] = [info["value"]]

        return output_cdf_filepath

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
            instrumentId = self.get_instrument_id()
            data_type = self.get_data_type()
            logical_source = f"{spacecraft_id}_{instrumentId}_{data_type}"
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
                f"{spacecraft_long_name}_{instrument_long_name}_{data_level_long_name}"
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
            # Get Parts
            mode = self.get_mode()
            dataLevel = self.get_data_level()
            odpd = self.get_optional_data_product_descriptor()
            # Build Derivation
            data_type = f"{mode}_{dataLevel}_{odpd}"
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
        _, version = version_str.split("v")
        return version

    def get_mode(self):
        """Function to get the mode attribute (TBS)"""
        # TODO
        return "mode"
