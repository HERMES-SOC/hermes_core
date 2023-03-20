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
                print(exc)

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

    def get_spacecraft_id(self):
        """Function to get Spacecraft ID from Source_name Global Attribute"""
        sc_id = self.target_dict["gAttrList"]["Source_name"]["value"]
        assert sc_id is not None
        # Wird Formatting
        if ">" in sc_id:
            short_name, _ = sc_id.split(">")
            sc_id = short_name.lower()  # Makse sure its all lowercase
        return sc_id

    def get_instrument_id(self):
        """
        Function to get Instrument ID from Descriptor Global Attribute

        Instrument of investigation identifier shortened to three
        letter acronym.
        """
        instr_id = self.target_dict["gAttrList"]["Descriptor"]["value"]
        assert instr_id is not None
        # Wird Formatting
        if ">" in instr_id:
            short_name, _ = instr_id.split(">")
            instr_id = short_name.lower()  # Makse sure its all lowercase
        return instr_id

    def get_data_level(self):
        """
        Function to get Data Level of CDF data

        The level to which the data product has been processed.
        """
        # TODO
        return "l0"

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
        return datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    def get_version(self):
        """
        Function to get the 3-part version number of the data product.
        """
        # TODO
        return self.target_dict["gAttrList"]["Data_version"]["value"].lower()

    def get_mode(self):
        """Function to get the mode attribute (TBS)"""
        # TODO
        return "mode"

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

        # SWITCH on attr_name
        if attr_name == "Logical_file_id":
            # Get Parts
            scID = self.get_spacecraft_id()
            instrumentId = self.get_instrument_id()
            dataLevel = self.get_data_level()
            startTime = self.get_start_time()
            version = self.get_version()
            # Build Derivation
            return f"{scID}_{instrumentId}_{dataLevel}_{startTime}_{version}"

        if attr_name == "Data_type":
            # Get Parts
            mode = self.get_mode()
            dataLevel = self.get_data_level()
            odpd = self.get_optional_data_product_descriptor()
            # Build Derivation
            return f"{mode}_{dataLevel}_{odpd}"

        return None

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
