from abc import ABC, abstractmethod
from pathlib import Path
import yaml
import spacepy
from spacepy.pycdf import CDF, _Hyperslice
from spacepy.pycdf.istp import FileChecks
import hermes_core
from hermes_core import log
from hermes_core.util import util
from hermes_core.util.exceptions import warn_user


class ScienceDataValidator(ABC):
    """
    Abstract base class for heliophysics data validators.
    """

    @abstractmethod
    def validate(self, file_path):
        """
        Validate the heliophysics data file.

        Parameters:
            file_path (str): The path to the data file.

        Returns:
            bool: True if the file is valid, False otherwise.
        """
        pass


class CDFValidator(ScienceDataValidator):
    """
    Validator for CDF files.
    """

    DEFAULT_GLOBAL_CDF_ATTRS_SCHEMA_FILE = "hermes_default_global_cdf_attrs_schema.yaml"
    DEFAULT_GLOBAL_CDF_ATTRS_FILE = "hermes_default_global_cdf_attrs.yaml"
    DEFAULT_VARIABLE_CDF_ATTRS_SCHEMA_FILE = "hermes_default_variable_cdf_attrs_schema.yaml"

    def __init__(self):
        super().__init__()

        # Data Validation, Complaiance, Derived Attributes
        self._global_attr_schema = self._load_default_global_attr_schema()

        # Data Validation and Compliance for Variable Data
        self._variable_attr_schema = self._load_default_variable_attr_schema()

    def _load_default_global_attr_schema(self) -> dict:
        # The Default Schema file is contained in the `hermes_core/data` directory
        default_schema_path = str(
            Path(hermes_core.__file__).parent
            / "data"
            / CDFValidator.DEFAULT_GLOBAL_CDF_ATTRS_SCHEMA_FILE
        )
        # Load the Schema
        return self._load_yaml_data(yaml_file_path=default_schema_path)

    def _load_default_variable_attr_schema(self) -> dict:
        # The Default Schema file is contained in the `hermes_core/data` directory
        default_schema_path = str(
            Path(hermes_core.__file__).parent
            / "data"
            / CDFValidator.DEFAULT_VARIABLE_CDF_ATTRS_SCHEMA_FILE
        )
        # Load the Schema
        return self._load_yaml_data(yaml_file_path=default_schema_path)

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

    def validate(self, file_path):
        """
        Validate the CDF file.

        Parameters:
            file_path (str): The path to the CDF file.

        Returns:
            bool: True if the file is a valid CDF, False otherwise.
        """
        # Initialize Validation Errrors
        validation_errors = []

        try:
            # Open CDF file with context manager
            with CDF(file_path) as cdf_file:
                # Verify that all `required` global attributes in the schema are present
                global_attr_validation_errors = self.validate_global_attr_schema(cdf_file=cdf_file)
                validation_errors.extend(global_attr_validation_errors)

                # Verify that all `required` variable attributes in the schema are present
                variable_attr_validation_errors = self.validate_variable_attr_schema(
                    cdf_file=cdf_file
                )
                validation_errors.extend(variable_attr_validation_errors)

                # Validate the CDF Using ISTP Module
                istp_validation_errors = FileChecks.all(f=cdf_file, catch=True)
                validation_errors.extend(istp_validation_errors)

        except IOError:
            validation_errors.append(f"Could not open CDF File at path: {file_path}")

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


class NetCDFValidator(ScienceDataValidator):
    """
    Validator for NetCDF files.
    """

    def validate(self, file_path):
        """
        Validate the NetCDF file.

        Parameters:
            file_path (str): The path to the NetCDF file.

        Returns:
            bool: True if the file is a valid NetCDF, False otherwise.
        """
        # Validation logic for NetCDF files
        pass


class FITSValidator(ScienceDataValidator):
    """
    Validator for FITS files.
    """

    def validate(self, file_path):
        """
        Validate the FITS file.

        Parameters:
            file_path (str): The path to the FITS file.

        Returns:
            bool: True if the file is a valid FITS, False otherwise.
        """
        # Validation logic for FITS files
        pass
