from abc import ABC, abstractmethod
from pathlib import Path
import yaml
import hermes_core
from hermes_core import log

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
        self._global_attr_schema = self._load_default_global_attr_schema()

        # Data Validation and Compliance for Variable Data
        self._variable_attr_schema = self._load_default_variable_attr_schema()

    @property
    def global_attribute_schema(self):
        return self._global_attr_schema

    @property
    def variable_attribute_schema(self):
        return self._variable_attr_schema

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


class NetCDFSchema(FileTypeSchema):
    """Schema for NetCDF files."""

    @property
    def global_attribute_schema(self):
        return {
            "attribute1": {"type": "string", "required": True},
            "attribute2": {"type": "int", "required": False},
            # Define more global attributes and their schemas
        }

    @property
    def variable_attribute_schema(self):
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


class FITSSchema(FileTypeSchema):
    """Schema for FITS files."""

    @property
    def global_attribute_schema(self):
        return {
            "attribute1": {"type": "string", "required": True},
            "attribute2": {"type": "int", "required": False},
            # Define more global attributes and their schemas
        }

    @property
    def variable_attribute_schema(self):
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
