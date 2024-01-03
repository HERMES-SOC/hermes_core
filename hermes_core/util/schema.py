"""
This module provides schema metadata derivations.
"""

from swxsoc_core.util.schema import SpaceWeatherDataSchema

__all__ = ["HermesDataSchema"]

DEFAULT_GLOBAL_CDF_ATTRS_SCHEMA_FILE = "hermes_default_global_cdf_attrs_schema.yaml"
DEFAULT_VARIABLE_CDF_ATTRS_SCHEMA_FILE = "hermes_default_variable_cdf_attrs_schema.yaml"


class HermesDataSchema(SpaceWeatherDataSchema):
    """Class representing the schema of a file type."""

    def __init__(self):
        super().__init__()
