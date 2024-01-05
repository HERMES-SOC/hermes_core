"""
This module provides schema metadata derivations.
"""
from pathlib import Path
from typing import Optional
from swxsoc_core.util.schema import SpaceWeatherDataSchema
import hermes_core

__all__ = ["HermesDataSchema"]

DEFAULT_GLOBAL_CDF_ATTRS_SCHEMA_FILE = "hermes_default_global_cdf_attrs_schema.yaml"
DEFAULT_VARIABLE_CDF_ATTRS_SCHEMA_FILE = "hermes_default_variable_cdf_attrs_schema.yaml"


class HermesDataSchema(SpaceWeatherDataSchema):
    """
    Class representing a schema for data requirements and formatting, specific to the
    HERMES Mission.

    There are two main componentes to the HERMES Data Schema, including both global and
    variable attribute information.

    Parameters
    ----------
    global_schema_layers :  `Optional[list[Path]]`
        Absolute file paths to global attribute schema files. These schema files are layered
        on top of one another in a latest-priority ordering. That is, the latest file that modifies
        a common schema attribute will take precedence over earlier values for a given attribute.
    variable_schema_layers :  `Optional[list[Path]]`
        Absolute file paths to variable attribute schema files. These schema files are layered
        on top of one another in a latest-priority ordering. That is, the latest file that modifies
        a common schema attribute will take precedence over earlier values for a given attribute.
    use_defaults: `Optional[bool]`
        Whether or not to load the default global and variable attribute schema files. These
        default schema files contain only the requirements for CDF ISTP validation.
    """

    def __init__(
        self,
        global_schema_layers: Optional[list[Path]] = None,
        variable_schema_layers: Optional[list[Path]] = None,
        use_defaults: Optional[bool] = True,
    ):
        # HERMES Default Global Schema
        global_schema_path = str(
            Path(hermes_core.__file__).parent
            / "data"
            / DEFAULT_GLOBAL_CDF_ATTRS_SCHEMA_FILE
        )
        # HERMES Default Variable Schema
        variable_schema_path = str(
            Path(hermes_core.__file__).parent
            / "data"
            / DEFAULT_VARIABLE_CDF_ATTRS_SCHEMA_FILE
        )

        # Seed Layers with Default
        if not use_defaults:
            _global_schema_layers = []
            _variable_schema_layers = []
        else:
            _global_schema_layers = [global_schema_path]
            _variable_schema_layers = [variable_schema_path]

        # Extend Custom Layers
        if global_schema_layers is not None and len(global_schema_layers) > 0:
            _global_schema_layers.extend(global_schema_layers)
        if variable_schema_layers is not None and len(variable_schema_layers) > 0:
            _variable_schema_layers.extend(variable_schema_layers)

        # Call SWxSOC Initialization to populate Schema
        super().__init__(
            global_schema_layers=_global_schema_layers,
            variable_schema_layers=_variable_schema_layers,
            use_defaults=use_defaults,
        )
