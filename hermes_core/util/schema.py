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

    There are two main components to the HERMES Data Schema, including both global and
    variable attribute information.

    Global schema information is loaded from YAML (dict-like) files in the following format:

    .. code-block:: yaml

        attribute_name:
            description: >
                Include a meaningful description of the attribute and context needed to understand
                its values.
            default: <string> # A default value for the attribute if needed/desired
            derived: <bool> # Whether or not the attribute's value can be derived using a python function
            derivation_fn: <string> # The name of a Python function to derive the value. Must be a function member of the schema class and match the signature below.
            required: <bool> # Whether the attribute is required
            validate: <bool> # Whether the attribute should be validated by the Validation module
            overwrite: <bool> # Whether an existing value for the attribute should be overwritten if a different value is derived.

    The signature for all functions to derive global attributes should follow the format below.
    The function takes in a parameter `data` which is a `HermesData` object, or that of an
    extended data class, and returns a single attribute value for the given attribute to be
    derived.

    .. code-block:: python

        def derivation_fn(self, data: HermesData):
            # ... do manipulations as needed from `data`
            return "attribute_value"

    Variable schema information is loaded from YAML (dict-like) files in the following format:

    .. code-block:: yaml

        attribute_key:
            attribute_name:
                description: >
                    Include a meaningful description of the attribute and context needed to understand
                    its values.
                derived: <bool> # Whether or not the attribute's value can be derived using a python function
                derivation_fn: <string> # The name of a Python function to derive the value. Must be a function member of the schema class and match the signature below.
                required: <bool> # Whether the attribute is required
                validate: <bool> # Whether the attribute should be validated by the Validation module
                overwrite: <bool> # Whether an existing value for the attribute should be overwritten if a different value is derived.
                valid_values: <list> # A list of valid values that the attribute can take. The value of the attribute is checked against the `valid_values` in the Validation module.
                alternate: <string> An additional attribute name that can be treated as an alternative of the given attribute.
        data:
            - attribute_name
            - ...
        support_data:
            - ...
        metadata:
            - ...

    The signature for all functions to derive variable attributes should follow the format below.
    The function takes in parameters `var_name`, `var_data`, and `guess_type`, where:

    - `var_name` is the variable name of the variable for which the attribute is being derived
    - `var_data` is the variable data of the variable for which the attribute is being derived
    - `guess_type` is the guessed CDF variable type of the data for which the attribute is being derived.

    The function must return a single attribute value for the given attribute to be derived.

    .. code-block:: python

        def derivation_fn(self, var_name: str, var_data: Union[Quantity, NDData, NDCube], guess_type: ctypes.c_long):
            # ... do manipulations as needed from data
            return "attribute_value"


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

    def _get_hermes_version(self, data):
        """Function to get the version of hermes_core used to generate the data"""
        attr_name = "hermes_version"
        if (attr_name not in data.meta) or (not data.meta[attr_name]):
            hermes_version = hermes_core.__version__
        else:
            hermes_version = data.meta[attr_name]
        return hermes_version
