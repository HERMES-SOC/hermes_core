import pytest
from collections import OrderedDict
from astropy.table import Table
from hermes_core.util.schema import HermesDataSchema

def test_hermes_data_schema():
    """Test Schema Template and Info Functions"""
    schema = HermesDataSchema()

    # Global Attribute Schema
    assert schema.global_attribute_schema is not None
    assert isinstance(schema.global_attribute_schema, dict)

    # Variable Attribute Schema
    assert schema.variable_attribute_schema is not None
    assert isinstance(schema.variable_attribute_schema, dict)

    # Default Global Attributes
    assert schema.default_global_attributes is not None
    assert isinstance(schema.default_global_attributes, dict)

    # Global Attribute Template
    assert schema.global_attribute_template() is not None
    assert isinstance(schema.global_attribute_template(), OrderedDict)

    # Measurement Attribute Template
    assert schema.measurement_attribute_template() is not None
    assert isinstance(schema.measurement_attribute_template(), OrderedDict)

    # Global Attribute Info
    assert schema.global_attribute_info() is not None
    assert isinstance(schema.global_attribute_info(), Table)
    assert isinstance(schema.global_attribute_info(attribute_name="Descriptor"), Table)
    with pytest.raises(KeyError):
        _ = schema.global_attribute_info(attribute_name="NotAnAttribute")

    # Measurement Attribute Info
    assert schema.measurement_attribute_info() is not None
    assert isinstance(schema.measurement_attribute_info(), Table)
    assert isinstance(
        schema.measurement_attribute_info(attribute_name="CATDESC"), Table
    )
    with pytest.raises(KeyError):
        _ = schema.measurement_attribute_info(attribute_name="NotAnAttribute")
