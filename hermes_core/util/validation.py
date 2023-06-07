from pathlib import Path
from abc import ABC, abstractmethod
from spacepy.pycdf import CDF
from spacepy.pycdf.istp import FileChecks, VariableChecks
from hermes_core.util.schema import CDFSchema


def validate(filepath):
    """
    Validate a data file such as a CDF.

    Parameters
    ----------
    filepath : `str`
        A fully specificed file path.

    Returns
    -------
    errors : `list[str]`
        A list of validation errors returned. A valid file will result in an emppty list being returned.
    """
    # Determine the file type
    file_extension = Path(filepath).suffix

    # Create the appropriate validator object based on file type
    if file_extension == ".cdf":
        validator = CDFValidator()
    elif file_extension == ".nc":
        validator = NetCDFValidator()
    else:
        raise ValueError(f"Unsupported file type: {file_extension}")

    # Call the validate method of the validator object
    return validator.validate(filepath)


class TimeDataValidator(ABC):
    """
    Abstract base class for heliophysics data validators.
    """

    @abstractmethod
    def validate(self, file_path):
        """
        Validate the heliophysics data file.

        Parameters
        ----------
        file_path : `str`
            The path to the data file.

        Returns
        -------
        errors : `list[str]`
            A list of validation errors returned. A valid file will result in an emppty list being returned.
        """
        pass


class CDFValidator(TimeDataValidator):
    """
    Validator for CDF files.
    """

    def __init__(self):
        super().__init__()

        # CDF Schema
        self.schema = CDFSchema()

    def validate(self, file_path):
        """
        Validate the CDF file.

        Parameters
        ----------
        file_path : `str`
            The path to the CDF file.

        Returns
        -------
        errors : `list[str]`
            A list of validation errors returned. A valid file will result in an emppty list being returned.
        """
        # Initialize Validation Errrors
        validation_errors = []

        try:
            # Open CDF file with context manager
            with CDF(file_path) as cdf_file:
                # Verify that all `required` global attributes in the schema are present
                global_attr_validation_errors = self._validate_global_attr_schema(cdf_file=cdf_file)
                validation_errors.extend(global_attr_validation_errors)

                # Verify that all `required` variable attributes in the schema are present
                variable_attr_validation_errors = self._validate_variable_attr_schema(
                    cdf_file=cdf_file
                )
                validation_errors.extend(variable_attr_validation_errors)

                # Validate the CDF Using ISTP Module `FileChecks` Class
                file_checks_errors = self._file_checks(cdf_file=cdf_file)
                validation_errors.extend(file_checks_errors)

        except IOError:
            validation_errors.append(f"Could not open CDF File at path: {file_path}")

        return validation_errors

    def _validate_global_attr_schema(self, cdf_file: CDF):
        """
        Function to ensure all required global attributes in the schema are present
        in the generated CDF File.
        """
        global_attr_validation_errors = []
        # Loop for each attribute in the schema
        for attr_name, attr_schema in self.schema.global_attribute_schema.items():
            # If it is a required attribute and not present
            if attr_schema["validate"] and (attr_name not in cdf_file.attrs):
                global_attr_validation_errors.append(
                    f"Required attribute ({attr_name}) not present in global attributes.",
                )
        return global_attr_validation_errors

    def _validate_variable_attr_schema(self, cdf_file: CDF):
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
        var_type_attrs = self.schema.variable_attribute_schema[var_type]

        # Get the `Var()` Class for the Variable
        var_data = cdf_file[var_name]

        # Loop for each Variable Attribute in the schema
        for attr_name in var_type_attrs:
            attr_schema = self.schema.variable_attribute_schema["attribute_key"][attr_name]
            # If it is a required attribute and not present
            if attr_schema["required"] and attr_name not in var_data.attrs:
                # Check to see if there is an "alternate" attribute
                if "alternate" not in attr_schema:
                    variable_errors.append(f"Variable: {var_name} missing '{attr_name}' attribute.")
                # If there is an alternate, and the alternate is not in the metadata
                if "alternate" in attr_schema and attr_schema["alternate"] not in var_data.attrs:
                    variable_errors.append(
                        f"Variable: {var_name} missing '{attr_name}' attribute. Alternative: {attr_schema['alternate']} not found."
                    )
            # Assume that the Attribue is Present in the metadata for the Variable
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
        # Validate `FORMAT` Attribute on the Variable
        if "FORMAT" in var_type_attrs and "FORMAT" in var_data.meta:
            format_errors = self._validate_format(cdf_file=cdf_file, var_name=var_name)
            if format_errors:
                variable_errors.append(format_errors)

        # Validate Variable using ISTP Module `VariableChecks` class
        variable_checks_errors = self._variable_checks(cdf_file=cdf_file, var_name=var_name)
        variable_errors.extend(variable_checks_errors)

        return variable_errors

    def _validate_format(self, cdf_file: CDF, var_name: str):
        # Save the Current Format
        variable_format = cdf_file[var_name].meta["FORMAT"]
        # Get the target Format for the Variable
        target_format = self.schema._format_helper(
            data=cdf_file, var_name=var_name, cdftype=cdf_file[var_name].type()
        )

        if variable_format != target_format:
            return f"Variable: {var_name} Attribute 'FORMAT' value '{variable_format}' does not match derrived format '{target_format}'"
        else:
            return None

    def _file_checks(self, cdf_file: CDF):
        """
        Function to call individual pieces of the `spacepy.pycdf.istp.FileChecks` Class.
        We do not want to run all validation checks from this class using the `all()` function
        so we break up the individual function calls here.
        """
        file_checks_errors = []

        check_fns = [
            FileChecks.empty_entry,
            FileChecks.filename,
            FileChecks.time_monoton,
            FileChecks.times,
        ]

        # Loop through the Functions we want to check
        for func in check_fns:
            # Try to call the given function and report errors
            try:
                file_checks_errors.extend(func(cdf_file))
            # If the function errors out or does not complete, report this an an error itself.
            except:
                file_checks_errors.append("Test {} did not complete.".format(func.__name__))

        return file_checks_errors

    def _variable_checks(self, cdf_file: CDF, var_name: str):
        """
        Function to call individual pieces of the `spacepy.pycdf.istp.VariableChecks` Class.
        We do not want to run all validation checks from this class using the `all()` function
        so we break up the individual function calls here.
        """
        variable_checks_errors = []

        check_fns = [
            # This function makes incorrect asumptions about the UNITS that must be placed on
            # DELTA_PLUS_VAR and DELTA_MINUS var metadata attributes.
            # VariableChecks.deltas,
            VariableChecks.depends,
            VariableChecks.depsize,
            VariableChecks.empty_entry,
            # This function makes incorrect assumptions that the variable name must exactly
            # match the FILEDNAM metadata attribute.
            # VariableChecks.fieldnam,
            VariableChecks.fillval,
            VariableChecks.recordcount,
            # This function makes incorrect assumptions about the valid DISPLAY_TYPE options
            # based on the shape of the variable data.
            # VariableChecks.validdisplaytype,
            VariableChecks.validrange,
            VariableChecks.validscale,
        ]

        # Loop through the Functions we want to check
        for func in check_fns:
            # Try to call the given function and report errors
            try:
                variable_checks_errors.extend(
                    ("{}: {}".format(var_name, e) for e in func(cdf_file[var_name]))
                )
            # If the function errors out or does not complete, report this an an error itself.
            except:
                variable_checks_errors.append(
                    "{}: Test {} did not complete.".format(var_name, func.__name__)
                )

        # for v in f:
        #     errors.extend(('{}: {}'.format(v, e)
        #                    for e in VariableChecks.all(f[v], catch=catch)))

        return variable_checks_errors


class NetCDFValidator(TimeDataValidator):
    """
    Validator for NetCDF files.
    """

    def validate(self, file_path):
        """
        Validate the NetCDF file.

        Parameters
        ----------
        file_path : `str`
            The path to the NetCDF file.

        Returns
        -------
        errors : `list[str]`
            A list of validation errors returned. A valid file will result in an emppty list being returned.
        """
        # Validation logic for NetCDF files
        pass
