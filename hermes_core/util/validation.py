from pathlib import Path
from abc import ABC, abstractmethod
import numpy as np
from spacepy.pycdf import CDF, CDFError
from spacepy.pycdf.istp import FileChecks, VariableChecks
from hermes_core.util.schema import HermesDataSchema

__all__ = ["validate", "CDFValidator"]


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
    else:
        raise ValueError(f"Unsupported file type: {file_extension}")

    # Call the validate method of the validator object
    return validator.validate(filepath)


class HermesDataValidator(ABC):
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


class CDFValidator(HermesDataValidator):
    """
    Validator for CDF files.
    """

    def __init__(self):
        super().__init__()

        # CDF Schema
        self.schema = HermesDataSchema()

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
            with CDF(file_path, readonly=True) as cdf_file:
                # Verify that all `required` global attributes in the schema are present
                global_attr_validation_errors = self._validate_global_attr_schema(
                    cdf_file=cdf_file
                )
                validation_errors.extend(global_attr_validation_errors)

                # Verify that all `required` variable attributes in the schema are present
                variable_attr_validation_errors = self._validate_variable_attr_schema(
                    cdf_file=cdf_file
                )
                validation_errors.extend(variable_attr_validation_errors)

                # Validate the CDF Using ISTP Module `FileChecks` Class
                file_checks_errors = self._file_checks(cdf_file=cdf_file)
                validation_errors.extend(file_checks_errors)

        except CDFError:
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
            # If it is a required attribute but null
            if (
                attr_schema["validate"]
                and (attr_name in cdf_file.attrs)
                and (
                    (cdf_file.attrs[attr_name][0] == "")
                    or (cdf_file.attrs[attr_name][0] is None)
                )
            ):
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
            attr_schema = self.schema.variable_attribute_schema["attribute_key"][
                attr_name
            ]
            # If it is a required attribute and not present
            if attr_schema["required"] and attr_name not in var_data.attrs:
                # Check to see if there is an "alternate" attribute
                if attr_schema["alternate"] is None:
                    variable_errors.append(
                        f"Variable: {var_name} missing '{attr_name}' attribute."
                    )
                # If there is an alternate, and the alternate is not in the metadata
                if (
                    "alternate" in attr_schema
                    and attr_schema["alternate"] is not None
                    and attr_schema["alternate"] not in var_data.attrs
                ):
                    variable_errors.append(
                        f"Variable: {var_name} missing '{attr_name}' attribute. Alternative: {attr_schema['alternate']} not found."
                    )
            # Assume that the Attribue is Present in the metadata for the Variable
            else:
                # If the Var Data can be Validated
                if (
                    "valid_values" in attr_schema
                    and attr_schema["valid_values"] is not None
                ):
                    attr_valid_values = attr_schema["valid_values"]
                    attr_value = var_data.attrs[attr_name]
                    if attr_value not in attr_valid_values:
                        variable_errors.append(
                            (
                                f"Variable: {var_name} Attribute '{attr_name}' not one of valid options.",
                                f"Was {attr_value}, expected one of {attr_valid_values}",
                            )
                        )

        # Validate Variable using ISTP Module `VariableChecks` class
        variable_checks_errors = self._variable_checks(
            cdf_file=cdf_file, var_name=var_name
        )
        variable_errors.extend(variable_checks_errors)

        return variable_errors

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
            except:  # noqa: E722
                file_checks_errors.append(
                    "Test {} did not complete.".format(func.__name__)
                )

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
            # This function makes incorrect assumtions that the FILLVAL must be derived from
            # the CDF data type of the variable. A FILLVAL should be allowed to be set as needed by
            # instrument team developers.
            # VariableChecks.fillval,
            VariableChecks.recordcount,
            # This function makes incorrect assumptions about the valid DISPLAY_TYPE options
            # based on the shape of the variable data.
            # VariableChecks.validdisplaytype,
            # This Function makes inforrect assumptions that the VLIDMIN and VLIDIMAX must be
            # derived from the CDF data type of the variable. A VALIDMIN and VALIDMAX should be
            # allowed to be set as needed by instrument team developers.
            self._validrange,
            self._validscale,
        ]

        # Loop through the Functions we want to check
        for func in check_fns:
            # Try to call the given function and report errors
            try:
                variable_checks_errors.extend(
                    ("{}: {}".format(var_name, e) for e in func(cdf_file[var_name]))
                )
            # If the function errors out or does not complete, report this an an error itself.
            except:  # noqa: E722
                variable_checks_errors.append(
                    "{}: Test {} did not complete.".format(var_name, func.__name__)
                )

        return variable_checks_errors

    def _validrange(self, v):
        """Check that all values are within VALIDMIN/VALIDMAX, or FILLVAL

        Compare all values of this variable to `VALIDMIN
        <https://spdf.gsfc.nasa.gov/istp_guide/vattributes.html#VALIDMIN>`_
        and ``VALIDMAX``; fails validation if any values are below
        VALIDMIN or above ``VALIDMAX`` unless equal to `FILLVAL
        <https://spdf.gsfc.nasa.gov/istp_guide/vattributes.html#FILLVAL>`_.

        Parameters
        ----------
        v : :class:`~spacepy.pycdf.Var`
            Variable to check

        Returns
        -------
        list of str
            Description of each validation failure.

        """
        return self._validhelper(v)

    def _validscale(self, v):
        """Check SCALEMIN<=SCALEMAX, and both in range for CDF datatype.

        Compares `SCALEMIN
        <https://spdf.gsfc.nasa.gov/istp_guide/vattributes.html#SCALEMIN>`_
        to ``SCALEMAX`` to make sure it isn't larger and both are
        within range of the variable CDF datatype.

        Parameters
        ----------
        v : :class:`~spacepy.pycdf.Var`
            Variable to check

        Returns
        -------
        list of str
            Description of each validation failure.

        """
        return self._validhelper(v, False)

    def _validhelper(self, v, rng=True):
        """Helper function for checking SCALEMIN/MAX, VALIDMIN/MAX

        Parameters
        ----------
        v : :class:`~spacepy.pycdf.Var`
            Variable to check

        rng : bool
            Do range check (True, default) or scale check (False)

        Returns
        -------
        list of str
            Description of each validation failure.
        """
        validscale = "VALID" if rng else "SCALE"
        whichmin, whichmax = (
            ("VALIDMIN", "VALIDMAX") if rng else ("SCALEMIN", "SCALEMAX")
        )
        errs = []
        vshape = v.shape
        minval, maxval = self.schema._get_minmax(v.type())
        if rng:
            data = v[...]
            is_fill = False
            if "FILLVAL" in v.attrs:
                filldtype = self.schema.numpytypedict.get(
                    v.attrs.type("FILLVAL"), object
                )
                if np.issubdtype(v.dtype, np.floating) and np.issubdtype(
                    filldtype, np.floating
                ):
                    is_fill = np.isclose(data, v.attrs["FILLVAL"])
                elif np.can_cast(np.asanyarray(v.attrs["FILLVAL"]), v.dtype):
                    is_fill = data == v.attrs["FILLVAL"]
        for which in (whichmin, whichmax):
            if which not in v.attrs:
                continue
            attrval = v.attrs[which]
            multidim = bool(np.shape(attrval))  # multi-dimensional
            if multidim:  # Compare shapes, require only 1D var
                # Match attribute dim to first non-record var dim
                firstdim = int(v.rv())
                if vshape[firstdim] != np.shape(attrval)[0]:
                    errs.append(
                        (
                            "{} element count {} does not match first data"
                            " dimension size {}."
                        ).format(which, np.shape(attrval)[0], v.shape[firstdim])
                    )
                    continue
                if len(vshape) != firstdim + 1:  # only one non-record dim
                    errs.append(
                        "Multi-element {} only valid with 1D variable.".format(which)
                    )
                    continue
                if firstdim:  # Add pseudo-record dim
                    attrval = np.reshape(attrval, (1, -1))
            # min, max, variable data all same dtype
            if not np.can_cast(np.asanyarray(attrval), np.asanyarray(minval).dtype):
                errs.append(
                    "{} type {} not comparable to variable type {}.".format(
                        which,
                        self.schema.cdftypenames[v.attrs.type(which)],
                        self.schema.cdftypenames[v.type()],
                    )
                )
                continue  # Cannot do comparisons
            if np.any((minval > attrval)) or np.any((maxval < attrval)):
                errs.append(
                    "{} ({}) outside valid data range ({},{}).".format(
                        which, attrval[0, :] if multidim else attrval, minval, maxval
                    )
                )
            if not rng or not len(v):  # nothing to compare
                continue
            # Always put numpy array on the left so knows to do element compare
            idx = (data < attrval) if which == whichmin else (data > attrval)
            idx = np.logical_and(idx, np.logical_not(is_fill))
            if idx.any():
                direction = "under" if which == whichmin else "over"
                if len(vshape) == 0:  # Scalar
                    errs.append(
                        "Value {} {} {} {}.".format(
                            data,
                            direction,
                            which,
                            attrval[0, :] if multidim else attrval,
                        )
                    )
                    continue
                badidx = np.nonzero(idx)
                badvals = data[badidx]
                if len(badidx) > 1:  # Multi-dimensional data
                    badidx = np.transpose(badidx)  # Group by value not axis
                else:
                    badidx = badidx[0]  # Just recover the index value
                if len(badvals) < 10:
                    badvalstr = ", ".join(str(d) for d in badvals)
                    badidxstr = ", ".join(str(d) for d in badidx)
                    errs.append(
                        "Value {} at index {} {} {} {}.".format(
                            badvalstr,
                            badidxstr,
                            direction,
                            which,
                            attrval[0, :] if multidim else attrval,
                        )
                    )
                else:
                    errs.append(
                        "{} values {} {} {}".format(
                            len(badvals),
                            direction,
                            which,
                            attrval[0, :] if multidim else attrval,
                        )
                    )
        if (whichmin in v.attrs) and (whichmax in v.attrs):
            if np.any(v.attrs[whichmin] > v.attrs[whichmax]):
                errs.append("{} > {}.".format(whichmin, whichmax))
        return errs
