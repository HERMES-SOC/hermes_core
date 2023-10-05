from abc import ABC, abstractmethod
from pathlib import Path
from collections import OrderedDict
from astropy.timeseries import TimeSeries
from astropy.time import Time
from astropy.nddata import NDData
from astropy.wcs import WCS
import astropy.units as u
from ndcube import NDCollection
from ndcube import NDCube
from hermes_core.util.exceptions import warn_user
from hermes_core.util.schema import HermesDataSchema

__all__ = ["CDFHandler"]

# ================================================================================================
#                                   ABSTRACT HANDLER
# ================================================================================================


class HermesDataIOHandler(ABC):
    """
    Abstract base class for handling input/output operations of heliophysics data.
    """

    @abstractmethod
    def load_data(self, file_path):
        """
        Load data from a file.

        Parameters
        ----------
        file_path : `str`
            A fully specified file path.

        Returns
        -------
        data : `~astropy.time.TimeSeries`
            An instance of `TimeSeries` containing the loaded data.
        """
        pass

    @abstractmethod
    def save_data(self, data, file_path):
        """
        Save data to a file.

        Parameters
        ----------
        data : `hermes_core.timedata.HermesData`
            An instance of `HermesData` containing the data to be saved.
        file_path : `str`
            The fully specified file path to save into.
        """
        pass


# ================================================================================================
#                                   CDF HANDLER
# ================================================================================================


class CDFHandler(HermesDataIOHandler):
    """
    A concrete implementation of HermesDataIOHandler for handling heliophysics data in CDF format.

    This class provides methods to load and save heliophysics data from/to a CDF file.
    """

    def __init__(self):
        super().__init__()

        # CDF Schema
        self.schema = HermesDataSchema()

    def load_data(self, file_path):
        """
        Load heliophysics data from a CDF file.

        Parameters
        ----------
        file_path : `str`
            The path to the CDF file.

        Returns
        -------
        data : `~astropy.time.TimeSeries`
            An instance of `TimeSeries` containing the loaded data.
        support : `dict`
            Non-record-varying data contained in the file
        """
        from spacepy.pycdf import CDF

        if not Path(file_path).exists():
            raise FileNotFoundError(f"CDF Could not be loaded from path: {file_path}")

        # Create a new TimeSeries
        ts = TimeSeries()
        # Create a Data Structure for Non-record Varying Data
        support = {}
        # Intermediate Type
        spectra = []

        # Open CDF file with context manager
        with CDF(file_path) as input_file:
            # Add Global Attributes from the CDF file to TimeSeries
            input_global_attrs = {}
            for attr_name in input_file.attrs:
                if len(input_file.attrs[attr_name]) == 0:
                    # gAttr is not set
                    input_global_attrs[attr_name] = ""
                elif len(input_file.attrs[attr_name]) > 1:
                    # gAttr is a List
                    input_global_attrs[attr_name] = input_file.attrs[attr_name][:]
                else:
                    # gAttr is a single value
                    input_global_attrs[attr_name] = input_file.attrs[attr_name][0]
            ts.meta.update(input_global_attrs)

            # First Variable we need to add is time/Epoch
            if "Epoch" in input_file:
                time_data = Time(input_file["Epoch"][:].copy())
                time_attrs = {}
                for attr_name in input_file["Epoch"].attrs:
                    time_attrs[attr_name] = input_file["Epoch"].attrs[attr_name]
                # Create the Time object
                ts["time"] = time_data
                # Create the Metadata
                ts["time"].meta = OrderedDict()
                ts["time"].meta.update(time_attrs)

            # Get all the Keys for Measurement Variable Data
            # These are Keys where the underlying object is a `dict` that contains
            # additional data, and is not the `EPOCH` variable
            variable_keys = filter(lambda key: key != "Epoch", list(input_file.keys()))
            # Add Variable Attributtes from the CDF file to TimeSeries
            for var_name in variable_keys:
                # Extract the Variable's Metadata
                var_attrs = {}
                for attr_name in input_file[var_name].attrs:
                    var_attrs[attr_name] = input_file[var_name].attrs[attr_name]

                # Extract the Variable's Data
                var_data = input_file[var_name][...]
                if input_file[var_name].rv():
                    # See if it is record-varying data with Units
                    if "UNITS" in var_attrs and len(var_data) == len(ts["time"]):
                        # Check if the variable is multi-dimensional
                        if len(var_data.shape) > 1:
                            try:
                                # Create an NDCube Object for the data
                                self._load_spectra_variable(
                                    spectra, var_name, var_data, var_attrs, ts.time
                                )
                            except ValueError:
                                warn_user(
                                    f"Cannot create NDCube for Spectra {var_name} with UNITS {var_attrs['UNITS']}. Creating Quantity with UNITS 'dimensionless_unscaled'."
                                )
                                # Swap Units
                                var_attrs["UNITS_DESC"] = var_attrs["UNITS"]
                                var_attrs[
                                    "UNITS"
                                ] = u.dimensionless_unscaled.to_string()
                                self._load_spectra_variable(
                                    spectra, var_name, var_data, var_attrs, ts.time
                                )
                        else:
                            # Load as Record-Varying `data`
                            try:
                                self._load_timeseries_variable(
                                    ts, var_name, var_data, var_attrs
                                )
                            except ValueError:
                                warn_user(
                                    f"Cannot create Quantity for Variable {var_name} with UNITS {var_attrs['UNITS']}. Creating Quantity with UNITS 'dimensionless_unscaled'."
                                )
                                # Swap Units
                                var_attrs["UNITS_DESC"] = var_attrs["UNITS"]
                                var_attrs[
                                    "UNITS"
                                ] = u.dimensionless_unscaled.to_string()
                                self._load_timeseries_variable(
                                    ts, var_name, var_data, var_attrs
                                )
                    else:
                        # Load as `support`
                        self._load_support_variable(
                            support, var_name, var_data, var_attrs
                        )
                else:
                    # Load Non-Record-Varying Data as `support`
                    self._load_support_variable(support, var_name, var_data, var_attrs)

        # Create a NDCollection
        if len(spectra) > 0:
            # Implement assertion that all spectra are aligned along time-varying dimension
            aligned_axes = tuple(0 for _ in spectra)
            spectra = NDCollection(spectra, aligned_axes=aligned_axes)
        else:
            spectra = NDCollection(spectra)

        # Return the given TimeSeries, NRV Data
        return ts, support, spectra

    def _load_timeseries_variable(self, ts, var_name, var_data, var_attrs):
        # Create the Quantity object
        var_data = u.Quantity(value=var_data, unit=var_attrs["UNITS"], copy=False)
        ts[var_name] = var_data
        # Create the Metadata
        ts[var_name].meta = OrderedDict()
        ts[var_name].meta.update(var_attrs)

    def _load_support_variable(self, support, var_name, var_data, var_attrs):
        # Create a NDData entry for the variable
        support[var_name] = NDData(data=var_data, meta=var_attrs)

    def _get_tensor_attribute(
        self, var_attrs, naxis, attribute_name, default_attribute
    ):
        """
        Function to get the `attribute_name` for each dimension of a multi-dimensional variable.

        For example if we have variable 'des_dist_brst' and we want to get the `.cunit` member
        for the WCS corresponding to the 'CUNIT' Keyword Attribute:
        - 'CUNIT1': 'eV'    (DEPEND_3: 'mms1_des_energy_brst')
        - 'CUNIT2': 'deg'   (DEPEND_2: 'mms1_des_theta_brst')
        - 'CUNIT3': 'deg'   (DEPEND_1: 'mms1_des_phi_brst' )
        - 'CUNIT4': 'ns'    (DEPEND_0: 'Epoch')

        We want to return a list of these units:
        ['eV', 'deg', 'deg', 'ns']
        """
        # Get `attribute_name` for each of the dimensions
        attr_values = []
        for dimension_i in range(naxis):
            dimension_attr_name = (
                f"{attribute_name}{dimension_i+1}"  # KeynameName Indexed 1-4 vs 0-3
            )
            if dimension_attr_name in var_attrs:
                attr_values.append(var_attrs[dimension_attr_name])
            else:
                attr_values.append(default_attribute)

        return attr_values

    def _get_world_coords(self, var_data, var_attrs, time):
        # Define WCS transformations in an astropy WCS object.

        # Get the N in var_attrs:
        if "WCSAXES" in var_attrs:
            # NOTE We have to cast this to an INT because spacepy does not let us directly set a
            # zAttr type when writing a variable attribute to a CDF. It tries to guess the type
            # of the attribute based on they type of the data.
            naxis = int(var_attrs["WCSAXES"])
        else:
            naxis = len(var_data.shape)
        wcs = WCS(naxis=naxis)

        for keyword, prop, default in self.schema.wcs_keyword_to_astropy_property:
            prop_value = self._get_tensor_attribute(
                var_attrs=var_attrs,
                naxis=naxis,
                attribute_name=keyword,
                default_attribute=default,
            )
            setattr(wcs.wcs, prop, prop_value)

        # wcs.wcs.ctype = 'WAVE', 'HPLT-TAN', 'HPLN-TAN'
        # wcs.wcs.cunit = 'keV', 'deg', 'deg'
        # wcs.wcs.cdelt = 0, 0, 0
        # wcs.wcs.crpix = 0, 0, 01
        # wcs.wcs.crval = 0, 0, 0
        # wcs.wcs.cname = 'wavelength', 'HPC lat', 'HPC lon'

        # TIME ATTRIBUTES
        wcs.wcs.timesys = "UTC"
        # Set the MJDREF (Modified  Julian Date Reference) to the start of the TimeSeries
        # An unexpected (feature?) of the WCS API is that MJDREF is an vector
        # attribute rather than a scalar attribute
        wcs.wcs.mjdref = [time[0].mjd, 0]
        wcs.wcs.timeunit = "ns"
        time_delta = time[1] - time[0]
        wcs.wcs.timedel = time_delta.to("ns").value

        return wcs

    def _load_spectra_variable(self, spectra, var_name, var_data, var_attrs, time):
        # Create a World Cordinate System for the Tensor
        var_wcs = self._get_world_coords(var_data, var_attrs, time)
        # Create a Cube
        var_cube = NDCube(
            data=var_data, wcs=var_wcs, meta=var_attrs, unit=var_attrs["UNITS"]
        )
        # Add to Spectra
        spectra.append((var_name, var_cube))

    def save_data(self, data, file_path):
        """
        Save heliophysics data to a CDF file.

        Parameters
        ----------
        data : `hermes_core.timedata.HermesData`
            An instance of `HermesData` containing the data to be saved.
        file_path : `str`
            The path to save the CDF file.

        Returns
        -------
        path : `str`
            A path to the saved file.
        """
        from spacepy.pycdf import CDF

        # Initialize a new CDF
        cdf_filename = f"{data.meta['Logical_file_id']}.cdf"
        output_cdf_filepath = str(Path(file_path) / cdf_filename)
        with CDF(output_cdf_filepath, masterpath="") as cdf_file:
            # Add Global Attriubtes to the CDF File
            self._convert_global_attributes_to_cdf(data, cdf_file)

            # Add zAttributes
            self._convert_variable_attributes_to_cdf(data, cdf_file)
        return output_cdf_filepath

    def _convert_global_attributes_to_cdf(self, data, cdf_file):
        # Loop though Global Attributes in target_dict
        for attr_name, attr_value in data.meta.items():
            # Make sure the Value is not None
            # We cannot add None Values to the CDF Global Attrs
            if attr_value is None:
                cdf_file.attrs[attr_name] = ""
            else:
                # Add the Attribute to the CDF File
                cdf_file.attrs[attr_name] = attr_value

    def _convert_variable_attributes_to_cdf(self, data, cdf_file):
        # Loop through Variable Attributes
        for var_name in data.timeseries.colnames:
            var_data = data.timeseries[var_name]
            if var_name == "time":
                # Add 'time' in the TimeSeries as 'Epoch' within the CDF
                cdf_file["Epoch"] = var_data.to_datetime()
                # Add the Variable Attributes
                for var_attr_name, var_attr_val in var_data.meta.items():
                    cdf_file["Epoch"].attrs[var_attr_name] = var_attr_val
            else:
                # Add the Variable to the CDF File
                cdf_file[var_name] = var_data.value
                # Add the Variable Attributes
                for var_attr_name, var_attr_val in var_data.meta.items():
                    if var_attr_val is None:
                        raise ValueError(
                            f"Variable {var_name}: Cannot Add vAttr: {var_attr_name}. Value was {str(var_attr_val)}"
                        )
                    else:
                        # Add the Attribute to the CDF File
                        cdf_file[var_name].attrs[var_attr_name] = var_attr_val

        # Loop through Non-Record-Varying Data
        for var_name, var_data in data.support.items():
            # Guess the data type to store
            # Documented in https://github.com/spacepy/spacepy/issues/707
            _, var_data_types, _ = self.schema._types(var_data.data)
            # Add the Variable to the CDF File
            cdf_file.new(
                name=var_name,
                data=var_data.data,
                type=var_data_types[0],
                recVary=False,
            )

            # Add the Variable Attributes
            for var_attr_name, var_attr_val in var_data.meta.items():
                if var_attr_val is None:
                    raise ValueError(
                        f"Variable {var_name}: Cannot Add vAttr: {var_attr_name}. Value was {str(var_attr_val)}"
                    )
                else:
                    # Add the Attribute to the CDF File
                    cdf_file[var_name].attrs[var_attr_name] = var_attr_val

        # Loop through High-Dimensional/Spectra Variables
        for var_name in data.spectra:
            var_data = data.spectra[var_name]
            # Add the Variable to the CDF File
            cdf_file[var_name] = var_data.data
            # Add the Variable Attributes
            for var_attr_name, var_attr_val in var_data.meta.items():
                if var_attr_val is None:
                    raise ValueError(
                        f"Variable {var_name}: Cannot Add vAttr: {var_attr_name}. Value was {str(var_attr_val)}"
                    )
                else:
                    # Add the Attribute to the CDF File
                    cdf_file[var_name].attrs[var_attr_name] = var_attr_val
