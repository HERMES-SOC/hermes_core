"""
This module provides a generic file reader.
"""
from astropy.timeseries import TimeSeries
from astropy import units as u

from hermes_core.util import parse_science_filename

__all__ = ["read_file"]


def read_file(data_filename: str) -> TimeSeries:
    """
    Read a level 1 or greater data file.

    Parameters
    ----------
    data_filename: str
        Fully specificied filename of an input file.

    Returns
    -------
    data: astropy.timeseries.TimeSeries
        A data structure.

    Examples
    --------
    >>> from hermes_core.io import read_file
    >>> data = read_file('hermes_sunsensor_l1_20220501_000000_v1.0.0.cdf')
    """

    file_metadata = parse_science_filename(data_filename)
    if file_metadata['level'] == 'l0':
        raise ValueError("This function does not support level 0 files.")

    # open the file and extract the data and put it into a TimeSeries if possible


    # example of a data structure
    data = TimeSeries(time_start='2016-03-22T12:30:31', time_delta=3 * u.s, n_samples=5)
    data['electron flux'] = range(5) * u.electron / (u.s * u.cm ** 2 * u.sr)

    return data