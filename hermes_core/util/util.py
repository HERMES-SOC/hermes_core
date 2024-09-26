"""
This module provides general utility functions.
"""

from swxsoc.util import util


__all__ = ["create_science_filename", "parse_science_filename"]


def create_science_filename(
    instrument: str,
    time: str,
    level: str,
    version: str,
    mode: str = "",
    descriptor: str = "",
    test: bool = False,
):
    """Return a compliant filename. The format is defined as

    hermes_{inst}_{mode}_{level}{test}_{descriptor}_{time}_v{version}.cdf

    This format is only appropriate for data level >= 1.

    Parameters
    ----------
    instrument : `str`
        The instrument name. Must be one of the following "eea", "nemesis", "merit", "spani"
    time : `str` (in isot format) or ~astropy.time
        The time
    level : `str`
        The data level. Must be one of the following "l0", "l1", "l2", "l3", "l4", "ql"
    version : `str`
        The file version which must be given as X.Y.Z
    descriptor : `str`
        An optional file descriptor.
    mode : `str`
        An optional instrument mode.
    test : bool
        Selects whether the file is a test file.

    Returns
    -------
    filename : `str`
        A CDF file name including the given parameters that matches the HERMES file naming conventions

    Raises
    ------
    ValueError: If the instrument is not recognized as one of the HERMES instruments
    ValueError: If the data level is not recognized as one of the HERMES valid data levels
    ValueError: If the data version does not match the HERMES data version formatting conventions
    ValueError: If the data product descriptor or instrument mode do not match the HERMES formatting conventions
    """
    return util.create_science_filename(
        instrument, time, level, version, mode, descriptor, test
    )


def parse_science_filename(filepath: str) -> dict:
    """
    Parses a science filename into its consitutient properties (instrument, mode, test, time, level, version, descriptor).

    Parameters
    ----------
    filepath: `str`
        Fully specificied filepath of an input file

    Returns
    -------
    result : `dict`
        A dictionary with each property.

    Raises
    ------
    ValueError: If the file's mission name is not "HERMES"
    ValueError: If the file's instreument name is not one of the HERMES instruments
    ValueError: If the data level >0 for packet files
    ValueError: If not a CDF File
    """
    return util.parse_science_filename(filepath)
