"""
This module provides general utility functions.
"""
from astropy.time import Time

__all__ = ["create_science_filename"]


def create_science_filename(
    instrument, time: str, mode="", level="l0", version="", descriptor="", test=False
):
    """Return a compliant filename root (without extension). The format is defined as

    hermes_{inst}_{mode}_{level}{test}_{descriptor}_{time}_v{version}

    Parameters
    ----------
    instrument : `str`
        The instrument name. Must be one of the following "eea", "nemesis", "merit", "spani"
    time : `str` (in isot format) or ~astropy.time
        The time
    level : str
        The data level. Must be one of the following
    version : str
        The file version which must be given as X.Y.Z
    descriptor : str
        An optional file descriptor.
    test : bool
        Selects whether the file is a test file.

    """
    time_format = "%Y%m%d_%H%M%S"
    test_str = ""
    instrument_shortnames = {"nemisis": "nms", "eea": "eea", "merit": "mrt", "spani": "spn"}
    valid_data_levels = ["l0", "l1", "ql", "l2", "l3", "l4"]

    if isinstance(time, str):
        time_str = Time(time, format="isot").strftime(time_format)
    else:
        time_str = time.strftime(time_format)

    if not instrument in instrument_shortnames.keys():
        raise ValueError("Instrument, {inst}, is not recognized. Must be ".format(inst=instrument))
    if level is valid_data_levels:
        raise ValueError("Level, {level}, is not recognized. Must be".format(level=level))
    if len(version.split(".")) != 3:
        raise ValueError(
            "Version, {version}, is not formatted correctly. Should be X.Y.Z".format(
                version=version
            )
        )
    if test is True:
        test_str = "test"

    filename = "hermes_{inst}_{mode}_{level}{test}_{descriptor}_{time}_v{version}".format(
        inst=instrument_shortnames[instrument],
        mode=mode,
        level=level,
        test=test_str,
        descriptor=descriptor,
        time=time_str,
        version=version,
    )
    filename = filename.replace("__", "_")  # reformat if mode or descriptor not given

    return filename
