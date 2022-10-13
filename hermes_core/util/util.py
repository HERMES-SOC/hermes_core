"""
This module provides general utility functions.
"""
import os

from astropy.time import Time

import hermes_core


__all__ = ["create_science_filename"]

TIME_FORMAT_L0 = "%Y%j-%H%M%S"
TIME_FORMAT = "%Y%m%d_%H%M%S"
VALID_DATA_LEVELS = ["l0", "l1", "ql", "l2", "l3", "l4"]


def create_science_filename(
    instrument, time, level, version, mode="", descriptor="", test=False
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
    """
    test_str = ""

    if isinstance(time, str):
        time_str = Time(time, format="isot").strftime(TIME_FORMAT)
    else:
        time_str = time.strftime(TIME_FORMAT)

    if instrument not in hermes_core.INST_NAMES:
        raise ValueError(
            f"Instrument, {instrument}, is not recognized. Must be one of {hermes_core.INST_NAMES}.")
    if level not in VALID_DATA_LEVELS[1:]:
        raise ValueError(
            f"Level, {level}, is not recognized. Must be one of {VALID_DATA_LEVELS[1:]}."
        )
    # check that version is in the right format with three parts
    if len(version.split(".")) != 3:
        raise ValueError(
            f"Version, {version}, is not formatted correctly. Should be X.Y.Z"
        )
    # check that version has integers in each part
    for item in version.split("."):
        try:
            int_value = int(item)
        except ValueError:
            raise ValueError(f"Version, {version}, is not all integers.")

    if test is True:
        test_str = "test"

    # the parse_science_filename function depends on _ not being present elsewhere
    if ("_" in mode) or ("_" in descriptor):
        raise ValueError(
            "The underscore symbol _ is not allowed in mode or descriptor."
        )

    filename = f"hermes_{hermes_core.INST_TO_SHORTNAME[instrument]}_{mode}_{level}{test_str}_{descriptor}_{time_str}_v{version}"
    filename = filename.replace("__", "_")  # reformat if mode or descriptor not given

    return filename + FILENAME_EXTENSION


def parse_science_filename(filename):
    """ "
    Parses a science filename into its consitutient properties (instrument, mode, test, time, level, version, descriptor).

    Parameters
    ----------
    filename: `str`
        The filename.

    Returns
    -------
    result : `dict`
        A dictionary with each property.
    """

    result = {
        "instrument": None,
        "mode": None,
        "test": False,
        "time": None,
        "level": None,
        "version": None,
        "descriptor": None,
    }

    file_name, file_ext = os.path.splitext(filename)

    filename_components = file_name.split("_")

    if filename_components[0] != hermes_core.MISSION_NAME:
        raise ValueError(
            f"File {filename} not recognized. Not a valid mission name."
        )

    if file_ext == ".bin":
        if filename_components[1] not in hermes_core.INST_TARGETNAMES:
            raise ValueError(
                f"File {filename} not recognized. Not a valid target name."
            )
        if filename_components[2] != VALID_DATA_LEVELS[0]:
            raise ValueError(
                f"Data level {filename_components[2]} is not correct for this file extension."
            )
        else:
            result["level"] = filename_components[2]
        #  reverse the dictionary to look up instrument name from the short name
        from_shortname = {v: k for k, v in hermes_core.INST_TO_TARGETNAME.items()}

        result["time"] = Time.strptime(filename_components[3], TIME_FORMAT_L0)

    elif file_ext == ".cdf":
        if filename_components[1] not in hermes_core.INST_SHORTNAMES:
            raise ValueError(
                "File {filename} not recognized. Not a valid instrument name."
            )

        #  reverse the dictionary to look up instrument name from the short name
        from_shortname = {v: k for k, v in hermes_core.INST_TO_SHORTNAME.items()}

        result["time"] = Time.strptime(
            filename_components[-3] + "_" + filename_components[-2], TIME_FORMAT
        )

        # mode and descriptor are optional so need to figure out if one or both or none is included
        if filename_components[2][0:2] not in VALID_DATA_LEVELS:
            # if the first component is not data level then it is mode and the following is data level
            result["mode"] = filename_components[2]
            result["level"] = filename_components[3].replace("test", "")
            if "test" in filename_components[3]:
                result["test"] = True
            if len(filename_components) == 8:
                result["descriptor"] = filename_components[4]
        else:
            result["level"] = filename_components[2].replace("test", "")
            if "test" in filename_components[2]:
                result["test"] = True
            if len(filename_components) == 7:
                result["descriptor"] = filename_components[3]
    else:
        raise ValueError(f"File extension {file_ext} not recognized.")

    result["instrument"] = from_shortname[filename_components[1]]
    result["version"] = filename_components[-1][1:]  # remove the v

    return result
