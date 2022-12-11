"""
A module for all things calibration.
"""
import numpy as np

import ccsdspy
from ccsdspy import PacketField

from hermes_core.util.util import create_science_filename, parse_science_filename
from hermes_core.io import read_file

__all__ = ["calibrate_file", "get_calibration_file", "read_calibration_file", "process_file", "plot_file"]

# the number of sun sensor vector readings in a vector packet
SUNSENSOR_NUM_VECTORS = 20


def process_file(data_filename: str) -> list:
    """
    This is the entry point for the pipeline processing. It runs all of the various processing steps required.

    Parameters
    ----------
    data_filename: str
        Fully specificied filename of an input file
    
    Returns
    -------
    output_filenames: list
        Fully specificied filenames the output filenames
    """
    output_files = []

    calibrated_file = calibrate_file(data_filename)
    #data_plot_files = plot_file(data_filename)
    #calib_plot_files = plot_file(calibrated_file)

    # add other tasks below
    return output_files


def calibrate_file(data_filename):
    """
    Given an input data file, raise it to the next level (e.g. level 0 to level 1, level 1 to quicklook) it and return a new file.

    Parameters
    ----------
    data_filename: str
        Fully specificied filename of the non-calibrated file (data level < 2)

    Returns
    -------
    output_filename: str
        Fully specificied filename of the non-calibrated file (data level < 2)

    Examples
    --------
    >>> from hermes_core.calibration import calibrate_file
    >>> level1_file = calibrate_file('hermes_SS_l0_2022239-000000_v0.bin')
    """
    # this function currently assumes that binary files have a single APID

    filedata = parse_science_filename(data_filename)

    if filedata["instrument"] == "sunsensor" and filedata["level"] == "l0":
        data = parse_sunsensor_vector_packets(data_filename)
        output_filename = sunsensor_vector_data_to_cdf(data, filedata)

    if filedata["instrument"] == "sunsensor" and filedata["level"] == "l1":
        data = read_file(data_filename)

    if filedata["level"] == "l2":
        calib_file = get_calibration_file(data_filename)
        if calib_file is None:
            raise ValueError("Calibration file for {} not found.".format(data_filename))
        else:
            calib_data = read_calibration_file(calib_file)

    return output_filename


def parse_sunsensor_vector_packets(data_filename):
    """
    Parse a level 0 sunsensor vector packet file.

    Parameters
    ----------
    data_filename: str
        Fully specificied filename

    Returns
    -------
    result: dict
        A dictionary of arrays which includes the ccsds header fields

    Examples
    --------
    >>> import hermes_core.calibration as calib
    >>> data_filename = "hermes_SS_l0_2022339-000000_v0.bin"
    >>> data_packets = calib.parse_sunsensor_vector_packets(data_filename)
    """

    # hardcoded packet fields because there are so many vector fields

    fields = {
        "SECONDS": ["uint", 32],
        "SUBSECONDS": ["uint", 16],
        "NUMVECTORSCONTAINED": ["int", 16],
        "VECTORPOLLINTERVAL": ["int", 16],
        "TEMPERATURE": ["uint", 16],
        "PROCESSINGTIME": ["uint", 32],
    }

    for i in range(SUNSENSOR_NUM_VECTORS):
        fields.update({f"X_INT_{i:02}": ["int", 16]})
        fields.update({f"Y_INT_{i:02}": ["int", 16]})
        fields.update({f"Z_INT_{i:02}": ["int", 16]})
        fields.update({f"FIT_QUALITY_{i:02}": ["int", 8]})
        fields.update({f"GEOMETRY_QUALITY_{i:02}": ["int", 8]})
        fields.update({f"VALIDITY_{i:02}": ["uint", 8]})
        fields.update({f"BYTEFILLER_{i:02}": ["uint", 8]})

    packet_fields = [
        PacketField(name=this_name, data_type=result[0], bit_length=result[1])
        for this_name, result in fields.items()
    ]

    pkt = ccsdspy.FixedLength(packet_fields)
    result = pkt.load(data_filename, include_primary_header=True)
    return result


def sunsensor_vector_data_to_cdf(data, metadata):
    """
    Write level 0 sunsensor vector data to a cdf file.

    Parameters
    ----------
    data: dict
        A dictionary of arrays which includes the ccsds header fields
    metadata: dict
        A metadata dictionary from the originating binary file

    Returns
    -------
    output_filename: str
        Fully specificied filename of cdf file

    Examples
    --------
    >>> from hermes_core.util.util import parse_science_filename
    >>> import hermes_core.calibration as calib
    >>> data_filename = "hermes_SS_l0_2022339-000000_v0.bin"
    >>> metadata = parse_science_filename(data_filename)
    >>> data_packets = calib.parse_sunsensor_vector_packets(data_filename)
    >>> cdf_filename = calib.sunsensor_vector_data_to_cdf(data_packets, metadata)
    """

    num_packets = len(data["SECONDS"])
    time = (
        data["SECONDS"] + data["SUBSECONDS"] / 1000.0
    )  # assume that subseconds are milliseconds

    vector_data = np.zeros((3, SUNSENSOR_NUM_VECTORS, num_packets), dtype="uint16")
    metadata_int = np.zeros((2, SUNSENSOR_NUM_VECTORS, num_packets), dtype="int8")
    metadata_uint = np.zeros((2, SUNSENSOR_NUM_VECTORS, num_packets), dtype="uint8")

    for i in range(num_packets):
        vector_data[0, :, i] = [
            data[f"X_INT_{j:02}"][i] for j in range(SUNSENSOR_NUM_VECTORS)
        ]
        vector_data[1, :, i] = [
            data[f"Y_INT_{j:02}"][i] for j in range(SUNSENSOR_NUM_VECTORS)
        ]
        vector_data[2, :, i] = [
            data[f"Z_INT_{j:02}"][i] for j in range(SUNSENSOR_NUM_VECTORS)
        ]
        metadata_int[0, :, i] = [
            data[f"FIT_QUALITY_{j:02}"][i] for j in range(SUNSENSOR_NUM_VECTORS)
        ]
        metadata_int[1, :, i] = [
            data[f"GEOMETRY_QUALITY_{j:02}"][i] for j in range(SUNSENSOR_NUM_VECTORS)
        ]
        metadata_uint[0, :, i] = [
            data[f"VALIDITY_{j:02}"][i] for j in range(SUNSENSOR_NUM_VECTORS)
        ]
        metadata_uint[1, :, i] = [
            data[f"BYTEFILLER_{j:02}"][i] for j in range(SUNSENSOR_NUM_VECTORS)
        ]

    # now open a cdf file and write these data into it
    cdf_filename = ""  # cdf_filename = create_science_filename()

    return cdf_filename


def get_calibration_file(data_filename: str, time=None) -> str:
    """
    Given a time, return the appropriate calibration file.

    Parameters
    ----------
    data_filename: str
        Fully specificied filename of the non-calibrated file (data level < 2)
    time: ~astropy.time.Time

    Returns
    -------
    calib_filename: str
        Fully specificied filename for the appropriate calibration file.

    Examples
    --------
    """
    return ''


def read_calibration_file(calib_filename: str) -> :
    """
    Given a calibration, return the calibration structure.

    Parameters
    ----------
    calib_filename: str
        Fully specificied filename of the non-calibrated file (data level < 2)

    Returns
    -------
    output_filename: str
        Fully specificied filename of the appropriate calibration file.

    Examples
    --------
    """

    # if can't read the file

    return None
