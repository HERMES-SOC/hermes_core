"""Tests for util.py"""
import pytest
import tempfile

from astropy.time import Time
from hermes_core.util import util


def test_science_filename_output():
    """Test expected output"""
    time = "2024-04-06T12:06:21"
    time_formatted = "20240406_120621"

    assert util.create_science_filename(
        "eea", time, level="l0", version="1.2.3"
    ) == "hermes_eea_l0_{}_v1.2.3".format(time_formatted)
    # merit
    assert util.create_science_filename(
        "merit", time, level="l2", version="2.4.5"
    ) == "hermes_mrt_l2_{}_v2.4.5".format(time_formatted)
    # nemisis and version
    assert util.create_science_filename(
        "nemisis", time, level="l2", version="1.3.5"
    ) == "hermes_nms_l2_{}_v1.3.5".format(time_formatted)
    # spani and level
    assert util.create_science_filename(
        "spani", time, level="l3", version="2.4.5"
    ) == "hermes_spn_l3_{}_v2.4.5".format(time_formatted)
    # mode
    assert util.create_science_filename(
        "spani", time, level="l3", mode="2s", version="2.4.5"
    ) == "hermes_spn_2s_l3_{}_v2.4.5".format(time_formatted)
    # test
    assert util.create_science_filename(
        "spani", time, level="l0", version="2.4.5", test=True
    ) == "hermes_spn_l0test_{}_v2.4.5".format(time_formatted)
    # all options
    assert util.create_science_filename(
        "spani",
        time,
        level="l3",
        mode="2s",
        descriptor="burst",
        version="2.4.5",
        test=True,
    ) == "hermes_spn_2s_l3test_burst_{}_v2.4.5".format(time_formatted)
    # Time object instead of str
    assert util.create_science_filename(
        "spani",
        Time(time),
        level="l3",
        mode="2s",
        descriptor="burst",
        version="2.4.5",
        test=True,
    ) == "hermes_spn_2s_l3test_burst_{}_v2.4.5".format(time_formatted)
    # Time object but created differently
    assert util.create_science_filename(
        "spani",
        Time(2460407.004409722, format="jd"),
        level="l3",
        mode="2s",
        descriptor="burst",
        version="2.4.5",
        test=True,
    ) == "hermes_spn_2s_l3test_burst_{}_v2.4.5".format(time_formatted)


def test_science_filename_exceptions():
    """Test for errors"""
    good_time = "2025-06-02T12:04:01"
    good_instrument = "eea"
    good_level = "l0"
    good_version = "1.3.4"
    with pytest.raises(ValueError):
        # not enough depth to version number
        util.create_science_filename(
            good_instrument, good_time, level=good_level, version="1.3"
        )
        util.create_science_filename(
            good_instrument, good_time, level=good_level, version="1"
        )
        util.create_science_filename(
            good_instrument, good_time, level=good_level, version="1.5.6.7"
        )
        util.create_science_filename(
            good_instrument, good_time, level=good_level, version="1.."
        )
        # a letter in version number
        util.create_science_filename(
            good_instrument, good_time, level=good_level, version="a.5.6"
        )

        # wrong level specification
        util.create_science_filename(
            good_instrument, good_time, level="la", version=good_version
        )
        util.create_science_filename(
            good_instrument, good_time, level="squirrel", version=good_version
        )

        # wrong instrument name
        util.create_science_filename(
            "eeb", good_time, level=good_level, version=good_version
        )
        util.create_science_filename(
            "fpi", good_time, level=good_level, version=good_version
        )
        util.create_science_filename(
            "potato", good_time, level=good_level, version=good_version
        )

        # bad time string
        # non-existent time
        util.create_science_filename(
            good_instrument,
            "2023-13-04T12:06:21",
            level=good_level,
            version=good_version,
        )
        # not isot format
        util.create_science_filename(
            "eeb", "2023/13/04 12:06:21", level=good_level, version=good_version
        )
        # not valid input for time
        util.create_science_filename(
            "eeb", time=12345345, level=good_level, version=good_version
        )
        # _ character in mode
        util.create_science_filename(
            "eeb", time=12345345, level=good_level, version=good_version, mode="o_o"
        )
        # _ character in descriptor
        util.create_science_filename(
            "eeb",
            time=12345345,
            level=good_level,
            version=good_version,
            descriptor="blue_green",
        )


def test_parse_science_filename_output():
    """Test for known outputs"""
    # all parameters
    input = {
        "instrument": "spani",
        "mode": "2s",
        "level": "l3",
        "test": False,
        "descriptor": "burst",
        "version": "2.4.5",
        "time": Time("2024-04-06T12:06:21"),
    }

    f = util.create_science_filename(
        input["instrument"],
        input["time"],
        input["level"],
        input["version"],
        test=input["test"],
        descriptor=input["descriptor"],
        mode=input["mode"],
    )
    assert util.parse_science_filename(f) == input

    # test only
    input = {
        "instrument": "nemisis",
        "level": "l3",
        "test": True,
        "version": "2.4.5",
        "time": Time("2024-04-06T12:06:21"),
        "mode": None,
        "descriptor": None,
    }

    f = util.create_science_filename(
        input["instrument"],
        input["time"],
        input["level"],
        input["version"],
        test=input["test"],
    )
    assert util.parse_science_filename(f) == input

    # descriptor only
    input = {
        "instrument": "spani",
        "mode": None,
        "level": "l3",
        "test": False,
        "descriptor": "burst",
        "version": "2.4.5",
        "time": Time("2024-04-06T12:06:21"),
    }

    f = util.create_science_filename(
        input["instrument"],
        input["time"],
        input["level"],
        input["version"],
        descriptor=input["descriptor"],
    )
    assert util.parse_science_filename(f) == input

    # mode only
    input = {
        "instrument": "nemisis",
        "mode": "2s",
        "level": "l2",
        "test": False,
        "descriptor": None,
        "version": "2.7.9",
        "time": Time("2024-04-06T12:06:21"),
    }

    f = util.create_science_filename(
        input["instrument"],
        input["time"],
        input["level"],
        input["version"],
        mode=input["mode"],
    )
    assert util.parse_science_filename(f) == input


def test_parse_science_filename_errors():
    """Test for errors"""
    with pytest.raises(ValueError):
        # wrong mission name
        f = "veeger_spn_2s_l3test_burst_20240406_120621_v2.4.5"
        util.parse_science_filename(f)

        # wrong instrument name
        f = "hermes_www_2s_l3test_burst_20240406_120621_v2.4.5"
        util.parse_science_filename(f)


testdata = [
    (b'Hello world', '3e25960a79dbc69b674cd4ec67a72c62'),
    (b'Hermes is the best!', 'fc31437934081bf3c7f5552908e3a440'),
    (b'to the moon!', 'fe83120ed5dbbf8ddc1ae0aef898efc7'),
]

@pytest.mark.parametrize("text,md5", testdata)
def test_hash_file(text, md5):
    """Test with known inputs and outputs"""
    f = tempfile.NamedTemporaryFile(delete=False)
    f.write(text)
    f.close()
    assert util.hash_file(f.name) == md5
