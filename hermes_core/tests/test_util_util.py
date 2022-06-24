"""Tests for util.py"""
from hermes_core.util import util
import pytest


def test_science_filename_exceptions():
    """Test for errors"""
    good_time = "2025-06-02T12:04:01"
    good_instrument = "eea"
    good_level = "l0"
    good_version = "1.3.4"
    with pytest.raises(ValueError):
        # not enough depth to version number
        util.create_science_filename(good_instrument, good_time, level=good_level, version="1.3")
        util.create_science_filename(good_instrument, good_time, level=good_level, version="1")

        # wrong level specification
        util.create_science_filename(good_instrument, good_time, level="la", version=good_version)
        util.create_science_filename(good_instrument, good_time, level="squirrel", version=good_version)

        # wrong instrument name
        util.create_science_filename("eeb", good_time, level=good_level, version=good_version)
        util.create_science_filename("fpi", good_time, level=good_level, version=good_version)
        util.create_science_filename("potato", good_time, level=good_level, version=good_version)

