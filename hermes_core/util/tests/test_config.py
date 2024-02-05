"""
Tests for the config module
"""

import io
import os
from pathlib import Path
from contextlib import redirect_stdout

import pytest

import hermes_core
from hermes_core import config
from hermes_core.util import HERMESUserWarning
from hermes_core.util.config import (
    CONFIG_DIR,
    _find_config_files,
    _get_user_configdir,
    _is_writable_dir,
    copy_default_config,
    dirs,
    get_and_create_download_dir,
    get_and_create_sample_dir,
    print_config,
)

USER = os.path.expanduser("~")


def test_is_writable_dir(tmpdir, tmp_path):
    assert _is_writable_dir(tmpdir)
    tmp_file = tmpdir.join("hello.txt")
    # Have to write to the file otherwise its seen as a directory(?!)
    tmp_file.write("content")
    # Checks directory with a file
    assert _is_writable_dir(tmpdir)
    # Checks a filepath instead of directory
    assert not _is_writable_dir(tmp_file)
