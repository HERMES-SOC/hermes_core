# see license/LICENSE.rst
try:
    from ._version import version as __version__
    from ._version import version_tuple
except ImportError:
    __version__ = "unknown version"
    version_tuple = (0, 0, "unknown version")

from hermes_core.util.config import load_config, print_config
from hermes_core.util.logger import _init_log

# Load user configuration
config = load_config()

log = _init_log(config=config)

# Then you can be explicit to control what ends up in the namespace,
__all__ = ["config", "print_config"]

MISSION_NAME = "hermes"
INST_NAMES = ["eea", "nemisis", "merit", "spani"]
INST_SHORTNAMES = ["eea", "nms", "mrt", "spn"]
INST_TARGETNAMES = ["EEA", "MAG", "MERIT", "SPANI"]
INST_TO_SHORTNAME = dict(zip(INST_NAMES, INST_SHORTNAMES))
INST_TO_TARGETNAME = dict(zip(INST_NAMES, INST_TARGETNAMES))

# log.info(f"hermes_core version: {__version__}")
