# see license/LICENSE.rst
try:
    from ._version import version as __version__
    from ._version import version_tuple
except ImportError:
    __version__ = "unknown version"
    version_tuple = (0, 0, "unknown version")

from hermes_core.util.config import load_config, print_config
from swxsoc.util.logger import _init_log
import swxsoc

# Load user configuration
config = load_config()
# Overwrite SWxSOC Config with HERMES Config
swxsoc.config = config

log = _init_log(config=config)

# Then you can be explicit to control what ends up in the namespace,
__all__ = ["config", "print_config"]

# log.info(f"hermes_core version: {__version__}")
