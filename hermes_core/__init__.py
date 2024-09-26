# see license/LICENSE.rst
try:
    from ._version import version as __version__
    from ._version import version_tuple
except ImportError:
    __version__ = "unknown version"
    version_tuple = (0, 0, "unknown version")

import os
import swxsoc

# Set the Mission Environment Variable to load the correct configuration
os.environ["SWXSOC_MISSION"] = "hermes"
swxsoc._reconfigure()
# Get the Updated Configuration
config = swxsoc.config

# Initialize HERMES Logger
log = swxsoc.util.logger._init_log(config)

# Then you can be explicit to control what ends up in the namespace,
__all__ = ["config"]

# log.info(f"hermes_core version: {__version__}")
