"""
This module provides errors/exceptions and warnings of general use.

Exceptions that are specific to a given package should **not** be here,
but rather in the particular package.

This code is based on that provided by SunPy see
    licenses/SUNPY.rst 
"""
import warnings

__all__ = [
    "HERMESWarning",
    "HERMESUserWarning",
    "HERMESDeprecationWarning",
    "HERMESPendingDeprecationWarning",
    "warn_user",
    "warn_deprecated",
    "warn_metadata",
]


class HERMESWarning(Warning):
    """
    The base warning class from which all Sunpy warnings should inherit.

    Any warning inheriting from this class is handled by the Sunpy
    logger. This warning should not be issued in normal code. Use
    "SunpyUserWarning" instead or a specific sub-class.
    """


class HERMESUserWarning(UserWarning, HERMESWarning):
    """
    The primary warning class for Sunpy.

    Use this if you do not need a specific type of warning.
    """


class HERMESDeprecationWarning(FutureWarning, HERMESWarning):
    """
    A warning class to indicate a deprecated feature.
    """


class HERMESPendingDeprecationWarning(PendingDeprecationWarning, HERMESWarning):
    """
    A warning class to indicate a soon-to-be deprecated feature.
    """


def warn_metadata(msg, stacklevel=1):
    """
    Raise a `SunpyMetadataWarning`.

    Parameters
    ----------
    msg : str
        Warning message.
    stacklevel : int
        This is interpreted relative to the call to this function,
        e.g. ``stacklevel=1`` (the default) sets the stack level in the
        code that calls this function.
    """
    warnings.warn(msg, HERMESMetadataWarning, stacklevel + 1)


def warn_user(msg, stacklevel=1):
    """
    Raise a `HERMESUserWarning`.

    Parameters
    ----------
    msg : str
        Warning message.
    stacklevel : int
        This is interpreted relative to the call to this function,
        e.g. ``stacklevel=1`` (the default) sets the stack level in the
        code that calls this function.
    """
    warnings.warn(msg, HERMESUserWarning, stacklevel + 1)


def warn_deprecated(msg, stacklevel=1):
    """
    Raise a `HERMESDeprecationWarning`.

    Parameters
    ----------
    msg : str
        Warning message.
    stacklevel : int
        This is interpreted relative to the call to this function,
        e.g. ``stacklevel=1`` (the default) sets the stack level in the
        code that calls this function.
    """
    warnings.warn(msg, HERMESDeprecationWarning, stacklevel + 1)
