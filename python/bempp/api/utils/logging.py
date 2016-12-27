"""This module contains functions to initialize the Python logger"""
#pylint: disable=import-self
#pylint: disable=no-member
#pylint: disable=redefined-builtin
#pylint: disable=invalid-name

from __future__ import absolute_import
import logging as _logging
import time as _time

# Logging levels

DEBUG = _logging.DEBUG
INFO = _logging.INFO
WARNING = _logging.WARNING
ERROR = _logging.ERROR
CRITICAL = _logging.CRITICAL

DEFAULT_FORMAT = '%(asctime)s:%(name)s:%(levelname)s: %(message)s'

def _init_logger():
    """Initialize the BEM++ logger."""

    logger = _logging.getLogger('BEMPP')
    logger.setLevel(_logging.DEBUG)
    logger.addHandler(_logging.NullHandler())
    return logger


def enable_console_logging(level=DEBUG, format=DEFAULT_FORMAT):
    """Enable console logging and return the console handler."""

    from bempp.api import LOGGER
    ch = _logging.StreamHandler()
    ch.setLevel(level)
    ch.setFormatter(_logging.Formatter(format))
    LOGGER.addHandler(ch)
    return ch


def enable_file_logging(file_name, level=DEBUG, format=DEFAULT_FORMAT):
    """Enable logging to a specific file."""

    from bempp.api import LOGGER
    fh = _logging.FileHandler(file_name)
    fh.setLevel(level)
    fh.setFormatter(_logging.Formatter(format))
    LOGGER.addHandler(fh)
    return fh


def set_logging_level(level):
    """Set the logging level."""

    from bempp.api import LOGGER
    LOGGER.setLevel(level)


def timeit(message):
    """Decorator to time a method in BEM++"""
    from bempp.api import LOGGER
    from bempp.api import global_parameters

    def timeit_impl(fun):
        """Implementation of timeit."""
        def timed_fun(*args, **kwargs):
            """The actual timer function."""
            if not global_parameters.verbosity.extended_verbosity:
                return fun(*args, **kwargs)

            st = _time.time()
            res = fun(*args, **kwargs)
            et = _time.time()
            LOGGER.info(message + " : {0:.3e}s".format(et-st))
            return res

        return timed_fun

    return timeit_impl

#pylint: disable=too-few-public-methods
class Timer:
    """Context manager to measure time in BEM++."""

    def __init__(self):
        """Constructor."""
        self.start = 0
        self.end = 0
        self.interval = 0

    def __enter__(self):
        self.start = _time.time()
        return self


    def __exit__(self, *args):
        self.end = _time.time()
        self.interval = self.end - self.start
