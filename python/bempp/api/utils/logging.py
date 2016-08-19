"""This module contains functions to initialize the Python logger"""
from __future__ import absolute_import
import logging as _logging

# Logging levels

DEBUG = _logging.DEBUG
INFO = _logging.INFO
WARNING = _logging.WARNING
ERROR = _logging.ERROR
CRITICAL = _logging.CRITICAL

DEFAULT_FORMAT = '%(asctime)s:%(name)s:%(levelname)s: %(message)s'

def _init_logger():

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
    import time
    from bempp.api import LOGGER

    def timeit_impl(fun):

        def timed_fun(*args, **kwargs):

            st = time.time()
            res = fun(*args, **kwargs)
            et = time.time()
            LOGGER.info(message + " : {0:.3e}s".format(et-st))
            return res

        return timed_fun

    return timeit_impl
