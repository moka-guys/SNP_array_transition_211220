#!/usr/bin/env python3
"""logger.py

Log messages using the python standard library logging module.
"""
import sys
import logging


class Logger(object):
    """
    Simple logging class

    Methods
        shutdown_logs()
            To prevent duplicate filehandlers and system handlers close and
            remove all handlers for all log files with a python logging object
        _get_file_handler()
            Returns the FileHandler associated with the logging object
        _get_stream_handler()
            Returns the StreamHandler associated with the logging object
        get_logger()
            Return a Python logging object
    """

    _formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s"
    )  # Log string format

    def __init__(self, logfile_path: str):
        """
        Constructor for the Logger class
            :param logfile_path (str): Logfile path
        """
        self.logger = self.get_logger("logger", logfile_path)

    def shutdown_logs(self):
        """
        To prevent duplicate filehandlers and system handlers close and
        remove all handlers for all log files that have a python logging object
        """
        for handler in self.logger.handlers[:]:
            self.logger.removeHandler(handler)
            handler.close()

    def _get_file_handler(self, filepath: str) -> logging.FileHandler:
        """
        Returns the FileHandler associated with the logging object
            :return file_handler (obj):   FileHandler object
        """
        file_handler = logging.FileHandler(filepath, mode="a", delay=True)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(self._formatter)
        return file_handler

    def _get_stream_handler(self) -> logging.StreamHandler:
        """
        Returns the StreamHandler associated with the logging object
            :return stream_handler (obj):   StreamHandler object
        """
        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.setLevel(logging.DEBUG)
        stream_handler.setFormatter(self._formatter)
        return stream_handler

    def get_logger(self, name: str, filepath: str) -> logging.Logger:
        """
        Return a Python logging object
            :param name (str):       Logger name
            :param filepath (str):   Logfile path
            :return logger (obj):   Python logging object
        """
        logger = logging.getLogger(name)
        logger.filepath = filepath
        logger.setLevel(logging.DEBUG)
        logger.addHandler(self._get_file_handler(filepath))
        logger.addHandler(self._get_stream_handler())
        return logger
