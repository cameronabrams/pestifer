# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
Defines the :class:`LogParser` class and its subclasses for parsing log files from various applications used by Pestifer, such as NAMD and Packmol.
"""
from __future__ import annotations
import logging
import time

from pathlib import Path

from ..util.progress import PestiferProgress
from ..util.stringthings import ByteCollector

logger = logging.getLogger(__name__)

def get_toflag(flag: str, bytes: str):
    """
    Get the substring of bytes up to the first occurrence of flag.
    If the flag is not found, return None and log a debug message.
    
    Parameters
    ----------
    flag : str
        The flag to search for in the bytes.
    bytes : bytes
        The bytes to search within.

    Returns
    -------
    bytes
        The substring of bytes up to the first occurrence of flag, or None if not found.
    """
    try:
        idx = bytes.index(flag)
    except ValueError:
        logger.debug(f'get_toflag: {flag} not found in {bytes}')
        return None
    substr = bytes[:idx]
    return substr.strip()

def get_toeol(flag: str, bytes: str):
    """
    Get the substring of bytes from the first occurrence of flag to the end of the line.
    If the flag is not found or there is no newline character after the flag, return None and log a debug message.
    
    Parameters
    ----------
    flag : str
        The flag to search for in the bytes.
    bytes : bytes
        The bytes to search within.

    Returns
    -------
    bytes
        The substring of bytes from the first occurrence of flag to the end of the line, or None if not found.
    """
    try:
        idx = bytes.index(flag) + len(flag)
    except ValueError:
        logger.debug(f'get_toeol: {flag} not found in {bytes}')
        return None
    if '\n' not in bytes[idx:]:
        logger.debug(f'get_toeol: no eol in {bytes}')
        return None
    eol = bytes[idx:].index('\n') + idx
    substr = bytes[idx:eol]
    return substr.strip()

def get_tokens(flag: str, bytes: str):
    """
    Get the tokens from the bytes that follow the first occurrence of flag.
    If the flag is not found or there are no tokens, return an empty list and log a debug message.

    Parameters
    ----------
    flag : str
        The flag to search for in the bytes.
    bytes : bytes
        The bytes to search within.

    Returns
    -------
    list
        A list of tokens found in the bytes, or an empty list if none are found.
    """
    # logger.debug(f'get_single: {flag} in {bytes}')
    toel = get_toeol(flag, bytes)
    if toel is None:
        return []
    return [x.strip() for x in toel.split()]

def get_single(flag: str, bytes: str):
    """
    Get the first token that follows the first occurrence of flag in the bytes.
    If the flag is not found or there are no tokens, return None and log a debug message.
    
    Parameters
    ----------
    flag : str
        The flag to search for in the bytes.
    bytes : bytes
        The bytes to search within.

    Returns
    -------
    bytes
        The first token found in the bytes, or None if not found.
    """
    # logger.debug(f'get_single: {flag} in {bytes}')
    vals = get_tokens(flag, bytes)
    if len(vals) == 0:
        return None
    return vals[0]

def get_values(flag: str, bytes: str, dtype=float):
    """
    Get the values from the bytes that follow the first occurrence of flag.
    If the flag is not found or there are no values, return an empty list and log a debug message.

    Parameters
    ----------
    flag : str
        The flag to search for in the bytes.
    bytes : bytes
        The bytes to search within.
    dtype : type
        The data type to convert the values to (default is float).

    Returns
    -------
    list
        A list of values found in the bytes, or an empty list if none are found.
    """
    # logger.debug(f'get_single: {flag} in {bytes}')
    vals = get_tokens(flag, bytes)
    if len(vals) == 0:
        return []
    rvals = []
    for v in vals:
        try:
            rvals.append(dtype(v))
        except ValueError:
            logger.debug(f'get_single: {flag} in {bytes} failed to convert {v} to {dtype}')
            return []
    return rvals

class LogParser(ByteCollector):
    """
    A base class for parsing log files from various applications used by Pestifer. This class is a subclass of :class:`ByteCollector <pestifer.util.stringthings.ByteCollector>` and provides methods for reading, updating, and dumping log data.
    """
    def __init__(self):
        super().__init__()
        self.progress = 0.0
        self.progress_bar = None

    def static(self, filename: str):
        """
        Initialize the LogParser from an existing, static file.
        """
        logger.debug(f'Initiating {self.__class__.__name__} from {filename}')
        with open(filename, 'r') as f:
            rawlines = f.read()
        self.update(rawlines)

    def update(self, bytes: str):
        """
        Update the LogParser with new bytes of data.
        """
        self.write(bytes)
        # logger.debug(f'Updating {self.__class__.__name__} with {len(bytes)} bytes -> {len(self.byte_collector)} bytes collected')

    def dump(self, basename: str = 'logparser'):
        """
        Dump the collected log data to a file with the specified basename.
        The file will be named `<basename>.log` and will contain the collected log data.
        
        Parameters
        ----------
        basename : str
            The base name for the log file. The final file will be named `<basename>.log`.
        """
        logger.debug(f'Dumping {self.__class__.__name__} to {basename}')
        self.write_file(f'{basename}.log')

    def measure_progress(self):
        """
        Measure the progress of the log parsing. This method is intended to be overridden by subclasses to provide specific progress measurement logic.
        """
        return self.progress

    def enable_progress_bar(self, PS: PestiferProgress = None):
        """
        Enable a progress bar for the log parser. If a progress bar instance is provided, it will be used to track the progress of the log parsing.
        
        Parameters
        ----------
        PS : PestiferProgress
            An instance of a progress bar to track the progress of the log parsing.
        """
        if PS is None:
            self.progress_bar = None
        else:
            self.progress_bar = PS
            self.progress_bar.register_update_function(self.measure_progress)

    def update_progress_bar(self):
        """
        Update the progress bar for the log parser.
        This method calls the `go` method of the progress bar instance, if it exists.
        """
        if self.progress_bar is not None:
            self.progress_bar.go()

    def success(self):
        """
        Check if the log parsing was successful. This method is intended to be overridden by subclasses to provide specific success criteria.
        
        Returns
        -------
        bool
            True if the log parsing was successful, False otherwise.
        """
        return False

    def follow(self, filename: str, sleep_interval: float = 0.5, timeout_intervals: int = 120):
        """
        Follow a log file, reading new lines as they are added.
        
        Parameters
        ----------
        filename : str
            The path to the log file to follow.
        sleep_interval : float
            The time to sleep between checks for new lines (in seconds).
        timeout_intervals : int
            The maximum number of intervals to wait before timing out.
        """
        logger.debug(f'Following {filename}')
        self.ntimeout = 0
        with open(filename, 'r') as f:
            while True:
                line = f.readline()
                if not line:
                    time.sleep(sleep_interval)
                    self.ntimeout += 1
                    if self.ntimeout > timeout_intervals:
                        logger.debug(f'Follow timed out after {self.ntimeout} intervals')
                        break
                    continue
                self.ntimeout = 0
                self.update(line)
                self.update_progress_bar()
                if self.success():
                    logger.debug(f'Follow succeeded after {self.ntimeout} intervals')
                    break

class VMDLogParser(LogParser):
    """
    A class for parsing VMD log files. This class is a subclass of :class:`LogParser <pestifer.logparsers.logparser.LogParser>` and provides methods for reading, updating, and dumping VMD log data.
    It also includes methods for processing specific lines in the log file, such as those containing information about the simulation, TCL commands, and energy calculations.
    
    Parameters
    ----------
    basename : str
        The base name for the log parser.
    """
    def __init__(self, basename: str = 'vmd-logparser'):
        super().__init__()
        self.basename = basename

    @classmethod
    def from_file(cls, filename: str | Path) -> VMDLogParser:
        """
        Create a VMDLogParser instance from a log file.
        """
        instance = cls()
        instance.static(filename)
        return instance

    def collect_validation_results(self):
        """
        Collect validation results from the log data.
        """
        results = []
        for line in self.byte_collector.splitlines():
            if "PASS" in line or "FAIL" in line:
                results.append(line)
        return results
