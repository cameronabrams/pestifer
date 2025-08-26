# Author Cameron F. Abrams, <cfa22@drexel.edu>

"""
psfgen log parsing utility
"""

import logging
import os
import re

from .logparser import LogParser, get_single

logger = logging.getLogger(__name__)

class PsfgenLogParser(LogParser):
    """
    A class for parsing Psfgen log files. This class is a subclass of :class:`LogParser <pestifer.logparsers.logparser.LogParser>` and provides methods for reading, updating, and dumping Psfgen log data.
    
    Parameters
    ----------
    basename : str
        The base name for the log parser. This is used to name the output log file.
    """
    info_keys = ['Info) ', 'Info: ']
    """
    A list of prefixes used to identify information lines in the Psfgen log file.
    """
    psfgen_key = 'psfgen) '
    """
    The prefix used to identify Psfgen-specific lines in the log file.
    """
    def __init__(self, basename='psfgen-logparser'):
        super().__init__()
        self.line_idx = [0]  # byte offsets of lines
        self.processed_line_idx = []
        self.metadata = {}
        self.basename = basename

    def update(self, bytes: str):
        """
        Update the Psfgen log parser with new bytes of data. This method appends the new bytes to the byte collector and processes the lines in the log file.
        
        Parameters
        ----------
        bytes : bytes
            The bytes to update the log parser with. This can be a string or bytes object containing the log data.
        """
        super().update(bytes)
        last_line_idx = self.line_idx[-1]
        addl_line_idx = [m.start()+1+last_line_idx for m in re.finditer(os.linesep,self.byte_collector[last_line_idx:])]
        scan_ldx = [last_line_idx]  # in case last line was incomplete in a previous pass
        if len(addl_line_idx) > 1:
            scan_ldx.extend(addl_line_idx)
        for i, j in zip(scan_ldx[:-1], scan_ldx[1:]):
            if i not in self.processed_line_idx:
                line = self.byte_collector[i:j]
                self.process_line(line)
                self.processed_line_idx.append(i)
        if len(addl_line_idx) > 1:
            self.line_idx.extend(addl_line_idx)
            last_line = self.byte_collector[addl_line_idx[-1]:]
            if last_line.endswith(os.linesep):
                self.process_line(last_line)
                self.processed_line_idx.append(addl_line_idx[-1])

    def process_line(self, line: str):
        """
        Process a line from the Psfgen log file. This method identifies the type of line (information, Psfgen command) and processes it accordingly.
        
        Parameters
        ----------
        line : str
            A line from the Psfgen log file to be processed. This can be a string or bytes object containing the log data.
        """
        assert line.endswith(os.linesep), f'process_line: {line} does not end with os.linesep'
        if line.startswith(self.info_keys[0]):
            o = len(self.info_keys[0])
            self.process_info_line(line[o:])
        elif line.startswith(self.info_keys[1]):
            o = len(self.info_keys[1])
            self.process_info_line(line[o:])
        elif line.startswith(self.psfgen_key):
            o = len(self.psfgen_key)
            self.process_psfgen_line(line[o:])

    def process_info_line(self, line: str):
        """ 
        Process an information line from the Psfgen log file. This method extracts metadata from the line and updates the `metadata` attribute accordingly.
        
        Parameters
        ----------
        line : str
            A line from the Psfgen log file that contains information to be processed. This can be a string or bytes object containing the log data.
        """
        if 'Exiting normally' in line:
            self.metadata['success'] = True
        elif line.startswith('VMD'):
            tokens = [x.strip() for x in line.split()]
            if len(tokens) > 1:
                if 'platform' not in self.metadata:
                    self.metadata['platform'] = tokens[2]
                if 'version' not in self.metadata:
                    self.metadata['version'] = tokens[4]

    def success(self):
        """
        Check if the Psfgen execution was successful. This method checks if the metadata contains the `success` key, which indicates that the execution was successful.
        
        Returns
        -------
        bool
            True if the log parsing was successful, False otherwise.
        """
        if 'success' in self.metadata:
            return self.metadata['success']
        return False

    def process_psfgen_line(self, line: str):
        """
        Process a Psfgen-specific line from the log file. This method extracts metadata from the line and updates the `metadata` attribute accordingly.
        
        Parameters
        ----------
        line : str
            A line from the Psfgen log file that contains Psfgen-specific information to be processed. This can be a string or bytes object containing the log data.
        """
        if line.startswith('total of'):
            if line.endswith('atoms'):
                self.metadata['number_of_atoms'] = int(get_single('total of', line))
            elif line.endswith('bonds'):
                self.metadata['number_of_bonds'] = int(get_single('total of', line))
            elif line.endswith('angles'):
                self.metadata['number_of_angles'] = int(get_single('total of', line))
            elif line.endswith('dihedrals'):
                self.metadata['number_of_dihedrals'] = int(get_single('total of', line))
            elif line.endswith('impropers'):
                self.metadata['number_of_impropers'] = int(get_single('total of', line))
            elif line.endswith('cross-terms'):
                self.metadata['number_of_cross_terms'] = int(get_single('total of', line))
            elif line.endswith('explicit exclusions'):
                self.metadata['number_of_exclusions'] = int(get_single('total of', line))
