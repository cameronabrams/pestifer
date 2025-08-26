# Author: Cameron F. Abrams, <cfa22@drexel.edu>

"""
PDB2PQR log parsing utility
"""

import logging
import os
import re

import pandas as pd

from .logparser import LogParser

logger = logging.getLogger(__name__)

class PDB2PQRLogParser(LogParser):
    """
    A class for parsing PDB2PQR log files. This class is a subclass of :class:`LogParser <pestifer.logparsers.logparser.LogParser>` and provides methods for reading, updating, and dumping PDB2PQR log data.
    
    Parameters
    ----------
    basename : str
        The base name for the log parser. This is used to name the output log file.
    """
    section_separator = '-'*80 # file is divided into sections
    """
    The separator used to identify sections in the PDB2PQR log file.
    """
    def __init__(self, basename='pdb2pqr-logparser'):
        super().__init__()
        self.processed_separator_idx = []
        self.processed_sections = []
        self.metadata = {}
        self.basename = basename
        self.progress = 0.0

    def update(self, bytes: str):
        """
        Update the PDB2PQR log parser with new bytes of data. This method appends the new bytes to the byte collector and processes the sections in the log file.
        
        Parameters
        ----------
        bytes : bytes
            The bytes to update the log parser with.
        """
        super().update(bytes)
        separator_idx = [0] + [m.start() for m in re.finditer(self.section_separator, self.byte_collector)]
        # logger.debug(f'update: found {len(separator_idx)-1} sections in {self.basename}')
        for i, j in zip(separator_idx[:-1], separator_idx[1:]):
            if i not in self.processed_separator_idx:
                self.processed_separator_idx.append(i)
                self.process_section(self.byte_collector[i+int(i>0)*len(self.section_separator):j])

    def finalize(self):
        """
        Finalize the PDB2PQR log parser by saving the metadata to a YAML file. This method writes the metadata to a YAML file with the base name of the log parser.
        """
        separator_idx=[0]+[m.start() for m in re.finditer(self.section_separator,self.byte_collector)]
        # logger.debug(f'update: found {len(separator_idx)-1} sections in {self.basename}')
        for i,j in zip(separator_idx[:-1],separator_idx[1:]):
            if i not in self.processed_separator_idx:
                self.processed_separator_idx.append(i)
                self.process_section(self.byte_collector[i+len(self.section_separator):j])

    def process_section(self, bytes: str):
        """
        Process a section of the PDB2PQR log file. This method identifies the type of section and processes it accordingly.
        """
        # logger.debug(f'process_section: {bytes[:50]}...')
        if 'SUMMARY OF THIS PREDICTION' in bytes:
            self.process_summary(bytes)

    def process_summary(self, bytes: str):
        """
        Process the summary section of the PDB2PQR log file. This method extracts information from the summary section and updates the `metadata` attribute with relevant data.
        
        Parameters
        ----------
        bytes : bytes
            The bytes representing the summary section to be processed. This can be a string or bytes object
            containing the summary data. 
        """
        # logger.debug(f'process_summary: {bytes[:50]}...')
        lines = bytes.split(os.linesep)
        expected_table = lines[3:]
        logger.debug(f'process_summary: expected_table begins with: {expected_table[0]}')
        table_lines = []
        for line in expected_table:
            if len(line.strip()) == 0:
                logger.debug('process_summary: empty line signals end of table')
                break
            tokens = [x.strip() for x in line.split() if x.strip()]
            resname = tokens[0]
            resnum = int(tokens[1])
            reschain = tokens[2]
            respka = float(tokens[3])
            resmodelpka = float(tokens[4])
            if len(tokens) > 5:
                resatomtype = tokens[5]
            else:
                resatomtype = None
            table_lines.append(dict(resname=resname, resnum=resnum, reschain=reschain, respka=respka, resmodelpka=resmodelpka, resatomtype=resatomtype))
        if len(table_lines) > 0:
            self.metadata['pka_table'] = pd.DataFrame(table_lines)
