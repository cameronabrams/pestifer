# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
Defines the :class:`LogParser` class and its subclasses for parsing log files from various applications used by Pestifer, such as NAMD and Packmol.
"""

import logging
import os
import re
import time
import yaml

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from enum import Enum

from ..core.stringthings import ByteCollector
from ..util.progress import PackmolProgress, NAMDProgress

logger=logging.getLogger(__name__)

def get_toflag(flag,bytes):
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
        idx=bytes.index(flag)
    except ValueError:
        logger.debug(f'get_toflag: {flag} not found in {bytes}')
        return None
    substr=bytes[:idx]
    return substr.strip()

def get_toeol(flag,bytes):
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
        idx=bytes.index(flag)+len(flag)
    except ValueError:
        logger.debug(f'get_toeol: {flag} not found in {bytes}')
        return None
    if '\n' not in bytes[idx:]:
        logger.debug(f'get_toeol: no eol in {bytes}')
        return None
    eol=bytes[idx:].index('\n')+idx
    substr=bytes[idx:eol]
    return substr.strip()

def get_tokens(flag,bytes):
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
    toel=get_toeol(flag,bytes)
    if toel is None:
        return []
    return [x.strip() for x in toel.split()]

def get_single(flag,bytes):
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
    vals=get_tokens(flag,bytes)
    if len(vals)==0:
        return None
    return vals[0]

def get_values(flag,bytes,dtype=float):
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
    vals=get_tokens(flag,bytes)
    if len(vals)==0:
        return []
    rvals=[]
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
        self.progress=0.0
        self.progress_bar=None

    def static(self,filename):
        """
        Initialize the LogParser from an existing, static file.
        """
        logger.debug(f'Initiating {self.__class__.__name__} from {filename}')
        with open(filename,'r') as f:
            rawlines=f.read()
        self.update(rawlines)

    def update(self,bytes):
        """
        Update the LogParser with new bytes of data.
        """
        self.write(bytes)
        # logger.debug(f'Updating {self.__class__.__name__} with {len(bytes)} bytes -> {len(self.byte_collector)} bytes collected')

    def dump(self,basename='logparser'):
        """
        Dump the collected log data to a file with the specified basename.
        The file will be named `<basename>.log` and will contain the collected log data.
        
        Parameters
        ----------
        basename : str
            The base name for the log file. The final file will be named `<basename>.log`.
        """
        logger.debug(f'Dumping {self.__class__.__name__} to {basename}')
        self.byte_collector.write_file(f'{basename}.log')

    def measure_progress(self):
        """
        Measure the progress of the log parsing. This method is intended to be overridden by subclasses to provide specific progress measurement logic.
        """
        return self.progress
    
    def enable_progress_bar(self,PS=None):
        """
        Enable a progress bar for the log parser. If a progress bar instance is provided, it will be used to track the progress of the log parsing.
        
        Parameters
        ----------
        PS : ProgressBar
            An instance of a progress bar to track the progress of the log parsing.
        """
        if PS is None:
            self.progress_bar=None
        else:
            self.progress_bar=PS
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

    def follow(self,filename,sleep_interval=0.5,timeout_intervals=120):
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
        self.ntimeout=0
        with open(filename,'r') as f:
            while True:
                line=f.readline()
                if not line:
                    time.sleep(sleep_interval)
                    self.ntimeout+=1
                    if self.ntimeout>timeout_intervals:
                        logger.debug(f'Follow timed out after {self.ntimeout} intervals')
                        break
                    continue
                self.ntimeout=0
                self.update(line)
                self.update_progress_bar()
                if self.success():
                    logger.debug(f'Follow succeeded after {self.ntimeout} intervals')
                    break

class VMDLog(LogParser):
    """
    A class for parsing VMD log files. This class is a subclass of :class:`LogParser <pestifer.util.logparsers.LogParser>` and provides methods for reading, updating, and dumping VMD log data.
    It also includes methods for processing specific lines in the log file, such as those containing information about the simulation, TCL commands, and energy calculations.
    
    Parameters
    ----------
    basename : str
        The base name for the log parser.
    """
    def __init__(self,basename='vmd-logparser'):
        super().__init__()
        self.basename=basename

class NAMDxst():
    """ 
    A class for parsing NAMD xst files, which contain information about the simulation cell dimensions.

    Parameters
    ----------
    filename : str
        The path to the NAMD xst file to parse.
    """
    def __init__(self,filename):
        celldf=pd.read_csv(filename,skiprows=2,header=None,sep=r'\s+',index_col=None)
        col='step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z s_x s_y s_z s_u s_v s_w'.split()[:len(celldf.columns)]
        celldf.columns=col
        self.df=celldf

    def add_file(self,filename):
        """
        Add another NAMD xst file to the existing data.
        
        Parameters
        ----------
        filename : str
            The path to the NAMD xst file to add.
        """
        celldf=pd.read_csv(filename,skiprows=2,header=None,sep=r'\s+',index_col=None)
        col='step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z s_x s_y s_z s_u s_v s_w'.split()[:len(celldf.columns)]
        celldf.columns=col
        self.df=pd.concat((self.df,celldf))

    def concat(self,other):
        """
        Concatenate another NAMDxst instance's data with the current instance's data.
        
        Parameters
        ----------
        other : NAMDxst
            Another instance of NAMDxst whose data will be concatenated with the current instance's data.
        """
        self.df=pd.concat((self.df,other.df))
        
class NAMDLog(LogParser):
    """
    A class for parsing NAMD log files. This class is a subclass of :class:`LogParser <pestifer.util.logparsers.LogParser>` and provides methods for reading, updating, and dumping NAMD log data.
    It also includes methods for processing specific lines in the log file, such as those containing information about the simulation, energy calculations, and pressure profiles.
    
    Parameters
    ----------
    basename : str
        The base name for the log parser. This is used to name the output log file.
    """

    info_key='Info: '
    """
    The key used to identify information lines in the NAMD log file.
    """
    tcl_key='TCL: '
    """
    The key used to identify TCL command lines in the NAMD log file.
    """
    energy_key='ENERGY: '
    """
    The key used to identify energy lines in the NAMD log file.
    """
    etitle_key='ETITLE: '
    """
    The key used to identify energy title lines in the NAMD log file.
    """
    pressureprofile_key='PRESSUREPROFILE: '
    """
    The key used to identify pressure profile lines in the NAMD log file.
    """
    struct_sep='*****************************'
    """
    The key used to identify the structure summary in the NAMD log file.
    """
    wallclock_key='WallClock: '
    """
    The key used to identify wall clock time lines in the NAMD log file.
    """
    restart_key='WRITING COORDINATES TO RESTART FILE AT STEP '
    """
    The key used to identify lines indicating that coordinates are being written to a restart file in the NAMD log file.
    """
    performance_key='PERFORMANCE: '
    """
    The key used to identify performance lines in the NAMD log file, e.g.,
    PERFORMANCE: 34000  averaging 23.789 ns/day, 0.00726387 sec/step with standard deviation 3.76697e-05
    """
    timing_key='TIMING: '
    """
    The key used to identify timing lines in the NAMD log file, e.g.,
    TIMING: 34000  CPU: 100.266, 0.00728793/step  Wall: 100.267, 0.00728689/step, 1.9962 hours remaining, 0.000000 MB of memory in use.
    """
    default_etitle_npt='TS,BOND,ANGLE,DIHED,IMPRP,ELECT,VDW,BOUNDARY,MISC,KINETIC,TOTAL,TEMP,POTENTIAL,TOTAL3,TEMPAVG,PRESSURE,GPRESSURE,VOLUME,PRESSAVG,GPRESSAVG'.split(',')
    """
    The default energy titles for NPT ensembles in NAMD log files.
    """
    default_etitle_nvt='TS,BOND,ANGLE,DIHED,IMPRP,ELECT,VDW,BOUNDARY,MISC,KINETIC,TOTAL,TEMP,POTENTIAL,TOTAL3,TEMPAVG'.split(',')  
    """
    The default energy titles for NVT ensembles in NAMD log files.
    """
    bail_out_key='Stack Traceback:'
    """
    The key used to identify lines indicating a stack traceback in the NAMD log file, which may indicate an error or crash.
    """

    def __init__(self,basename='namd-logparser'):
        super().__init__()
        self.line_idx=[0] # byte offsets of lines
        self.processed_line_idx=[]
        self.time_series_data={}
        self.metadata={}
        self.dataframes={}
        self.reading_structure_summary=False
        self.basename=basename
        self._line_processors={
            self.info_key: self.process_info_line,
            self.tcl_key: self.process_tcl_line,
            self.energy_key: self.process_energy_line,
            self.pressureprofile_key: self.process_pressureprofile_line,
            self.wallclock_key: self.process_wallclock_line,
            self.restart_key: self.process_restart_line,
            self.performance_key: self.process_performance_line,
            self.timing_key: self.process_timing_line
        }
    
    @classmethod
    def from_file(cls,filename,passfilter=[]):
        """
        Create a NAMDLog instance from an existing NAMD log file.
        
        Parameters
        ----------
        filename : str
            The path to the NAMD log file to read.
        
        Returns
        -------
        NAMDLog
            An instance of NAMDLog with the data from the specified file.
        """
        logger.debug(f'Creating {cls.__name__} from {filename}')
        instance=cls()
        instance.static(filename,passfilter=passfilter)
        return instance
    
    def static(self,filename,passfilter=[]):
        """
        Initialize the NAMDLog from an existing, static file.
        
        Parameters
        ----------
        filename : str
            The path to the NAMD log file to read.
        filter : list
            A list of strings to filter the lines in the log file. Only lines containing these strings will be processed.
        """
        logger.debug(f'Initiating {self.__class__.__name__} from {filename}')
        with open(filename,'r',encoding='utf-8') as f:
            raw=f.read()
        self.write(raw)
        rawlines=raw.splitlines(keepends=True)
        for line in rawlines:
            if not passfilter or any(f in line for f in passfilter):
                self.process_line(line)
        self.finalize()

    def process_struct_summ_datum(self,line):
        """
        Process a line from the structure summary section of the NAMD log file.
        """
        o=len(self.info_key)
        if line.endswith('FIXED ATOMS\n'):
            self.metadata['number_of_fixed_atoms']=int(get_toflag('FIXED ATOMS',line))
        elif line.endswith('ATOMS\n'):
            self.metadata['number_of_atoms']=int(get_toflag('ATOMS',line))
        elif line.endswith('BONDS\n') and not 'RIGID' in line:
            self.metadata['number_of_bonds']=int(get_toflag('BONDS',line))
        elif line.endswith('ANGLES\n'):
            self.metadata['number_of_angles']=int(get_toflag('ANGLES',line))
        elif line.endswith('DIHEDRALS\n'):
            self.metadata['number_of_dihedrals']=int(get_toflag('DIHEDRALS',line))
        elif line.endswith('IMPROPERS\n'):
            self.metadata['number_of_impropers']=int(get_toflag('IMPROPERS',line))
        elif line.endswith('CROSSTERMS\n'):
            self.metadata['number_of_crossterms']=int(get_toflag('CROSSTERMS',line))
        elif line.endswith('EXCLUSIONS\n'):
            self.metadata['number_of_exclusions']=int(get_toflag('EXCLUSIONS',line))
        elif line.endswith('RIGID BONDS\n'):
            self.metadata['number_of_rigid_bonds']=int(get_toflag('RIGID BONDS',line))
        elif line.endswith('ATOMS IN LARGEST HYDROGEN GROUP\n'):
            self.metadata['atoms_in_largest_hydrogen_group']=int(get_toflag('ATOMS IN LARGEST HYDROGEN GROUP',line))
        elif line.startswith('ATOM DENSITY ='):
            self.metadata['atom_density']=float(get_single('ATOM DENSITY =',line))
        elif 'TOTAL MASS =' in line:
            self.metadata['total_mass']=float(get_single('TOTAL MASS =',line))
        elif 'TOTAL CHARGE =' in line:
            self.metadata['total_charge']=float(get_single('TOTAL CHARGE =',line))
        elif 'MASS DENSITY =' in line:
            self.metadata['mass_density']=float(get_single('MASS DENSITY =',line))

    def process_info_line(self,line):
        """
        Process a line from the information section of the NAMD log file.
        """
        # logger.debug(f'process_info: {line}')
        if line.startswith('STRUCTURE SUMMARY:'):
            self.reading_structure_summary=True
        if line.startswith(self.struct_sep):
            self.reading_structure_summary=False
        if self.reading_structure_summary:
            self.process_struct_summ_datum(line)
            return
        if line.startswith(f'TIMESTEP'):
            self.metadata['timestep']=float(get_single('TIMESTEP',line))
        elif 'FIRST TIMESTEP' in line:
            self.metadata['first_timestep']=int(get_single('FIRST TIMESTEP',line))
        elif 'NUMBER OF STEPS' in line:
            self.metadata['number_of_steps']=int(get_single('NUMBER OF STEPS',line))
        elif 'RANDOM NUMBER SEED' in line:
            self.metadata['random_number_seed']=int(get_single('RANDOM NUMBER SEED',line))
        elif 'RESTART FILENAME' in line:
            self.metadata['restart_filename']=get_single('RESTART FILENAME',line)
        elif 'RESTART FREQUENCY' in line:
            self.metadata['restart_frequency']=int(get_single('RESTART FREQUENCY',line))
        elif 'OUTPUT FILENAME' in line:
            self.metadata['output_filename']=get_single('OUTPUT FILENAME',line)
        elif 'ENERGY OUTPUT STEPS' in line:
            self.metadata['energy_output_steps']=int(get_single('ENERGY OUTPUT STEPS',line))
        elif 'LANGEVIN DYNAMICS ACTIVE' in line:
            self.metadata['ensemble']='NVT'
        elif 'LANGEVIN PISTON PRESSURE CONTROL ACTIVE' in line:
            self.metadata['ensemble']='NPT'
        elif 'PERIODIC CELL BASIS 1' in line:
            self.metadata['periodic_cell_basis_1']=get_values('PERIODIC CELL BASIS 1',line,dtype=float)
        elif 'PERIODIC CELL BASIS 2' in line:
            self.metadata['periodic_cell_basis_2']=get_values('PERIODIC CELL BASIS 2',line,dtype=float)
        elif 'PERIODIC CELL BASIS 3' in line:
            self.metadata['periodic_cell_basis_3']=get_values('PERIODIC CELL BASIS 3',line,dtype=float)

    def process_tcl_line(self,line):
        """
        Process a line from the TCL command section of the NAMD log file.
        """
        if line.startswith('Running for'):
            self.metadata['running_for']=int(get_single('Running for',line))
        elif line.startswith('Minimizing for'):
            self.metadata['minimizing_for']=int(get_single('Minimizing for',line))
            self.metadata['ensemble']='minimize'

    def process_energy_line(self,line):
        """
        Process a line from the energy section of the NAMD log file.

        Parameters
        ----------
        line : str
            A line from the NAMD log file that contains energy data.
        """
        tokens=[x.strip() for x in line.split()]
        if 'etitle' not in self.metadata:
            if len(tokens)==len(self.default_etitle_npt):
                logger.debug(f'process_energy: {len(tokens)} tokens found, using default etitle for NPT')
                self.metadata['etitle']=self.default_etitle_npt
            elif len(tokens)==len(self.default_etitle_nvt):
                logger.debug(f'process_energy: {len(tokens)} tokens found, using default etitle for NVT')
                self.metadata['etitle']=self.default_etitle_nvt
            else:
                logger.debug(f'process_energy: {len(tokens)} tokens found, but either {len(self.default_etitle_nvt)} or {len(self.default_etitle_npt)} expected')
                return
        else:
            if len(tokens)!=len(self.metadata['etitle']):
                logger.debug(f'process_energy: {len(tokens)} tokens found, but {len(self.metadata["etitle"])} expected')
                return
        tokens[0]=int(tokens[0])
        for i in range(1,len(tokens)):
            tokens[i]=float(tokens[i])
        new_line={k:v for k,v in zip(self.metadata['etitle'],tokens)}
        if 'VOLUME' in new_line:
            new_line['DENSITY']=self.metadata['total_mass']/new_line['VOLUME']
        if not 'energy' in self.time_series_data:
            self.time_series_data['energy']=[]
        self.time_series_data['energy'].append(new_line)
        # if self.energy_df.empty:
        #     self.energy_df=pd.DataFrame(new_line,index=[0])
        # else:
        #     self.energy_df=pd.concat([self.energy_df,pd.DataFrame(new_line,index=[0])],ignore_index=True)

    def process_energy_title(self,line):
        """
        Process a line from the energy title section of the NAMD log file.

        Parameters
        ----------
        line : str
            A line from the NAMD log file that contains the energy titles.
        """
        tokens=[x.strip() for x in line.split()]
        if 'etitle' not in self.metadata:
            self.metadata['etitle']=tokens
        else:
            if len(tokens)!=len(self.metadata['etitle']):
                logger.debug(f'process_energy_title: {len(tokens)} tokens found, but {len(self.metadata["etitle"])} expected')

    def process_restart_line(self,line):
        """
        Process a line from the restart section of the NAMD log file.

        Parameters
        ----------
        line : str
            A line from the NAMD log file that contains restart information.
        """
        # the only thing in line should be the time step at which the last restart was written
        if not 'restart' in self.time_series_data:
            self.time_series_data['restart']=[]
        ts=int(line.strip())
        self.time_series_data['restart'].append(ts)

    def process_pressureprofile_line(self,line):
        """
        Process a line from the pressure profile section of the NAMD log file.

        Parameters
        ----------
        line : str
            A line from the NAMD log file that contains pressure profile data.
        """
        tokens=[x.strip() for x in line.split()]
        TS=int(tokens[0])
        logger.debug(f'process_pressureprofile_line: TS {tokens[0]}')
        for i in range(1,len(tokens)):
            tokens[i]=float(tokens[i])
        this_col=tokens[1:]
        new_line={'TS':TS}
        new_line.update({k:v for k,v in zip([str(x) for x in np.arange(len(this_col))],this_col)})
        if not 'number_of_pressure_slabs' in self.metadata:
            self.metadata['number_of_pressure_slabs']=len(this_col)
        if not 'pressureprofile' in self.time_series_data:
            self.time_series_data['pressureprofile']=[]
        self.time_series_data['pressureprofile'].append(new_line)
        # if self.pressureprofile_df.empty:
        #     # we need to set the indices to be the slab indices
        #     self.pressureprofile_df=pd.DataFrame(this_col,index=np.arange(self.metadata['number_of_pressure_slabs']),columns=[str(tokens[0])])
        # else:
        #     # append this col to the dataframe
        #     new_col=pd.DataFrame(this_col,index=np.arange(self.metadata['number_of_pressure_slabs']),columns=[str(tokens[0])])
        #     self.pressureprofile_df=pd.concat([self.pressureprofile_df,new_col],axis=1)
        #         # if self.pressureprofile_df.shape[1]==2:
        #             logger.debug(f'{self.pressureprofile_df.head().to_string()}')

    def process_wallclock_line(self,line):
        """
        Process a line from the wall clock time section of the NAMD log file.

        Parameters
        ----------
        line : str
            A line from the NAMD log file that contains the wall clock time.
        """
        tokens=[x.strip() for x in line.split()]
        self.metadata['wallclock_time']=float(tokens[0])

    def process_performance_line(self,line):
        """
        Process a line from the performance section of the NAMD log file, e.g.,
        PERFORMANCE: 34000  averaging 23.789 ns/day, 0.00726387 sec/step with standard deviation 3.76697e-05

        Parameters
        ----------
        line : str
            A line from the NAMD log file that contains performance data.
        """
        tokens=[x.strip() for x in line.split()]
        if len(tokens)<5:
            logger.debug(f'process_performance_line: {line} does not have enough tokens')
            return
        if not 'performance' in self.time_series_data:
            self.time_series_data['performance']=[]
        self.time_series_data['performance'].append({
            'steps': int(tokens[0]),
            'ns_per_day': float(tokens[2]),
            'sec_per_step': float(tokens[4]),
            'std_dev': float(tokens[-1]),
        })

    def process_timing_line(self,line):
        """
        Process a line from the timing section of the NAMD log file, e.g.,
        TIMING: 34000  CPU: 100.266, 0.00728793/step  Wall: 100.267, 0.00728689/step, 1.9962 hours remaining, 0.000000 MB of memory in use.

        Parameters
        ----------
        line : str
            A line from the NAMD log file that contains timing data.
        """
        tokens=[x.strip(' :,').replace('/step','').replace(',','') for x in line.split()]
        if len(tokens)<11:
            logger.debug(f'process_timing_line: {line} does not have enough tokens')
            return
        if not 'timing' in self.time_series_data:
            self.time_series_data['timing']=[]
        self.time_series_data['timing'].append({
            'steps': int(tokens[0]),
            'cpu_time': float(tokens[2]),
            'cpu_per_step': float(tokens[3]),
            'wall_time': float(tokens[5]),
            'wall_per_step': float(tokens[6]),
            'hours_remaining': float(tokens[7]),
            'memory_in_use': float(tokens[10]),
        })

    def process_line(self,line):
        """
        Process a line from the NAMD log file. This method identifies the type of line (information, TCL command, energy, pressure profile, or wall clock time) and processes it accordingly.
        
        Parameters
        ----------
        line : str
            A line from the NAMD log file to be processed.
        """
        assert line.endswith(os.linesep),f'process_line: {line} does not end with os.linesep'
        if self.bail_out_key in line:
            logger.debug(f'process_line: {line} contains bail out key, stopping processing')
            return -1
        if line.startswith(self.etitle_key):
            if 'etitle' in self.metadata:
                # logger.debug(f'process_energy_title: {line} already has etitle')
                return 0
            self.process_energy_title(line[len(self.etitle_key):])
            return 0
        for key,func in self._line_processors.items():
            if line.startswith(key):
                func(line[len(key):])
                return 0
    
    def update(self,bytes):
        """
        Update the NAMD log parser with new bytes of data. This method appends the new bytes to the byte collector and processes the lines in the log file. It identifies the end of each line and processes each line based on its content.  This is best used on a log file that is being written to, such as a live NAMD simulation log file.  For static files, use the :meth:`static <pestifer.util.logparsers.NAMDLog.static>` method instead.

        Parameters
        ----------
        bytes : bytes
            The bytes to update the log parser with. This can be a string or bytes object containing the log data.
        """
        super().update(bytes) # this just appends the bytes to the byte_collector
        last_line_idx=self.line_idx[-1] # recall byte index of last line processed
        addl_line_idx=[m.start()+1+last_line_idx for m in re.finditer(os.linesep,self.byte_collector[last_line_idx:])]
        scan_ldx=[last_line_idx] # in case last line was incomplete in a previous pass
        if len(addl_line_idx)>1:
            scan_ldx.extend(addl_line_idx)
        for i,j in zip(scan_ldx[:-1],scan_ldx[1:]):
            if i not in self.processed_line_idx:
                line=self.byte_collector[i:j]
                result=self.process_line(line)
                if result==-1:  # bail out key found, stop processing
                    return
                self.processed_line_idx.append(i)
        if len(addl_line_idx)>1:
            self.line_idx.extend(addl_line_idx)
            last_line=self.byte_collector[addl_line_idx[-1]:]
            if last_line.endswith(os.linesep):
                self.process_line(last_line)
                self.processed_line_idx.append(addl_line_idx[-1])

    def measure_progress(self):
        if 'number_of_steps' not in self.metadata:
            # logger.debug('measure_progress: number_of_steps not in metadata')
            return 0.0
        number_of_steps=self.metadata['number_of_steps']
        if 'first_timestep' not in self.metadata:
            # logger.debug('measure_progress: first_timestep not in metadata')
            return 0.0
        first_time_step=self.metadata['first_timestep']
        # if self.energy_df.empty:
        #     # logger.debug('measure_progress: energy_df is empty')
        #     return 0.0
        if 'running_for' in self.metadata:
            running_for=self.metadata['running_for']
            if running_for>0:
                number_of_steps=running_for
        elif 'minimizing_for' in self.metadata:
            minimizing_for=self.metadata['minimizing_for']
            if minimizing_for>0:
                number_of_steps=minimizing_for
        last_row=self.time_series_data['energy'][-1] if 'energy' in self.time_series_data and len(self.time_series_data['energy'])>0 else None
        if last_row is None:
            logger.debug('measure_progress: last_row is None')
            return 0.0
        if 'TS' not in last_row:
            logger.debug('measure_progress: TS not in last_row')
            return 0.0
        most_recent_time_step=last_row['TS']
        complete_steps=most_recent_time_step-first_time_step
        return complete_steps/number_of_steps
    
    def success(self):
        """ 
        Check if the NAMD log parsing was successful. This method checks if the metadata contains the ``wallclock_time`` key, which indicates that the log file has been processed successfully. 
        """
        return 'wallclock_time' in self.metadata

    def finalize(self):
        """
        Finalize the log parsing by creating dataframes for each time series.
        """
        for key in self.time_series_data:
            self.dataframes[key]=pd.DataFrame(self.time_series_data[key])
        if 'first_timestep' not in self.metadata:
            if 'energy' in self.dataframes:
                self.metadata['first_timestep']=int(self.dataframes['energy'].iloc[0]['TS'])
    
    def write_csv(self):
        """
        Write the parsed data to CSV files. This method creates a CSV file for each dataframe in the `dataframes` attribute, using the basename provided during initialization.
        The files will be named `<basename>-<key>.csv`, where `<key>` is the key of the dataframe in the `dataframes` dictionary.
        """
        for key in self.dataframes:
            self.dataframes[key].to_csv(f'{self.basename}-{key}.csv',index=False)

class PackmolLog(LogParser):
    """
    A class for parsing Packmol log files. This class is a subclass of :class:`LogParser <pestifer.util.logparsers.LogParser>` and provides methods for reading, updating, and dumping Packmol log data.

    Parameters
    ----------
    basename : str
        The base name for the log parser. This is used to name the output log file.
    """
    banner_separator='#'*80  # banners can occur within sections
    """
    The separator used to identify banners in the Packmol log file.
    """
    section_separator='-'*80 # file is divided into sections
    """
    The separator used to identify sections in the Packmol log file.
    """
    def __init__(self,basename='packmol-logparser'):
        super().__init__()
        self.processed_separator_idx=[]
        self.processed_sections=[]
        self.processed_banners=[]
        self.metadata={}
        self.gencan={}
        self.molmoves=[]
        self.basename=basename
        self.progress=0.0

    def update(self,bytes):
        """
        Update the Packmol log parser with new bytes of data. This method appends the new bytes to the byte collector and processes the sections in the log file.
        """
        super().update(bytes)
        separator_idx=[0]+[m.start() for m in re.finditer(self.section_separator,self.byte_collector)]
        for i,j in zip(separator_idx[:-1],separator_idx[1:]):
            if i not in self.processed_separator_idx:
                self.processed_separator_idx.append(i)
                self.process_section(self.byte_collector[i:j])

    def measure_progress(self):
        """
        Measure the progress of the Packmol log parsing. This method calculates the progress based on the number of GENCAN loops completed and the maximum number of GENCAN loops.
        It updates the `progress` attribute with the calculated progress value.
        """
        self.progress=self.metadata['current_total_gencan_loops']/self.metadata['max_total_gencan_loops'] if 'current_total_gencan_loops' in self.metadata and 'max_total_gencan_loops' in self.metadata else 0.0
        return super().measure_progress()
    
    def process_banner(self,bytes):
        """
        Process a banner in the Packmol log file. This method identifies the type of banner and updates the metadata accordingly.
        
        Parameters
        ----------
        bytes : bytes
            The bytes representing the banner to be processed. This can be a string or bytes object containing the banner data.
        """
        if 'Version' in bytes:
            vidx=bytes.index('Version')+len('Version')+1
            vstr=bytes[vidx:].split()[0]
            self.metadata['version']=vstr
            self.phase='initialization'
        elif 'Building initial approximation' in bytes:
            self.metadata['initial_approximation']={}
            self.phase='initial_approximation'
        elif 'Objective function at initial point' in bytes:
            idx=bytes.index('Objective function at initial point')+len('Objective function at initial point')+1
            ofstr=bytes[idx:].split()[0]
            self.metadata['initial_objective_function']=float(ofstr)
            self.phase='packing_molecules'
        elif 'Packing all molecules together' in bytes:
            self.phase='packing_all_molecules_together'
        elif 'Packing molecules of type' in bytes:
            idx=bytes.index('Packing molecules of type')+len('Packing molecules of type')+1
            mtype=bytes[idx:].split()[0]
            self.phase=f'packing_molecules_of_type_{mtype}'
            if not 'molecule_types_packed' in self.metadata:
                self.metadata['molecule_types_packed']=[]
            self.metadata['molecule_types_packed'].append(mtype)
        self.processed_banners.append(bytes)

    def process_header(self,bytes):
        """
        Process the header section of the Packmol log file. This method extracts metadata from the header and updates the `metadata` attribute with relevant information.
        
        Parameters
        ----------
        bytes : bytes
            The bytes representing the header section to be processed. This can be a string or bytes object
            containing the header data."""
        self.metadata['url']=get_single('Userguide at:',bytes)
        flag='Types of coordinate files specified:'
        self.metadata['coordinate_file_types']=bytes[bytes.index(flag)+len(flag):].split()[0]
        flag='Periodic boundary condition activated:'
        self.metadata['pbc_boxsize']=[float(bytes[bytes.index(flag)+len(flag):].split()[x].strip()) for x in range(0,3)]
        flag='PBC Reference box:'
        self.metadata['pbc_reference_box']=[float(bytes[bytes.index(flag)+len(flag):].split()[x].strip()) for x in range(0,6)]
        flag='Seed for random number generator:'
        self.metadata['seed']=int(bytes[bytes.index(flag)+len(flag):].split()[0].strip())
        flag='Output file:'
        self.metadata['output_file']=bytes[bytes.index(flag)+len(flag):].split()[0].strip()
        cf_idx=[m.start() for m in re.finditer('Reading coordinate file:',bytes)]
        self.metadata['coordinate_files']=[]
        for i in cf_idx:
            idx=i+len('Reading coordinate file:')
            self.metadata['coordinate_files'].append(bytes[idx:].split()[0].strip())
        flag='Number of independent structures:'
        self.metadata['number_of_independent_structures']=int(bytes[bytes.index(flag)+len(flag):].split()[0].strip())
        str_idx=[m.start() for m in re.finditer('Structure',bytes)]
        # logger.debug(f'Found {len(str_idx)} structures')
        self.metadata['structures']=[]
        for i in str_idx:
            idx=i+len('Structure')
            eol=bytes[idx:].index('\n')+idx
            substr=bytes[idx:eol]
            # logger.debug(f'Found structure: {substr} {idx} {eol}')
            tokens=substr.replace('(','').replace(')','').replace(':','').split()
            idx=int(tokens[0])
            fl=tokens[1]
            na=int(tokens[2])
            self.metadata['structures'].append(dict(idx=idx,file=fl,natoms=na))
        flag='Maximum number of GENCAN loops for all molecule packing:'
        self.metadata['max_gencan_loops_all']=int(bytes[bytes.index(flag)+len(flag):].split()[0].strip())
        flag='Maximum number of GENCAN loops for type:'
        mgl_idx=[m.start() for m in re.finditer(flag,bytes)]
        for i in mgl_idx:
            idx=i+len(flag)
            eol=bytes[idx:].index('\n')+idx
            substr=bytes[idx:eol]
            tokens=substr.replace(':','').split()
            mtype=int(tokens[0])
            mloops=int(tokens[1])
            match=next((d for d in self.metadata['structures'] if d['idx'] == mtype), None)
            if match:
                match['max_gencan_loops']=mloops
                match['current_gencan_loops']=0
            else:
                logger.error(f'process_header: No match for mtype {mtype} in structures {self.metadata["structures"]}')
        self.metadata['max_total_gencan_loops']=self.metadata['max_gencan_loops_all']+sum([d['max_gencan_loops'] for d in self.metadata['structures'] if 'max_gencan_loops' in d])
        self.metadata['current_total_gencan_loops']=0
        flag='Distance tolerance:'
        self.metadata['distance_tolerance']=float(bytes[bytes.index(flag)+len(flag):].split()[0].strip())
        flag='Number of molecules of type'
        nm_idx=[m.start() for m in re.finditer(flag,bytes)]
        for i in nm_idx:
            idx=i+len(flag)
            eol=bytes[idx:].index('\n')+idx
            substr=bytes[idx:eol]
            tokens=substr.replace(':','').split()
            mtype=int(tokens[0])
            n=int(tokens[1])
            match=next((d for d in self.metadata['structures'] if d['idx'] == mtype), None)
            if match:
                match['number_of_molecules']=n
        flag='Total number of atoms:'
        self.metadata['total_number_atoms']=int(bytes[bytes.index(flag)+len(flag):].split()[0].strip())
        flag='Total number of molecules:'
        self.metadata['total_number_molecules']=int(bytes[bytes.index(flag)+len(flag):].split()[0].strip())
        flag='Number of fixed molecules:'
        self.metadata['number_of_fixed_molecules']=int(bytes[bytes.index(flag)+len(flag):].split()[0].strip())
        flag='Number of free molecules:'
        self.metadata['number_of_free_molecules']=int(bytes[bytes.index(flag)+len(flag):].split()[0].strip())
        flag='Number of variables:'
        self.metadata['number_of_variables']=int(bytes[bytes.index(flag)+len(flag):].split()[0].strip())
        flag='Total number of fixed atoms:'
        self.metadata['total_number_fixed_atoms']=int(bytes[bytes.index(flag)+len(flag):].split()[0].strip())
        flag='Maximum internal distance of type'
        mid_idx=[m.start() for m in re.finditer(flag,bytes)]
        for i in mid_idx:
            idx=i+len(flag)
            eol=bytes[idx:].index('\n')+idx
            substr=bytes[idx:eol]
            tokens=substr.replace(':','').split()
            mtype=int(tokens[0])
            mid=float(tokens[1])
            match=next((d for d in self.metadata['structures'] if d['idx'] == mtype), None)
            if match:
                match['maximum_internal_distance']=mid

    def process_molecule_report(self,bytes):
        """
        Process a report about molecules in the Packmol log file. This method extracts information about the molecules, such as their type, function values before and after moving, and restraint violations.
        
        Parameters
        ----------
        bytes : bytes
            The bytes representing the molecule report to be processed. This can be a string or bytes object
            containing the report data.
        """
        flag='Molecules of type'
        idx=bytes.index(flag)+len(flag)
        eol=bytes[idx:].index('\n')+idx
        substr=bytes[idx:eol]
        tokens=substr.replace(':','').split()
        mtype=int(tokens[0])
        match=next((d for d in self.metadata['structures'] if d['idx'] == mtype), None)
        subflag='before_adjusting'
        if 'Adjusting random positions to fit the constraints.' in bytes:
            subflag='after_adjusting'
        # logger.debug(f'subflag: {subflag}')
        if match:
            match[subflag]={}
            flag='Function value before moving molecules:'
            match[subflag]['function_value_before_moving_molecules']=[]
            fbm_idx=[m.start() for m in re.finditer(flag,bytes)]
            for i in fbm_idx:
                idx=i+len(flag)
                eol=bytes[idx:].index('\n')+idx
                substr=bytes[idx:eol]
                tokens=substr.replace(':','').split()
                val=float(tokens[0])
                match[subflag]['function_value_before_moving_molecules'].append(val)
            flag='Function value after moving molecules:'
            match[subflag]['function_value_after_moving_molecules']=[]
            fam_idx=[m.start() for m in re.finditer(flag,bytes)]
            for i in fam_idx:
                idx=i+len(flag)
                eol=bytes[idx:].index('\n')+idx
                substr=bytes[idx:eol]
                tokens=substr.replace(':','').split()
                val=float(tokens[0])
                match[subflag]['function_value_after_moving_molecules'].append(val)
            flag='Restraint-only function value:'
            match[subflag]['restraint_only_function_value']=float(bytes[bytes.index(flag)+len(flag):].split()[0].strip())
            flag='Maximum violation of the restraints:'
            match[subflag]['maximum_violation_of_the_restraints']=float(bytes[bytes.index(flag)+len(flag):].split()[0].strip())

    def process_gencan_report(self,bytes):
        """
        Process a report about the GENCAN (GENeric CANonical) packing algorithm in the Packmol log file. This method extracts information about the GENCAN loops, such as the iteration number, radius scaling, function values, and violations of target distances and constraints.
        It updates the ``gencan`` attribute with the extracted information.
        
        Parameters
        ----------
        bytes : bytes
            The bytes representing the GENCAN report to be processed. This can be a string or bytes object
            containing the report data."""
        G=None
        ptok=self.phase.split('_')
        if ptok[1]=='all':
            mtype='all'
        else:
            mtype=int(ptok[4])
        # logger.debug(f'process_gencan_report: {mtype}')
        if mtype not in self.gencan:
            self.gencan[mtype]=[]
        G=self.gencan[mtype]
        if G is None:
            return
        iternum=int(get_single('Starting GENCAN loop:',bytes))
        radscal=float(get_single('Scaling radii by:',bytes))
        fval_last=float(get_single('Function value from last GENCAN loop: f =',bytes))
        fval_best=float(get_single('Best function value before: f =',bytes))
        max_viol_dist=float(get_single('Maximum violation of target distance:',bytes))
        max_viol_constr=float(get_single('Maximum violation of the constraints:',bytes))
        dsofar=dict(iteration=iternum,radius_scaling=radscal,function_value_last=fval_last,function_value_best=fval_best,max_viol_dist=max_viol_dist,max_viol_constr=max_viol_constr)
        if 'All-type function value:' in bytes:
            all_type_f=float(get_single('All-type function value:',bytes))
            dsofar['all_type_function_value']=all_type_f
        G.append(dsofar)
        self.metadata['current_total_gencan_loops']+=1
        match=next((d for d in self.metadata['structures'] if d['idx'] == mtype), None)
        if match:
            match['current_gencan_loops']+=1

    def process_gencan_success(self,bytes):
        """
        Process a report about the successful completion of the GENCAN packing algorithm in the Packmol log file. This method extracts information about the molecule type, objective function value, maximum violation of target distance, and maximum constraint violation.
        It updates the metadata of the corresponding molecule type in the ``metadata['structures']`` attribute
        with the extracted information.
        
        Parameters
        ----------
        bytes : bytes
            The bytes representing the GENCAN success report to be processed. This can be a string or bytes object
            containing the report data. 
        """
        mtype=int(get_single('Packing solved for molecules of type',bytes))
        obj_func_val=float(get_single('Objective function value:',bytes))
        max_viol_target_dist=float(get_single('Maximum violation of target distance:',bytes))
        max_viol_constr=float(get_single('Max. constraint violation:',bytes))
        match=next((d for d in self.metadata['structures'] if d['idx'] == mtype), None)
        if match:
            match['gencan_success']=dict(objective_function_value=obj_func_val,max_viol_target_distance=max_viol_target_dist,max_viol_constr=max_viol_constr)
            unneeded_loops=match.get('max_gencan_loops',0)-match.get('current_gencan_loops',0)
            self.metadata['max_total_gencan_loops']-=(unneeded_loops-1)

    def process_moving_worst_molecules(self,bytes):
        """
        Process a report about moving the worst molecules in the Packmol log file. This method extracts information about the function values before and after moving molecules, the types of molecules moved, and their percentages.
        
        Parameters
        ----------
        bytes : bytes
            The bytes representing the report about moving worst molecules to be processed. This can be a string or bytes object
            containing the report data.
        """
        top_idx=[m.start() for m in re.finditer('Moving',bytes)]
        cont_idx=top_idx[0]+len('Moving worst molecules ...')
        mbytes=bytes[cont_idx:]
        fvbefore=float(get_single('Function value before moving molecules:',mbytes))
        flag='Moving '
        moltypes=[]
        pct=[]
        fam_idx=[m.start() for m in re.finditer(flag,mbytes)]
        for i in fam_idx:
            idx=i+len(flag)
            eol=mbytes[idx:].index('\n')+idx
            substr=mbytes[idx:eol]
            tokens=substr.replace(':','').split()
            val=int(tokens[0])
            cnt=int(tokens[-1])
            moltypes.append((val,cnt))
        flag='Type'
        fam_idx=[m.start() for m in re.finditer(flag,mbytes)]
        for i in fam_idx:
            idx=i+len(flag)
            eol=mbytes[idx:].index('\n')+idx
            substr=mbytes[idx:eol]
            tokens=substr.replace(':','').split()
            val=int(tokens[0])
            cnt=float(tokens[-1].replace('%',''))
            pct.append((val,cnt))
        fvafter=float(get_single('Function value after moving molecules:',mbytes))
        result=dict(function_value_before_moving_molecules=fvbefore,
                    function_value_after_moving_molecules=fvafter,
                    moltypes=moltypes,
                    pcts=pct)
        self.molmoves.append(result)

    def process_section(self,bytes):
        """
        Process a section of the Packmol log file. This method identifies the type of section and processes it accordingly.
        
        Parameters
        ----------
        bytes : bytes
            The bytes representing the section to be processed. This can be a string or bytes object containing
            the section data.
        """
        banner_idx=[m.start() for m in re.finditer(self.banner_separator,bytes)]
        nbsep=len(banner_idx)
        banners=[]
        if nbsep>1:
            for i,j in zip(banner_idx[:-1:2],banner_idx[1::2]):
                banners.append(bytes[i:j])
            if nbsep%2==1:
                banners.append(bytes[banner_idx[-2]:banner_idx[-1]])
        for b in banners:
            self.process_banner(b)

        if 'Reading input file...' in bytes:
            self.process_header(bytes)
        elif 'Molecules of type:' in bytes:
            self.process_molecule_report(bytes)
        elif 'Starting GENCAN loop:' in bytes:
            self.process_gencan_report(bytes)
        elif 'Packing solved for molecules of type' in bytes:
            self.process_gencan_success(bytes)
        elif 'Moving worst molecules' in bytes:
            self.process_moving_worst_molecules(bytes)

    def finalize(self):
        """
        Finalize the Packmol log parser by saving the metadata and GENCAN data to files. This method writes the metadata to a YAML file and the GENCAN data to CSV files.
        It also generates plots of the GENCAN function values over iterations for each molecule type and saves them as PNG files.
        
        Returns
        -------
        str
            The filename of the generated PNG file containing the GENCAN function value plots.  
        """
        with open(f'{self.basename}_packmol-results.yaml','w') as f:
            yaml.dump(self.metadata,f,default_flow_style=False)
        self.gencan_df={}
        for k,v in self.gencan.items():
            self.gencan_df[k]=pd.DataFrame(v)
            self.gencan_df[k].to_csv(f'{self.basename}_{k}_packmol.csv',index=False)
        fig,ax=plt.subplots(1,len(self.gencan_df),figsize=(4*len(self.gencan_df),4))
        for i,(k,v) in enumerate(self.gencan_df.items()):
            if len(v)==0:
                continue
            ax[i].plot(v['iteration'],v['function_value_last'],label='Last function value')
            ax[i].plot(v['iteration'],v['function_value_best'],label='Best function value')
            ax[i].set_title(f'Molecule type {k}')
            ax[i].set_xlabel('Iteration')
            ax[i].set_ylabel('Function value')
            ax[i].legend()
            ax[i].set_yscale('log')
            ax[i].grid(True)
        plt.tight_layout()
        plt.savefig(f'{self.basename}_packmol.png')
        plt.close()
        return f'{self.basename}_packmol.png'

class PsfgenLog(LogParser):
    """
    A class for parsing Psfgen log files. This class is a subclass of :class:`LogParser <pestifer.util.logparsers.LogParser>` and provides methods for reading, updating, and dumping Psfgen log data.
    
    Parameters
    ----------
    basename : str
        The base name for the log parser. This is used to name the output log file.
    """
    info_keys=['Info) ', 'Info: ']
    """
    A list of prefixes used to identify information lines in the Psfgen log file.
    """
    psfgen_key='psfgen) '
    """
    The prefix used to identify Psfgen-specific lines in the log file.
    """
    def __init__(self,basename='psfgen-logparser'):
        super().__init__()
        self.line_idx=[0] # byte offsets of lines
        self.processed_line_idx=[]
        self.metadata={}
        self.basename=basename

    def update(self,bytes):
        """
        Update the Psfgen log parser with new bytes of data. This method appends the new bytes to the byte collector and processes the lines in the log file.
        
        Parameters
        ----------
        bytes : bytes
            The bytes to update the log parser with. This can be a string or bytes object containing the log data.
        """
        super().update(bytes)
        last_line_idx=self.line_idx[-1]
        addl_line_idx=[m.start()+1+last_line_idx for m in re.finditer(os.linesep,self.byte_collector[last_line_idx:])]
        scan_ldx=[last_line_idx] # in case last line was incomplete in a previous pass
        if len(addl_line_idx)>1:
            scan_ldx.extend(addl_line_idx)
        for i,j in zip(scan_ldx[:-1],scan_ldx[1:]):
            if i not in self.processed_line_idx:
                line=self.byte_collector[i:j]
                self.process_line(line)
                self.processed_line_idx.append(i)
        if len(addl_line_idx)>1:
            self.line_idx.extend(addl_line_idx)
            last_line=self.byte_collector[addl_line_idx[-1]:]
            if last_line.endswith(os.linesep):
                self.process_line(last_line)
                self.processed_line_idx.append(addl_line_idx[-1])

    def process_line(self,line):
        """
        Process a line from the Psfgen log file. This method identifies the type of line (information, Psfgen command) and processes it accordingly.
        
        Parameters
        ----------
        line : str
            A line from the Psfgen log file to be processed. This can be a string or bytes object containing the log data.
        """
        assert line.endswith(os.linesep),f'process_line: {line} does not end with os.linesep'
        if line.startswith(self.info_keys[0]):
            o=len(self.info_keys[0])
            self.process_info_line(line[o:])
        elif line.startswith(self.info_keys[1]):
            o=len(self.info_keys[1])
            self.process_info_line(line[o:])
        elif line.startswith(self.psfgen_key):
            o=len(self.psfgen_key)
            self.process_psfgen_line(line[o:])

    def process_info_line(self,line):
        """ 
        Process an information line from the Psfgen log file. This method extracts metadata from the line and updates the `metadata` attribute accordingly.
        
        Parameters
        ----------
        line : str
            A line from the Psfgen log file that contains information to be processed. This can be a string or bytes object containing the log data.
        """
        if 'Exiting normally' in line:
            self.metadata['success']=True
        elif line.startswith('VMD'):
            tokens=[x.strip() for x in line.split()]
            if len(tokens)>1:
                if 'platform' not in self.metadata:
                    self.metadata['platform']=tokens[2]
                if 'version' not in self.metadata:
                    self.metadata['version']=tokens[4]

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

    def process_psfgen_line(self,line):
        """
        Process a Psfgen-specific line from the log file. This method extracts metadata from the line and updates the `metadata` attribute accordingly.
        
        Parameters
        ----------
        line : str
            A line from the Psfgen log file that contains Psfgen-specific information to be processed. This can be a string or bytes object containing the log data.
        """
        if line.startswith('total of'):
            if line.endswith('atoms'):
                self.metadata['number_of_atoms']=int(get_single('total of',line))
            elif line.endswith('bonds'):
                self.metadata['number_of_bonds']=int(get_single('total of',line))
            elif line.endswith('angles'):
                self.metadata['number_of_angles']=int(get_single('total of',line))
            elif line.endswith('dihedrals'):
                self.metadata['number_of_dihedrals']=int(get_single('total of',line))
            elif line.endswith('impropers'):
                self.metadata['number_of_impropers']=int(get_single('total of',line))
            elif line.endswith('cross-terms'):
                self.metadata['number_of_cross_terms']=int(get_single('total of',line))
            elif line.endswith('explicit exclusions'):
                self.metadata['number_of_exclusions']=int(get_single('total of',line))

class PDB2PQRLog(LogParser):
    """
    A class for parsing PDB2PQR log files. This class is a subclass of :class:`LogParser <pestifer.util.logparsers.LogParser>` and provides methods for reading, updating, and dumping PDB2PQR log data.
    
    Parameters
    ----------
    basename : str
        The base name for the log parser. This is used to name the output log file.
    """
    section_separator='-'*80 # file is divided into sections
    """
    The separator used to identify sections in the PDB2PQR log file.
    """
    def __init__(self,basename='pdb2pqr-logparser'):
        super().__init__()
        self.processed_separator_idx=[]
        self.processed_sections=[]
        self.metadata={}
        self.basename=basename
        self.progress=0.0

    def update(self,bytes):
        """
        Update the PDB2PQR log parser with new bytes of data. This method appends the new bytes to the byte collector and processes the sections in the log file.
        
        Parameters
        ----------
        bytes : bytes
            The bytes to update the log parser with.
        """
        super().update(bytes)
        separator_idx=[0]+[m.start() for m in re.finditer(self.section_separator,self.byte_collector)]
        # logger.debug(f'update: found {len(separator_idx)-1} sections in {self.basename}')
        for i,j in zip(separator_idx[:-1],separator_idx[1:]):
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

    def process_section(self,bytes):
        """
        Process a section of the PDB2PQR log file. This method identifies the type of section and processes it accordingly.
        """
        # logger.debug(f'process_section: {bytes[:50]}...')
        if 'SUMMARY OF THIS PREDICTION' in bytes:
            self.process_summary(bytes)

    def process_summary(self,bytes):
        """
        Process the summary section of the PDB2PQR log file. This method extracts information from the summary section and updates the `metadata` attribute with relevant data.
        
        Parameters
        ----------
        bytes : bytes
            The bytes representing the summary section to be processed. This can be a string or bytes object
            containing the summary data. 
        """
        # logger.debug(f'process_summary: {bytes[:50]}...')
        lines=bytes.split(os.linesep)
        expected_table=lines[3:]
        logger.debug(f'process_summary: expected_table begins with: {expected_table[0]}')
        table_lines=[]
        for line in expected_table:
            if len(line.strip())==0:
                logger.debug('process_summary: empty line signals end of table')
                break
            tokens=[x.strip() for x in line.split() if x.strip()]
            resname=tokens[0]
            resnum=int(tokens[1])
            reschain=tokens[2]
            respka=float(tokens[3])
            resmodelpka=float(tokens[4])
            if len(tokens)>5:
                resatomtype=tokens[5]
            else:
                resatomtype=None
            table_lines.append(dict(resname=resname,resnum=resnum,reschain=reschain,respka=respka,resmodelpka=resmodelpka,resatomtype=resatomtype))
        if len(table_lines)>0:
            self.metadata['pka_table']=pd.DataFrame(table_lines)

def subcommand_follow_namd_log(filename,basename=None):
    """ 
    Follow a NAMD log file and parse it
    """
    if not os.path.exists(filename):
        logger.debug(f'File {filename} does not exist')
        return -1
    if basename is None:
        basename=os.path.splitext(os.path.basename(filename))[0]
    namd_log=NAMDLog(basename=basename)
    PS=NAMDProgress()
    namd_log.enable_progress_bar(PS)
    namd_log.follow(filename)
    if not namd_log.success():
        logger.debug('NAMD log file did not complete successfully')
        return -2
    namd_log.finalize()
    return 0
    