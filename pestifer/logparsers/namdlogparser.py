# Author: Cameron F. Abrams, <cfa22@drexel.edu>

"""
NAMD log parsing utility
"""

import logging
import os
import re 

import numpy as np
import pandas as pd

from pathlib import Path

from .logparser import LogParser, get_single, get_toflag, get_values

from ..util.progress import NAMDProgress
from ..util.stringthings import my_logger

logger = logging.getLogger(__name__)

class NAMDxstParser(LogParser):
    """ 
    A class for parsing NAMD xst files, which contain information about the simulation cell dimensions.

    Parameters
    ----------
    filename : str
        The path to the NAMD xst file to parse.
    """
    def __init__(self, basename: str = 'namd-xstparser'):
        self.basename = basename
        self.filename = f'{basename}.xst'
        self.dataframe: pd.DataFrame | None = None
        super().__init__()

    @classmethod
    def from_file(cls, basename: str = 'namd-xstparser'):
        """
        Generate a NAMDxst instance from an existing NAMD xst file.
        """
        logger.debug(f'Creating {cls.__name__} from {basename}')
        instance = cls(basename)
        if not os.path.exists(instance.filename):
            # throw a warning and return None
            logger.debug(f'FYI: No {instance.filename} exists for this run.')
            return None
        instance.dataframe = pd.read_csv(instance.filename, skiprows=2, header=None, sep=r'\s+', index_col=None)
        col = 'TS a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z s_x s_y s_z s_u s_v s_w'.split()[:len(instance.dataframe.columns)]
        instance.dataframe.columns = col
        return instance

class NAMDLogParser(LogParser):
    """
    A class for parsing NAMD log files. This class is a subclass of :class:`LogParser <pestifer.logparsers.logparser.LogParser>` and provides methods for reading, updating, and dumping NAMD log data.
    It also includes methods for processing specific lines in the log file, such as those containing information about the simulation, energy calculations, and pressure profiles.
    
    Parameters
    ----------
    basename : str
        The base name for the log parser. This is used to name the output log file.
    """

    info_key = 'Info: '
    """
    The key used to identify information lines in the NAMD log file.
    """
    tcl_key = 'TCL: '
    """
    The key used to identify TCL command lines in the NAMD log file.
    """
    energy_key = 'ENERGY: '
    """
    The key used to identify energy lines in the NAMD log file.
    """
    etitle_key = 'ETITLE: '
    """
    The key used to identify energy title lines in the NAMD log file.
    """
    pressureprofile_key = 'PRESSUREPROFILE: '
    """
    The key used to identify pressure profile lines in the NAMD log file.
    """
    struct_sep = '*****************************'
    """
    The key used to identify the structure summary in the NAMD log file.
    """
    wallclock_key = 'WallClock: '
    """
    The key used to identify wall clock time lines in the NAMD log file.
    """
    restart_key = 'WRITING COORDINATES TO RESTART FILE AT STEP '
    """
    The key used to identify lines indicating that coordinates are being written to a restart file in the NAMD log file.
    """
    performance_key = 'PERFORMANCE: '
    """
    The key used to identify performance lines in the NAMD log file, e.g.,
    PERFORMANCE: 34000  averaging 23.789 ns/day, 0.00726387 sec/step with standard deviation 3.76697e-05
    """
    timing_key = 'TIMING: '
    """
    The key used to identify timing lines in the NAMD log file, e.g.,
    TIMING: 34000  CPU: 100.266, 0.00728793/step  Wall: 100.267, 0.00728689/step, 1.9962 hours remaining, 0.000000 MB of memory in use.
    """
    default_etitle_npt = 'TS,BOND,ANGLE,DIHED,IMPRP,ELECT,VDW,BOUNDARY,MISC,KINETIC,TOTAL,TEMP,POTENTIAL,TOTAL3,TEMPAVG,PRESSURE,GPRESSURE,VOLUME,PRESSAVG,GPRESSAVG'.split(',')
    """
    The default energy titles for NPT ensembles in NAMD log files.
    """
    default_etitle_nvt = 'TS,BOND,ANGLE,DIHED,IMPRP,ELECT,VDW,BOUNDARY,MISC,KINETIC,TOTAL,TEMP,POTENTIAL,TOTAL3,TEMPAVG'.split(',')
    """
    The default energy titles for NVT ensembles in NAMD log files.
    """
    bail_out_key = 'Stack Traceback:'
    """
    The key used to identify lines indicating a stack traceback in the NAMD log file, which may indicate an error or crash.
    """

    def __init__(self, basename='namd-logparser'):
        super().__init__()
        self.line_idx = [0]  # byte offsets of lines
        self.processed_line_idx = []
        self.time_series_data = {}
        self.metadata = {}
        self.dataframes = {}
        self.reading_structure_summary = False
        self.basename = basename
        self._line_processors = {
            self.energy_key: self.process_energy_line,
            self.pressureprofile_key: self.process_pressureprofile_line,
            self.restart_key: self.process_restart_line,
            self.performance_key: self.process_performance_line,
            self.timing_key: self.process_timing_line,  
            self.info_key: self.process_info_line,
            self.tcl_key: self.process_tcl_line,
            self.wallclock_key: self.process_wallclock_line
        }
        self.filename = f'{basename}.log'
    
    @classmethod
    def from_file(cls, filename: Path | str, passfilter: list[str] = []):
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
        instance = cls()
        instance.filename = filename
        instance.basename = os.path.splitext(os.path.basename(filename))[0]
        logger.debug(f'instance.basename: {instance.basename}  filename: {instance.filename}')
        instance.static(filename, passfilter=passfilter)
        return instance

    def static(self, filename: Path | str, passfilter: list[str] = []):
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
            raw = f.read()
        self.write(raw)
        rawlines = raw.splitlines(keepends=True)
        for line in rawlines:
            if not passfilter or any(field in line for field in passfilter):
                self.process_line(line)
        self.finalize()

    def process_struct_summ_datum(self, line: str):
        """
        Process a line from the structure summary section of the NAMD log file.
        """
        if line.endswith('FIXED ATOMS\n'):
            self.metadata['number_of_fixed_atoms'] = int(get_toflag('FIXED ATOMS', line))
        elif line.endswith('ATOMS\n'):
            self.metadata['number_of_atoms'] = int(get_toflag('ATOMS', line))
        elif line.endswith('BONDS\n') and not 'RIGID' in line:
            self.metadata['number_of_bonds'] = int(get_toflag('BONDS', line))
        elif line.endswith('ANGLES\n'):
            self.metadata['number_of_angles'] = int(get_toflag('ANGLES', line))
        elif line.endswith('DIHEDRALS\n'):
            self.metadata['number_of_dihedrals'] = int(get_toflag('DIHEDRALS', line))
        elif line.endswith('IMPROPERS\n'):
            self.metadata['number_of_impropers'] = int(get_toflag('IMPROPERS', line))
        elif line.endswith('CROSSTERMS\n'):
            self.metadata['number_of_crossterms'] = int(get_toflag('CROSSTERMS', line))
        elif line.endswith('EXCLUSIONS\n'):
            self.metadata['number_of_exclusions'] = int(get_toflag('EXCLUSIONS', line))
        elif line.endswith('RIGID BONDS\n'):
            self.metadata['number_of_rigid_bonds'] = int(get_toflag('RIGID BONDS', line))
        elif line.endswith('ATOMS IN LARGEST HYDROGEN GROUP\n'):
            self.metadata['atoms_in_largest_hydrogen_group'] = int(get_toflag('ATOMS IN LARGEST HYDROGEN GROUP', line))
        elif line.startswith('ATOM DENSITY ='):
            self.metadata['atom_density'] = float(get_single('ATOM DENSITY =', line))
        elif 'TOTAL MASS =' in line:
            self.metadata['total_mass'] = float(get_single('TOTAL MASS =', line))
        elif 'TOTAL CHARGE =' in line:
            self.metadata['total_charge'] = float(get_single('TOTAL CHARGE =', line))
        elif 'MASS DENSITY =' in line:
            self.metadata['mass_density'] = float(get_single('MASS DENSITY =', line))

    def process_info_line(self, line: str):
        """
        Process a line from the information section of the NAMD log file.
        """
        # logger.debug(f'process_info: {line}')
        if line.startswith('STRUCTURE SUMMARY:'):
            self.reading_structure_summary = True
        if line.startswith(self.struct_sep):
            self.reading_structure_summary = False
        if self.reading_structure_summary:
            self.process_struct_summ_datum(line)
            return
        if line.startswith(f'TIMESTEP'):
            self.metadata['timestep'] = float(get_single('TIMESTEP', line))
        elif 'FIRST TIMESTEP' in line:
            self.metadata['first_timestep'] = int(get_single('FIRST TIMESTEP', line))
        elif 'NUMBER OF STEPS' in line:
            self.metadata['number_of_steps'] = int(get_single('NUMBER OF STEPS', line))
        elif 'RANDOM NUMBER SEED' in line:
            self.metadata['random_number_seed'] = int(get_single('RANDOM NUMBER SEED', line))
        elif 'RESTART FILENAME' in line:
            self.metadata['restart_filename'] = get_single('RESTART FILENAME', line)
        elif 'RESTART FREQUENCY' in line:
            self.metadata['restart_frequency'] = int(get_single('RESTART FREQUENCY', line))
        elif 'OUTPUT FILENAME' in line:
            self.metadata['output_filename'] = get_single('OUTPUT FILENAME', line)
        elif 'ENERGY OUTPUT STEPS' in line:
            self.metadata['energy_output_steps'] = int(get_single('ENERGY OUTPUT STEPS', line))
        elif 'LANGEVIN DYNAMICS ACTIVE' in line:
            self.metadata['ensemble'] = 'NVT'
        elif 'LANGEVIN PISTON PRESSURE CONTROL ACTIVE' in line:
            self.metadata['ensemble'] = 'NPT'
        elif 'SHAPE OF CELL IS CONSTRAINED IN X-Y PLANE' in line:
            self.metadata['ensemble'] = 'NPAT'
        elif 'PERIODIC CELL BASIS 1' in line:
            self.metadata['periodic_cell_basis_1'] = get_values('PERIODIC CELL BASIS 1', line, dtype=float)
        elif 'PERIODIC CELL BASIS 2' in line:
            self.metadata['periodic_cell_basis_2'] = get_values('PERIODIC CELL BASIS 2', line, dtype=float)
        elif 'PERIODIC CELL BASIS 3' in line:
            self.metadata['periodic_cell_basis_3'] = get_values('PERIODIC CELL BASIS 3', line, dtype=float)
        elif 'SLAB THICKNESS:' in line:
            self.metadata['slab_thickness'] = float(get_single('SLAB THICKNESS:', line))
        elif 'NUMBER OF SLABS:' in line:
            self.metadata['number_of_pressure_slabs'] = int(get_single('NUMBER OF SLABS:', line))

    def process_tcl_line(self, line: str):
        """
        Process a line from the TCL command section of the NAMD log file.
        """
        if line.startswith('Running for'):
            self.metadata['running_for'] = int(get_single('Running for', line))
        elif line.startswith('Minimizing for'):
            self.metadata['minimizing_for'] = int(get_single('Minimizing for', line))
            self.metadata['ensemble'] = 'minimize'

    def process_energy_line(self, line: str):
        """
        Process a line from the energy section of the NAMD log file.

        Parameters
        ----------
        line : str
            A line from the NAMD log file that contains energy data.
        """
        tokens = [x.strip() for x in line.split()]
        if 'etitle' not in self.metadata:
            if len(tokens) == len(self.default_etitle_npt):
                # logger.debug(f'process_energy: {len(tokens)} tokens found, using default etitle for NPT')
                self.metadata['etitle'] = self.default_etitle_npt
            elif len(tokens) == len(self.default_etitle_nvt):
                # logger.debug(f'process_energy: {len(tokens)} tokens found, using default etitle for NVT')
                self.metadata['etitle'] = self.default_etitle_nvt
            else:
                # logger.debug(f'process_energy: {len(tokens)} tokens found, but either {len(self.default_etitle_nvt)} or {len(self.default_etitle_npt)} expected')
                return
        else:
            if len(tokens) != len(self.metadata['etitle']):
                # logger.debug(f'process_energy: {len(tokens)} tokens found, but {len(self.metadata["etitle"])} expected')
                return
        tokens[0] = int(tokens[0])
        for i in range(1, len(tokens)):
            tokens[i] = float(tokens[i])
        new_line = {k: v for k, v in zip(self.metadata['etitle'], tokens)}
        if 'energy' not in self.time_series_data:
            self.time_series_data['energy'] = []
            if 'first_timestep' not in self.metadata:
                self.metadata['first_timestep'] = new_line['TS']
        self.time_series_data['energy'].append(new_line)

    def process_energy_title(self, line: str):
        """
        Process a line from the energy title section of the NAMD log file.

        Parameters
        ----------
        line : str
            A line from the NAMD log file that contains the energy titles.
        """
        tokens = [x.strip() for x in line.split()]
        if 'etitle' not in self.metadata:
            self.metadata['etitle'] = tokens
        else:
            if len(tokens) != len(self.metadata['etitle']):
                logger.debug(f'process_energy_title: {len(tokens)} tokens found, but {len(self.metadata["etitle"])} expected')

    def process_restart_line(self, line: str):
        """
        Process a line from the restart section of the NAMD log file.

        Parameters
        ----------
        line : str
            A line from the NAMD log file that contains restart information.
        """
        # the only thing in line should be the time step at which the last restart was written
        if 'restart' not in self.time_series_data:
            self.time_series_data['restart'] = []
        ts = int(line.strip())
        self.time_series_data['restart'].append(ts)

    def process_pressureprofile_line(self, line: str):
        """
        Process a line from the pressure profile section of the NAMD log file.

        Parameters
        ----------
        line : str
            A line from the NAMD log file that contains pressure profile data.
        """
        slab_thickness = self.metadata.get('slab_thickness', None)
        if slab_thickness is None:
            logger.debug('process_pressureprofile_line: slab_thickness is not defined in metadata')
            return
        tokens = [x.strip() for x in line.split()]
        TS = int(tokens[0])
        # logger.debug(f'process_pressureprofile_line: TS {tokens[0]}')
        for i in range(1, len(tokens)):
            tokens[i] = float(tokens[i])
        this_col = tokens[1:]
        pressure_series = {'TS': TS}
        # depth_series = {'TS': TS}
        pressure_series.update({k: v for k, v in zip([f"{'xyz'[i % 3]}_{(i//3)}" for i in range(len(this_col))], this_col)})
        # depth_series.update(   {k: v for k, v in zip(range(len(this_col)//3), [(i+0.5)*slab_thickness for i in range(len(this_col)//3)])})
        if 'number_of_pressure_slabs' not in self.metadata:
            self.metadata['number_of_pressure_slabs'] = len(this_col)
        if 'pressureprofile' not in self.time_series_data:
            self.time_series_data['pressureprofile'] = []
        # if 'depthprofile' not in self.time_series_data:
        #     self.time_series_data['depthprofile'] = []
        # self.time_series_data['depthprofile'].append(depth_series)
        self.time_series_data['pressureprofile'].append(pressure_series)

    def process_wallclock_line(self, line: str):
        """
        Process a line from the wall clock time section of the NAMD log file.

        Parameters
        ----------
        line : str
            A line from the NAMD log file that contains the wall clock time.
        """
        tokens = [x.strip() for x in line.split()]
        self.metadata['wallclock_time'] = float(tokens[0])

    def process_performance_line(self, line: str):
        """
        Process a line from the performance section of the NAMD log file, e.g.,
        PERFORMANCE: 34000  averaging 23.789 ns/day, 0.00726387 sec/step with standard deviation 3.76697e-05

        Parameters
        ----------
        line : str
            A line from the NAMD log file that contains performance data.
        """
        tokens = [x.strip() for x in line.split()]
        if len(tokens) < 5:
            logger.debug(f'process_performance_line: {line} does not have enough tokens')
            return
        if 'performance' not in self.time_series_data:
            self.time_series_data['performance'] = []
        self.time_series_data['performance'].append({
            'steps': int(tokens[0]),
            'ns_per_day': float(tokens[2]),
            'sec_per_step': float(tokens[4]),
            'std_dev': float(tokens[-1]),
        })

    def process_timing_line(self, line: str):
        """
        Process a line from the timing section of the NAMD log file, e.g.,
        TIMING: 34000  CPU: 100.266, 0.00728793/step  Wall: 100.267, 0.00728689/step, 1.9962 hours remaining, 0.000000 MB of memory in use.
        or
        TIMING: 19000  CPU: 130.551, 0.025328/step  Wall: 140.3    , 0.0272845/step , 6.33326 ns/days       , 0.0621481 hours remaining, 0.000000 MB of memory in use.

        Parameters
        ----------
        line : str
            A line from the NAMD log file that contains timing data.
        """
        elements = line.split(',')
        if len(elements) < 5:
            # logger.debug(f'process_timing_line: {line} does not have enough elements')
            return
        if len(elements) == 5:
            hasnsperday = False
        elif len(elements) == 6:
            hasnsperday = True
        else:
            logger.debug(f'process_timing_line: {line} has too many elements ({len(elements)})')
            return

        # element[0] contributes TS and cpu_time
        # element[1] contributes cpu_per_step and wall_time
        # element[2] contributes wall_per_step
        # if hasnsperday:
        #    element[3] contributes ns_per_day
        #    element[4] contributes hours_remaining
        #    element[5] contributes memory_in_use
        # else 
        #     element[3] contributes hours_remaining
        #     element[4] contributes memory_in_use        

        # process element[0] "19000  CPU: 130.551"
        tokens = [x.strip() for x in elements[0].split()]
        if len(tokens) != 3:
            logger.debug(f'process_timing_line: {line} does not have enough tokens in first element')
            return
        if not tokens[0].isdigit():
            logger.debug(f'process_timing_line: {line} first token is not a digit: {tokens[0]}')
            return
        TS = int(tokens[0])
        cpu_time = float(tokens[2])

        # process element[1] "0.025328/step  Wall: 140.3    "
        tokens = [x.strip() for x in elements[1].replace('/step', '').split()]  # [0.025328, Wall:, 140.3]
        if len(tokens) != 3:
            logger.debug(f'process_timing_line: {line} does not have enough tokens in second element')
            return
        if not tokens[0].replace('.', '', 1).isdigit():
            logger.debug(f'process_timing_line: {line} first token is not a digit: {tokens[0]}')
            return
        cpu_per_step = float(tokens[0])
        wall_time = float(tokens[2])

        # process element[2] "0.00728689/step"
        tokens = [x.strip() for x in elements[2].replace('/step', '').split()]  # [0.00728689]
        if len(tokens) != 1:
            logger.debug(f'process_timing_line: {line} does not have enough tokens in third element')
            return
        if not tokens[0].replace('.', '', 1).isdigit():
            logger.debug(f'process_timing_line: {line} first token is not a digit: {tokens[0]}')
            return
        wall_per_step = float(tokens[0])
        if hasnsperday:
        #     process element[3] "6.33326 ns/days       "
            tokens = [x.strip() for x in elements[3].replace('ns/days', '').split()]
            if len(tokens) != 1:
                logger.debug(f'process_timing_line: {line} does not have enough tokens ({tokens})in fourth element')
                return
            if not tokens[0].replace('.', '', 1).isdigit():
                logger.debug(f'process_timing_line: {line} first token is not a digit: {tokens[0]}')
                return
            ns_per_day = float(tokens[0])  # orphaned
        #     process element[4] "0.0621481 hours remaining"
            tokens = [x.strip() for x in elements[4].replace('hours remaining', '').split()]  # [0.0621481]
            if len(tokens) != 1:
                logger.debug(f'process_timing_line: {line} does not have enough tokens ({tokens}) in fifth element')
                return
            if not tokens[0].replace('.', '', 1).isdigit():
                logger.debug(f'process_timing_line: {line} first token is not a digit: {tokens[0]}')
                return
            hours_remaining = float(tokens[0])
        #     process element[5] "0.000000 MB of memory in use."
            tokens = [x.strip() for x in elements[5].replace('MB of memory in use.', '').split()]  # [0.000000]
            if len(tokens) != 1:
                logger.debug(f'process_timing_line: {line} does not have enough tokens in sixth element')
                return
            if not tokens[0].replace('.', '', 1).isdigit():
                logger.debug(f'process_timing_line: {line} first token is not a digit: {tokens[0]}')
                return
            memory_in_use = float(tokens[0])
        else:
        #     process element[3] "0.0621481 hours remaining"
            tokens = [x.strip() for x in elements[3].replace('hours remaining', '').split()]  # [0.0621481]
            if len(tokens) != 1:
                logger.debug(f'process_timing_line: {line} does not have enough tokens in fourth element {tokens}')
                return
            if not tokens[0].replace('.', '', 1).isdigit():
                logger.debug(f'process_timing_line: {line} first token is not a digit: {tokens[0]}')
                return
            hours_remaining = float(tokens[0])
        #     process element[4] "0.000000 MB of memory in use."
            tokens = [x.strip() for x in elements[4].replace('MB of memory in use.', '').split()]  # [0.000000]
            if len(tokens) != 1:
                logger.debug(f'process_timing_line: {line} does not have enough tokens in fifth element {tokens}')
                return
            if not tokens[0].replace('.', '', 1).isdigit():
                logger.debug(f'process_timing_line: {line} first token is not a digit: {tokens[0]}')
                return
            memory_in_use = float(tokens[0])
        if not 'timing' in self.time_series_data:
            self.time_series_data['timing'] = []
        self.time_series_data['timing'].append({
            'steps': TS,
            'cpu_time': cpu_time,
            'cpu_per_step': cpu_per_step,
            'wall_time': wall_time,
            'wall_per_step': wall_per_step,
            'hours_remaining': hours_remaining,
            'memory_in_use': memory_in_use,
        })

    def process_line(self, line: str):
        """
        Process a line from the NAMD log file. This method identifies the type of line (information, TCL command, energy, pressure profile, or wall clock time) and processes it accordingly.
        
        Parameters
        ----------
        line : str
            A line from the NAMD log file to be processed.
        """
        assert line.endswith(os.linesep), f'process_line: {line} does not end with os.linesep'
        if self.bail_out_key in line:
            logger.debug(f'process_line: {line} contains bail out key, stopping processing')
            return -1
        if line.startswith(self.etitle_key):
            if 'etitle' in self.metadata:
                # logger.debug(f'process_energy_title: {line} already has etitle')
                return 0
            self.process_energy_title(line[len(self.etitle_key):])
            return 0
        for key, func in self._line_processors.items():
            if line.startswith(key):
                func(line[len(key):])
                return 0

    def update(self, bytes: str):
        """
        Update the NAMD log parser with new bytes of data. This method appends the new bytes to the byte collector and processes the lines in the log file. It identifies the end of each line and processes each line based on its content.  This is best used on a log file that is being written to, such as a live NAMD simulation log file.  For static files, use the :meth:`static <pestifer.logparsers.namdlogparser.NAMDLog.static>` method instead.

        Parameters
        ----------
        bytes : bytes
            The bytes to update the log parser with. This can be a string or bytes object containing the log data.
        """
        super().update(bytes) # this just appends the bytes to the byte_collector
        last_line_idx = self.line_idx[-1] # recall byte index of last line processed
        addl_line_idx = [m.start()+1+last_line_idx for m in re.finditer(os.linesep,self.byte_collector[last_line_idx:])]
        scan_ldx = [last_line_idx] # in case last line was incomplete in a previous pass
        if len(addl_line_idx) > 1:
            scan_ldx.extend(addl_line_idx)
        for i, j in zip(scan_ldx[:-1], scan_ldx[1:]):
            if i not in self.processed_line_idx:
                line = self.byte_collector[i:j]
                result = self.process_line(line)
                if result == -1:  # bail out key found, stop processing
                    return
                self.processed_line_idx.append(i)
        if len(addl_line_idx) > 1:
            self.line_idx.extend(addl_line_idx)
            last_line = self.byte_collector[addl_line_idx[-1]:]
            if last_line.endswith(os.linesep):
                self.process_line(last_line)
                self.processed_line_idx.append(addl_line_idx[-1])

    def measure_progress(self):
        """
        Measure the progress of the NAMD simulation based on the metadata and time series data. This method calculates the fraction of completed steps relative to the total number of steps, using the first time step and the last recorded time step in the energy data.
        """
        if 'number_of_steps' not in self.metadata:
            # logger.debug('measure_progress: number_of_steps not in metadata')
            return 0.0
        number_of_steps = self.metadata['number_of_steps'] # this will be zero for a minimization
        if 'first_timestep' not in self.metadata:
            # logger.debug('measure_progress: first_timestep not in metadata')
            return 0.0
        first_time_step = self.metadata['first_timestep']
        if 'running_for' in self.metadata:
            running_for = self.metadata['running_for']
            if running_for > 0:
                number_of_steps = running_for
        elif 'minimizing_for' in self.metadata:
            minimizing_for = self.metadata['minimizing_for']
            if minimizing_for > 0:
                number_of_steps = minimizing_for
        last_row = self.time_series_data['energy'][-1] if 'energy' in self.time_series_data and len(self.time_series_data['energy']) > 0 else None
        if last_row is None:
            # logger.debug('measure_progress: last_row is None')
            return 0.0
        if 'TS' not in last_row:
            # logger.debug('measure_progress: TS not in last_row')
            return 0.0
        most_recent_time_step = last_row['TS']
        complete_steps = most_recent_time_step - first_time_step
        return complete_steps / number_of_steps

    def success(self):
        """ 
        Check if the NAMD log parsing was successful. This method checks if the metadata contains the ``wallclock_time`` key, which indicates that the log file has been processed successfully. 
        """
        return 'wallclock_time' in self.metadata

    def finalize(self):
        """
        Finalize the log parsing by creating dataframes for each time series.
        """
        # parse the XST file
        logger.debug('finalize namdlog parser metadata:')
        my_logger(self.metadata, logger.debug)
        self.auxlogparser = NAMDxstParser.from_file(basename=os.path.splitext(self.filename)[0])
        for key in self.time_series_data:
            self.dataframes[key] = pd.DataFrame(self.time_series_data[key])
        if self.auxlogparser:
            self.dataframes['xst'] = self.auxlogparser.dataframe
        # add a 'DENSITY' column to the energy dataframe
        if 'energy' in self.dataframes:
            if 'total_mass' in self.metadata and 'VOLUME' in self.dataframes['energy'].columns:
                self.dataframes['energy']['DENSITY'] = self.metadata['total_mass'] / self.dataframes['energy']['VOLUME']
        # convert the integer slab index column headings in the pressure profile dataframe to floating point z-coordinates using metadata['slab_thickness']
        # if 'pressureprofile' in self.dataframes:
        #     if 'number_of_pressure_slabs' in self.metadata and 'slab_thickness' in self.metadata:
        #         slab_thickness = self.metadata['slab_thickness']
        #         number_of_pressure_slabs = self.metadata['number_of_pressure_slabs']
        #         z_coords = [(i + 0.5) * slab_thickness for i in range(number_of_pressure_slabs)]
        # If we did not find the first time step, infer it from the energy log
        if 'first_timestep' not in self.metadata:
            if 'energy' in self.dataframes:
                self.metadata['first_timestep'] = int(self.dataframes['energy'].iloc[0]['TS'])
        return self
    
    def write_csv(self):
        """
        Write the parsed data to CSV files. This method creates a CSV file for each dataframe in the `dataframes` attribute, using the basename provided during initialization.
        The files will be named `<basename>-<key>.csv`, where `<key>` is the key of the dataframe in the `dataframes` dictionary.
        """
        self.csvfilenames = {}
        for key in self.dataframes:
            self.dataframes[key].to_csv(f'{self.basename}-{key}.csv', index=False)
            self.csvfilenames[key] = f'{self.basename}-{key}.csv'
        return self.csvfilenames

def subcommand_follow_namd_log(filename, basename: str | None = None):
    """ 
    Follow a NAMD log file and parse it
    """
    if not os.path.exists(filename):
        logger.debug(f'File {filename} does not exist')
        return -1
    if basename is None:
        basename = os.path.splitext(os.path.basename(filename))[0]
    namd_log = NAMDLogParser(basename=basename)
    PS = NAMDProgress()
    namd_log.enable_progress_bar(PS)
    namd_log.follow(filename)
    if not namd_log.success():
        logger.debug('NAMD log file did not complete successfully')
        return -2
    namd_log.finalize()
    return 0