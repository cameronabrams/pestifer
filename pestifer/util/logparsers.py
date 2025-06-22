# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" Defines the PackmolLog, PsfgenLog, and NAMDLog classes for parsing
    log files
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

from ..stringthings import ByteCollector
from ..util.progress import PackmolProgress, NAMDProgress

logger=logging.getLogger(__name__)

def get_toflag(flag,bytes):
    try:
        idx=bytes.index(flag)
    except ValueError:
        logger.debug(f'get_toflag: {flag} not found in {bytes}')
    substr=bytes[:idx]
    return substr.strip()

def get_toeol(flag,bytes):
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
    # logger.debug(f'get_single: {flag} in {bytes}')
    toel=get_toeol(flag,bytes)
    if toel is None:
        return []
    return [x.strip() for x in toel.split()]

def get_single(flag,bytes):
    # logger.debug(f'get_single: {flag} in {bytes}')
    vals=get_tokens(flag,bytes)
    if len(vals)==0:
        return None
    return vals[0]

def get_values(flag,bytes,dtype=float):
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
    def __init__(self):
        super().__init__()
        self.progress=0.0
        self.progress_bar=None

    def static(self,filename):
        logger.debug(f'Initiating {self.__class__.__name__} from {filename}')
        with open(filename,'r') as f:
            rawlines=f.read()
        self.update(rawlines)

    def update(self,bytes):
        self.write(bytes)
        # logger.debug(f'Updating {self.__class__.__name__} with {len(bytes)} bytes -> {len(self.byte_collector)} bytes collected')

    def dump(self,basename='logparser'):
        logger.debug(f'Dumping {self.__class__.__name__} to {basename}')
        self.byte_collector.write_file(f'{basename}.log')

    def measure_progress(self):
        return self.progress
    
    def enable_progress_bar(self,PS=None):
        if PS is None:
            self.progress_bar=None
        else:
            self.progress_bar=PS
            self.progress_bar.register_update_function(self.measure_progress)

    def update_progress_bar(self):
        if self.progress_bar is not None:
            self.progress_bar.go()

    def success(self):
        return False

    def follow(self,filename,sleep_interval=0.5,timeout_intervals=120):
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
    def __init__(self,basename='vmd-logparser'):
        super().__init__()
        self.basename=basename

class NAMDxst():
    def __init__(self,filename):
        celldf=pd.read_csv(filename,skiprows=2,header=None,sep=r'\s+',index_col=None)
        col='step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z s_x s_y s_z s_u s_v s_w'.split()[:len(celldf.columns)]
        celldf.columns=col
        self.df=celldf
    def add_file(self,filename):
        celldf=pd.read_csv(filename,skiprows=2,header=None,sep=r'\s+',index_col=None)
        col='step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z s_x s_y s_z s_u s_v s_w'.split()[:len(celldf.columns)]
        celldf.columns=col
        self.df=pd.concat((self.df,celldf))
    def concat(self,other):
        self.df=pd.concat((self.df,other.df))
        
class NAMDLog(LogParser):
    info_key='Info: '
    tcl_key='TCL: '
    energy_key='ENERGY: '
    etitle_key='ETITLE: '
    pressureprofile_key='PRESSUREPROFILE: '
    struct_sep='*****************************'
    wallclock_key='WallClock: '
    default_etitle_npt='TS,BOND,ANGLE,DIHED,IMPRP,ELECT,VDW,BOUNDARY,MISC,KINETIC,TOTAL,TEMP,POTENTIAL,TOTAL3,TEMPAVG,PRESSURE,GPRESSURE,VOLUME,PRESSAVG,GPRESSAVG'.split(',')
    default_etitle_nvt='TS,BOND,ANGLE,DIHED,IMPRP,ELECT,VDW,BOUNDARY,MISC,KINETIC,TOTAL,TEMP,POTENTIAL,TOTAL3,TEMPAVG'.split(',')  
    def __init__(self,basename='namd-logparser'):
        super().__init__()
        self.line_idx=[0] # byte offsets of lines
        self.processed_line_idx=[]
        self.energy_df=pd.DataFrame()
        self.pressureprofile_df=pd.DataFrame() # each column is a pressure profile at a given timestep
        self.metadata={}
        self.reading_structure_summary=False
        self.basename=basename

    def process_struct_summ_datum(self,line):
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
            logger.debug(f'process_info: {line}')
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
        if line.startswith('Running for'):
            self.metadata['running_for']=int(get_single('Running for',line))
        elif line.startswith('Minimizing for'):
            self.metadata['minimizing_for']=int(get_single('Minimizing for',line))
            self.metadata['ensemble']='minimize'

    def process_energy_line(self,line):
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
        if self.energy_df.empty:
            self.energy_df=pd.DataFrame(new_line,index=[0])
        else:
            self.energy_df=pd.concat([self.energy_df,pd.DataFrame(new_line,index=[0])],ignore_index=True)

    def process_energy_title(self,line):
        tokens=[x.strip() for x in line.split()]
        if 'etitle' not in self.metadata:
            self.metadata['etitle']=tokens
        else:
            if len(tokens)!=len(self.metadata['etitle']):
                logger.debug(f'process_energy_title: {len(tokens)} tokens found, but {len(self.metadata["etitle"])} expected')


    def update(self,bytes):
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

    def process_pressureprofile_line(self,line):
        tokens=[x.strip() for x in line.split()]
        tokens[0]=int(tokens[0])
        logger.debug(f'process_pressureprofile_line: TS {tokens[0]}')
        for i in range(1,len(tokens)):
            tokens[i]=float(tokens[i])
        this_col=tokens[1:]
        if not 'number_of_pressure_slabs' in self.metadata:
            self.metadata['number_of_pressure_slabs']=len(this_col)
        if self.pressureprofile_df.empty:
            # we need to set the indices to be the slab indices
            self.pressureprofile_df=pd.DataFrame(this_col,index=np.arange(self.metadata['number_of_pressure_slabs']),columns=[str(tokens[0])])
        else:
            # append this col to the dataframe
            new_col=pd.DataFrame(this_col,index=np.arange(self.metadata['number_of_pressure_slabs']),columns=[str(tokens[0])])
            self.pressureprofile_df=pd.concat([self.pressureprofile_df,new_col],axis=1)
                # if self.pressureprofile_df.shape[1]==2:
                #     logger.debug(f'{self.pressureprofile_df.head().to_string()}')

    def process_wallclock_line(self,line):
        tokens=[x.strip() for x in line.split()]
        self.metadata['wallclock_time']=float(tokens[0])

    def process_line(self,line):
        assert line.endswith(os.linesep),f'process_line: {line} does not end with os.linesep'
        if line.startswith(self.info_key):
            o=len(self.info_key)
            self.process_info_line(line[o:])
        elif line.startswith(self.tcl_key):
            o=len(self.tcl_key)
            self.process_tcl_line(line[o:])
        elif line.startswith(self.etitle_key):
            if 'etitle' in self.metadata:
                # logger.debug(f'process_energy_title: {line} already has etitle')
                return
            o=len(self.etitle_key)
            self.process_energy_title(line[o:])
        elif line.startswith(self.energy_key):
            o=len(self.energy_key)
            self.process_energy_line(line[o:])
            if 'first_timestep' not in self.metadata:
                self.metadata['first_timestep']=int(self.energy_df.iloc[0]['TS'])
        elif line.startswith(self.pressureprofile_key):
            o=len(self.pressureprofile_key)
            self.process_pressureprofile_line(line[o:])
        elif line.startswith(self.wallclock_key):
            o=len(self.wallclock_key)
            self.process_wallclock_line(line[o:])

    def measure_progress(self):
        if 'number_of_steps' not in self.metadata:
            # logger.debug('measure_progress: number_of_steps not in metadata')
            return 0.0
        number_of_steps=self.metadata['number_of_steps']
        if 'first_timestep' not in self.metadata:
            # logger.debug('measure_progress: first_timestep not in metadata')
            return 0.0
        first_time_step=self.metadata['first_timestep']
        if self.energy_df.empty:
            # logger.debug('measure_progress: energy_df is empty')
            return 0.0
        if 'running_for' in self.metadata:
            running_for=self.metadata['running_for']
            if running_for>0:
                number_of_steps=running_for
        elif 'minimizing_for' in self.metadata:
            minimizing_for=self.metadata['minimizing_for']
            if minimizing_for>0:
                number_of_steps=minimizing_for
        last_row=self.energy_df.iloc[-1]
        most_recent_time_step=last_row['TS']
        complete_steps=most_recent_time_step-first_time_step
        return complete_steps/number_of_steps
    
    def success(self):
        return 'wallclock_time' in self.metadata

    def finalize(self):
        if not self.energy_df.empty:
            self.energy_df.to_csv(f'{self.basename}-energy.csv',index=False)
        if not self.pressureprofile_df.empty:
            self.pressureprofile_df.to_csv(f'{self.basename}-pressureprofile.csv',index=True)

class PackmolPhase(Enum):
    UNINITIALIZED = 0
    INITIALIZATION = 1

class PackmolLog(LogParser):
    banner_separator='#'*80  # banners can occur within sections
    section_separator='-'*80 # file is divided into sections
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
        super().update(bytes)
        separator_idx=[0]+[m.start() for m in re.finditer(self.section_separator,self.byte_collector)]
        for i,j in zip(separator_idx[:-1],separator_idx[1:]):
            if i not in self.processed_separator_idx:
                self.processed_separator_idx.append(i)
                self.process_section(self.byte_collector[i:j])

    def measure_progress(self):
        self.progress=self.metadata['current_total_gencan_loops']/self.metadata['max_total_gencan_loops'] if 'current_total_gencan_loops' in self.metadata and 'max_total_gencan_loops' in self.metadata else 0.0
        return super().measure_progress()
    
    def process_banner(self,bytes):
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
    info_keys=['Info) ', 'Info: ']
    psfgen_key='psfgen) '
    def __init__(self,basename='psfgen-logparser'):
        super().__init__()
        self.line_idx=[0] # byte offsets of lines
        self.processed_line_idx=[]
        self.metadata={}
        self.basename=basename

    def update(self,bytes):
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
        if 'success' in self.metadata:
            return self.metadata['success']
        return False

    def process_psfgen_line(self,line):
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

def subcommand_follow_namd_log(filename,basename=None):
    """ Follow a NAMD log file and parse it
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
    