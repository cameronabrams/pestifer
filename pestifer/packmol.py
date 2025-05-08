# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" Defines the PackmolInputWriter class for generating packmol input files 
    Defines the PackmolLog class for parsing packmol log files
"""

import datetime
import logging
import os
import re

import numpy as np
import pandas as pd

from .command import Command
from .util.progress import PackmolProgress
from .scriptwriters import Filewriter
from .stringthings import ByteCollector, FileCollector, striplist

logger=logging.getLogger(__name__)

def get_single(flag,bytes,default='0'):
    # logger.debug(f'get_single: {flag} in {bytes}')
    try:
        idx=bytes.index(flag)+len(flag)
    except ValueError:
        logger.debug(f'get_single: {flag} not found in {bytes}')
        return default
    eol=bytes[idx:].index('\n')+idx
    substr=bytes[idx:eol].split()[0]
    return substr.strip()

class PackmolInputWriter(Filewriter):
    def __init__(self,config):
        super().__init__(comment_char='#')
        self.indent=4*' '
        self.config=config
        self.progress=self.config.progress
        self.F=FileCollector()
        self.default_ext='.inp'
        self.default_script=f'packmol{self.default_ext}'
        self.scriptname=self.default_script
    
    def newscript(self,basename=None):
        timestampstr=datetime.datetime.today().ctime()
        if basename:
            self.basename=basename
        else:
            self.basename=os.path.splitext(self.default_script)[0]
        self.scriptname=f'{self.basename}{self.default_ext}'
        self.newfile(self.scriptname)
        self.banner(f'{__package__}: {self.basename}{self.default_ext}')
        self.banner(f'Created {timestampstr}')

    def writescript(self):
        self.writefile()

    def runscript(self,*args,**options):
        assert hasattr(self,'scriptname'),f'No scriptname set.'
        self.logname=f'{self.basename}.log'
        logger.debug(f'Log file: {self.logname}')
        cmd=Command(f'{self.config.shell_commands["packmol"]} < {self.scriptname}')
        progress_struct=None
        if self.progress:
            progress_struct=PackmolProgress()
        return cmd.run(ignore_codes=[173],logfile=self.logname,progress=progress_struct)

class PackmolLog(ByteCollector):
    banner_separator='#'*80
    section_separator='-'*80
    def __init__(self):
        super().__init__()
        self.processed_separator_idx=[]
        self.processed_sections=[]
        self.processed_banners=[]
        self.metadata={}
        self.gencan={}

    def update(self,bytes):
        self.write(bytes)
        separator_idx=[0]+[m.start() for m in re.finditer(self.section_separator,self.byte_collector)]
        for i,j in zip(separator_idx[:-1],separator_idx[1:]):
            if i not in self.processed_separator_idx:
                self.processed_separator_idx.append(i)
                self.process_section(self.byte_collector[i:j])
    
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

    def process_gencan_success(self,bytes):
        mtype=int(get_single('Packing solved for molecules of type',bytes))
        obj_func_val=float(get_single('Objective function value:',bytes))
        max_viol_target_dist=float(get_single('Maximum violation of target distance:',bytes))
        max_viol_constr=float(get_single('Max. constraint violation:',bytes))
        match=next((d for d in self.metadata['structures'] if d['idx'] == mtype), None)
        if match:
            match['gencan_success']=dict(objective_function_value=obj_func_val,max_viol_target_distance=max_viol_target_dist,max_viol_constr=max_viol_constr)

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

    def finalize(self):
        self.gencan_df={}
        for k,v in self.gencan.items():
            self.gencan_df[k]=pd.DataFrame(v)
 