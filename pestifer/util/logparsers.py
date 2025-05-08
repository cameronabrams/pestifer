# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" Defines the PackmolLog class and NAMDLog classes for parsing 
    log files
"""

import logging
import os
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ..stringthings import ByteCollector, my_logger

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

class PackmolLog(ByteCollector):
    banner_separator='#'*80  # banners can occur within sections
    section_separator='-'*80 # file is divided into sections
    def __init__(self):
        super().__init__()
        self.processed_separator_idx=[]
        self.processed_sections=[]
        self.processed_banners=[]
        self.metadata={}
        self.gencan={}

    def static(self,filename):
        logger.debug(f'Initiating PackmolLog from {filename}')
        with open(filename,'r') as f:
            rawlines=f.read()
        self.update(rawlines)

    def update(self,bytes):
        logger.debug(f'Updating PackmolLog with {len(bytes)} bytes')
        self.write(bytes)
        # check to see if we need to update the parsing
        separator_idx=[0]+[m.start() for m in re.finditer(self.section_separator,self.byte_collector)]
        for i,j in zip(separator_idx[:-1],separator_idx[1:]):
            if i not in self.processed_separator_idx:
                self.processed_separator_idx.append(i)
                self.process_section(self.byte_collector[i:j])
    
    def dump(self,basename='packmol'):
        logger.debug(f'Dumping PackmolLog to {basename}')
        self.byte_collector.write_file(f'{basename}.log')

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

    def finalize(self,basename='packmol-gencan'):
        self.gencan_df={}
        for k,v in self.gencan.items():
            self.gencan_df[k]=pd.DataFrame(v)
            self.gencan_df[k].to_csv(f'{basename}_{k}.csv',index=False)
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
        plt.savefig(f'{basename}.png')

def getinfo(name,line):
    parseline=line[len('Info:'):].split()
    while '=' in parseline:
        parseline.remove('=')
    toks=name.split()
    while len(toks)>0:
        parseline.remove(toks.pop(0))
    return parseline[0]

def gettclinfo(name,line):
    assert name in line
    namelen=len(name.split())
    tok=line.split()[namelen+1]
    return int(tok)

class NAMDxst:
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
        
class NAMDLog(ByteCollector):
    groupnames=['ETITLE:','ENERGY:','charmrun>','Charm++>','Info:','TIMING:','WallClock:','WRITING','TCL:']
    infonames=['TIMESTEP','FIRST TIMESTEP','NUMBER OF STEPS','RANDOM NUMBER SEED','TOTAL MASS','TOTAL CHARGE','RESTART FILENAME','RESTART FREQUENCY','OUTPUT FILENAME']
    countnames=['ATOMS','BONDS','ANGLES','DIHEDRALS','IMPROPERS','CROSSTERMS']
    def __init__(self,filename,inherited_etitles=[]):
        self.successful_run=False
        logger.debug(f'Initiating NAMDLog from {filename}')
        if inherited_etitles:
            logger.debug(f'  Explicit etitles: {inherited_etitles}')
        self.filename=filename
        self.inherited_etitles=inherited_etitles
        if not os.path.exists(filename):
            raise FileNotFoundError(f'{filename} not found')
        with open(filename,'r') as f:
            rawlines=f.read().split('\n')
        self.numlines=len(rawlines)
        logger.debug(f'{filename}: {self.numlines} total lines')
        self.groups={}
        self.info={}
        self.counts={}
        for l in rawlines:
            if len(l)>0:
                tok=l.split()[0]
                if tok in self.groupnames:
                    if not tok in self.groups:
                        self.groups[tok]=[]
                    self.groups[tok].append(l)
        info_records=self.groups.get('Info:',[])
        logger.debug(f'Number of Info: records: {len(info_records)}')
        if len(info_records)==0:
            return None
        check_tokens=[]
        for _ in range(10):
            tokens=info_records[_].split()
            if len(tokens)>1:
                check_tokens.append(tokens[1])
        logger.debug(f'check_tokens {check_tokens}')
        if not 'NAMD' in check_tokens:
            logger.debug(f'{filename} is evidently not a NAMD log')
            return None
        for info in info_records:
            for il in self.infonames:
                if info[6:].startswith(il):
                    logger.debug(f'Getting data for Info: record name {il} from record {info}')
                    self.info[il]=getinfo(il,info)
            for cn in self.countnames:
                if info.endswith(cn):
                    tokens=info.split()
                    if len(tokens)==3 and tokens[1].isdigit():
                        self.counts[cn]=int(tokens[1])
        tcl_records=self.groups.get('TCL:',[])
        for tclr in tcl_records:
            if tclr[5:].startswith('Running for'):
                tokens=tclr.split()
                self.requested_timesteps=int(tokens[3])
                
        for k,v in self.groups.items():
            logger.debug(f'{filename}: {k} {len(v)} lines')

        self.output_timestep=None
        self.restart_timestep=None
        for rmsg in self.groups['WRITING']:
            tok=rmsg.split()
            if tok[1]=='COORDINATES' and tok[3]=='OUTPUT':
                self.output_timestep=int(tok[-1])
            if tok[1]=='COORDINATES' and tok[3]=='RESTART':
                self.restart_timestep=int(tok[-1])
    
    def energy(self):
        logger.debug(f'energy() call on log from {self.filename}')
        self.etitles=[]
        if 'ETITLE:' not in self.groups or len(self.groups['ETITLE:'])==0:
            logger.debug(f'No ETITLE records found in {self.filename}')
        else:
            self.etitles=self.groups['ETITLE:'][0].split()[1:]
        logger.debug(f'{self.etitles}')
        if not self.etitles:
            if self.inherited_etitles:
                self.etitles=self.inherited_etitles
            else:
                logger.debug(f'{self.filename} has no ETITLE: records and we are not inheriting any from a previous log.')
                return self
        edata={}
        for e in self.etitles:
            edata[e]=[]
        for e in self.groups['ENERGY:']:
            dtype=float
            if e=='TS':
                dtype=int
            ed=e.split()[1:]
            for f,d in zip(self.etitles,ed):
                edata[f].append(dtype(d))
        self.edata=pd.DataFrame(edata)
        if 'VOLUME' in self.edata:
            self.edata['DENSITY']=np.reciprocal(self.edata['VOLUME'])*float(self.info['TOTAL MASS'])
        my_logger(self.edata.head(),logger.debug)
        return self
    
    def success(self):
        return len(self.groups['WallClock:'])==1
    


            