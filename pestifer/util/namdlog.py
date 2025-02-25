# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" A class for parsing NAMD log files and xst files
"""

import pandas as pd
import numpy as np
import os
import logging
logger=logging.getLogger(__name__)
from ..stringthings import my_logger

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
        
class NAMDLog:
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
    


            