# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import pandas as pd
import numpy as np
import os
import logging
logger=logging.getLogger(__name__)
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


class NAMDLog:
    groupnames=['ETITLE:','ENERGY:','charmrun>','Charm++>','Info:','TIMING:','WallClock:']
    infonames=['TIMESTEP','NUMBER OF STEPS','RANDOM NUMBER SEED','TOTAL MASS']
    def __init__(self,filename,inherited_etitles=[]):
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
        for l in rawlines:
            if len(l)>0:
                tok=l.split()[0]
                if tok in self.groupnames:
                    if not tok in self.groups:
                        self.groups[tok]=[]
                    self.groups[tok].append(l)
        logger.debug(f'Number of Info: records: {len(self.groups.get("Info:",[]))}')
        for info in self.groups.get('Info:',[]):
            for il in self.infonames:
                if il in info:
                    logger.debug(f'Getting data for Info: record name {il} from record {info}')
                    self.info[il]=getinfo(il,info)

        for k,v in self.groups.items():
            logger.debug(f'{filename}: {k} {len(v)} lines')
    
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
        return self
    
    def success(self):
        return len(self.groups['WallClock:'])==1
    
class NAMDProgress:
    groupnames=['ENERGY:','Info:','TCL:']
    infonames=['FIRST TIMESTEP']
    tclnames=['Running for','Minimizing for']
    def __init__(self):
        self.groups={}
        self.info={}
        self.tcl={}
    def parse_f(self,logstring):
        self.groups={}
        self.info={}
        self.tcl={}
        loglines=logstring.split('\n')
        for l in loglines:
            # print(l)
            if len(l)>0:
                tok=l.split()[0]
                if tok in self.groupnames:
                    if not tok in self.groups:
                        self.groups[tok]=[]
                    # if tok=='TCL:': logger.debug(f'TCL appends {l} {"Minimizing for" in l}')
                    self.groups[tok].append(l)
                # else:
                #     logger.debug(f'Ignoring {tok}')
        for info in self.groups.get('Info:',[]):
            for il in self.infonames:
                if il in info:
                    self.info[il]=getinfo(il,info)
        for tclline in self.groups.get('TCL:',[]):
            for tl in self.tclnames:
                if tl in tclline:
                    # logger.debug(f'parsing {tl} {tclline}')
                    self.tcl[tl]=gettclinfo(tl,tclline)
    def init_f(self,logstring):
        if not logstring:
            return False
        self.parse_f(logstring)
        first_step=self.info.get('FIRST TIMESTEP',None)
        # logger.debug(f'tcl lines: {self.tcl}')
        num_steps=self.tcl.get('Running for',None)
        nummin_steps=self.tcl.get('Minimizing for',None)
        # logger.debug(f'nummin_steps {nummin_steps}')
        if nummin_steps:
            result={'first_step':0}
            result['num_steps']=nummin_steps
            # logger.debug(f'result {result}')
            return result
        elif num_steps:
            if first_step:
                result={'first_step':int(first_step)}
                result['num_steps']=num_steps
                return result
        return False
    def meas_f(self,logstring):
        if not logstring:
            return False
        self.parse_f(logstring)
        if 'ENERGY:' in self.groups and len(self.groups['ENERGY:'])>0:
            last_line=self.groups['ENERGY:'][-1]
            if len(last_line)<20:
                return False
            tok=last_line.split()
            if len(tok)>1:
                current_time_step=int(tok[1])
                # print(f'\ntook {current_time_step}\n')
                return dict(current_time_step=current_time_step)
        return False
    def comp_f(self,init_obj,meas_obj):
        if not init_obj or not meas_obj:
            return 0
        fac=(meas_obj['current_time_step']-init_obj['first_step'])/init_obj['num_steps']
        assert fac<=1.0,f'error: {meas_obj} {init_obj}'
        return fac

            