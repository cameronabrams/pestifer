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