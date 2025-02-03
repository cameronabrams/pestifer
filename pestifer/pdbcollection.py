# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#
# Manages the collection of PDB files used as inputs for packmol
#
# Collections are subdivided by charmmff 'streams'
# 
import os
import glob
import logging
import shutil
import yaml
logger=logging.getLogger(__name__)

class PDBInput:
    def __init__(self,dir=None,solo=None,noh=False):
        if dir:
            if noh:
                self.conformers=glob.glob(os.path.join(dir,f'*-noh.pdb'))
            else:
                self.conformers=[x for x in glob.glob(os.path.join(dir,f'*.pdb')) if not 'noh' in x]
            self.conformers.sort()
            info=os.path.join(dir,'info.yaml')
            with open(info,'r') as f:
                self.info=yaml.safe_load(f)
        elif solo:
            self.conformers=[solo]
            self.info={}
    
    def checkout(self,index=0):
        basename=os.path.basename(self.conformers[index])
        if not os.path.exists(basename):
            shutil.copy(self.conformers[index],basename)
        return basename

    def get_charge(self):
        return self.info.get('charge',None)
    
    def get_ref_length(self,index=0):
        conformers=self.info.get('conformers',[])
        basename=os.path.basename(self.conformers[index])
        for c in conformers:
            if c['pdb']==basename:
                return c['head-tail-length']
        return 0.0

    def get_ref_atoms(self):
        return self.info.get('reference-atoms',{})
    
    def get_parameters(self):
        return self.info.get('parameters',[])

class PDBCollection:
    _user_count=0
    def __init__(self,basepath=''):
        self.collections={}
        self.registercollection(basepath,'base')
        # logger.debug(str(self.collections['base']))
    
    def registercollection(self,pathname,pathkey):
        self.collections[pathkey]={}
        self.collections[pathkey]['path']=pathname
        if os.path.isdir(pathname):
            streams_full=[x for x in glob.glob(os.path.join(pathname,'*')) if os.path.isdir(x)]
        else:
            logger.warning(f'{pathname} is not found.')
            return
        self.collections[pathkey]['streams']={os.path.basename(x):[] for x in streams_full}
        for stream in self.collections[pathkey]['streams']:
            streamdir=os.path.join(pathname,stream)
            pdbdirs_full=[x for x in glob.glob(os.path.join(streamdir,'*')) if os.path.isdir(x)]
            raw_pdbs=[os.path.splitext(os.path.basename(x))[0] for x in glob.glob(os.path.join(streamdir,'*')) if x.endswith('.pdb')]
            raw_pdbs.sort()
            pdbdirs=[os.path.basename(x) for x in pdbdirs_full]
            pdbdirs.sort()
            self.collections[pathkey]['streams'][stream]=pdbdirs if len(pdbdirs)>0 else raw_pdbs

    def add_usercollection(self,userpath=''):
        self._user_count+=1
        self.registercollection(userpath,f'user{self._user_count}')

    def search_collection(self,collection_name,name,noh=False):
        c=self.collections.get(collection_name,{})
        s=c.get('streams',{})
        for k,v in s.items():
            if name in v:
                fullname=os.path.join(c['path'],k,name)
                if os.path.isdir(fullname):
                    return PDBInput(dir=fullname,noh=noh)
                else:
                    return PDBInput(solo=os.path.join(f'{fullname}.pdb'))
        return None

    def get_pdb(self,name,noh=False):
        for colname in self.collections.keys():
            res=self.search_collection(colname,name,noh=noh)
            if res:
                return res
        return None


