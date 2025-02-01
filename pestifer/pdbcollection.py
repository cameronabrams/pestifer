# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#
# Manages the collection of PDB files used as inputs for packmol
#
# Collections are subdivided by charmmff 'streams'
# 
import os
import glob
import logging
logger=logging.getLogger(__name__)

class PDBCollection:
    _user_count=0
    def __init__(self,basepath=''):
        self.collections={}
        self.registercollection(basepath,'base')
        logger.debug(str(self.collections['base']))
    
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

    def search_collection(self,collection_name,name):
        c=self.collections.get(collection_name,{})
        s=c.get('streams',{})
        for k,v in s.items():
            if name in v:
                fullname=os.path.join(c['path'],k,name)
                if os.path.isdir(fullname):
                    return glob.glob(os.path.join(fullname,f'{name}-??.pdb'))
                else:
                    return os.path.join(f'{fullname}.pdb')
        return None

    def get_pdb_path(self,name):
        for colname in self.collections.keys():
            res=self.search_collection(colname,name)
            if res:
                return res
        return None


