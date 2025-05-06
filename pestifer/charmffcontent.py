# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import logging
import os
from .pdbcollection import PDBCollection

logger=logging.getLogger(__name__)

comment_these_out=['set','if','WRNLEV','BOMLEV']
def copy_charmm_top(sysfile,localfile=''):
    if localfile=='':
        localfile=os.path.split(sysfile)[1]
    with open(sysfile,'r') as f:
        lines=f.read().split('\n')
    with open(localfile,'w') as f:
        for l in lines:
            is_comment=any([l.startswith(x) for x in comment_these_out])
            if not is_comment:
                f.write(l+'\n')

class CHARMFFContent:
    """ A class for handling all CHARMM force field content
    """
    def __init__(self,abs_dir):
        self.abs_dir=abs_dir
        self.resource_dir,self.charmmff_dir=os.path.split(abs_dir)
        assert self.charmmff_dir=='charmmff','expected charmmff directory'
        cwd=os.getcwd()
        os.chdir(self.resource_dir)
        self.dirtree={a:[b,c] for a,b,c in os.walk('charmmff')}
        os.chdir(cwd)
        self.charmm_elements=self.dirtree['charmmff'][0]
        logger.debug(f'Charmmff elements: {self.charmm_elements}')
        self.pdb_collection=None
        if 'pdb' in self.charmm_elements:
            self.pdb_collection=PDBCollection(os.path.join(self.resource_dir,'charmmff','pdb'))
        self.custom_files=self.dirtree['charmmff/custom'][1]

        self.par=[x for x in self.dirtree['charmmff/toppar'][1] if x.startswith('par')]
        self.top=[x for x in self.dirtree['charmmff/toppar'][1] if x.startswith('top')]
        self.toppar=[x for x in self.dirtree['charmmff/toppar'][1] if x.startswith('toppar')]
        self.streams=self.dirtree['charmmff/toppar/stream'][0]

    def copy_charmmfile_local(self,basename):
        """ Given a basename for any charmmff file, copy from the existing unaltered charmmff installation to the cwd, and modify it for use directly in NAMD if necessary"""
        if os.path.exists(basename):
            logger.debug(f'{basename} already exists in {os.getcwd()}')
            return basename
        if os.sep in basename:
            # this is a path
            logger.debug(f'expected a basename and got a path {basename}')
            basename=os.path.split(basename)[1]
            logger.debug(f'truncated to basename {basename}')
        if basename.startswith('par_') or basename.startswith('top_'):
            if basename in self.par or basename in self.top: 
                # this is a toplevel parameter or topology file
                src=os.path.join(self.resource_dir,'charmmff/toppar',basename)
                copy_charmm_top(src)
                return basename
        elif basename.startswith('toppar'):
            logger.debug(f'basename {basename} is a stream file')
            if basename in self.toppar:
                # this is a toplevel toppar file
                src=os.path.join(self.resource_dir,'charmmff/toppar',basename)
                copy_charmm_top(src)
                return basename
            else: # this is a stream file
                for stream in self.streams:
                    searchlist=self.dirtree[f'charmmff/toppar/stream/{stream}'][1]
                    logger.debug(f'searching for {basename} in {searchlist}')
                    if basename in searchlist:
                        src=os.path.join(self.resource_dir,'charmmff/toppar/stream',stream,basename)
                        copy_charmm_top(src)
                        return basename
        logger.debug(f'searching for {basename} in {self.custom_files}')
        if basename in self.custom_files:
            src=os.path.join(self.resource_dir,'charmmff/custom',basename)
            copy_charmm_top(src)
            return basename
        logger.warning(f'{basename} not found')
        return None
    
    def clean_local_charmmff_files(self):
        """ remove all local charmmff files
        """
        for f in os.listdir('.'):
            if f.startswith('par_') or f.startswith('top_') or f.startswith('toppar'):
                os.remove(f)