# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import logging
import os
import tarfile
from .pdbcollection import PDBCollection

logger=logging.getLogger(__name__)

comment_these_out=['set','if','WRNLEV','BOMLEV']

class CHARMMFFContent:
    """ A class for handling all CHARMM force field content.  The CHARMM force field is
    stored in a tarbal downloaded directly from the MacKerell lab at the University of Michigan.
    https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/toppar_c36_jul24.tgz
    """
    def __init__(self,abs_dir):
        self.abs_dir=abs_dir
        self.resource_dir,self.charmmff_dir=os.path.split(abs_dir)
        assert self.charmmff_dir=='charmmff','expected charmmff directory'
        cwd=os.getcwd()
        os.chdir(self.resource_dir)
        self.tarfile=tarfile.open(os.path.join('charmmff','toppar_c36_jul24.tgz'),'r:gz')
        self.tarmembers=self.tarfile.getmembers()
        self.toplevel_par={os.path.basename(x.name):x.name for x in self.tarmembers if x.isfile() and x.name.startswith('toppar/par_') and 'history' not in x.name}
        self.toplevel_top={os.path.basename(x.name):x.name for x in self.tarmembers if x.isfile() and x.name.startswith('toppar/top_') and 'history' not in x.name}
        self.toplevel_toppar={os.path.basename(x.name):x.name for x in self.tarmembers if x.isfile() and x.name.startswith('toppar/') and 'history' not in x.name}
        self.filenamemap={**self.toplevel_par,**self.toplevel_top,**self.toplevel_toppar}
        self.streams=[os.path.basename(x.name) for x in self.tarmembers if x.isdir() and x.name.startswith('toppar/stream/') and 'history' not in x.name]
        self.streamfiles={}
        for stream in self.streams:
            self.streamfiles[stream]={os.path.basename(x.name):x.name 
                                      for x in self.tarmembers if x.isfile() 
                                      and x.name.startswith(f'toppar/stream/{stream}/') and 'history' not in x.name}
            self.filenamemap.update(self.streamfiles[stream])
        check_basenames=list(self.filenamemap.keys())
        assert len(check_basenames)==len(set(check_basenames)),f'found duplicate basenames in charmmff tarball: {check_basenames}'
        self.dirtree={a:[b,c] for a,b,c in os.walk('charmmff')}
        os.chdir(cwd)
        self.charmm_elements=self.dirtree['charmmff'][0]
        logger.debug(f'Charmmff elements: {self.charmm_elements}')
        self.pdb_collection=None
        if 'pdb' in self.charmm_elements:
            self.pdb_collection=PDBCollection(os.path.join(self.resource_dir,self.charmmff_dir,'pdb'))
        self.custom_files=self.dirtree[f'{self.charmmff_dir}/custom'][1]
        for f in self.custom_files:
            assert f not in self.filenamemap, f'custom file {f} already exists in filenamemap'
            self.filenamemap[f]=os.path.join(self.abs_dir,'custom',f)

    def __del__(self):
        self.tarfile.close()

    def copy_charmmfile_local(self,basename):
        """ Given a basename for any charmmff file, extract from the existing unaltered charmmff installation
            and modify it for use directly in NAMD if necessary"""
        if os.path.exists(basename):
            logger.debug(f'{basename} already exists in {os.getcwd()}')
            return basename
        if os.sep in basename:
            # this is a path
            logger.debug(f'expected a basename and got a path {basename}')
            basename=os.path.split(basename)[1]
            logger.debug(f'truncated to basename {basename}')
        if basename in self.custom_files:
            with open(self.filenamemap[basename]) as file:
                lines=file.read().split('\n')
                with open(basename,'w') as f:
                    for l in lines:
                        is_comment=any([l.startswith(x) for x in comment_these_out])
                        if not is_comment:
                            f.write(l+'\n')
                        else:
                            f.write('! commented out by pestifer:\n')
                            f.write(f'! {l}\n')
        # logger.debug(f'searching for {basename} in {self.custom_files}')
        elif basename in self.filenamemap:
            with self.tarfile.extractfile(self.filenamemap[basename]) as file:
                lines=file.read().decode('utf-8').split('\n')
                with open(basename,'w') as f:
                    for l in lines:
                        is_comment=any([l.startswith(x) for x in comment_these_out])
                        if not is_comment:
                            f.write(l+'\n')
                        else:
                            f.write('! commented out by pestifer:\n')
                            f.write(f'! {l}\n')
        else:
            logger.warning(f'{basename} not found in charmmff')
        return basename
    
    def clean_local_charmmff_files(self):
        """ remove all local charmmff files
        """
        for f in os.listdir('.'):
            if f.startswith('par_') or f.startswith('top_') or f.startswith('toppar'):
                os.remove(f)