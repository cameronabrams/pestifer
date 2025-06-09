# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import logging
import os
import tarfile
from .pdbrepository import PDBRepository


logger=logging.getLogger(__name__)

comment_these_out=['set','if','WRNLEV','BOMLEV']

class CHARMMFFContent:
    """ A class for handling all CHARMM force field content.  The CHARMM force field is
    stored in a tarball downloaded directly from the MacKerell lab at the University of Michigan.
    https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/toppar_c36_jul24.tgz
    """
    def __init__(self,charmmff_path='',tarfilename='toppar_c36_jul24.tgz'):
        self.tarfile=None
        self.tarfilename=tarfilename
        if not charmmff_path:
            charmmff_path='.'
        if not os.path.isdir(charmmff_path):
            raise NotADirectoryError(f'Expected a directory at {charmmff_path}, but it is not a directory')
        self.charmmff_path=os.path.abspath(charmmff_path)
        self.basename=os.path.basename(self.charmmff_path)
        self.parent_path=os.path.dirname(self.charmmff_path)
        cwd=os.getcwd()
        os.chdir(self.parent_path)
        self.dirtree={a:[b,c] for a,b,c in os.walk(self.basename)}
        os.chdir(cwd)
        self.charmm_elements=self.dirtree[self.basename][0]
        logger.debug(f'Charmmff elements: {self.charmm_elements}')
        self.pdb_repository=None
        if 'pdbrepository' in self.charmm_elements:
            self.pdb_repository=PDBRepository(os.path.join(self.charmmff_path,'pdbrepository'))
        self.load_charmmff(tarfilename)
        self.custom_files=self.dirtree[f'{self.basename}/custom'][1]
        for f in self.custom_files:
            assert f not in self.filenamemap, f'custom file {f} already exists in filenamemap'
            self.filenamemap[f]=os.path.join(self.charmmff_path,'custom',f)

    def __del__(self):
        self.tarfile.close()

    def load_charmmff(self,tarfilename='toppar_c36_jul24.tgz'):
        def okfilename(name):
            """ Check if a filename is ok to use in the tarfile """
            return not 'history' in name and (name.endswith('.str') or name.endswith('.prm') or name.endswith('.rtf'))

        """ Load the CHARMM force field tarball from the specified path """
        if not os.path.exists(os.path.join(self.charmmff_path,tarfilename)):
            raise FileNotFoundError(f'CHARMM force field tarball {tarfilename} not found in {self.charmmff_path}')
        if self.tarfile is not None:
            self.tarfile.close()
        logger.debug(f'Loading CHARMM force field tarball {tarfilename} from {self.charmmff_path}')
        self.tarfile=tarfile.open(os.path.join(self.charmmff_path,tarfilename),'r:gz')
        self.tarmembers=self.tarfile.getmembers()
        self.toplevel_par={os.path.basename(x.name):x.name for x in self.tarmembers if x.isfile() and x.name.startswith('toppar/par_') and okfilename(x.name)}
        self.toplevel_top={os.path.basename(x.name):x.name for x in self.tarmembers if x.isfile() and x.name.startswith('toppar/top_') and okfilename(x.name)}
        self.toplevel_toppar={os.path.basename(x.name):x.name for x in self.tarmembers if x.isfile() and x.name.startswith('toppar/toppar_') and okfilename(x.name)}
        self.filenamemap={**self.toplevel_par,**self.toplevel_top,**self.toplevel_toppar}
        self.streams=[os.path.basename(x.name) for x in self.tarmembers if x.isdir() and x.name.startswith('toppar/stream/') and 'history' not in x.name]
        self.streamfiles={}
        for stream in self.streams:
            self.streamfiles[stream]={os.path.basename(x.name):x.name 
                                      for x in self.tarmembers if x.isfile() 
                                      and x.name.startswith(f'toppar/stream/{stream}/') and okfilename(x.name)}
            self.filenamemap.update(self.streamfiles[stream])
        check_basenames=list(self.filenamemap.keys())
        assert len(check_basenames)==len(set(check_basenames)),f'found duplicate basenames in charmmff tarball: {check_basenames}'

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
                lines=file.readlines()
                logger.debug(f'found {len(lines)} lines in {basename} in custom files')
                with open(basename,'w') as f:
                    for l in lines:
                        is_comment=any([l.startswith(x) for x in comment_these_out])
                        if not is_comment:
                            f.write(l+'\n')
                        else:
                            f.write('! commented out by pestifer:\n')
                            f.write(f'! {l}\n')
        elif basename in self.filenamemap:
            logger.debug(f'found {basename} in at {self.filenamemap[basename]} in tarball')
            with self.tarfile.extractfile(self.filenamemap[basename]) as file:
                logger.debug(f' file has type {type(file)}')
                lines=file.read().decode().splitlines()
                logger.debug(f'type of lines is {type(lines)}')
                logger.debug(f'found {len(lines)} lines in {basename} in tarfile')
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
            if f.startswith('par') or f.startswith('top_') or f.startswith('toppar') or f.startswith('charmm') or f.endswith('.str') or f.endswith('.prm') or f.endswith('.rtf'):
                os.remove(f)