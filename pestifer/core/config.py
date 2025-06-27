# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" The Config class:  parsing the config file and creating the Config object 

    Pestifer's user-configuration input uses ycleptic, an enhanced, YAML-based
    configuration file manager.  The Config object is a descendent of the
    Yclept class.  It also houses the ResourceManager object, which manages
    access to the contents of Pestifer's "resources" directory.

"""

import os
import logging
import shutil
from GPUtil import getGPUs
from importlib.metadata import version
from ycleptic.yclept import Yclept

from .resourcemanager import ResourceManager
from .stringthings import my_logger

logger=logging.getLogger(__name__)
segtype_of_resname={}
charmm_resname_of_pdb_resname={}
res_321={}
res_123={}

class Config(Yclept):
    def __init__(self,userfile='',userdict={},quiet=False):
        vrep=f'ycleptic v. {version("ycleptic")}\npidibble v. {version("pidibble")}'
        self.RM=ResourceManager()
        logger.debug(f'config quiet? {quiet}')
        if not quiet:
            my_logger(vrep,logger.info,just='<',frame='*',fill='',no_indent=True)
            my_logger(str(self.RM),logger.debug,just='<',frame='*',fill='',no_indent=True)
        # resolve full pathname of YCleptic base config for this application
        basefile=self.RM.get_ycleptic_config()
        assert os.path.exists(basefile)
        # ycleptic's init:
        super().__init__(basefile,userfile=userfile,userdict=userdict)
        processor_info=self.processor_info()
        if not quiet:
            my_logger(processor_info,logger.info,just='<',frame='*',fill='')
        self._set_internal_shortcuts()
        self._set_shell_commands(verify_access=(userfile!=''))
        self.RM.update_charmmff(self['user']['charmmff'].get('tarball',''))
        self.RM.update_pdbrepository(self['user']['charmmff'].get('pdbcollections',[]))

    def processor_info(self):
        self.slurmvars={k:os.environ[k] for k in os.environ if 'SLURM' in k}
        self.local_ncpus=os.cpu_count()
        self.gpus_allocated=''
        self.ngpus=0
        self.gpu_devices=''
        retstr=''
        if self.slurmvars and 'SLURM_NNODES' in self.slurmvars and 'SLURM_NTASKS_PER_NODE' in self.slurmvars:
            # we are in a batch execution managed by slurm
            nnodes=int(self.slurmvars['SLURM_NNODES'])
            ntaskspernode=int(self.slurmvars['SLURM_NTASKS_PER_NODE'])
            ncpus=nnodes*ntaskspernode
            retstr+=f'SLURM: {nnodes} nodes; {ncpus} cpus'
            if 'SLURM_JOB_GPUS' in self.slurmvars:
                self.gpu_devices=self.slurmvars['SLURM_JOB_GPUS']
                self.ngpus=len(self.gpus_allocated.split(','))
                ess='s' if self.ngpus>1 else ''
                retstr+f'; {self.ngpus} gpu{ess}'
        else:
            retstr+=f'Local: {self.local_ncpus} cpus'
            ncpus=self.local_ncpus
            gpus=getGPUs()
            if len(gpus)>0:
                self.ngpus=len(gpus)
                self.gpu_devices=','.join([str(x.id) for x in gpus])
                ess='s' if self.ngpus>1 else ''
                retstr+=f'; {self.ngpus} gpu{ess}'

        self.ncpus=ncpus
        return retstr

    def _set_shell_commands(self,verify_access=True):
        required_commands=['charmrun','namd3','vmd','catdcd','packmol']
        command_alternates={'namd3':'namd2'}
        self.shell_commands={}
        for rq in required_commands:
            self.shell_commands[rq]=self['user']['paths'][rq]
            fullpath=shutil.which(self.shell_commands[rq])
            if not fullpath:
                logger.warning(f'{self.shell_commands[rq]}: not found.')
                if rq in command_alternates:
                    rqalt=command_alternates[rq]
                    self.shell_commands[rq]=self['user']['paths'][rqalt]
                    if verify_access:
                        assert shutil.which(self.shell_commands[rq]),f'Alternate command {self.shell_commands[rq]} not found.'
            if verify_access:
                assert os.access(fullpath,os.X_OK),f'You do not have permission to execute {fullpath}'
        self.namd_type=self['user']['namd']['processor-type']
        cn='namd3gpu'
        self.shell_commands[cn]=self['user']['paths'][cn]
        fullpath=shutil.which(self.shell_commands[cn])
        if not fullpath:
            logger.warning(f'{self.shell_commands[cn]}: not found')
        if self.shell_commands['namd3gpu']==self.shell_commands['namd3']:
            logger.warning(f'CPU and GPU namd3 have the same path: {self.shell_commands["namd3gpu"]}')
        if self.namd_type=='gpu': # make sure we can run the GPUResident namd3
            if verify_access:
                assert os.access(fullpath,os.X_OK),f'You do not have permission to execute {fullpath}'
        self.namd_deprecates=self['user']['namd']['deprecated3']

    def _set_internal_shortcuts(self):
        self.progress=len(self.slurmvars)==0
        RM=self.RM
        self.tcl_root=RM.get_tcldir()
        assert os.path.exists(self.tcl_root)
        self.tcl_pkg_path=RM.get_tcl_pkgdir()
        assert os.path.exists(self.tcl_pkg_path)
        self.tcl_script_path=RM.get_tcl_scriptsdir()
        assert os.path.exists(self.tcl_script_path)
        self.vmd_startup_script=os.path.join(self.tcl_root,'vmdrc.tcl')
        assert os.path.exists(self.vmd_startup_script)
        # self.charmmff_toppar_path=RM.get_charmmff_toppardir()
        # assert os.path.exists(self.charmmff_toppar_path)
        self.charmmff_custom_path=RM.get_charmmff_customdir()
        assert os.path.exists(self.charmmff_custom_path)
        self.charmmff_pdbrepository=RM.charmmff_content.pdbrepository
        self.user_charmmff_toppar_path=''
        if hasattr(self,'user'):
            self.user_charmmff_toppar_path=os.path.join(self['user']['charmff'],'toppar')
            assert os.path.exists(self.user_charmmff_toppar_path)
        self.namd_config_defaults=self['user']['namd']
        self.segtypes=self['user']['psfgen']['segtypes']
        for stn,stspec in self.segtypes.items():
            if stspec and 'resnames' in stspec:
                initresnames=stspec['resnames']
                for r in initresnames:
                    if len(r)>4 and r.isalnum():
                        rabbrv=r[:4]
                        # logger.debug(r,rabbrv)
                        if not rabbrv in stspec['resnames']:
                            stspec['resnames'].append(rabbrv)
                            # logger.debug(f'Adding abbreviated resname {rabbrv} to resnames for segtype {stn}')
        self.pdb_to_charmm_resnames={}
        for alias in self['user']['psfgen']['aliases']:
            tok=alias.split()
            if tok[0]=='residue':
                self.pdb_to_charmm_resnames[tok[1]]=tok[2]
        self.segtype_of_resname={}
        for st in self['user']['psfgen']['segtypes']:
            res=self['user']['psfgen']['segtypes'][st].get('resnames',[])
            for r in res:
                self.segtype_of_resname[r]=st

        # globals
        charmm_resname_of_pdb_resname.update(self.pdb_to_charmm_resnames)
        for p,c in charmm_resname_of_pdb_resname.items():
            self.segtype_of_resname[c]=self.segtype_of_resname[p]
        segtype_of_resname.update(self.segtype_of_resname)
        res_123.update(self['user']['psfgen']['segtypes']['protein']['invrescodes'])
        res_321.update(self['user']['psfgen']['segtypes']['protein']['rescodes'])

