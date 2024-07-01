# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#
# Class for handling running external commands
#
import logging
import subprocess
import os
from io import StringIO
import json
import shutil
from glob import glob

logger=logging.getLogger(__name__)

class CondaCheck:
    """Class for interfacing with conda environments"""
    def __init__(self):
        self.conda_exe=os.environ.get('CONDA_EXE',None)
        assert os.access(self.conda_exe,os.X_OK)
        if not self.conda_exe:
            logger.info(f'Shell environment variable CONDA_EXE is not set')
            logger.info(f'No conda found')
        else:
            logger.debug(f'conda executable: {self.conda_exe}')
            check_result=subprocess.run(f'{self.conda_exe} info --json',
                                        shell=True, 
                                        executable='/bin/bash',
                                        check=True,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)
            stdout=check_result.stdout.decode('utf-8')
            with StringIO(stdout) as f:
                self.conda_info=json.load(f)

            self.conda_root=self.conda_info['root_prefix']
            self.conda_root_writable=self.conda_info['root_writable']
            self.active_env=self.conda_info.get('active_prefix_name',None)
            if not self.active_env:
                logger.debug(f'You are not running {__package__} in a conda environment, but conda was detected on your system.')

            conda_envs=self.conda_info.get('envs',[])
            self.conda_envs=['base']
            if len(conda_envs)>1:
                for env in conda_envs[1:]:
                    self.conda_envs.append(os.path.split(env)[-1])
            self.init_shell=os.path.join(f'{self.conda_root}','etc','profile.d','conda.sh')
            assert os.path.exists(self.init_shell)

    def env_lib_dir(self,env):
        if env=='base':
            edir=os.path.join(self.conda_root,'lib')
        else:
            edir=os.path.join(self.conda_root,'envs',env,'lib')
        return edir

    def env_python_dir(self,env):
        libdir=self.env_lib_dir(env)
        trial=glob(os.path.join(libdir,'python3*'))
        for t in trial:
            if not os.path.islink(t):
                break
        return t        

    def env_python_site_packages_dir(self,env):
        pdir=self.env_python_dir(env)
        pspdir=os.path.join(pdir,'site-packages')
        assert os.path.isdir(pspdir),f'No dir {pspdir} found'
        return pspdir

    # def check_packmol_memgen_pdbremix(self,env):
    #     pspdir=self.env_python_site_packages_dir(env)
    #     pmlib=os.path.join(pspdir,'packmol_memgen','lib','pdbremix')
    #     if os.path.exists(pmlib):
    #         the_bad_file=os.path.join(pmlib,'v3numpy.py')
    #         assert os.path.exists(the_bad_file),f'{the_bad_file}: not found.'
    #         with open(the_bad_file,'r') as f:
    #             file_contents=f.read()
    #         pos_result=file_contents.find('dtype=np.float')
    #         neg_result=file_contents.find('dtype=np.float64')
    #         if pos_result!=-1 and neg_result==-1:
    #             logger.debug(f'venv {env} has an unpatched version of packmol-memgen/pdbremix')
    #             logger.debug(f'You should update this environment to support ambertools 23.6 or better')
    #             return 'UNPATCHED'
    #         return 'PATCHED'
    #     return 'UNFOUND'

    # def patch_packmol_memgen_pdbremix(self,env):
    #     if self.check_packmol_memgen_pdbremix(env)=='UNPATCHED':
    #         pspdir=self.env_python_site_packages_dir(env)
    #         pmlib=os.path.join(pspdir,'packmol_memgen','lib','pdbremix')
    #         CWD=os.getcwd()
    #         os.chdir(pmlib)
    #         the_bad_file='v3numpy.py'
    #         with open(the_bad_file,'r') as f:
    #             file_contents=f.read()
    #         patched_result=file_contents.replace('dtype=np.float','dtype=np.float64')
    #         shutil.copy(the_bad_file,'v3numpy.py_orig')
    #         with open(the_bad_file,'w') as f:
    #             f.write('# patched by pestifer\n')
    #             f.write(patched_result)
    #         os.chdir(CWD)
    #         return True
    #     return False

    # def restore_packmol_memgen_pdbremix(self,env):
    #     if self.check_packmol_memgen_pdbremix(env)=='PATCHED':
    #         pspdir=self.env_python_site_packages_dir(env)
    #         pmlib=os.path.join(pspdir,'packmol_memgen','lib','pdbremix')
    #         CWD=os.getcwd()
    #         os.chdir(pmlib)
    #         assert os.path.exists('v3numpy.py_orig')
    #         shutil.move('v3numpy.py_orig','v3numpy.py')
    #         os.chdir(CWD)

    def check_envs(self):
        result={'good':[],'bad':[],'none':[]}
        for e in self.conda_envs:
            check_result=self.check_packmol_memgen_pdbremix(e)
            if check_result=='UNPATCHED':
                result['bad'].append(e)
            elif check_result=='PATCHED':
                result['good'].append(e)
            else:
                result['none'].append(e)
        return result
    
    def info(self):
        if not self.conda_root:
            return f'No conda environments available.  Some functionality may not be available.'
        return f'Conda root: {self.conda_root}\nActive env: {self.active_env}\nEnvs: {self.conda_envs}\nInit shell: {self.init_shell}'

    def env_exists(self,envname):
        return envname in self.conda_envs
    
    def get_package_version(self,envname,pkgname,from_list=False):
        try:
            if envname!=self.active_env:
                logging.debug(f'Checking version of {pkgname} in non-active environment {envname}')
                if not from_list:
                    pcs=subprocess.Popen(f"source {self.init_shell}\nconda activate {envname}\npython -c 'import {pkgname}; print({pkgname}.__version__)'",
                    shell=True, executable='/bin/bash',stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
                    stdout,stderr=pcs.communicate()
                    results=stdout.split('\n')
                    return results[0]
                else:
                    c=f'source {self.init_shell}\nconda activate {envname}\nconda list {pkgname}'
                    logger.debug(f'Issuing {c}')
                    check_result=subprocess.run(c,shell=True,executable='/bin/bash',check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                    results=check_result.stdout.decode('utf-8').split('\n')
                    err=check_result.stderr.decode('utf-8').split('\n')
                    if results[-2][0]!='#':
                        version=results[-2].split()[1]
                        return version
                    else:
                        logger.debug(f'so {results}')
                        logger.debug(f'se {err}')
                        logger.debug(f'Could not determine version of {pkgname} in env {envname} via "conda list"')
                        return None
            else:
                if not from_list:
                    from importlib.metadata import version
                    return version(pkgname)
                else:
                    check_result=subprocess.run(f'conda list {pkgname}',
                    shell=True, executable='/bin/bash', check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                    results=check_result.stdout.decode('utf-8').split('\n')
                    if results[-2][0]!='#':
                        version=results[-2].split()[1]
                        return version
                    else:
                        logger.debug(f'Could not determine version of {pkgname} in env {envname} via "conda list"')
                        return None                    
        except Exception as e:
            logger.debug(f'{e}')
            logger.debug(f'Could not determine version of {pkgname} in env {envname}')
            return None

    def condafy(self,command,env=None):
        if not env or not self.env_exists(env) or env==self.active_env:
            return command
        commands=command.split('\n')
        if len(commands)>1:
            if commands[0]==f'source {self.init_shell}':
                logger.debug(f'You are trying to condafy a command list that is already condafied.')
                if commands[1]==f'conda activate {env}':
                    logger.debug(f'You have also specified the environment that is already specified in this condafied command.')
                else:
                    commands[1]=f'conda activate {env}'
        else:
            commands=[f'source {self.init_shell}',f'conda activate {env}']+commands 
        return '\n'.join(commands)
    
class Command:
    divider_line_length=55
    def __init__(self,command:str,*args,**options):
        self.command=command
        self.args=args
        self.options=options
        self.c=f'{self.command} '+' '.join(args)+' '.join([f'-{k} {v}' for k,v in self.options.items()])
        self.stdout=''
        self.stderr=''
        self.is_condafied=False

    def condarun(self,env=None,Condaspec=None,**kwargs):
        if Condaspec!=None:
            self.c=Condaspec.condafy(self.c,env=env)
            self.is_condafied=True
            logger.debug(f'Condafied to {self.c}')
        return self.run(override=kwargs.get('override',()),ignore_codes=kwargs.get('ignore_codes',[]),quiet=kwargs.get('quiet',True),progress=kwargs.get('progress',False))

    def run(self,override=(),ignore_codes=[],quiet=True,progress=False):
        if not quiet:
            logger.debug(f'{self.c}')
        # This runs the command
        process=subprocess.Popen(self.c,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        # Here is where I can institute a progress indicator
        # if progress:
        #     pass # not yet implemented
        # else:
        # # This holds execution of Python until the command exits
        self.stdout,self.stderr=process.communicate()
        if process.returncode!=0 and not process.returncode in ignore_codes:
            logger.error(f'Returncode: {process.returncode}')
            if len(self.stdout)>0:
                logger.error('stdout buffer follows\n'+'*'*self.divider_line_length+'\n'+self.stdout+'\n'+'*'*self.divider_line_length)
            if len(self.stderr)>0:
                logger.error('stderr buffer follows\n'+'*'*self.divider_line_length+'\n'+self.stderr+'\n'+'*'*self.divider_line_length)
            raise subprocess.SubprocessError(f'Command "{self.c}" failed with returncode {process.returncode}')
        else:
            # logger.info(f'Returncode: {process.returncode}.')
            if len(override)==2:
                needle,msg=override
                if needle in self.stdout or needle in self.stderr:
                    logger.info(f'Returncode: {process.returncode}, but another error was detected:')
                    logger.error(msg)
                    if len(self.stdout)>0:
                        logger.error('stdout buffer follows\n'+'*'*self.divider_line_length+'\n'+self.stdout+'\n'+'*'*self.divider_line_length)
                    if len(self.stderr)>0:
                        logger.error('stderr buffer follows\n'+'*'*self.divider_line_length+'\n'+self.stderr+'\n'+'*'*self.divider_line_length)