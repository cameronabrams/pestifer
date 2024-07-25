# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#
# Class for handling running external commands
#
import logging
import subprocess
import os
import atexit
import signal
from io import StringIO
import json
import shutil
from glob import glob
# from .progress import PestiferProgress
import time
from .stringthings import ByteCollector

logger=logging.getLogger(__name__)

class CondaCheck:
    """Class for interfacing with conda environments"""
    def __init__(self):
        self.conda_env={k:os.environ[k] for k in os.environ if "CONDA" in k}
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
    
    def info(self):
        if not self.conda_root:
            return f'No conda environments available.  Some functionality may not be available.'
        return f'Conda root: {self.conda_root}\nActive env: {self.active_env}'
        # return f'Conda root: {self.conda_root}\nActive env: {self.active_env}\nEnvs: {self.conda_envs}\nInit shell: {self.init_shell}'

    def env_exists(self,envname):
        return envname in self.conda_envs
    
    def get_package_version(self,pkgname,env=None,from_list=False):
        if not env: env=self.active_env
        try:
            if env!=self.active_env:
                logging.debug(f'Checking version of {pkgname} in non-active environment {env}')
                if not from_list:
                    pcs=subprocess.Popen(f"source {self.init_shell}\nconda activate {env}\npython -c 'import {pkgname}; print({pkgname}.__version__)'",
                    shell=True, executable='/bin/bash',stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
                    stdout,stderr=pcs.communicate()
                    results=stdout.split('\n')
                    return results[0]
                else:
                    c=f'source {self.init_shell}\nconda activate {env}\nconda list {pkgname}'
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
                        logger.debug(f'Could not determine version of {pkgname} in env {env} via "conda list"')
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
                        logger.debug(f'Could not determine version of {pkgname} in env {env} via "conda list"')
                        return None                    
        except Exception as e:
            logger.debug(f'{e}')
            logger.debug(f'Could not determine version of {pkgname} in env {env}')
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
        return self.run(override=kwargs.get('override',()),ignore_codes=kwargs.get('ignore_codes',[]),quiet=kwargs.get('quiet',True),progress=kwargs.get('progress',None),logfile=kwargs.get('logfile',None))

    def run(self,logfile=None,override=(),ignore_codes=[],quiet=True,progress=None,**kwargs):
        """ runs this Command instance
        
        Parameters:
        -----------

        logfile : str
            name of log file to write process' stdout and stderr. If None (default), stdout and stderr are retained
            only in this instance's stdout and stderr attributes.
        override : tuple
            tuple composed of a "needle" and a "message".  If the needle is found in the stdout or stderr of
            the process, the message is displayed and an error is thrown, halting the program.
        ignore_codes : list
            a list of integer exit codes that are ignored in addition to 0
        quiet : boolean
        progress: PestiferProgress instance (or descendant) 
            used for progress bar/elapsed time displays

        """
        if not quiet:
            logger.debug(f'{self.c}')
        log=None
        if not logfile:
            logger.debug(f'no logfile specified for {self.c}')
        else:
            if os.path.exists(logfile) and not kwargs.get('overwrite_logs',False):
                nlogs=len(glob(f'%{logfile}'))
                shutil.move(logfile,f'%{logfile}-{nlogs+1}%')
            log=open(logfile,'w')
        if progress and not progress.unmeasured:
            bytes=ByteCollector()
        process=subprocess.Popen(self.c,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        global pid
        pid=process.pid
        def kill_child():
            if pid is None:
                pass
            else:
                try:
                    os.kill(pid, signal.SIGTERM)
                except:
                    pass
        atexit.register(kill_child)
        otally=0
        ttally=0
        t2tally=0
        ncalls=0
        self.stdout=''
        self.stderr=''
        while True:
            b=time.process_time_ns()
            output=process.stdout.readline()
            self.stdout+=output
            e=time.process_time_ns()
            ncalls+=1
            ttally+=(e-b)/1.e6
            otally+=len(output)
            if log: log.write(output)
            if output=='' and process.poll() is not None:
                break
            if progress:
                if not progress.unmeasured:
                    bytes.write(output)
                    progress.go(bytes.byte_collector)
                else:
                    b=time.process_time_ns()
                    progress.go()
                    e=time.process_time_ns()
                    t2tally+=(e-b)/1.e6
                    # logger.debug(f'{len(output):>5d} B {otally:>9,d} B {ttally/ncalls:>12.6f} ms/readline call {t2tally/ncalls:>12.6f} ms/pbar call')
        if progress:
            print()
        if logfile:
            logger.debug(f'Log written to {logfile}')
        remaining_stdout,self.stderr=process.communicate()
        # logger.debug(f'stdout [{self.stdout}]')
        # logger.debug(f'remaining stdout [{remaining_stdout}]')
        # logger.debug(f'stderr [{self.stderr}]')
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
                    if len(self.stdout)>0 and needle in self.stdout:
                        logger.error('stdout buffer follows\n'+'*'*self.divider_line_length+'\n'+self.stdout+'\n'+'*'*self.divider_line_length)
                    if len(self.stderr)>0 and needle in self.stderr:
                        logger.error('stderr buffer follows\n'+'*'*self.divider_line_length+'\n'+self.stderr+'\n'+'*'*self.divider_line_length)