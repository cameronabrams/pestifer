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
                                        # executable='/bin/bash',
                                        # check=True,
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
                logger.debug(f'Checking version of {pkgname} in non-active environment {env}')
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
                    shell=True,check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
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
