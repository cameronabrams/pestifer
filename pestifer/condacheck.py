# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#
# Checks for conda environments
import subprocess
import os
import logging
logger=logging.getLogger(__name__)

class CondaCheck:
    def __init__(self):
        self.has_conda=os.environ.get('CONDA_EXE',None)
        self.current_conda_env=None
        if not self.has_conda:
            logger.debug(f'No conda detected')
        else:
            self.conda_prefix=os.environ.get('CONDA_PREFIX',None)
            if not self.conda_prefix:
                self.conda_prefix='/'.join(self.has_conda.split('/')[:-2])
                logger.debug('You are not running pestifer in a conda environment, but conda was detected on your system.')
            else:
                if 'envs' in self.conda_prefix:
                    self.current_conda_env=self.conda_prefix.split('/')[-1]
                else:
                    self.current_conda_env='base'
                logger.debug(f'Pestifer is running in conda environment "{self.current_conda_env}"')

            self.base_conda_prefix=os.environ.get('CONDA_PREFIX_1',self.conda_prefix)
            self.init_shell=f'{self.base_conda_prefix}/etc/profile.d/conda.sh'
            check_result=subprocess.run(f"""source {self.init_shell}
            conda env list""",
            shell=True, executable='/bin/bash', check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            lines=check_result.stdout.decode('utf-8').split('\n')
            self.conda_envs=[]
            for l in lines:
                if len(l)>0 and l[0]!='#':
                    tok=l.split()
                    self.conda_envs.append(tok)
    
    def env_exists(self,envname):
        return envname in [x[0] for x in self.conda_envs]
    
    def get_package_version(self,envname,pkgname,from_list=False):
        try:
            if envname!=self.current_conda_env:
                if not from_list:
                    check_result=subprocess.run(f"""source {self.init_shell}
                    conda activate {envname}
                    python -c 'import {pkgname}; print({pkgname}.__version__)'""",
                    shell=True, executable='/bin/bash', check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                    results=check_result.stdout.decode('utf-8').split('\n')
                    return results[0]
                else:
                    check_result=subprocess.run(f"""source {self.init_shell}
                    conda activate {envname}
                    conda list {pkgname}""",
                    shell=True, executable='/bin/bash', check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                    results=check_result.stdout.decode('utf-8').split('\n')
                    if results[-2][0]!='#':
                        version=results[-2].split()[1]
                        return version
                    else:
                        logger.debug(f'Could not determine version of {pkgname} in env {envname} via "conda list"')
                        return None
            else:
                if not from_list:
                    from importlib.metadata import version
                    return version(pkgname)
                else:
                    check_result=subprocess.run(f"""conda list {pkgname}""",
                    shell=True, executable='/bin/bash', check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                    results=check_result.stdout.decode('utf-8').split('\n')
                    if results[-2][0]!='#':
                        version=results[-2].split()[1]
                        return version
                    else:
                        logger.debug(f'Could not determine version of {pkgname} in env {envname} via "conda list"')
                        return None                    
        except:
            logger.debug(f'Could not determine version of {pkgname} in env {envname}')
            return None