# just so we can maybe do a dynamic install, to compile stuff in c
# not currently used

from setuptools import setup
from distutils.errors import DistutilsSetupError
from distutils import log as distutils_logger
import os
import subprocess

def compile_TcL_modules(sources):
    if sources is None or not isinstance(sources, (list, tuple)):
        raise DistutilsSetupError(f'{sources} must be a list or tuple of *.c files')
    sources = list(sources)
    if len(sources)>1:
        sources_path = os.path.commonpath(sources)
    else:
        sources_path = os.path.dirname(sources[0])
    sources_path = os.path.realpath(sources_path)
    if not sources_path.endswith(os.path.sep):
        sources_path+= os.path.sep

    if not os.path.exists(sources_path) or not os.path.isdir(sources_path):
        raise DistutilsSetupError()

    output_dir = os.path.realpath(os.path.join(sources_path,'..','lib'))
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_lib = 'bondstruct.so'

    distutils_logger.info('Will execute the following command in with subprocess.Popen: \n{0}'.format(
            'make static && mv {0} {1}'.format(output_lib, os.path.join(output_dir, output_lib))))


    make_process = subprocess.Popen(f'make {output_lib} && mv {output_lib} {os.path.join(output_dir, output_lib)}',
                                    cwd=sources_path,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    shell=True)
    stdout, stderr = make_process.communicate()
    distutils_logger.debug(stdout)
    if stderr:
        raise DistutilsSetupError('An ERROR occured while running the '
                                    'Makefile for the {0} resource. '
                                'Error status: {1}'.format(output_lib, stderr))

TcL_module_sources=['pestifer/Resources/tcl/modules/src/bondstruct.c']

setup(name='pestifer')

compile_TcL_modules(TcL_module_sources)