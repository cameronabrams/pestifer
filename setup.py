# just so we can maybe do a dynamics install, to compile stuff in c

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.errors import DistutilsSetupError
from distutils import log as distutils_logger
import os
import subprocess

myso=Extension('pestifier',['tcl/modules/src/bondstruct.c'])

# overloading build_ext to compile the tcl module "bondstruct"
class specialized_build_ext(build_ext, object):
    special_extension = myso.name
    def build_extension(self, ext):
        if ext.name!=self.special_extension:
            # Handle unspecial extensions with the parent class' method
            super(specialized_build_ext, self).build_extension(ext)
        else:
            # Handle special extension
            sources = ext.sources
            if sources is None or not isinstance(sources, (list, tuple)):
                raise DistutilsSetupError(
                       "in 'ext_modules' option (extension '%s'), "
                       "'sources' must be present and must be "
                       "a list of source filenames" % ext.name)
            sources = list(sources)

            if len(sources)>1:
                sources_path = os.path.commonpath(sources)
            else:
                sources_path = os.path.dirname(sources[0])
            sources_path = os.path.realpath(sources_path)
            if not sources_path.endswith(os.path.sep):
                sources_path+= os.path.sep

            if not os.path.exists(sources_path) or not os.path.isdir(sources_path):
                raise DistutilsSetupError(
                       "in 'extensions' option (extension '%s'), "
                       "the supplied 'sources' base dir "
                       "must exist" % ext.name)

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
                                          'Makefile for the {0} library. '
                                        'Error status: {1}'.format(output_lib, stderr))

setup(ext_modules=[myso],cmdclass={'build_ext':specialized_build_ext})