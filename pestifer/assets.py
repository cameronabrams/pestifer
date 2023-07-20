### Makes contents of application Library and Scripts directories accessible

import logging
import os
import importlib.resources
logger=logging.getLogger(__name__)

class Assets:
    def __init__(self,libpackage='Library',scriptpackage='Scripts'):
        self.LibraryPath=''
        self.ScriptsPath=''
        try:
            with importlib.resourses.path(libpackage,'__init__.py') as f:
                self.LibraryPath=os.path.split(os.path.abspath(f))[0]
        except:
            raise ImportError(f'Could not find package {libpackage} in your installation of Pestifer.')
        try:
            with importlib.resourses.path(scriptpackage,'__init__.py') as f:
                self.ScriptsPath=os.path.split(os.path.abspath(f))[0]
        except:
            raise ImportError(f'Could not find package {scriptpackage} in your installation of Pestifer.')


    def info(self):
        msg=f'Library is at {self.LibraryPath}\n'
        msg+=f'Scripts are at {self.ScriptsPath}\n'
        return msg
