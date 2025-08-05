# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import datetime
import logging
import os

from .filewriter import Filewriter
from ..core.command import Command
from ..core.stringthings import ByteCollector, FileCollector
from ..util.logparsers import PackmolLog
from ..util.progress import PackmolProgress

logger = logging.getLogger(__name__)
class PackmolScripter(Filewriter):
    """
    This class extends the Filewriter class to provide functionality for creating and managing Packmol scripts.
    
    Parameters
    ----------
    config : Config
        The configuration object containing settings for the script.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(comment_char=kwargs.get('comment_char', '#'))
        self.indent = 4 * ' '
        self.progress = kwargs.get('progress', PackmolProgress())
        self.packmol = kwargs.get('packmol')
        self.F = FileCollector()
        self.default_ext = '.inp'
        self.default_script = f'packmol{self.default_ext}'
        self.scriptname = self.default_script

    def newscript(self,basename=None):
        """
        Create a new Packmol input script.

        Parameters
        ----------
        basename : str, optional
            The base name for the script file (without extension). If not provided, a default name will be used.
        """
        timestampstr=datetime.datetime.today().ctime()
        if basename:
            self.basename=basename
        else:
            self.basename=os.path.splitext(self.default_script)[0]
        self.scriptname=f'{self.basename}{self.default_ext}'
        self.newfile(self.scriptname)
        self.banner(f'{__package__}: {self.basename}{self.default_ext}')
        self.banner(f'Created {timestampstr}')

    def writescript(self):
        """
        Write the Packmol input script to the file specified in the ``scriptname`` attribute.
        """
        self.writefile()

    def runscript(self,*args,**options):
        """
        Run the Packmol script using the specified shell command.
        This method constructs a command to execute Packmol with the specified script and options.
        """
        assert hasattr(self,'scriptname'),f'No scriptname set.'
        self.logname=f'{self.basename}.log'
        self.logparser=PackmolLog(basename=self.basename)
        logger.debug(f'Log file: {self.logname}')
        cmd=Command(f'{self.packmol} < {self.scriptname}')
        progress_struct=None
        if self.progress:
            progress_struct=PackmolProgress()
            self.logparser.enable_progress_bar(progress_struct)
        return cmd.run(ignore_codes=[173],logfile=self.logname,logparser=self.logparser)
    