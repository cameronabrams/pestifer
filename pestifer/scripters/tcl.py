# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import logging
import os

from .generic import GenericScripter
from ..util.provenance import stamp as provenance_stamp
from ..util.stringthings import FileCollector

logger = logging.getLogger(__name__)

class TcLScripter(GenericScripter):
    """
    This class extends the GenericScripter class to provide functionality for creating and managing Tcl scripts.

    Parameters
    ----------
    config : Config
        The configuration object containing settings for the script.
    """
    def __init__(self, *args, **kwargs):
        # forward the rest: GenericScripter picks the build seed out of namd_config, and
        # dropping kwargs here would silently unstamp every .tcl and .namd script
        super().__init__(*args, **kwargs)
        self.progress = kwargs.get('progress', True)
        self.F = FileCollector()
        self.default_ext = '.tcl'
        self.default_script = f'pestifer-script{self.default_ext}'
        self.scriptname = self.default_script

    def newscript(self, basename=None):
        """
        Initialize a new Tcl script with a specified basename.
        If no basename is provided, a default script name is used.

        Parameters
        ----------
        basename : str, optional
            The base name for the script file. If not provided, a default name is used.
        """
        if basename:
            self.basename = basename
        else:
            self.basename = os.path.splitext(self.default_script)[0]
        self.scriptname = f'{self.basename}{self.default_ext}'
        self.newfile(self.scriptname)
        msg = f'{__package__}: {self.basename}{self.default_ext}'
        self.comment(msg)
        self.banner(provenance_stamp(getattr(self, 'build_seed', None)))

    def writescript(self):
        """
        Writes the Tcl script to the file specified by the ``scriptname`` attribute.
        """
        self.writefile()

    def addfile(self, filename):
        """
        Add a file to the list of files associated with the script; these might be inputs or outputs (mostly outputs).
        This method appends the specified filename to the internal ``FileCollector``.

        Parameters
        ----------
        filename : str
            The name of the file to be added.
        """
        self.F.append(filename)


