# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import datetime
import logging
import os
import re
import shutil
import subprocess

from .generic import GenericScripter
from ..core.command import Command
from ..util.stringthings import FileCollector
from ..logparsers.packmollogparser import PackmolLogParser
from ..util.progress import PackmolProgress

logger = logging.getLogger(__name__)

MIN_PACKMOL_VERSION = (20, 15, 1)

_PACKMOL_VERSION_RE = re.compile(r'Version\s+(\d+)\.(\d+)\.(\d+)')


def get_packmol_version(executable):
    """
    Probe ``executable`` by sending it an empty stdin and parse the version banner.
    Returns a ``(major, minor, patch)`` tuple, or ``None`` if the banner cannot be parsed.
    """
    try:
        proc = subprocess.run(
            [executable], input='', capture_output=True, text=True, timeout=10,
        )
    except (OSError, subprocess.TimeoutExpired) as e:
        logger.debug(f'Could not probe {executable!r} for version: {e}')
        return None
    m = _PACKMOL_VERSION_RE.search(proc.stdout)
    if not m:
        return None
    return tuple(int(x) for x in m.groups())


def _is_conda_install(executable_path):
    """
    Return True if ``executable_path`` sits inside a conda environment
    (i.e. its env root contains a ``conda-meta`` directory).
    """
    if not executable_path:
        return False
    env_root = os.path.dirname(os.path.dirname(executable_path))
    return os.path.isdir(os.path.join(env_root, 'conda-meta'))


_AMBERTOOLS_HINT = (
    ' The packmol bundled inside "packmol-memgen" (part of AmberTools, typically '
    'installed via conda) is too old for pestifer. Download and compile packmol '
    '>= 20.15.1 from https://github.com/m3g/packmol and put it ahead of the conda '
    'version on your PATH.'
)


def check_packmol_version(executable, min_version=MIN_PACKMOL_VERSION):
    """
    Raise ``RuntimeError`` if ``executable``'s reported version is below ``min_version``.
    Returns the detected version tuple on success.
    """
    found = get_packmol_version(executable)
    min_str = '.'.join(str(x) for x in min_version)
    resolved = shutil.which(executable) or executable
    hint = _AMBERTOOLS_HINT if _is_conda_install(resolved) else ''
    if found is None:
        raise RuntimeError(
            f'Could not determine version of packmol at {resolved!r}; '
            f'pestifer requires packmol >= {min_str}.' + hint
        )
    if found < min_version:
        found_str = '.'.join(str(x) for x in found)
        raise RuntimeError(
            f'packmol at {resolved!r} is version {found_str}; '
            f'pestifer requires >= {min_str}.' + hint
        )
    return found

class PackmolScripter(GenericScripter):
    """
    This class extends the GenericScripter class to provide functionality for creating and managing Packmol scripts.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(comment_char=kwargs.get('comment_char', '#'))
        self.indent = 4 * ' '
        self.progress = kwargs.get('progress', True)
        self.packmol = kwargs.get('packmol')
        self.F = FileCollector()
        self.default_ext = '.inp'
        self.default_script = f'packmol{self.default_ext}'
        self.scriptname = self.default_script

    def newscript(self, basename=None):
        """
        Create a new Packmol input script.

        Parameters
        ----------
        basename : str, optional
            The base name for the script file (without extension). If not provided, a default name will be used.
        """
        timestampstr = datetime.datetime.today().ctime()
        if basename:
            self.basename = basename
        else:
            self.basename = os.path.splitext(self.default_script)[0]
        self.scriptname = f'{self.basename}{self.default_ext}'
        self.newfile(self.scriptname)
        self.banner(f'{__package__}: {self.basename}{self.default_ext}')
        self.banner(f'Created {timestampstr}')

    def writescript(self):
        """
        Write the Packmol input script to the file specified in the ``scriptname`` attribute.
        """
        self.writefile()

    log_error_patterns = (
        re.compile(r'\bERROR\b'),
        re.compile(r'Could not open file'),
        re.compile(r'Could not find'),
    )

    def runscript(self, *args, **options):
        """
        Run the Packmol script using the specified shell command.
        This method constructs a command to execute Packmol with the specified script and options.
        After packmol exits, the log is scanned for error markers; if any are found a non-zero
        code is returned so callers raise.
        """
        assert hasattr(self, 'scriptname'), f'No scriptname set.'
        self.logname = f'{self.basename}.log'
        self.logparser = PackmolLogParser(basename=self.basename)
        logger.debug(f'Log file: {self.logname}')
        cmd = Command(f'{self.packmol} < {self.scriptname}')
        progress_struct = None
        if self.progress:
            progress_struct = PackmolProgress()
            self.logparser.enable_progress_bar(progress_struct)
        rc = cmd.run(ignore_codes=[173], logfile=self.logname, logparser=self.logparser)
        log_errors = self._scan_log_for_errors()
        if log_errors:
            logger.error(f'Packmol log {self.logname} contains {len(log_errors)} error line(s):')
            for line in log_errors:
                logger.error(f'  {line}')
            if rc == 0:
                rc = 1
        return rc

    def _scan_log_for_errors(self):
        """
        Scan ``self.logname`` for packmol error markers. Returns a list of offending lines
        (empty if the log is clean or absent).
        """
        if not os.path.exists(self.logname):
            return []
        hits = []
        with open(self.logname) as f:
            for line in f:
                if any(p.search(line) for p in self.log_error_patterns):
                    hits.append(line.rstrip())
        return hits
