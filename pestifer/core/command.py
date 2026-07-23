# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#
"""
Class for handling and running external commands and managing their output.
"""

import atexit
import logging
import os
import signal
import shutil
import subprocess
import sys
import threading
import time

from glob import glob

from ..logparsers import LogParser
from ..util.util import running_under_pytest

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Clean shutdown of external child processes (NAMD/VMD).
#
# Each external command is launched in its own session (``start_new_session``),
# so the child is a process-group leader and its whole tree -- e.g.
# ``charmrun`` -> ``namd3`` -> N PE processes -- can be taken down with a single
# ``killpg``.  Without this, interrupting pestifer (Ctrl-C, or ``kill`` from a
# scheduler) left orphaned NAMD jobs pinning every core, because the tracked pid
# was the shell/charmrun wrapper and SIGTERM to it never reached ``namd3``.
# ---------------------------------------------------------------------------

# pid -> was it launched in its own session?  Session-detached children (NAMD via
# charmrun) are torn down with ``killpg`` on their whole tree; children left in
# pestifer's own session (VMD -- see below) must be signaled by pid only, never by
# process group, or we would kill pestifer itself.
_active_children: dict[int, bool] = {}
_active_child_lock = threading.Lock()
_signal_handlers_installed = False


def _register_child(pid: int, new_session: bool = True) -> None:
    with _active_child_lock:
        _active_children[pid] = new_session


def _unregister_child(pid: int) -> None:
    with _active_child_lock:
        _active_children.pop(pid, None)


def _terminate_children(sig: int) -> int:
    """Send *sig* to every live child.  Session-detached children are signaled by
    process group (``killpg``); others by pid.  Returns the count signaled."""
    with _active_child_lock:
        children = list(_active_children.items())
    n = 0
    for pid, new_session in children:
        try:
            if new_session:
                os.killpg(os.getpgid(pid), sig)
            else:
                os.kill(pid, sig)
            n += 1
        except (ProcessLookupError, OSError):
            pass
    return n


def _cleanup_children_atexit() -> None:
    """Backstop for normal interpreter exit: make sure nothing was left running."""
    if _terminate_children(signal.SIGTERM):
        time.sleep(0.3)
        _terminate_children(signal.SIGKILL)


atexit.register(_cleanup_children_atexit)


def _signal_handler(signum, frame):
    name = signal.Signals(signum).name
    with _active_child_lock:
        n = len(_active_children)
    logger.warning(f'Received {name}: shutting down and terminating '
                   f'{n} running child process(es)...')
    _terminate_children(signal.SIGTERM)
    time.sleep(0.5)
    _terminate_children(signal.SIGKILL)
    logger.warning(f'pestifer run INTERRUPTED by {name} before completion. '
                   f'Any partially built files remain in {os.getcwd()}; the build is '
                   f'incomplete and should be restarted.')
    # exit through the normal Python path so logging/atexit flush cleanly, with an
    # exit status that reflects the signal (128 + signum, the shell convention)
    sys.exit(128 + signum)


def install_signal_handlers() -> None:
    """Install SIGINT/SIGTERM handlers that tear down child processes cleanly.

    Idempotent, and a no-op off the main thread (where ``signal.signal`` is
    illegal) or under pytest (so the test runner keeps its own signal handling)."""
    global _signal_handlers_installed
    if _signal_handlers_installed:
        return
    if threading.current_thread() is not threading.main_thread() or running_under_pytest():
        return
    for s in (signal.SIGINT, signal.SIGTERM):
        try:
            signal.signal(s, _signal_handler)
        except (ValueError, OSError):
            pass
    _signal_handlers_installed = True

class Command:
    """ 
    Class for running external commands in a subprocess.
    This class allows you to create a command with its arguments and options, and then run it while capturing its output.
    The command is run in a shell, and you can specify a logfile to write the output to.
    If the command returns a non-zero exit code, an error is logged, and the stdout and stderr buffers are printed.
    You can also specify a tuple of (needle, message) to override the default behavior and log a custom message if the needle is found in the stdout or stderr.
    The command can be run with a progress bar and elapsed time display using a LogParser instance.
    The command can also be run quietly, suppressing the output to the console.
    """

    divider_line_length = 55
    """
    The length of the divider line used in logging output to separate sections of the log.
    """
 
    def __init__(self, command: str, *args, **options):
        """ 
        Initializes a Command instance with a command, its arguments, and options.

        Parameters
        ----------
        command : str
            The command to be executed.
        args : tuple
            A tuple of arguments to be passed to the command.
        options : dict
            A dictionary of options to be passed to the command, where keys are option names and values are option values.
        """
        self.command = command
        self.args = args
        self.options = options
        self.c = f'{self.command} ' + ' '.join(args) + ' '.join([f'-{k} {v}' for k, v in self.options.items()])
        self.stdout = ''
        self.stderr = ''

    def run(self, logfile=None, override=(), ignore_codes=[], quiet=True, logparser: LogParser = None, log_stderr=False, new_session=True, stdin=None, **kwargs):
        """
        Runs this Command instance
        
        Parameters
        ----------
        logfile : str
            name of log file to write process' stdout and stderr. If None (default), stdout and stderr are retained
            only in this instance's stdout and stderr attributes.
        override : tuple
            tuple composed of a "needle" and a "message".  If the needle is found in the stdout or stderr of
            the process, the message is displayed and an error is thrown, halting the program.
        ignore_codes : list
            a list of integer exit codes that are ignored in addition to 0
        quiet : bool
            if True, suppresses all output from the command
        logparser: LogParser
            used for progress bar/elapsed time displays and log parsing
        """
        _pytest = running_under_pytest()
        if not quiet:
            logger.debug(f'{self.c}')
        log = None
        if not logfile:
            if not quiet:
                logger.debug(f'No logfile specified for {self.c}')
        else:
            if os.path.exists(logfile) and not kwargs.get('overwrite_logs', False):
                nlogs = len(glob(f'%{logfile}'))
                shutil.move(logfile, f'%{logfile}-{nlogs+1}%')
                logger.debug(f'Rotating {logfile} to %{logfile}-{nlogs+1}%')
            log = open(logfile, 'w')
            logger.debug(f'Opened {logfile} for writing')

        if log_stderr:
            stderr_redirect = subprocess.STDOUT
        else:
            stderr_redirect = subprocess.PIPE
        # By default each external command runs in its own session so its whole process
        # tree (e.g. charmrun -> namd3 -> PEs) can be torn down with one killpg on
        # interrupt.  VMD passes new_session=False: a new session strips the controlling
        # terminal, and a VMD launcher that wraps the binary in ``rlwrap`` (common in
        # site VMD 2.x builds) then exits immediately without running the script -- so
        # VMD stays in pestifer's session and is signaled by pid instead.
        install_signal_handlers()
        process = subprocess.Popen(self.c, shell=True, stdout=subprocess.PIPE, stderr=stderr_redirect,
                                   text=True, start_new_session=new_session, stdin=stdin)
        _register_child(process.pid, new_session=new_session)
        self.stdout = ''
        self.stderr = ''
        output = ''

        try:
            while True:
                output = process.stdout.readline()
                self.stdout += output
                if logparser:
                    logparser.update(output)
                    if not _pytest:
                        logparser.update_progress_bar()
                if log:
                    log.write(output)
                    log.flush()
                if output == '' and process.poll() is not None:
                    break
            if hasattr(logparser, 'progress_bar') and logparser.progress_bar is not None:
                if not _pytest:
                    logparser.progress_bar.finish()
                else:
                    print()
            if logfile:
                logger.debug(f'Log written to {logfile}')
                log.close()
            remaining_stdout, self.stderr = process.communicate()
            self.stdout += remaining_stdout
        finally:
            _unregister_child(process.pid)
        if logparser:
            logparser.update(remaining_stdout)
            if not _pytest and not (hasattr(logparser, 'progress_bar') and logparser.progress_bar is not None):
                logparser.update_progress_bar()
            if hasattr(logparser, 'finalize'):
                logparser.finalize()
            if hasattr(logparser, 'write_csv'):
                logparser.write_csv()
        if process.returncode != 0 and not process.returncode in ignore_codes:
            logger.error(f'Returncode: {process.returncode}')
            if len(self.stdout) > 0:
                logger.error('stdout buffer follows\n' + '*' * self.divider_line_length + '\n' + self.stdout + '\n' + '*' * self.divider_line_length)
            if len(self.stderr) > 0:
                logger.error('stderr buffer follows\n' + '*' * self.divider_line_length + '\n' + self.stderr + '\n' + '*' * self.divider_line_length)
            return process.returncode
        if len(override) == 2:
            needle, msg = override
            if needle in self.stdout or needle in self.stderr:
                logger.info(f'Returncode: {process.returncode}, but another error was detected:')
                logger.error(msg)
                if len(self.stdout) > 0 and needle in self.stdout:
                    logger.error('stdout buffer follows\n' + '*' * self.divider_line_length + '\n' + self.stdout + '\n' + '*' * self.divider_line_length)
                if len(self.stderr) > 0 and needle in self.stderr:
                    logger.error('stderr buffer follows\n' + '*' * self.divider_line_length + '\n' + self.stderr + '\n' + '*' * self.divider_line_length)
        return 0