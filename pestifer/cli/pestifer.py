# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
Defines the command-line interface for pestifer
"""
import argparse as ap
import importlib.metadata
import logging
import os
import shutil
import sys

__pestifer_version__ = importlib.metadata.version("pestifer")
from ..util.stringthings import banner
from ..subcommands import _subcommands
from ..core.errors import PestiferError

logger = logging.getLogger(__name__)

logging.getLogger("pidibble").setLevel(logging.WARNING)
logging.getLogger("ycleptic").setLevel(logging.WARNING)
logging.getLogger("matplotlib").setLevel(logging.WARNING)
logging.getLogger("filelock").setLevel(logging.WARNING)

class NiceHelpFormatter(ap.HelpFormatter):
    """
    A custom help formatter that provides a more readable help output.
    It formats the help text with a maximum help position and width.
    """
    def __init__(self, prog):
        super().__init__(prog, max_help_position=40, width=100)

def subcommand_exit_code(result) -> int:
    """Exit status implied by whatever a subcommand returned.

    A subcommand that ran tasks reports failure either as an int or via an ``exit_code``
    attribute on the object it returns (the Controller); without honoring that, the process
    exited 0 for an aborted build and any caller keying off the exit status recorded it as a
    success.

    Most subcommands, though, return ``True`` to mean "done" -- and ``bool`` is a subclass of
    ``int``, so treating a bare int as a status turned every one of those into exit 1.  Booleans
    are therefore success, and only a genuine int (or a genuine int ``exit_code``) is a status.
    """
    if result is None or isinstance(result, bool):
        return 0
    if isinstance(result, int):
        return result
    code = getattr(result, 'exit_code', 0)
    if isinstance(code, bool) or not isinstance(code, int):
        return 0
    return code


def cli():
    """
    Command-line interface for pestifer.
    """

    parser = ap.ArgumentParser(formatter_class=NiceHelpFormatter)
    parser.add_argument(
        "--banner",
        default=True,
        action=ap.BooleanOptionalAction,
        help="enable or disable the banner"
    )
    parser.add_argument(
        '--kick-ass',
        default=False,
        action=ap.BooleanOptionalAction,
        help=ap.SUPPRESS)
    parser.add_argument(
        '--log-level',
        type=str,
        default='debug',
        choices=['info', 'debug', 'warning'],
        help='logging level (default: %(default)s)')
    parser.add_argument(
        '--log-file',
        type=str,
        default=None,
        help="log file (default: 'pestifer_diagnostics.log'; build/run derive it from the config name)")
    parser.add_argument(
        '--version',
        action='version',
        version=f'%(prog)s {__pestifer_version__}'
    )
    subparsers = parser.add_subparsers(
        title="Available commands (use \"pestifer <command> --help\" for help with any command)",
        dest="command",
        metavar="",
        required=False
    )

    for subcommand in _subcommands:
        subcommand.add_subparser(subparsers)

    args = parser.parse_args()
    loglevel_numeric = getattr(logging, args.log_level.upper())
    default_log_file_func = getattr(args, 'default_log_file_func', None)
    log_file = args.log_file or (default_log_file_func(args) if default_log_file_func else None)
    if log_file:
        if os.path.exists(log_file):
            shutil.copyfile(log_file, log_file + '.bak')
        logging.basicConfig(filename=log_file, filemode='w', format='%(asctime)s %(name)s %(message)s', level=loglevel_numeric)
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(logging.Formatter('%(levelname)s> %(message)s'))
    logging.getLogger('').addHandler(console)
    # Emit the banner on the same stream the console log handler uses (stderr), not stdout.
    # Under `pestifer ... > run.log 2>&1` the two streams share one file, but a redirected
    # stdout is block-buffered and only flushes at exit, so a banner printed to stdout lands
    # at the *end* of the log -- or mid-file once its buffer fills -- even though it is written
    # first.  stderr stays line-buffered when redirected, so writing there keeps the banner at
    # the top where it belongs.
    banner(lambda line: print(line, file=console.stream, flush=True), args)
    try:
        result = args.func(args)
    except PestiferError as e:
        logger.error(str(e))
        sys.exit(1)
    exit_code = subcommand_exit_code(result)
    if exit_code:
        logger.error(f'pestifer exiting with code {exit_code}')
        sys.exit(exit_code)
