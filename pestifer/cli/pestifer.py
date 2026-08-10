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
        args.func(args)
    except PestiferError as e:
        logger.error(str(e))
        sys.exit(1)
