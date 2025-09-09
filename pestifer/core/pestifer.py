# Author: Cameron F. Abrams <cfa22@drexel.edu>
""" 
Defines the command-line interface for pestifer
"""
import argparse as ap
import importlib.metadata
import logging

__pestifer_version__ = importlib.metadata.version("pestifer")
from ..util.stringthings import banner
from ..subcommands import _subcommands

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
        help="Enable or disable the banner"
    )
    parser.add_argument(
        '--kick-ass', 
        default=False, 
        action=ap.BooleanOptionalAction, 
        help=ap.SUPPRESS)
    parser.add_argument(
        '--log-level', 
        type=str, 
        default=None, 
        choices=['info', 'debug', 'warning'], 
        help='Logging level (default: %(default)s)')
    parser.add_argument(
        '--log-file', 
        type=str, 
        default='pestifer_diagnostics.log', 
        help='Log file (default: %(default)s)')
    subparsers = parser.add_subparsers(
        title="Available commands (use \"pestifer <command> --help\" for help with any command)",
        dest="command",
        metavar="",
        required=False
    )

    for subcommand in _subcommands:
        subcommand.add_subparser(subparsers)

    args = parser.parse_args()
    banner(print, args)
    args.func(args)
