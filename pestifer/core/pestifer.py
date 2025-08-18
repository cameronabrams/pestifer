# Author: Cameron F. Abrams <cfa22@drexel.edu>
""" 
Defines all subcommands of the pestifer command
"""
import argparse as ap
import importlib.metadata
import logging
import os
import textwrap

from pathlib import Path

__pestifer_version__ = importlib.metadata.version("pestifer")
from ..util.stringthings import banner, _banner_message, _enhanced_banner_message, oxford
from ..subcommands import _subcommands

logger=logging.getLogger(__name__)

logging.getLogger("pidibble").setLevel(logging.WARNING)
logging.getLogger("ycleptic").setLevel(logging.WARNING)
logging.getLogger("matplotlib").setLevel(logging.WARNING)


class NiceHelpFormatter(ap.HelpFormatter):
    """
    A custom help formatter that provides a more readable help output.
    It formats the help text with a maximum help position and width.
    """
    def __init__(self, prog):
        super().__init__(prog, max_help_position=40, width=100)

class BanneredHelpFormatter(NiceHelpFormatter):
    """
    A custom help formatter that includes a banner at the top of the help output.
    It inherits from NiceHelpFormatter and overrides the format_help method to include a banner message.
    """
    def format_help(self):
        # Return the banner followed by the regular help
        return textwrap.dedent(_banner_message) + '\n\n' + super().format_help()

def cli():
    """
    Command-line interface for pestifer.
    This function sets up the argument parser, defines commands, and executes the appropriate command based on the user's input.
    It includes commands for running simulations, fetching examples, showing resources, and more.
    """
    package_path = Path(__file__).resolve().parent.parent.parent
    is_source_package_with_git = os.path.isdir(os.path.join(package_path, '.git'))

    parser = ap.ArgumentParser(formatter_class=NiceHelpFormatter)
    parser.add_argument('--no-banner', default=False, action='store_true', help='turn off the banner')
    parser.add_argument('--kick-ass-banner', default=False, action='store_true', help=ap.SUPPRESS)
    parser.add_argument('--log-level', type=str, default=None, choices=['info', 'debug', 'warning'], help='Logging level (default: %(default)s)')
    parser.add_argument('--log-file', type=str, default='pestifer_diagnostics.log', help='Log file (default: %(default)s)')
    subparsers = parser.add_subparsers(
        title="available commands (use \"pestifer <command> --help\" for help with any command)",
        dest="command",
        metavar="",
        required=False
    )

    for subcommand in _subcommands:
        subcommand.add_subparser(subparsers)

    args = parser.parse_args()

    if not args.no_banner:
        banner(print)
    if args.kick_ass_banner:
        print(_enhanced_banner_message)

    if hasattr(args, 'func'):
        args.func(args)
    else:
        my_list = oxford(subcommand.name for subcommand in _subcommands)
        print(f'No subcommand found. Expected one of {my_list}')