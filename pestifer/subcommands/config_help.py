# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
The config-help subcommand.  This subcommand exposes Ycleptic's built-in interactive help system for
exploring Pestifer's base configuration.  The user can specify the entry point using the
``directives`` argument, and can disable interactivity using the ``--no-interactive`` flag.
"""
import argparse as ap

from dataclasses import dataclass

from . import Subcommand

from ..core.config import Config

@dataclass
class ConfigHelpSubcommand(Subcommand):
    name: str = 'config-help'
    short_help: str = "show help for configuration options"
    long_help: str = "Display help information for all available configuration options."

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        config = Config().configure_new()
        directives = args.directives
        iprompt='pestifer-help: ' if args.interactive else ''
        config.console_help(directives, interactive_prompt=iprompt, exit=True)
        return True

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('directives', type=str, nargs='*', default=[], help='specific directives at which the help system begins (default is at the root)')
        self.parser.add_argument('--interactive', default=True,
            action=ap.BooleanOptionalAction, help='enable/disable interactive mode to query config options (default: %(default)s)')
