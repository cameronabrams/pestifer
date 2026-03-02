# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
The config-default subcommand.  Echoes all or part of a configuration script with default
values, depending on the ``directives`` entry point.
"""
import yaml
import argparse as ap

from dataclasses import dataclass

from . import Subcommand

from ..core.config import Config

@dataclass
class ConfigDefaultSubcommand(Subcommand):
    name: str = 'config-default'
    short_help: str = "show default configuration options"
    long_help: str = "Display default values for all available configuration options."

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        config = Config().configure_new()
        directives = args.directives
        specs=config.make_default_specs(*directives)
        print(yaml.dump(specs))
        return True

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('directives', type=str, nargs='*', default=[], help='specific directives to get default values for (default: all)')
