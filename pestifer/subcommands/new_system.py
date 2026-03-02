# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
The new-system subcommand.  Generates a minimal or full system configuration script for a given
PDB or UniProt ID.
"""
import os
import argparse as ap

from dataclasses import dataclass

from . import Subcommand

from ..core.resourcemanager import ResourceManager

@dataclass
class NewSystemSubcommand(Subcommand):
    name: str = 'new-system'
    short_help: str = 'create a new system script from scratch'
    long_help: str = 'Generate a new system script with a basic template.'

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        r = ResourceManager()
        build_type = 'full' if args.full else 'minimal'
        outputfilename = args.output if args.output else f'{args.id}.yaml'
        title = args.title if args.title else ''
        if os.path.exists(outputfilename):
            raise FileExistsError(f'Output file {outputfilename} already exists; please remove it or choose a different name')
        r.example_manager.new_example_yaml(db_id=args.id, build_type=build_type, outputfilename=outputfilename, title=title)
        return True

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('id', type=str, default=None, help='PDB/UniProt ID  of the new system to create')
        self.parser.add_argument('--full', default=False, action=ap.BooleanOptionalAction, help='enable/disable full set of tasks in the new system configuration (default: %(default)s)')
        self.parser.add_argument('--output', type=str, default=None, help='output filename for the new system configuration (default: <id>.yaml)')
        self.parser.add_argument('--title', type=str, default=None, help='title for the new system configuration (default: none)')
        return self.parser
