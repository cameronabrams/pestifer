# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
The show-resources subcommand.  Shows available resources for the current configuration,
including TcL scripts, examples, and CHARMM force field data.
"""
import argparse as ap

from dataclasses import dataclass

from . import Subcommand

from ..core.resourcemanager import ResourceManager

@dataclass
class ShowResourcesSubcommand(Subcommand):
    name: str = 'show-resources'
    short_help: str = 'show available resources for the current configuration'
    long_help: str = 'Display a list of all available resources for the current configuration.'

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        charmmff_config = {}
        if args.user_pdbcollection:
            charmmff_config['pdbrepository'] = args.user_pdbcollection
        if args.charmmff_release:
            charmmff_config['release'] = args.charmmff_release
        r = ResourceManager(charmmff_config=charmmff_config)
        resource_type = args.resource_type
        specs = {resource_type: []}
        if hasattr(args, resource_type):
            for subresource in getattr(args, resource_type):
                specs[resource_type].append(subresource)
        r.show(out_stream=print, components=specs, fullnames=args.fullnames, missing_fullnames=r.labels.residue_fullnames)
        return True

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('resource_type', type=str, default='examples', help='type of resource to show; [tcl|examples|charmmff]')
        self.parser.add_argument('--charmmff', type=str, nargs='+', default=[], help='show sub-resources of charmmff resources (\'toppar\', \'custom\', \'pdb\')')
        self.parser.add_argument('--fullnames', default=False, action='store_true', help='show full names of any residues shown with --charmmff pdb')
        self.parser.add_argument('--user-pdbcollection', type=str, nargs='+', default=[], help='additional collections of PDB files outside pestifer installation')
        self.parser.add_argument('--charmmff-release', type=str, default='', help='CHARMMFF release to use (e.g. "February2026"); defaults to the newest available')
        return self.parser
