# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
The show-resources subcommand.  Shows available resources for the current configuration,
including TcL scripts, examples, and CHARMM force field data.
"""
import argparse as ap

from dataclasses import dataclass

from . import Subcommand

from ..core.resourcemanager import ResourceManager


def _report_resname(info: dict, out_stream=print):
    """Format one :meth:`ResourceManager.lookup_resname` result as a short block."""
    out_stream(info['resname'])
    if info['in_topology']:
        line = f"  topology: {info['kind']} defined in {info['topfile']}"
        if info['segtype']:
            line += f" (segtype: {info['segtype']})"
        out_stream(line)
        if info['charmm_alias']:
            out_stream(f"           (looked up under CHARMM resname {info['charmm_alias']})")
    else:
        out_stream('  topology: not found in any CHARMM topology or stream file')
    if info['in_pdbrepository']:
        conf = f"{info['nconformers']} conformer" + ('s' if info['nconformers'] != 1 else '')
        long = f' -- "{info["longname"]}"' if info['longname'] else ''
        out_stream(f"  PDB repo: coordinates available ({conf}){long}")
    else:
        out_stream('  PDB repo: no coordinates in the built-in PDB repository')
    out_stream('')


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
        if resource_type == 'resname':
            names = getattr(args, 'query', []) or []
            if not names:
                print('Usage: pestifer show-resources resname RESNAME [RESNAME ...]')
                return True
            for name in names:
                _report_resname(r.lookup_resname(name), out_stream=print)
            return True
        specs = {resource_type: []}
        if hasattr(args, resource_type):
            for subresource in getattr(args, resource_type):
                specs[resource_type].append(subresource)
        r.show(out_stream=print, components=specs, fullnames=args.fullnames, missing_fullnames=r.labels.residue_fullnames)
        return True

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('resource_type', type=str, default='examples',
                                 choices=['tcl', 'examples', 'charmmff', 'resname'],
                                 help='type of resource to show')
        self.parser.add_argument('query', type=str, nargs='*', default=[],
                                 help='for resource_type=resname: one or more residue names to look up '
                                      '(case-insensitive)')
        self.parser.add_argument('--charmmff', type=str, nargs='+', default=[], help='show sub-resources of charmmff resources (\'toppar\', \'custom\', \'pdb\')')
        self.parser.add_argument('--fullnames', default=False, action='store_true', help='show full names of any residues shown with --charmmff pdb')
        self.parser.add_argument('--user-pdbcollection', type=str, nargs='+', default=[], help='additional collections of PDB files outside pestifer installation')
        self.parser.add_argument('--charmmff-release', type=str, default='', help='CHARMMFF release to use (e.g. "February2026"); defaults to the newest available')
        return self.parser
