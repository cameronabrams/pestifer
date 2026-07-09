# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
The show-resources subcommand.  Shows available resources for the current configuration,
including TcL scripts, examples, and CHARMM force field data.
"""
import argparse as ap

from dataclasses import dataclass

from . import Subcommand

from ..core.resourcemanager import ResourceManager


_SOURCE_LABEL = {
    'standard': 'native to the CHARMM release',
    'custom': 'pestifer built-in custom file',
    'user': 'user-custom file (e.g. ~/.pestifer/toppar or a configured searchpath)',
}


def _report_resname(info: dict, out_stream=print):
    """Format one :meth:`ResourceManager.lookup_resname` result as a short block."""
    out_stream(info['resname'])
    if info['in_topology']:
        line = f"  topology: {info['kind']} defined in {info['topfile']}"
        if info['segtype']:
            line += f" (segtype: {info['segtype']})"
        out_stream(line)
        if info.get('charmm_synonym'):
            out_stream(f'            "{info["charmm_synonym"]}"')
        if info.get('source'):
            out_stream(f"  source:   {_SOURCE_LABEL.get(info['source'], info['source'])}")
        if info['charmm_alias']:
            out_stream(f"           (looked up under CHARMM resname {info['charmm_alias']})")
    else:
        out_stream('  topology: not found in any CHARMM topology or stream file')
    if info.get('is_patch'):
        out_stream('  PDB repo: not applicable (a patch modifies a residue; it has no coordinates of its own)')
    elif info['in_pdbrepository']:
        conf = f"{info['nconformers']} conformer" + ('s' if info['nconformers'] != 1 else '')
        long = f' -- "{info["longname"]}"' if info['longname'] else ''
        out_stream(f"  PDB repo: coordinates available ({conf}){long}")
        detail = []
        if info.get('charge') is not None:
            detail.append(f"charge {info['charge']:+g}")
        if info.get('head_tail_length'):
            detail.append(f"head-to-tail {info['head_tail_length']:.1f} Å")
        if detail:
            out_stream(f"            ({', '.join(detail)})")
    else:
        out_stream('  PDB repo: no coordinates in the built-in PDB repository')
    out_stream('')


_SOURCE_TAG = {'standard': 'std', 'custom': 'custom', 'user': 'user'}


def _report_resname_compact(info: dict, out_stream=print):
    """One aligned line summarizing a resname (used for substring search results)."""
    if info['in_topology']:
        kind = 'RESI' if info['kind'] and 'RESI' in info['kind'] else 'PRES'
        where = f"{kind} {info['segtype'] or ''}".strip()
    else:
        where = '-'
    src = _SOURCE_TAG.get(info.get('source'), '-')
    if info.get('is_patch'):
        pdb = 'n/a'
    elif info['in_pdbrepository']:
        pdb = str(info['nconformers'])
    else:
        pdb = '-'
    syn = f"  {info['charmm_synonym']}" if info.get('charmm_synonym') else ''
    out_stream(f"  {info['resname']:<8s} {where:<14s} {src:<7s} pdb: {pdb:<3s}{syn}")


def _overview(out_stream=print):
    """Print a one-line summary of each resource type (the no-argument default)."""
    out_stream('pestifer packaged resources -- use `show-resources <type>`:')
    out_stream('  examples   bundled example builds     (show-resources examples)')
    out_stream('  pdb-repo   PDB-repository residues     (show-resources pdb-repo [--fullnames])')
    out_stream('  tcl        packaged Tcl scripts        (show-resources tcl)')
    out_stream('  resname    look up / search residues   (show-resources resname RESN ... | resname --contains STR)')


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
        if resource_type is None:
            _overview(out_stream=print)
            return True
        if resource_type == 'resname':
            contains = getattr(args, 'contains', None)
            names = getattr(args, 'query', []) or []
            if contains:
                matches = r.search_resnames(contains)
                print(f'{len(matches)} residue name(s) contain "{contains}":')
                for name in matches:
                    _report_resname_compact(r.lookup_resname(name), out_stream=print)
                return True
            if not names:
                print('Usage: pestifer show-resources resname RESNAME [RESNAME ...]  (or --contains SUBSTRING)')
                return True
            for name in names:
                _report_resname(r.lookup_resname(name), out_stream=print)
            return True
        if resource_type == 'pdb-repo':
            r.charmmff_content.provision_pdbrepository()
            r.charmmff_content.pdbrepository.show(print, fullnames=args.fullnames,
                                                  missing_fullnames=r.labels.residue_fullnames)
            return True
        r.show(out_stream=print, components={resource_type: []}, fullnames=args.fullnames,
               missing_fullnames=r.labels.residue_fullnames)
        return True

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('resource_type', type=str, nargs='?', default=None,
                                 choices=['tcl', 'examples', 'pdb-repo', 'resname'],
                                 help='type of resource to show; omit for an overview of all types')
        self.parser.add_argument('query', type=str, nargs='*', default=[],
                                 help='for resource_type=resname: one or more residue names to look up '
                                      '(case-insensitive)')
        self.parser.add_argument('--contains', type=str, default=None,
                                 help='with resource_type=resname: list all known residue names '
                                      'containing this substring (case-insensitive)')
        self.parser.add_argument('--fullnames', default=False, action='store_true', help='with pdb-repo: show the full descriptive name of each residue')
        self.parser.add_argument('--user-pdbcollection', type=str, nargs='+', default=[], help='additional collections of PDB files outside pestifer installation')
        self.parser.add_argument('--charmmff-release', type=str, default='', help='CHARMMFF release to use (e.g. "February2026"); defaults to the newest available')
        return self.parser
