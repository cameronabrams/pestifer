# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
The new-system subcommand.  Generates a minimal or full system configuration script for a given
PDB or UniProt ID.
"""
import logging
import os
import argparse as ap

from dataclasses import dataclass

from . import Subcommand

from ..core.errors import PestiferError
from ..core.resourcemanager import ResourceManager

logger = logging.getLogger(__name__)


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
            raise PestiferError(f'Output file {outputfilename} already exists; please remove it or choose a different name')
        findings = None
        if getattr(args, 'inspect', False):
            findings = _inspect(args.id)
        r.example_manager.new_example_yaml(db_id=args.id, build_type=build_type, outputfilename=outputfilename, title=title, findings=findings)
        return True

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('id', type=str, default=None, help='PDB/UniProt ID  of the new system to create')
        self.parser.add_argument('--full', default=False, action=ap.BooleanOptionalAction, help='enable/disable full set of tasks in the new system configuration (default: %(default)s)')
        self.parser.add_argument('--inspect', default=False, action='store_true', help='fetch and inspect the structure for missing residues (REMARK 465) and engineered mutations/tags (SEQADV), and annotate the generated config with the corresponding mod stubs')
        self.parser.add_argument('--output', type=str, default=None, help='output filename for the new system configuration (default: <id>.yaml)')
        self.parser.add_argument('--title', type=str, default=None, help='title for the new system configuration (default: none)')
        return self.parser


def _inspect(db_id: str):
    """Fetch+inspect ``db_id`` for the ``--inspect`` annotation, trying PDB then mmCIF. Returns the
    :class:`~pestifer.core.system_inspector.Findings`, or ``None`` if inspection fails entirely
    (the config is still generated, just without annotation)."""
    from ..core.system_inspector import inspect_structure
    source_db = 'alphafold' if (len(db_id) != 4 and db_id[:1] in ('P', 'O')) else 'rcsb'
    last_err = None
    for fmt in ('pdb', 'cif'):
        try:
            return inspect_structure(db_id, source_format=fmt, source_db=source_db)
        except Exception as e:            # network/parse/format problems -- try the other format
            last_err = e
            logger.debug(f'--inspect: {fmt} inspection of {db_id} failed: {e}')
    logger.warning(f'--inspect: could not fetch/inspect {db_id} as PDB or mmCIF ({last_err}); '
                   f'generating the config without inspection annotation.')
    return None
