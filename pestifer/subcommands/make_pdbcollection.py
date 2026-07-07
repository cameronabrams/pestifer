# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
The make-pdb-collection subcommand.  This allows a user to generate their own collection of
sample PDB files for residues defined in the CHARMM force field (or in files that have
CHARMM-format RESI blocks). These sample PDB's are used as inputs to the grid packer when
building membranes.
"""
import argparse as ap

from dataclasses import dataclass

from . import Subcommand
from ..charmmff.make_pdb_collection import make_pdb_collection
from ..charmmff.make_solvent_box import build_solvent_entry

@dataclass
class MakePDBCollectionSubcommand(Subcommand):
    name: str = 'make-pdb-collection'
    short_help: str = "create a PDB collection from a set of input files"
    long_help: str = "Generate a PDB collection from a list of PDB files or IDs."

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        """
        Create a PDB collection.

        The default (flag-based) mode samples conformers of RESIs from a CHARMM stream or
        topology file (single-molecule entries for the grid membrane packer).  The
        ``solvent`` mode instead builds a pre-equilibrated periodic *box* of one solvent
        (a ``kind: box`` entry for VMD's ``solvate`` plugin).
        """
        if getattr(args, 'mode', None) == 'solvent':
            build_solvent_entry(args)
        else:
            make_pdb_collection(args)
        return True

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)

        # ---- solvent box-build mode (a curated single solvent, not a stream scan) -------
        # optional verb; when absent, the flag-based conformer/stream pipeline below runs
        modes = self.parser.add_subparsers(dest='mode', metavar='MODE',
                                           help="'solvent' to build a pre-equilibrated solvent box; "
                                                "omit for the default conformer/stream sampler")
        sol = modes.add_parser('solvent', help='build a pre-equilibrated periodic solvent box (kind: box entry)')
        sol.add_argument('--resname', type=str, required=True, help='RESI name of the solvent to build a box for')
        sol.add_argument('--nmol', type=int, default=216, help='number of solvent molecules in the box (default: %(default)s)')
        sol.add_argument('--density', type=float, default=1.0, help='target bulk density in g/cc used to size the initial box (default: %(default)s)')
        sol.add_argument('--temperature', type=float, default=300.0, help='equilibration temperature in K (default: %(default)s)')
        sol.add_argument('--pressure', type=float, default=1.0, help='equilibration pressure in bar (default: %(default)s)')
        sol.add_argument('--minimize-steps', dest='minimize_steps', type=int, default=1000, help='minimization steps before NPT (default: %(default)s)')
        sol.add_argument('--npt-steps', dest='npt_steps', type=int, default=50000, help='NPT equilibration steps (default: %(default)s)')
        sol.add_argument('--seed', type=int, default=None, help='random seed for the packing orientations (default: none)')
        sol.add_argument('--key-atom', dest='key_atom', type=str, default=None, help="key atom name for solvate -ks (default: the residue's first atom)")
        sol.add_argument('--refic-idx', dest='refic_idx', type=int, default=0, help='reference IC index for building the single molecule (default: %(default)s)')
        sol.add_argument('--output-dir', dest='output_dir', type=str, default=None, help="output directory for the entry; the entry lands in <output-dir>/<RESN>/ (default: 'solvent')")
        sol.add_argument('--charmmff-release', dest='charmmff_release', type=str, default='', help='CHARMMFF release to use (default: newest available)')
        sol.add_argument('--no-cleanup', dest='cleanup', action='store_false', help='keep the NAMD scratch working directory')

        self.parser.add_argument('--streamID', type=str, default=None, help='charmmff stream to scan')
        self.parser.add_argument('--substreamID', type=str, default='', help='charmmff substream to scan')
        self.parser.add_argument('--topfile', type=str, default=None, help='charmmff-format topology file to scan')
        self.parser.add_argument('--force', default=False, action='store_true', help='force overwrite of any existing molecules in the database')
        self.parser.add_argument('--cleanup', default=True, action=ap.BooleanOptionalAction, help='clean up all working files (default: %(default)s)')
        self.parser.add_argument('--resname', type=str, default='', help='single resname to generate')
        self.parser.add_argument('--take-ic-from', type=str, default='', help='alternate resname to take ICs from if this resname has bad ICs')
        self.parser.add_argument('--force-constant', type=float, default=1.0, help='harmonic force constant used in non-equilibrium MD to stretch a molecule (default: %(default)s)')
        self.parser.add_argument('--lenfac', type=float, default=1.4, help='this factor times topological distance is the cartesian distance to which you want to stretch a molecule (default: %(default)s)')
        self.parser.add_argument('--minimize-steps', type=int, default=500, help='number of minimization steps immediately after each build (default: %(default)s)')
        self.parser.add_argument('--sample-steps', type=int, default=5000, help='number of sample steps (default: %(default)s)')
        self.parser.add_argument('--nsamples', type=int, default=10, help='number of samples (default: %(default)s)')
        self.parser.add_argument('--sample-temperature', type=float, default=300.0, help='temperature for sampling (default: %(default)s)')
        self.parser.add_argument('--output-dir', type=str, default=None, help='name of output directory relative to CWD; defaults to streamID')
        self.parser.add_argument('--fail-dir', type=str, default='fails', help='name of output directory for failed runs relative to CWD (default: %(default)s)')
        self.parser.add_argument('--refic-idx', type=int, default=0, help='index of reference IC to use to build a single molecule (default: %(default)s)')
        self.parser.add_argument('--charmmff-release', type=str, default='', help='CHARMMFF release to use (e.g. "February2026"); defaults to the newest available')
        return self.parser