# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
The make-pdbcollection subcommand.  This allows a user to generate their own collection of 
sample PDB files for residues defined in the CHARMM force field (or in files that have 
CHARMM-format RESI blocks). These sample PDB's are specifically used as inputs to 
packmol.
"""
from dataclasses import dataclass

from . import Subcommand
from ..charmmff.make_pdb_collection import make_pdb_collection
import argparse as ap

@dataclass
class MakePDBCollectionSubcommand(Subcommand):
    name: str = 'make-pdbcollection'
    short_help: str = "Create a PDB collection from a set of input files"
    long_help: str = "Generate a PDB collection from a list of PDB files or IDs."

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        """
        Create a PDB collection from a set of input files.
        """
        make_pdb_collection(args)
        return True

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('--streamID', type=str, default=None, help='charmmff stream to scan')
        self.parser.add_argument('--substreamID', type=str, default='', help='charmmff substream to scan')
        self.parser.add_argument('--topfile', type=str, default=None, help='charmmff-format topology file to scan')
        self.parser.add_argument('--force', default=False, action='store_true', help='force overwrite of any existing molecules in the database')
        self.parser.add_argument('--cleanup', default=True, action=ap.BooleanOptionalAction, help='clean up all working files (default: %(default)s)')
        self.parser.add_argument('--resname', type=str, default='', help='single resname to generate')
        self.parser.add_argument('--take-ic-from', type=str, default='', help='alternate resname to take ICs from if this resname has bad ICs')
        self.parser.add_argument('--log-level', type=str, default='debug', choices=['info', 'debug', 'warning'], help='Logging level (default: %(default)s)')
        self.parser.add_argument('--force-constant', type=float, default=1.0, help='harmonic force constant used in non-equilibrium MD to stretch a molecule (default: %(default)s)')
        self.parser.add_argument('--lenfac', type=float, default=1.4, help='this factor times topological distance is the cartesian distance to which you want to stretch a molecule (default: %(default)s)')
        self.parser.add_argument('--minimize-steps', type=int, default=500, help='number of minimization steps immediately after each build (default: %(default)s)')
        self.parser.add_argument('--sample-steps', type=int, default=5000, help='number of sample steps (default: %(default)s)')
        self.parser.add_argument('--nsamples', type=int, default=10, help='number of samples (default: %(default)s)')
        self.parser.add_argument('--sample-temperature', type=float, default=300.0, help='temperature for sampling (default: %(default)s)')
        self.parser.add_argument('--output-dir', type=str, default=None, help='name of output directory relative to CWD; defaults to streamID')
        self.parser.add_argument('--fail-dir', type=str, default='fails', help='name of output directory for failed runs relative to CWD (default: %(default)s)')
        self.parser.add_argument('--refic-idx', type=int, default=0, help='index of reference IC to use to build a single molecule (default: %(default)s)')
        return self.parser