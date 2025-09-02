# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
The desolvate subcommand.  This subcommand allows a user to specify a system state 
plus a trajectory which will then be converted to a solvent-free state and trajectory.
The state is represented minimally by a PDB, PSF, and XSC file, and the trajectory can 
be a list of one or more DCD files, ordered chronologically.  The user can specify the 
exact collection atoms to keep using the --keepatsel option; by default, this 
is "protein or glycan or lipid".

Example:
--------

.. code-block:: bash

   $ pestifer desolvate --psf input.psf --pdb input.pdb --dcd-infiles input1.dcd input2.dcd

This will create the new PSF file ``dry.psf`` and a single DCD file concatentating the inputs, 
called ``dry.dcd``, by default.  The output file names can be specified using the --psf-outfile
and --dcd-outfile options.

``desolvate`` invokes both VMD to create an atom index file and generate the stripped psf file,
and ``catdcd`` to process the trajectory files.  The atom index file name will be ``dry.idx``
by default, but can be specified using the --idx-outfile option.

"""

import argparse as ap

from dataclasses import dataclass

from . import Subcommand
from ..core.controller import Controller
from ..core.config import Config

@dataclass
class DesolvateSubcommand(Subcommand):
    name: str = 'desolvate'
    short_help: str = "Desolvate a system"
    long_help: str = "Remove solvent molecules from a solvated system."
    func_returns_type: type = dict

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        config = Config().configure_new()
        C = Controller(
            config, userspecs={
                'tasks': [{
                    'desolvate': {
                        'basename': 'desolvate',
                        'keepatselstr': args.keepatselstr,
                        'psf':args.psf,
                        'pdb':args.pdb,
                        'dcd_infiles':args.dcd_infiles,
                        'dcd_outfile':args.dcd_outfile,
                        'psf_outfile':args.psf_outfile,
                        'idx_outfile':args.idx_outfile,
                        'dcd_stride':args.dcd_stride
                    }
                }]
            }
        )
        return C.do_tasks()

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('--psf', type=str, help='name of input PSF file')
        self.parser.add_argument('--pdb', type=str, help='name of input PDB file (optional)')
        self.parser.add_argument('--dcd-infiles', type=str, nargs='+', default=[], help='list of input dcd files in chronological order')
        self.parser.add_argument('--keepatselstr', type=str, default='protein or glycan or lipid', help='VMD atomsel string for atoms you want to keep (default: "%(default)s")')
        self.parser.add_argument('--psf-outfile', type=str, default='dry.psf', help='name of output PSF file to create (default: %(default)s)')
        self.parser.add_argument('--dcd-outfile', type=str, default='dry.dcd', help='name of DCD output file to create (default: %(default)s)')
        self.parser.add_argument('--idx-outfile', type=str, default='dry.idx', help='name of index file for catdcd to create  (default: %(default)s)')
        self.parser.add_argument('--dcd-stride', type=int, default=1, help='stride in number of frames for catdcd (default: %(default)s)')
        return self.parser