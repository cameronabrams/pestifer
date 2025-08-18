# Author: Cameron F. Abrams <cfa22@drexel.edu>

from dataclasses import dataclass

from . import Subcommand
from ..core.controller import Controller
from ..core.config import Config

@dataclass
class DesolvateSubcommand(Subcommand):
    name: str = 'desolvate'
    short_help: str = "Desolvate a system"
    long_help: str = "Remove water molecules from a solvated system."
    return_type: type = dict

    @staticmethod
    def func(args, **kwargs):
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
        report = C.do_tasks()
        return report

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