# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
The density-profile subcommand.  Computes and plots species-resolved mass-density
profiles (water, lipid, protein, ions) along the bilayer normal z for a membrane
system, from a PSF, a single coordinate frame (PDB or NAMD binary ``.coor``), and
an XSC cell file.

The bilayer midplane is centered at z=0 and bulk solvent is made continuous through
the periodic boundary.  For a multicomponent bilayer, ``--lipid-components`` adds a
per-lipid-species profile in addition to the total lipid profile.

Example:
--------

.. code-block:: bash

   $ pestifer density-profile --basename my_system --lipid-components

This reads ``my_system.psf``, ``my_system.coor`` (or ``my_system.pdb``), and
``my_system.xsc`` and writes ``my_system-density-profile.png``.
"""
import logging
import os

import argparse as ap

from dataclasses import dataclass

from . import Subcommand
from ..util.densityprofile import DensityProfile

logger = logging.getLogger(__name__)


@dataclass
class DensityProfileSubcommand(Subcommand):
    name: str = 'density-profile'
    log_file: str = 'density-profile.log'
    short_help: str = "plot species-resolved density profiles along z"
    long_help: str = ("Compute and plot water/lipid/protein/ion mass-density profiles "
                      "along the bilayer normal from a PSF, a coordinate frame, and an XSC.")
    func_returns_type: type = bool

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        psf, coor, xsc = args.psf, args.coor, args.xsc
        if args.basename:
            psf = psf or f'{args.basename}.psf'
            xsc = xsc or f'{args.basename}.xsc'
            if not coor:
                for ext in ('.coor', '.pdb'):
                    if os.path.exists(f'{args.basename}{ext}'):
                        coor = f'{args.basename}{ext}'
                        break
        missing = [n for n, v in (('psf', psf), ('coor', coor), ('xsc', xsc)) if not v]
        if missing:
            raise ValueError(f'missing input(s): {", ".join(missing)} '
                             '(give --basename or --psf/--coor/--xsc)')
        out = args.out or (f'{args.basename}-density-profile.png' if args.basename
                           else 'density-profile.png')
        dp = DensityProfile(psf, coor, xsc)
        if args.lipid_components and len(dp.lipid_resnames) > 1:
            logger.info(f'lipid components: {", ".join(dp.lipid_resnames)}')
        dp.plot(out, title=args.title, dz=args.dz,
                lipid_components=args.lipid_components, figsize=tuple(args.figsize))
        logger.info(f'wrote {out}')
        print(out)
        return True

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('--basename', type=str, default=None,
                                 help='basename for <basename>.psf/.coor(.pdb)/.xsc and output')
        self.parser.add_argument('--psf', type=str, default=None, help='input PSF file')
        self.parser.add_argument('--coor', type=str, default=None,
                                 help='input coordinate frame (PDB or NAMD binary .coor)')
        self.parser.add_argument('--xsc', type=str, default=None, help='input XSC cell file')
        self.parser.add_argument('--out', type=str, default=None,
                                 help='output PNG (default: <basename>-density-profile.png)')
        self.parser.add_argument('--title', type=str, default='', help='plot title')
        self.parser.add_argument('--dz', type=float, default=1.0,
                                 help='z-slab thickness in Angstroms (default: %(default)s)')
        self.parser.add_argument('--lipid-components', action='store_true',
                                 help='also plot a profile for each lipid species '
                                      '(in addition to the total lipid profile)')
        self.parser.add_argument('--figsize', type=float, nargs=2, default=(6.4, 4.4),
                                 help='figure size in inches (default: %(default)s)')
        return self.parser
