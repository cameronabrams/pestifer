# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
The pressure-profile-ewald subcommand.  Reconstructs a *complete* NAMD pressure profile --
reciprocal-space term included -- from a production run's own profile and its trajectory.

A run with ``pressureProfile on`` omits the PME reciprocal-space contribution, and NAMD makes
``pressureProfileEwald`` mutually exclusive with the ordinary profile, so no single run reports the
whole thing.  This replays the trajectory for the reciprocal term alone and adds it, frame by frame
and matched by MD step, to the profile the production run already wrote.  The production run must
use ``fullElectFrequency 1``.  See :mod:`pestifer.util.ppewald` for the mechanism and the
measurements behind it.

Example:
--------

.. code-block:: bash

   $ pestifer pressure-profile-ewald --namd-config prod.namd --log prod.log --dcd prod.dcd

This writes ``ppreplay-pressure-profile.png`` and ``ppreplay-pressure-profile.csv``.
"""
import logging
import os

import argparse as ap

from dataclasses import dataclass

from . import Subcommand
from ..util.ppewald import (replay, write_csv, plot, check_against_tensor,
                             closure_warnings)
from ..util.provenance import stamp as provenance_stamp

logger = logging.getLogger(__name__)


@dataclass
class PressureProfileEwaldSubcommand(Subcommand):
    name: str = 'pressure-profile-ewald'
    group: str = 'After the run'
    log_file: str = 'pressure-profile-ewald.log'
    short_help: str = "reconstruct a complete pressure profile, reciprocal term included"
    long_help: str = ("Replay a production trajectory for the PME reciprocal-space pressure "
                      "profile, which NAMD will not report alongside the ordinary one, and add it "
                      "frame by frame to the production run's own profile.")
    func_returns_type: type = bool

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        for name, path in (('--namd-config', args.namd_config), ('--log', args.log),
                           ('--dcd', args.dcd)):
            if not os.path.exists(path):
                raise FileNotFoundError(f'{name} {path}: no such file')
        grid = tuple(args.ewald_grid) if args.ewald_grid else None
        result = replay(args.namd_config, args.dcd, args.log,
                        workdir=args.workdir, ewald_grid=grid, stride=args.stride,
                        temperature=args.temperature, namd=args.namd, nprocs=args.np,
                        prefix=args.prefix, accept_mts=args.accept_mts)
        if args.slabs is not None and args.slabs != result.nslabs:
            raise RuntimeError(f'--slabs {args.slabs}, but the production profile has '
                               f'{result.nslabs} slabs; the replay always uses the production '
                               'slab count')
        logger.info(f'reconstructed {result.nframes} frames x {result.nslabs} slabs')

        check = check_against_tensor(result)
        if check is None:
            logger.warning(
                'the production log carries no pressure tensor, so the reconstruction was NOT '
                'checked.  Set outputPressure equal to pressureProfileFreq in the production run '
                'to get one.')
        else:
            dev, real = check['deviation'], check['real_only_deviation']
            fmt = lambda v: ' / '.join(f'{x:+8.2f}' for x in v)
            logger.info('agreement with the production run\'s own pressure tensor (frame mean of '
                        f'slab average - PRESSAVG, xx / yy / zz, over {check["nframes"]} frames):')
            logger.info(f'    real space only : {fmt(real)} bar')
            logger.info(f'    complete        : {fmt(dev)} bar')
            for reason in closure_warnings(check):
                logger.warning(reason)

        base = os.path.join(args.workdir, f'{args.prefix}-pressure-profile')
        csv_path, png_path = f'{base}.csv', f'{base}.png'
        write_csv(result, csv_path)
        plot(result, png_path, title=args.title, figsize=tuple(args.figsize),
             stamp=provenance_stamp())
        logger.info(f'wrote {csv_path} and {png_path}')
        print(csv_path)
        print(png_path)
        return True

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('--namd-config', type=str, required=True,
                                 help='the production NAMD config the trajectory came from; its '
                                      'force field, cutoffs, PME and restraint settings are '
                                      'reused verbatim so the reciprocal term matches the run')
        self.parser.add_argument('--log', type=str, required=True,
                                 help='the production NAMD log: its PRESSUREPROFILE records are '
                                      'the real-space half, and its PRESSAVG tensor (written with '
                                      'outputPressure = pressureProfileFreq) is the self-check')
        self.parser.add_argument('--dcd', type=str, required=True,
                                 help='trajectory of that run (must carry unit cell data, i.e. it '
                                      'was written with DCDUnitCell yes, at a multiple of '
                                      'pressureProfileFreq)')
        self.parser.add_argument('--slabs', type=int, default=None,
                                 help='optional check: the slab count you expect; the replay '
                                      'always uses the production profile\'s own')
        self.parser.add_argument('--ewald-grid', type=int, nargs=3, default=None,
                                 metavar=('X', 'Y', 'Z'),
                                 help='pressureProfileEwald grid (default: sized from the cell at '
                                      'about 3 A; raise it and re-run to check convergence)')
        self.parser.add_argument('--stride', type=int, default=1,
                                 help='replay every Nth frame (default: %(default)s)')
        self.parser.add_argument('--temperature', type=float, default=300.0,
                                 help='temperature NAMD needs to start; it does not enter the '
                                      'reciprocal term (default: %(default)s)')
        self.parser.add_argument('--accept-mts', action='store_true',
                                 help='proceed even though the production run used '
                                      'fullElectFrequency > 1 (the result double-counts part of '
                                      'the reciprocal term)')
        self.parser.add_argument('--workdir', type=str, default='.',
                                 help='directory for the replay configs, logs and output')
        self.parser.add_argument('--prefix', type=str, default='ppreplay',
                                 help='basename for generated files (default: %(default)s)')
        self.parser.add_argument('--namd', type=str, default='namd3',
                                 help='NAMD executable (default: %(default)s)')
        self.parser.add_argument('--np', type=int, default=1,
                                 help='processors for NAMD (+pN) (default: %(default)s)')
        self.parser.add_argument('--title', type=str, default='', help='plot title')
        self.parser.add_argument('--figsize', type=float, nargs=2, default=(15, 5),
                                 help='figure size in inches (default: %(default)s)')
        return self.parser
