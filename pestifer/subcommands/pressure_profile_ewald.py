# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
The pressure-profile-ewald subcommand.  Reconstructs a *complete* NAMD pressure profile --
reciprocal-space term included -- by replaying a saved trajectory twice and summing the halves.

A run with ``pressureProfile on`` omits the PME reciprocal-space contribution, and NAMD makes
``pressureProfileEwald`` mutually exclusive with the ordinary profile, so no single run reports the
whole thing.  Computing the Ewald half inline is impractical (it evaluates every step, not every
sampled frame); replaying the trajectory costs the inline penalty divided by the sampling stride.
See :mod:`pestifer.util.ppewald` for the mechanism and its validation.

Example:
--------

.. code-block:: bash

   $ pestifer pressure-profile-ewald --namd-config prod.namd --dcd prod.dcd --slabs 20

This reads the production config for its force field, cutoffs and PME settings, emits two replay
configs, runs both over ``prod.dcd``, and writes ``ppreplay-pressure-profile.png`` and
``ppreplay-pressure-profile.csv``.
"""
import logging
import os

import argparse as ap

from dataclasses import dataclass

from . import Subcommand
from ..util.ppewald import replay, write_csv, plot, validate_against_total_pressure
from ..util.provenance import stamp as provenance_stamp

logger = logging.getLogger(__name__)


@dataclass
class PressureProfileEwaldSubcommand(Subcommand):
    name: str = 'pressure-profile-ewald'
    group: str = 'After the run'
    log_file: str = 'pressure-profile-ewald.log'
    short_help: str = "reconstruct a complete pressure profile, reciprocal term included"
    long_help: str = ("Replay a trajectory twice -- once for the real-space profile, once for the "
                      "PME reciprocal-space profile, which NAMD will not report together -- and "
                      "sum them frame by frame into the complete pressure profile.")
    func_returns_type: type = bool

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        for name, path in (('--namd-config', args.namd_config), ('--dcd', args.dcd)):
            if not os.path.exists(path):
                raise FileNotFoundError(f'{name} {path}: no such file')
        grid = tuple(args.ewald_grid)
        result = replay(args.namd_config, args.dcd,
                        workdir=args.workdir, slabs=args.slabs, ewald_grid=grid,
                        stride=args.stride, temperature=args.temperature,
                        namd=args.namd, nprocs=args.np, prefix=args.prefix)
        logger.info(f'reconstructed {result.nframes} frames x {result.nslabs} slabs')

        check = validate_against_total_pressure(result)
        if check is not None:
            # The real-space pass reports the total pressure of a complete force evaluation in its
            # ENERGY lines, so this is an independent check that summing the halves was right.
            logger.info('agreement with NAMD\'s own reported pressure '
                        f'(mean |slab-average - PRESSURE|, over {result.nframes} frames):')
            logger.info(f'    real space only : {check["real_only_deviation"]:10.2f} bar')
            logger.info(f'    reconstructed   : {check["reconstructed_deviation"]:10.2f} bar')
            if check['reconstructed_deviation'] > check['real_only_deviation']:
                logger.warning('the reconstruction agrees WORSE than the real-space half alone; '
                               'treat this profile as suspect')

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
                                      'reused verbatim so the profile matches the run')
        self.parser.add_argument('--dcd', type=str, required=True,
                                 help='trajectory to replay (must carry unit cell data, i.e. it '
                                      'was written with DCDUnitCell yes)')
        self.parser.add_argument('--slabs', type=int, default=20,
                                 help='pressureProfileSlabs for both passes (default: %(default)s)')
        self.parser.add_argument('--ewald-grid', type=int, nargs=3, default=(10, 10, 10),
                                 metavar=('X', 'Y', 'Z'),
                                 help='pressureProfileEwald grid (default: %(default)s, NAMD\'s '
                                      'own default; raise it and re-run to check convergence)')
        self.parser.add_argument('--stride', type=int, default=1,
                                 help='replay every Nth frame (default: %(default)s)')
        self.parser.add_argument('--temperature', type=float, default=300.0,
                                 help='temperature for velocity initialization; affects only the '
                                      'kinetic term (default: %(default)s)')
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
