# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
The mdplot subcommand.  This subcommand exposes the ``mdplot`` task standalone interface for generating plots of data extracted from NAMD log and xst files.
"""
import logging

import argparse as ap

from dataclasses import dataclass

from pestifer.core.artifacts import FileArtifactList

from . import Subcommand

from ..core.controller import Controller
from ..core.config import Config

logger = logging.getLogger(__name__)

@dataclass
class MDPlotSubcommand(Subcommand):
    name: str = 'mdplot'
    log_file: str = 'mdplot.log'
    short_help: str = "generate plots from MD simulation data"
    long_help: str = "Generates plots from all combined MD simulations"
    func_returns_type: type = FileArtifactList

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        config = Config().configure_new()
        timeseries = []
        if args.timeseries:
            timeseries += args.timeseries
        if args.timecoseries:
            timeseries += [args.timecoseries]
        C = Controller().configure(
            config, userspecs={
                'tasks': [{
                    'mdplot': {
                        'basename': args.basename,
                        'reprocess-logs': True,
                        'logs': args.logs,
                        'figsize': args.figsize,
                        'timeseries': timeseries,
                        'profiles': args.profiles,
                        'profiles-per-block': args.profiles_per_block,
                        'colormap': args.colormap,
                        'colormap-direction': args.colormap_direction,
                        'legend': True,
                        'grid': True,
                        # getattr: programmatic callers build a Namespace by hand and need not
                        # know about every flag the parser happens to define
                        'stage-markers': not getattr(args, 'no_stage_markers', False),
                        'block-average': getattr(args, 'block_average', -1),
                        'histograms': getattr(args, 'histograms', []) or [],
                        'histogram-tail': getattr(args, 'histogram_tail', 0.5),
                        'lipids-per-leaflet': getattr(args, 'lipids_per_leaflet', 0),
                        'units': {
                            'density': 'g/cc',
                            'a_x': 'Å',
                            'b_y': 'Å',
                            'c_z': 'Å',
                            'pressure': 'bar'
                        }
                    }
                }]
            },
            terminate=False, validate=False  # standalone plot utility, not a build pipeline
        )
        C.tasks[0].taskname = args.basename
        report = C.do_tasks()
        all_fileartifacts = C.pipeline.get_all_file_artifacts()
        logger.debug(f'{len(all_fileartifacts)} file artifacts found:')
        artifact_keys = set([art.key for art in all_fileartifacts])
        chk=0
        for key in artifact_keys:
            logger.debug(f'Artifacts with key {key}:')
            key_artifacts = all_fileartifacts.filter_by_key(key)
            for artifact in key_artifacts:
                chk+=1
                logger.debug(f'  {artifact.name}')
        logger.debug(f'Total of {chk} file artifacts found.')
        return all_fileartifacts

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)
        self.parser.add_argument('--logs', type=str, default=[], nargs='+', help='list of one more NAMD logs in chronological order')
        self.parser.add_argument('--basename', type=str, default='mdplot', help='basename of output files')
        self.parser.add_argument('--figsize', type=int, nargs=2, default=[9, 6],
                                 help='figure size in inches (default: %(default)s)')
        self.parser.add_argument('--timecoseries', type=str, default=[], nargs='+', help='timeseries to plot on same axes')
        self.parser.add_argument('--timeseries', type=str, default=['density'], nargs='+', help='timeseries to plot')
        self.parser.add_argument('--profiles', type=str, default=[], nargs='*', help='profiles (along z) to plot')
        self.parser.add_argument('--profiles-per-block', type=int, default=100, help='number of profiles to plot per block (default: %(default)s)')
        self.parser.add_argument('--colormap', type=str, default='viridis', help='matplotlib colormap for multiple traces on a single plot (default: %(default)s)')
        self.parser.add_argument('--colormap-direction', type=int, choices=[1,-1], default=1, help='direction of colormap (1 or -1) (default: %(default)s)')
        self.parser.add_argument('--block-average', type=int, default=-1, metavar='N',
                                 help='draw a rolling mean over N samples with a +/- 1 sigma band, '
                                      'with the raw trace faint behind it (default: off)')
        self.parser.add_argument('--no-stage-markers', default=False, action='store_true',
                                 help='do not mark the boundaries between simulation stages')
        self.parser.add_argument('--histograms', type=str, default=[], nargs='+', metavar='NAME',
                                 help='quantities whose distribution to plot, one figure each')
        self.parser.add_argument('--histogram-tail', type=float, default=0.5, metavar='FRAC',
                                 help='trailing fraction of each series a histogram uses '
                                      '(default: %(default)s)')
        self.parser.add_argument('--lipids-per-leaflet', type=int, default=0, metavar='N',
                                 help='lipids per leaflet, enabling the derived apl trace')
        return self.parser