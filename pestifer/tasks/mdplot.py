# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`MDPlotTask` class for making plots of energy-like quantities from NAMD runs.
This class is a descendant of the :class:`BaseTask <pestifer.tasks.basetask.BaseTask>` class and is used to extract energy-like data from NAMD log files,
pressure profiles, and XST files, and to generate plots based on this data.
It handles the collection of energy data from multiple NAMD runs, creates CSV files for energy and pressure profile data,
and generates plots for specified traces and profiles.
The plots can include energy traces, pressure profiles, and histograms, with options for units, legends, and grid lines.
It also supports the extraction of data from existing NAMD log files and XST files, allowing for flexible data visualization.

Usage as a task in a build workflow is described in the :ref:`config_ref tasks mdplot` documentation.  This module is also used in standalone form by the :ref:`subs_mdplot` command.

"""
import logging
import os
import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from matplotlib import colormaps as cmaps

from .basetask import BaseTask
from .mdtask import MDTask
from ..core.errors import PestiferBuildError
from ..util.density_convergence import total_atoms
from ..util.units import g_per_amu, A3_per_cm3
from ..logparsers import NAMDLogParser
from ..util.stringthings import to_latex_math
from ..core.artifacts import *

logger = logging.getLogger(__name__)

# beyond this many stages the labels overlap into noise; the rules alone still carry the story
_MAX_LABELED_STAGES = 8
# a label needs this fraction of the x range clear of the previous one to stay readable
_MIN_LABEL_GAP_FRAC = 0.05
# Observables whose fluctuation amplitude dwarfs the signal, so a raw trace hides the behavior
# being looked for.  NAMD's PRESSAVG/GPRESSAVG are already running averages and are left alone;
# TEMP is thermostatted and stays legible raw.
_AUTO_SMOOTH_COLUMNS = {'PRESSURE', 'GPRESSURE'}
_AUTO_SMOOTH_DIVISOR = 40      # window ~ 1/40 of the series
_AUTO_SMOOTH_MIN = 5
_AUTO_SMOOTH_MAX = 200


def auto_block_average(column: str, n_samples: int) -> int:
    """Window to smooth *column* with when `block-average` is left on automatic, or 0.

    Scaled to the series so a short warm-up is not flattened and a long production run is not
    left ragged."""
    if column.upper() not in _AUTO_SMOOTH_COLUMNS:
        return 0
    window = max(_AUTO_SMOOTH_MIN, min(n_samples // _AUTO_SMOOTH_DIVISOR, _AUTO_SMOOTH_MAX))
    return window if n_samples > 2 * window else 0
logging.getLogger("matplotlib").setLevel(logging.WARNING)

_NAMD_DEFAULT_UNITS = {
    'BOND': 'kcal/mol', 'ANGLE': 'kcal/mol', 'DIHED': 'kcal/mol',
    'IMPRP': 'kcal/mol', 'ELECT': 'kcal/mol', 'VDW': 'kcal/mol',
    'BOUNDARY': 'kcal/mol', 'MISC': 'kcal/mol', 'KINETIC': 'kcal/mol',
    'TOTAL': 'kcal/mol', 'POTENTIAL': 'kcal/mol', 'TOTAL3': 'kcal/mol',
    'TEMP': 'K', 'TEMPAVG': 'K',
    'PRESSURE': 'bar', 'GPRESSURE': 'bar', 'PRESSAVG': 'bar', 'GPRESSAVG': 'bar',
    'VOLUME': 'Å³',
}



# '<controller>-<task>-<run>_<label>' is the basename every NAMD run writes its files under, so a
# stage -- one task's contiguous series of runs -- is recoverable from the filename alone.  Doing
# it this way rather than from the recorded MD history is deliberate: the standalone re-plot path
# (`pestifer mdplot --logs ...`) has no pipeline and hence no history, and stage boundaries are
# exactly as meaningful there.
_STAGE_RE = re.compile(r'^(?P<controller>\d+)-(?P<task>\d+)-(?P<run>\d+)_(?P<label>.*)$')


def stage_of(csv_basename: str):
    """``(stage_key, label)`` for a run's basename, or ``(None, None)`` if it does not parse."""
    m = _STAGE_RE.match(csv_basename)
    if not m:
        return None, None
    return f"{m.group('controller')}-{m.group('task')}", m.group('label')


def select_commensurable_runs(entries: list, natoms: int | None) -> list:
    """The trailing run of MD history entries describing a system of ``natoms`` atoms.

    Walks backwards from the most recent run and stops at the first one describing a different
    system, which is what separates "dynamics on the system I am holding" from dynamics on an
    earlier, different one.  A `solvate` between two MD stages changes the atom count and so ends
    the span; a `validate` or a `ring_check` does not and so does not.

    Entries with no recorded count are kept rather than treated as a boundary: an unknown count is
    not evidence of a different system, and treating it as one would silently truncate the plot --
    the exact failure this replaces.  With ``natoms`` unknown, every entry is kept.
    """
    if natoms is None:
        return list(entries)
    keep = []
    for entry in reversed(entries):
        recorded = entry.get('natoms')
        if recorded is not None and recorded != natoms:
            break
        keep.append(entry)
    keep.reverse()
    return keep


class MDPlotTask(BaseTask):
    """ 
    A class for making plots of energy-like quantities from a series of one or more NAMD 
    runs.  Since NAMD runs are always invoked using a log-parser, a csv file is created that 
    contains all energy-like data from the run. 
    """
    _yaml_header = 'mdplot'
    """
    YAML header for the MDPlotTask, used to identify the task in configuration files as part of a ``tasks`` list.
    """

    @classmethod
    def pipeline_contract(cls, specs):
        from .pipeline_contract import TaskContract, MD_OUTPUT
        # plots molecular-dynamics output; needs an md task to have run
        return TaskContract(requires=(MD_OUTPUT,), provides=())





    def _resolve_units(self, name, df_of_column):
        """``(column_key, dataframe, multiplier, unit_label)`` for a requested quantity.

        Shared by the timeseries and histogram paths: a histogram that ignored the unit spec would
        label a density axis 0.62 (raw amu/A^3) where the timeseries beside it reads 1.03 g/cc, and
        the same number would appear in two different units in one figure set.
        """
        units = 1.0
        unitspec = self.specs.get('units', {}).get(name, '*')
        if unitspec == '*':
            unitspec = _NAMD_DEFAULT_UNITS.get(name.upper(), '*')
            units = 1.0
        elif name == 'density':
            if unitspec in ['g_per_cc', 'g/cc', 'g_per_cm3', 'g/cm3']:
                units = g_per_amu * A3_per_cm3
            else:
                logger.debug(f'Unitspec "{unitspec}" not recognized.')
                units = 1.0
        key = name.upper()
        df = df_of_column.get(key, None)
        if df is None:
            key = name.lower()
            df = df_of_column.get(key, None)
        if df is not None and unitspec == 'kcal/mol' and df[key].abs().max() > 1000:
            units = 1e-3
            unitspec = '1000 kcal/mol'
        return key, df, units, unitspec

    def _add_derived_columns(self):
        """Add cell-geometry columns that the xst series implies but does not carry.

        NAMD logs the three cell vectors; what a membrane build actually wants to look at is the
        lateral area, the volume, and -- given a lipid count, which no log records -- the area per
        lipid.  Computing them here makes them ordinary traces, so they work with every existing
        option (units, block averaging, stage markers) instead of needing bespoke plotting.

        Uses the exact vector forms, |a x b| and |a . (b x c)|, rather than a_x*b_y and its
        product: those agree only for an orthorhombic cell, and a triclinic one would be silently
        wrong.
        """
        df = self.dataframes.get('xst')
        if df is None or df.empty:
            return
        needed = ('a_x', 'a_y', 'a_z', 'b_x', 'b_y', 'b_z', 'c_x', 'c_y', 'c_z')
        if not all(c in df.columns for c in needed):
            logger.debug('xst series lacks full cell vectors; skipping derived cell columns')
            return
        a = df[['a_x', 'a_y', 'a_z']].to_numpy(dtype=float)
        b = df[['b_x', 'b_y', 'b_z']].to_numpy(dtype=float)
        c = df[['c_x', 'c_y', 'c_z']].to_numpy(dtype=float)
        axb = np.cross(a, b)
        df['area'] = np.linalg.norm(axb, axis=1)
        df['volume'] = np.abs(np.einsum('ij,ij->i', axb, c))
        with np.errstate(divide='ignore', invalid='ignore'):
            df['aspect'] = np.divide(np.linalg.norm(a, axis=1), np.linalg.norm(b, axis=1))
        added = ['area', 'volume', 'aspect']
        n_lipids = int(self.specs.get('lipids-per-leaflet', 0) or 0)
        if n_lipids > 0:
            df['apl'] = df['area'] / n_lipids
            added.append('apl')
        self.dataframes['xst'] = df
        logger.debug(f'{self.taskname}: derived cell columns available as traces: {added}')

    def _plot_histograms(self, histograms, df_of_column, output_dir):
        """One distribution figure per named quantity, built from the equilibrated tail.

        A timeseries shows whether a quantity settled; this shows what it settled at and how
        tightly.  The early part of a run is the approach to equilibrium, and including it makes
        the distribution bimodal and its mean meaningless, so only the trailing fraction is used.
        """
        tail_frac = float(self.specs.get('histogram-tail', 0.5))
        tail_frac = min(max(tail_frac, 0.0), 1.0)
        figsize = self.specs.get('figsize') or (9, 6)
        axis_labels = self.specs.get('axis-labels', {})
        for name in histograms:
            key, df, units, unitspec = self._resolve_units(name, df_of_column)
            if df is None:
                logger.debug(f'No data found for histogram {name}. Skipping...')
                continue
            values = (df[key] * units).dropna()
            if len(values) < 4:
                logger.debug(f'Too few samples ({len(values)}) to histogram {name}. Skipping...')
                continue
            tail = values.iloc[int(len(values) * (1.0 - tail_frac)):]
            fig, ax = plt.subplots(1, 1, figsize=figsize)
            nbins = max(10, min(int(np.sqrt(len(tail))), 60))
            ax.hist(tail, bins=nbins, color=self.colormap(0.5), edgecolor='white', linewidth=0.4)
            mean, sd = float(tail.mean()), float(tail.std(ddof=1)) if len(tail) > 1 else 0.0
            ax.axvline(mean, color='0.25', linestyle='-', linewidth=1.2)
            ax.axvspan(mean - sd, mean + sd, color='0.4', alpha=0.15, linewidth=0)
            label = to_latex_math(axis_labels.get(name, name))
            ax.set_xlabel(label + (f' ({unitspec})' if unitspec and unitspec != '*' else ''))
            ax.set_ylabel('count')
            ax.set_title(f'{label}: '
                         f'mean {mean:.4g} $\\pm$ {sd:.3g} over the trailing '
                         f'{tail_frac:.0%} ({len(tail)} samples)', fontsize=9)
            if self.specs.get('grid', False):
                ax.grid(True, axis='y')
            png_path = os.path.join(output_dir, f'{self.basename}-{name}-hist.png')
            plt.savefig(png_path, bbox_inches='tight')
            self.register(png_path, key=f'{name}-histogram-plot', artifact_type=PNGImageFileArtifact, keep=True)
            plt.clf()
            logger.debug(f'{self.taskname}: histogram of {name} -> {png_path}')

    def _draw_stage_markers(self, ax, key, df, x_values):
        """Draw a labeled vertical rule where each simulation stage begins.

        A plot normally spans several NAMD runs concatenated into one curve -- a minimization, a
        warm-up, a self-terminating equilibration -- and without the boundaries a reader cannot
        tell which part of the trace is which, or that a discontinuity is a change of ensemble
        rather than a physical event.

        Boundaries were recorded as row indices while the series was assembled, so they map onto
        whichever x axis is in use.  The first stage is skipped (its rule would sit on the axis)
        and labels are dropped when there are too many stages to read, which is the usual case for
        a chunked task -- there the boundaries themselves are the useful signal.
        """
        series_key = None
        for name in self.dataframes:
            if self.dataframes[name] is df:
                series_key = name
                break
        boundaries = self.stage_boundaries.get(series_key, [])
        if len(boundaries) < 2:
            return
        xv = list(x_values)
        label_them = len(boundaries) <= _MAX_LABELED_STAGES
        ymin, ymax = ax.get_ylim()
        # Stages are wildly uneven in length -- a 1 ps minimization sits beside a 160 ps
        # equilibration -- so labeling every boundary piles them on top of each other at the left
        # edge.  Label only where there is room; the rules still mark every boundary.
        x_span = (max(xv) - min(xv)) if len(xv) > 1 else 0.0
        min_gap = _MIN_LABEL_GAP_FRAC * x_span
        last_labeled_x = None
        n_marked = n_labeled = 0
        for row, label in boundaries[1:]:                 # the first stage starts at the axis
            if row <= 0 or row >= len(xv):
                continue
            x = xv[row]
            ax.axvline(x, color='0.55', linestyle='--', linewidth=0.8, zorder=0)
            n_marked += 1
            if not label_them:
                continue
            if last_labeled_x is not None and x_span > 0 and (x - last_labeled_x) < min_gap:
                continue
            ax.annotate(label, xy=(x, ymax), xytext=(2, -2), textcoords='offset points',
                        rotation=90, va='top', ha='left', fontsize=6, color='0.35', zorder=0)
            last_labeled_x = x
            n_labeled += 1
        logger.debug(f'{self.taskname}: marked {n_marked} stage boundary/boundaries on the '
                     f'{series_key} plot ({n_labeled} labeled)')

    def _collect_from_md_history(self) -> bool:
        """Fill ``self.csvartifacts`` from the pipeline's MD history.  Returns False if there is
        no history to read, leaving the caller to fall back to the old prior-task walk.

        The history is every NAMD run in order, so the question this answers is *where to start*:
        which runs describe the system we are holding now.  That is decided by atom count, not by
        task adjacency -- a `solvate` between two MD stages makes everything before it a different
        system, while a `validate` or a `ring_check` between them changes nothing.  Walking the
        history backwards while the atom count matches gets the first case right (the pre-solvation
        minimize is correctly dropped) and, unlike the adjacency rule it replaces, the second one
        too (an intervening non-MD task no longer truncates the plot).
        """
        history = self.get_current_artifact('md_history')
        if history is None or not history.data:
            return False
        entries = list(history.data)
        state: StateArtifacts = self.get_current_artifact('state')
        natoms = None
        if state is not None and getattr(state, 'psf', None) is not None:
            try:
                natoms = total_atoms(state.psf.path)
            except Exception as e:
                logger.debug(f'{self.taskname}: could not count atoms of the current state ({e})')
        keep = select_commensurable_runs(entries, natoms)
        if not keep:
            return False
        dropped = len(entries) - len(keep)
        stages = {(e['controller'], e['task_index']) for e in keep}
        logger.debug(f'MD history: {len(keep)} run(s) across {len(stages)} stage(s) match the '
                     f'current system ({natoms} atoms)'
                     + (f'; {dropped} earlier run(s) describe a different system' if dropped else ''))
        for entry in keep:
            for key in entry.get('csv_keys', []):
                artifact = CSVDataFileArtifact(f"{entry['csv_basename']}-{key}", key=f'{key}-csv')
                if artifact.exists():
                    self.csvartifacts.append(artifact)
                else:
                    logger.debug(f'MD history references a missing CSV: {artifact.name}')
        return len(self.csvartifacts) > 0

    def do(self):
        self.next_basename()
        my_logger(self.specs, logger.debug)
        output_dir = self.specs.get('output_dir', 'mdplots')
        os.makedirs(output_dir, exist_ok=True)
        self.reprocess_logs = self.specs.get('reprocess-logs', False)
        self.explicit_logs = self.specs.get('logs', [])
        self.running_sums = self.specs.get('running_sums', ['cpu_time', 'wall_time'])
        self.dataframes: dict[str, pd.DataFrame] = {}
        # series name -> [(row index where a stage starts, label), ...]
        self.stage_boundaries: dict[str, list] = {}
        self.colormapname = self.specs.get('colormap', 'tab10')
        self.colormap = cmaps.get(self.colormapname, cmaps['tab10'])
        self.colormap_direction = self.specs.get('colormap-direction', 1)
        # task list
        self.priortasklist: list[BaseTask] = []
        priortaskpointer = getattr(self, 'prior', None)
        while priortaskpointer is not None and isinstance(priortaskpointer, MDTask):
            self.priortasklist.append(priortaskpointer)
            priortaskpointer = priortaskpointer.prior
        self.priortasklist = self.priortasklist[::-1]
        logger.debug(f'Found {len(self.priortasklist)} prior MD tasks: {[pt.index for pt in self.priortasklist]}.')
        self.csvartifacts = FileArtifactList([])
        if self.reprocess_logs:
            logger.debug(f'Reprocessing logs: {self.reprocess_logs}')
            namdlog_objs = []
            if self.explicit_logs:
                for f in self.explicit_logs:
                    logger.debug(f'Extracting data from {f}')
                    the_log = NAMDLogParser.from_file(f)
                    csvs_generated = the_log.write_csv()
                    for key in csvs_generated:
                        artifact = self.register(csvs_generated[key], key=f'{key}-csv', artifact_type=CSVDataFileArtifact) 
                        if artifact is None:
                            raise PestiferBuildError(f'CSV file {csvs_generated[key]} does not exist.')
                        else:
                            self.csvartifacts.append(artifact)
                    namdlog_objs.append(the_log)
        elif self._collect_from_md_history():
            pass
        elif len(self.priortasklist) > 0:
            # fallback for a pipeline with no recorded history (e.g. artifacts carried in from an
            # older run): the original adjacency walk, which stops at the first non-MD task
            logger.debug(f'No MD history recorded; falling back to the prior-task walk: '
                         f'{[pt.index for pt in self.priortasklist]}')
            for pt in self.priortasklist:
                logger.debug(f'Extracting data from prior task {pt.taskname}-{pt.index}')
                artifactfile_collection = list(pt.get_my_artifactfile_collection().filter_by_artifact_type(CSVDataFileArtifact))
                for pt_artifact in artifactfile_collection:
                    self.csvartifacts.append(pt_artifact)
        else:
            raise PestiferBuildError('No CSV artifacts found in prior tasks.')

        if len(self.csvartifacts) == 0:
            raise PestiferBuildError('No CSV artifacts found.  Cannot extract time series data.')
        logger.debug(f'Found {len(self.csvartifacts)} CSV artifacts.')

        self.all_time_series_names = list(set([art.key.replace('-csv', '') for art in self.csvartifacts]))
        logger.debug(f'Found {len(self.all_time_series_names)} time series names: {self.all_time_series_names}')
    
        for tst in self.all_time_series_names:
            if not tst in self.dataframes:
                self.dataframes[tst] = pd.DataFrame()
            csv_artifact_collection = list(self.csvartifacts.filter_by_key(tst + '-csv'))
            if len(csv_artifact_collection) == 0:
                logger.debug(f'No CSV artifact found for {tst}. Skipping...')
                continue
            last_stage = None
            for csvartifact in csv_artifact_collection:
                logger.debug(f'Collecting data from CSV file {csvartifact.path.name}')
                csvname = csvartifact.path.name
                # note where each stage starts, as a row index into the series being assembled;
                # row indices survive the later conversion of the x axis to simulation time
                stage_key, stage_label = stage_of(csvartifact.name[:-len(f'-{tst}.csv')]
                                                  if csvartifact.name.endswith(f'-{tst}.csv')
                                                  else csvartifact.name)
                if stage_key is not None and stage_key != last_stage:
                    self.stage_boundaries.setdefault(tst, []).append(
                        (len(self.dataframes[tst]), stage_label))
                    last_stage = stage_key
                try:
                    newdf = pd.read_csv(csvname, header=0, index_col=None)
                except:
                    logger.debug(f'For some reason, I could not read this csv file')
                    logger.debug(f'{csvname} does not exist or is not a valid csv file')
                    continue
                # logger.debug(f'newdf shape: {newdf.shape}')
                # show the range of the first column
                # if not newdf.empty:
                #     logger.debug(f' -> first column range: {newdf.iloc[:,0].min()} - {newdf.iloc[:,0].max()}')
                if not self.dataframes[tst].empty and any(col in newdf.columns for col in self.running_sums):
                    # shift any columns designated as running sums by the final value of the previous dataframe
                    for col in self.running_sums:
                        if col in newdf.columns:
                            newdf[col] = newdf[col] + self.dataframes[tst][col].iloc[-1]
                if not self.dataframes[tst].empty and self.dataframes[tst].iloc[-1,0] == newdf.iloc[0,0]:
                    # logger.debug(f'Dropping first row of newdf to avoid duplicate time step {newdf.iloc[0,0]}')
                    # drop the first row of newdf to avoid duplicate time step
                    newdf = newdf.iloc[1:,:]
                self.dataframes[tst] = pd.concat([self.dataframes[tst], newdf], ignore_index=True)
                logger.debug(f'{tst} dataframe shape: {self.dataframes[tst].shape}')

        # save each new dataframe to a csv file
        for key,df in self.dataframes.items():
            if not df.empty:
                csvname=f'{self.basename}-{key}.csv'
                df.to_csv(csvname,index=False)
                self.register(f'{self.basename}-{key}', key=f'{key}-csv', artifact_type=CSVDataFileArtifact)

        # compute time_ps for the energy dataframe after all CSVs are concatenated;
        # using cumsum(dt_fs * delta_TS) correctly handles chained runs where dt may vary
        if 'energy' in self.dataframes and 'dt_fs' in self.dataframes['energy'].columns:
            df_e = self.dataframes['energy']
            ts_arr = df_e['TS'].values.astype(float)
            dt_arr = df_e['dt_fs'].values
            delta_ts = np.diff(ts_arr, prepend=ts_arr[0])
            df_e['time_ps'] = np.cumsum(dt_arr * delta_ts / 1000.0)
            self.dataframes['energy'] = df_e

        self._add_derived_columns()

        # build a dictionary of column headings:dataframe pairs
        df_of_column = {}
        for key, df in self.dataframes.items():
            for col in df.columns[1:]: # ignore first column, which is usually 'TS' or 'step'
                if col not in df_of_column:
                    df_of_column[col] = df
                else:
                    logger.debug(f'Column {col} found in multiple dataframes')

        timeseries = self.specs.get('timeseries', [])
        time_step_column_names = self.specs.get('time_step_column_names', ['TS', 'step', 'steps'])
        histograms = self.specs.get('histograms', [])
        profiles = self.specs.get('profiles', [])
        if len(profiles) > 0:
            # must be sure c_z is tracked in parallel with profiles
            has_cz = 'xst' in self.dataframes
            if not has_cz:
                logger.debug('not tracking box depth, so depth profiles will not be plotted')
                profiles = []
            else:
                profiles = [p for p in profiles if f'{p}profile' in self.dataframes]
                if not profiles:
                    logger.debug('no profile dataframes available; skipping profile plots')
                else:
                    if self.dataframes['xst'].shape[0] != self.dataframes[f'{profiles[0]}profile'].shape[0]:
                        logger.debug(f'xst: {self.dataframes["xst"].shape[0]} rows, {profiles[0]}profile: {self.dataframes[f"{profiles[0]}profile"].shape[0]} rows')
                    dp = self.dataframes['xst']
                    pp = self.dataframes[f'{profiles[0]}profile']
                    pp['c_z'] = dp[dp['TS'].isin(pp['TS'])]['c_z'].values
                    self.dataframes[f'{profiles[0]}profile'] = pp
        profiles_per_block = self.specs.get('profiles-per-block', 100)
        legend = self.specs.get('legend', True)
        stage_markers = self.specs.get('stage-markers', True)
        block_average = int(self.specs.get('block-average', -1))
        grid = self.specs.get('grid', False)
        if histograms:
            self._plot_histograms(histograms, df_of_column, output_dir)
        logger.debug(f'Timeseries to plot: {timeseries}')
        for trace in timeseries:
            unitspecs = []
            figsize = self.specs.get('figsize') or (9, 6)
            fig, ax = plt.subplots(1, 1, figsize=figsize)
            if type(trace) != list:
                tracelist = [trace]
            else:
                tracelist = trace

            # Determine time-axis settings once per figure using the first
            # trace whose dataframe carries time_ps.
            time_unit = None   # 'ps' or 'ns', or None → fall back to TS index
            time_scale = 1.0
            eff_dt_scaled = None  # dt in time_unit per TS step (for secondary axis)
            ts_offset = 0.0
            for t_i in tracelist:
                df_check = df_of_column.get(t_i.upper(), df_of_column.get(t_i.lower()))
                if df_check is not None and 'time_ps' in df_check.columns:
                    max_time_ps = float(df_check['time_ps'].max())
                    if max_time_ps >= 1000.0:
                        time_unit, time_scale = 'ns', 1e-3
                    else:
                        time_unit, time_scale = 'ps', 1.0
                    ts_arr = df_check['TS'].values.astype(float)
                    ts_offset = float(ts_arr[0])
                    ts_range = float(ts_arr[-1] - ts_arr[0])
                    if ts_range > 0:
                        eff_dt_ps = max_time_ps / ts_range
                    elif 'dt_fs' in df_check.columns:
                        eff_dt_ps = float(df_check['dt_fs'].iloc[0]) / 1000.0
                    else:
                        eff_dt_ps = None
                    if eff_dt_ps and eff_dt_ps > 0:
                        eff_dt_scaled = eff_dt_ps * time_scale
                    break

            boundary_series = None
            for idx, t_i in enumerate(tracelist):
                key, df, units, unitspec = self._resolve_units(t_i, df_of_column)
                unitspecs.append(unitspec)
                if df is None:
                    logger.debug(f'No data found for trace {t_i}. Skipping...')
                    continue
                if time_unit is not None and 'time_ps' in df.columns:
                    x_values = df['time_ps'] * time_scale
                else:
                    time_step_column = None
                    for tcn in time_step_column_names:
                        if tcn in df.columns:
                            time_step_column = tcn
                            break
                    if time_step_column is None:
                        logger.debug(f'No time step column found for trace {t_i}. Skipping...')
                        continue
                    x_values = df[time_step_column]
                color = self.colormap(idx / max(1, len(tracelist)-1))
                if self.colormap_direction == -1 and len(tracelist) > 1:
                    color = self.colormap(1.0 - idx / max(1, len(tracelist)-1))
                label = key
                if label.endswith('_time'):
                    label = label.replace('_time', ' time')
                plot_label = label.title() if '_' not in label else r'$'+label+r'$'
                y_values = df[key] * units
                window = block_average
                if window < 0:                      # automatic: smooth only the noisy observables
                    window = auto_block_average(key, len(y_values))
                    if window:
                        logger.debug(f'{self.taskname}: {key} is smoothed automatically over '
                                     f'{window} samples (raw series kept visible)')
                if window > 1 and len(y_values) > window:
                    # the raw series stays visible but recedes; the eye follows the mean, and the
                    # band shows the fluctuation the mean was taken over
                    rolled = y_values.rolling(window, center=True, min_periods=1)
                    mean, sd = rolled.mean(), rolled.std().fillna(0.0)
                    ax.plot(x_values, y_values, color=color, alpha=0.25, linewidth=0.8)
                    ax.fill_between(x_values, mean - sd, mean + sd, color=color, alpha=0.20,
                                    linewidth=0)
                    ax.plot(x_values, mean, label=plot_label, color=color, linewidth=1.6)
                else:
                    ax.plot(x_values, y_values, label=plot_label, color=color)
                if stage_markers and boundary_series is None:
                    boundary_series = (key, df, x_values)

            if stage_markers and boundary_series is not None:
                self._draw_stage_markers(ax, *boundary_series)

            if time_unit is not None:
                ax.set_xlabel(f'simulation time ({time_unit})')
                if eff_dt_scaled is not None:
                    _ts0, _dt = ts_offset, eff_dt_scaled
                    ax_top = ax.secondary_xaxis(
                        'top',
                        functions=(
                            lambda t, ts0=_ts0, dt=_dt: ts0 + t / dt,
                            lambda ts, ts0=_ts0, dt=_dt: (ts - ts0) * dt,
                        )
                    )
                    ax_top.set_xlabel('time step')
            else:
                ax.set_xlabel('time step')
            axis_labels = self.specs.get('axis-labels', {})
            ax.set_ylabel(', '.join([to_latex_math(axis_labels.get(n, n)) + ' (' + u + ')' for n, u in zip(tracelist, unitspecs)]))
            if legend and len(tracelist) > 1:
                ax.legend()
            if grid:
                ax.grid(True)
            tracename = '-'.join(tracelist)
            png_path = os.path.join(output_dir, f'{self.basename}-{tracename}.png')
            plt.savefig(png_path, bbox_inches='tight')
            self.register(png_path, key=f'{tracename}-timeseries-plot', artifact_type=PNGImageFileArtifact, keep=True)
            plt.clf()
        for profile in profiles:
            if profile == 'pressure':
                df = self.dataframes.get('pressureprofile', None)
                if df is None or dp is None:
                    logger.debug(f'No pressure profile data found. Skipping...')
                    continue
                if not df.empty:
                    unitspec = self.specs.get('units', {}).get('pressure', 'bar')
                    if unitspec == '*':
                        units = 1.0
                    else:
                        if unitspec in ['bar', 'atm']:
                            units = 1.0
                        elif unitspec in ['Pa', 'pascal']:
                            units = 1e-5
                        elif unitspec in ['kPa', 'kilopascal']:
                            units = 1e-2
                        elif unitspec in ['MPa', 'megapascal']:
                            units = 1.0
                        elif unitspec in ['GPa', 'gigapascal']:
                            units = 1e3
                        else:
                            logger.debug(f'Unitspec "{unitspec}" not recognized.')
                            units = 1.0
                    figsize = self.specs.get('figsize') or (21, 6)
                    fig, ax = plt.subplots(1, 3, figsize=figsize, sharey=True)
                    ax[0].set_xlabel(r'$\frac{1}{2}$($P_{xx} + P_{yy}$) '+f'({unitspec})')
                    ax[0].set_ylabel('z (Å)')
                    ax[1].set_xlabel(r'$P_{zz}$ '+f'({unitspec})')
                    ax[2].set_xlabel(r'$\frac{1}{3}(P_{xx} + P_{yy} + P_{zz})$ '+f'({unitspec})')
                    nprofiles = df.shape[0]
                    nblocks = nprofiles // profiles_per_block
                    if nprofiles % profiles_per_block > 0:
                        nblocks += 1
                    if nblocks == 0:
                        nblocks = 1
                    logger.debug(f'Number of profiles: {nprofiles}, profiles per block: {profiles_per_block}, nblocks: {nblocks}')
                    for i in range(nblocks):
                        start = i * profiles_per_block
                        end = start + profiles_per_block
                        if end > nprofiles:
                            end = nprofiles
                        if start < nprofiles:
                            pprofilex = df.iloc[start:end, 1:-1:3].mean(axis=0).to_numpy()
                            pprofiley = df.iloc[start:end, 2:-1:3].mean(axis=0).to_numpy()
                            pprofilexy = (pprofilex + pprofiley) / 2.0
                            pprofilez = df.iloc[start:end, 3:-1:3].mean(axis=0).to_numpy()
                            pprofile = (pprofilex + pprofiley + pprofilez) / 3.0
                            pprofilex_std = df.iloc[start:end, 1:-1:3].std(axis=0).to_numpy()
                            pprofiley_std = df.iloc[start:end, 2:-1:3].std(axis=0).to_numpy()
                            pprofilez_std = df.iloc[start:end, 3:-1:3].std(axis=0).to_numpy()
                            profilexy_std = (pprofilex_std + pprofiley_std) / 2.0
                            profile_std = (pprofilex_std + pprofiley_std + pprofilez_std) / 3.0
                            # std is NaN for a single-profile block (ddof=1); compute the
                            # range nan-safely and only set explicit limits when finite,
                            # otherwise let matplotlib autoscale (a NaN limit raises)
                            with np.errstate(invalid='ignore'):
                                lo = np.nanmin([pprofilex.min() - pprofilex_std.min(),
                                                pprofiley.min() - pprofiley_std.min(),
                                                pprofilez.min() - pprofilez_std.min()])
                                hi = np.nanmax([pprofilex.max() + pprofilex_std.max(),
                                                pprofiley.max() + pprofiley_std.max(),
                                                pprofilez.max() + pprofilez_std.max()])
                            xmin = (np.round(lo, 0) // 500 - 1) * 500 * units
                            xmax = (np.round(hi, 0) // 500 + 1) * 500 * units
                            if np.isfinite(xmin) and np.isfinite(xmax) and xmin < xmax:
                                ax[0].set_xlim(xmin, xmax)
                                ax[1].set_xlim(xmin, xmax)
                                ax[2].set_xlim(xmin, xmax)
                            # get an average box depth for this time interval
                            Lz = df['c_z'].iloc[start:end].mean(axis=0)
                            logger.debug(f'Lz {Lz}')
                            dprofile = np.linspace(0, Lz.max(), df.shape[1]//3)
                            # plot pressure and stdev along x-axis with shaded region and slab index on y axis in reverse numerical order
                            color = self.colormap(i / max(1, nblocks-1))
                            if self.colormap_direction == -1:
                                color = self.colormap(1.0 - i / max(1, nblocks-1))
                            ax[0].plot(pprofilexy * units, dprofile, label=f'TS {df.iloc[start, 0]}-{df.iloc[end - 1, 0]}', color=color)
                            ax[0].fill_betweenx(dprofile, 
                                             pprofilexy * units - profilexy_std * units, 
                                             pprofilexy * units + profilexy_std * units, 
                                             alpha=0.2,
                                             color=color)
                            ax[1].plot(pprofilez * units, dprofile, label=f'TS {df.iloc[start, 0]}-{df.iloc[end - 1, 0]}', color=color)
                            ax[1].fill_betweenx(dprofile, 
                                             pprofilez * units - pprofilez_std * units, 
                                             pprofilez * units + pprofilez_std * units, 
                                             alpha=0.2,
                                             color=color)
                            ax[2].plot(pprofile * units, dprofile, label=f'TS {df.iloc[start, 0]}-{df.iloc[end - 1, 0]}', color=color)
                            ax[2].fill_betweenx(dprofile, 
                                             pprofile * units - profile_std * units, 
                                             pprofile * units + profile_std * units, 
                                             alpha=0.2,
                                             color=color)
                    ax[0].axvline(x=0, color='k', linestyle='--')
                    ax[1].axvline(x=0, color='k', linestyle='--')
                    if legend:
                        # place a legend only in ax[2] but outside the plot area
                        ax[2].legend(loc='upper left', bbox_to_anchor=(-0.5, 1))
                    if grid:
                        ax[0].grid(True)
                        ax[1].grid(True)
                        ax[2].grid(True)
                    png_path = os.path.join(output_dir, f'{self.basename}-pressureprofile.png')
                    plt.savefig(png_path, bbox_inches='tight')
                    self.register(png_path, key='pressureprofile-plot', artifact_type=PNGImageFileArtifact, keep=True)
                    plt.clf()
                else:
                    logger.debug(f'No pressure profile data.  Skipping...')
                    continue
            else:
                logger.debug(f'Profile {profile} not recognized.  Skipping...')
                continue
        self.result = 0
        return self.result
    