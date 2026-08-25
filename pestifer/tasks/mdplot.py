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
from ..util.provenance import stamp_figure, stamp as provenance_stamp
from ..util.density_convergence import total_atoms
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
    'VOLUME': 'Å³', 'DENSITY': 'g/cc',
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
            # The DENSITY column is stored in g/cc by the log parser, so asking for g/cc is a no-op.
            # It used to be stored raw and converted here; converting again would double-apply the
            # factor for any config carrying this spec.
            if unitspec in ['g_per_cc', 'g/cc', 'g_per_cm3', 'g/cm3']:
                units = 1.0
                unitspec = 'g/cc'
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



    def _series_from_logs(self, logs: list) -> dict:
        """Parse NAMD logs into concatenated dataframes, without writing CSVs.

        Used for overlay sources, which are read only to be drawn: writing their CSVs would
        scatter files named after another run through this run's directory.
        """
        from ..logparsers import NAMDLogParser
        series = {}
        for f in logs:
            try:
                parsed = NAMDLogParser.from_file(f)
            except Exception as e:
                logger.debug(f'overlay: could not parse {f} ({e})')
                continue
            for key, df in getattr(parsed, 'dataframes', {}).items():
                if df is None or df.empty:
                    continue
                if key not in series:
                    series[key] = df.copy()
                    continue
                prev = series[key]
                if not prev.empty and not df.empty and prev.iloc[-1, 0] == df.iloc[0, 0]:
                    df = df.iloc[1:, :]        # same boundary timestep as the previous chunk
                series[key] = pd.concat([prev, df], ignore_index=True)
        for key, df in series.items():
            if 'dt_fs' in df.columns and 'TS' in df.columns:
                ts = df['TS'].values.astype(float)
                df['time_ps'] = np.cumsum(df['dt_fs'].values * np.diff(ts, prepend=ts[0]) / 1000.0)
                series[key] = df
        return series

    def _plot_overlays(self, overlays, quantities, df_of_column, output_dir, block_average, grid):
        """One figure per quantity, one curve per run.

        Stage markers are deliberately omitted here: the runs being compared have their own,
        differing stage structures, and drawing several sets of boundaries on one axes says less
        than it obscures.
        """
        sources = [('this run', df_of_column)]
        for entry in overlays:
            if not isinstance(entry, dict):
                logger.debug(f'overlay entry {entry!r} is not a dict; skipping')
                continue
            label, logs = entry.get('label'), entry.get('logs', [])
            if not logs:
                logger.debug(f'overlay {label!r} names no logs; skipping')
                continue
            series = self._series_from_logs(logs)
            if not series:
                logger.debug(f'overlay {label!r} yielded no data; skipping')
                continue
            columns = {}
            for df in series.values():
                for col in df.columns[1:]:
                    columns.setdefault(col, df)
            sources.append((label or 'overlay', columns))
        if len(sources) < 2:
            logger.debug('no usable overlay sources; skipping overlay figures')
            return
        axis_labels = self.specs.get('axis-labels', {})
        for name in quantities:
            figsize = self.specs.get('figsize') or (9, 6)
            fig, ax = plt.subplots(1, 1, figsize=figsize)
            drawn, unitspec_seen = 0, '*'
            for i, (label, columns) in enumerate(sources):
                key, df, units, unitspec = self._resolve_units(name, columns)
                if df is None:
                    logger.debug(f'overlay: {label} has no {name}')
                    continue
                x_values, x_desc = self._x_axis_for(df)
                if x_values is None:
                    continue
                y_values = df[key] * units
                window = block_average if block_average >= 0 else auto_block_average(key, len(y_values))
                color = self.colormap(i / max(1, len(sources) - 1))
                if window > 1 and len(y_values) > window:
                    y_values = y_values.rolling(window, center=True, min_periods=1).mean()
                ax.plot(x_values, y_values, label=label, color=color, linewidth=1.3)
                unitspec_seen = unitspec if unitspec != '*' else unitspec_seen
                drawn += 1
            if drawn < 2:
                plt.close(fig)
                logger.debug(f'overlay for {name}: fewer than two runs had data; skipping')
                continue
            ax.set_xlabel(x_desc or 'time step')
            lbl = to_latex_math(axis_labels.get(name, name))
            ax.set_ylabel(lbl + (f' ({unitspec_seen})' if unitspec_seen != '*' else ''))
            ax.legend()
            if grid:
                ax.grid(True)
            png_path = os.path.join(output_dir, f'{self.basename}-{name}-overlay.png')
            stamp_figure(fig, self.build_stamp())
            fig.savefig(png_path, bbox_inches='tight')
            self.register(png_path, key=f'{name}-overlay-plot', artifact_type=PNGImageFileArtifact, keep=True)
            plt.close(fig)
            logger.debug(f'{self.taskname}: overlay of {name} across {drawn} runs -> {png_path}')

    def _plot_timeseries(self, timeseries, df_of_column, output_dir, block_average,
                         stage_markers, legend, grid):
        """Draw one figure per requested trace group.

        Extracted from ``do`` so it matches its siblings (``_plot_panels``, ``_plot_overlays``,
        ``_plot_histograms``) and so its figure content can be asserted on without running a
        build: a caller can invoke it directly and inspect the axes it produced.
        """
        time_step_column_names = self.specs.get('time_step_column_names', ['TS', 'step', 'steps'])
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
            # Count what is actually drawn, not what was asked for.  A requested trace whose data
            # is absent is skipped below, and if *every* trace is skipped the figure is empty --
            # in which case a legend call warns, and saving would leave a blank PNG registered as
            # a build artifact with only a debug-level breadcrumb explaining why.
            plotted, missing = 0, []
            for idx, t_i in enumerate(tracelist):
                key, df, units, unitspec = self._resolve_units(t_i, df_of_column)
                unitspecs.append(unitspec)
                if df is None:
                    logger.debug(f'No data found for trace {t_i}. Skipping...')
                    missing.append(str(t_i))
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
                        missing.append(str(t_i))
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
                plotted += 1
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
            # Legend on what was drawn, not on what was requested: asking for a legend over
            # zero labelled artists is what emits matplotlib's "No artists with labels" warning.
            if legend and plotted > 1:
                ax.legend()
            if grid:
                ax.grid(True)
            tracename = '-'.join(tracelist)
            if plotted == 0:
                logger.warning(f'{self.taskname}: no data for {tracename} '
                               f'(missing: {", ".join(missing) or "unknown"}); no plot written')
                plt.clf()
                continue
            if missing:
                logger.warning(f'{self.taskname}: plotting {tracename} without '
                               f'{", ".join(missing)} -- no data for those traces')
            png_path = os.path.join(output_dir, f'{self.basename}-{tracename}.png')
            stamp_figure(plt.gcf(), self.build_stamp())
            plt.savefig(png_path, bbox_inches='tight')
            self.register(png_path, key=f'{tracename}-timeseries-plot', artifact_type=PNGImageFileArtifact, keep=True)
            plt.clf()
    def _plot_panels(self, panels, df_of_column, output_dir, block_average, stage_markers,
                     legend, grid):
        """Stacked panels sharing one x axis.

        Three quantities in three figures cannot be read together -- the eye cannot align three
        separate x axes -- so a shared-axis stack is what makes a claim like "the area kept
        relaxing after the density had settled" visible at all.  Each entry is one panel; a list
        entry overlays several traces in that panel.
        """
        groups = [p if isinstance(p, list) else [p] for p in panels]
        groups = [g for g in groups if g]
        if not groups:
            return
        figsize = self.specs.get('figsize') or (9, 2.6 * len(groups) + 1.2)
        fig, axes = plt.subplots(len(groups), 1, figsize=figsize, sharex=True, squeeze=False)
        axes = [a[0] for a in axes]
        axis_labels = self.specs.get('axis-labels', {})
        x_desc = None
        drawn = 0
        for ax, group in zip(axes, groups):
            unitspecs, boundary = [], None
            for idx, name in enumerate(group):
                key, df, units, unitspec = self._resolve_units(name, df_of_column)
                if df is None:
                    logger.debug(f'No data found for panel trace {name}. Skipping...')
                    unitspecs.append(unitspec)
                    continue
                x_values, this_desc = self._x_axis_for(df)
                if x_values is None:
                    continue
                x_desc = x_desc or this_desc
                y_values = df[key] * units
                window = block_average
                if window < 0:
                    window = auto_block_average(key, len(y_values))
                color = self.colormap(idx / max(1, len(group) - 1))
                if window > 1 and len(y_values) > window:
                    rolled = y_values.rolling(window, center=True, min_periods=1)
                    mean, sd = rolled.mean(), rolled.std().fillna(0.0)
                    ax.plot(x_values, y_values, color=color, alpha=0.25, linewidth=0.8)
                    ax.fill_between(x_values, mean - sd, mean + sd, color=color, alpha=0.20,
                                    linewidth=0)
                    ax.plot(x_values, mean, label=key, color=color, linewidth=1.6)
                else:
                    ax.plot(x_values, y_values, label=key, color=color)
                unitspecs.append(unitspec)
                boundary = boundary or (key, df, x_values)
                drawn += 1
            ylabel = ', '.join(to_latex_math(axis_labels.get(n, n)) +
                               (f' ({u})' if u and u != '*' else '')
                               for n, u in zip(group, unitspecs))
            ax.set_ylabel(ylabel, fontsize=9)
            if grid:
                ax.grid(True)
            if legend and len(group) > 1:
                ax.legend(fontsize=8)
            if stage_markers and boundary is not None:
                self._draw_stage_markers(ax, *boundary)
        if not drawn:
            plt.close(fig)
            logger.debug('no panel traces had data; skipping the panel figure')
            return
        axes[-1].set_xlabel(x_desc or 'time step')
        fig.align_ylabels(axes)
        fig.tight_layout()
        name = '_'.join('-'.join(g) for g in groups)
        png_path = os.path.join(output_dir, f'{self.basename}-panels-{name}.png')
        stamp_figure(fig, self.build_stamp())
        fig.savefig(png_path, bbox_inches='tight')
        self.register(png_path, key='panels-plot', artifact_type=PNGImageFileArtifact, keep=True)
        plt.close(fig)
        logger.debug(f'{self.taskname}: {len(groups)}-panel figure -> {png_path}')

    def _x_axis_for(self, df):
        """``(x values, axis label)`` for a dataframe: simulation time when it carries one."""
        if 'time_ps' in df.columns:
            max_ps = float(df['time_ps'].max())
            if max_ps >= 1000.0:
                return df['time_ps'] * 1e-3, 'simulation time (ns)'
            return df['time_ps'], 'simulation time (ps)'
        for tcn in self.specs.get('time_step_column_names', ['TS', 'step', 'steps']):
            if tcn in df.columns:
                return df[tcn], 'time step'
        return None, None

    def _write_summary_table(self, quantities, df_of_column, output_dir):
        """Per-stage summary: what a reader usually squints at a figure to extract.

        One row per simulation stage with its step range and, for each plotted quantity, the mean
        and standard deviation over that stage -- so "the density settled at 1.03 during the NPT
        stage" is a number that can be pasted or diffed rather than eyeballed off a curve.
        """
        rows = []
        for name in quantities:
            key, df, units, unitspec = self._resolve_units(name, df_of_column)
            if df is None:
                continue
            series_key = next((k for k, v in self.dataframes.items() if v is df), None)
            boundaries = self.stage_boundaries.get(series_key, [])
            if not boundaries:
                boundaries = [(0, series_key or 'all')]
            values = df[key] * units
            step_col = next((c for c in ('TS', 'step', 'steps') if c in df.columns), None)
            for i, (start, label) in enumerate(boundaries):
                end = boundaries[i + 1][0] if i + 1 < len(boundaries) else len(values)
                chunk = values.iloc[start:end].dropna()
                if chunk.empty:
                    continue
                row = dict(stage=label, quantity=name, units=unitspec,
                           samples=len(chunk), mean=float(chunk.mean()),
                           std=float(chunk.std(ddof=1)) if len(chunk) > 1 else 0.0)
                if step_col is not None:
                    row['first_step'] = int(df[step_col].iloc[start])
                    row['last_step'] = int(df[step_col].iloc[min(end, len(df)) - 1])
                rows.append(row)
        if not rows:
            logger.debug('no data for a summary table; skipping')
            return
        table = pd.DataFrame(rows)
        csvname = os.path.join(output_dir, f'{self.basename}-summary.csv')
        table.to_csv(csvname, index=False)
        self.register(csvname, key='summary-csv', artifact_type=CSVDataFileArtifact, keep=True)
        logger.info(f'{self.taskname}: per-stage summary ({len(rows)} rows) -> {csvname}')

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
            stamp_figure(plt.gcf(), self.build_stamp())
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

    def build_stamp(self) -> str:
        """The provenance mark for figures this task writes.

        Overridden because re-plotting is the one path where the running configuration is *not*
        the provenance of the data.  ``pestifer mdplot`` and ``reprocess-logs: true`` draw figures
        from logs some other run produced -- possibly on another machine, possibly years earlier --
        and stamping those with this configuration's seed would attribute the data to a run that
        never made it.  NAMD records the seed it used in its own log, so that is what is used when
        the logs agree on one; when they disagree (a build's chunks each get a derived seed) or
        record none, the mark falls back to the version alone rather than inventing a seed.
        """
        if getattr(self, 'reprocess_logs', False):
            seeds = {s for s in getattr(self, '_log_seeds', set()) if s}
            return provenance_stamp(seeds.pop() if len(seeds) == 1 else None)
        return super().build_stamp()

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
            # the seeds actually recorded by the runs being re-plotted -- see build_stamp
            self._log_seeds = set()
            if self.explicit_logs:
                for f in self.explicit_logs:
                    logger.debug(f'Extracting data from {f}')
                    the_log = NAMDLogParser.from_file(f)
                    self._log_seeds.add(the_log.metadata.get('random_number_seed'))
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
        self._plot_timeseries(timeseries, df_of_column, output_dir, block_average,
                              stage_markers, legend, grid)
        panels = self.specs.get('panels', [])
        if panels:
            self._plot_panels(panels, df_of_column, output_dir, block_average, stage_markers,
                              legend, grid)
        overlays = self.specs.get('overlay', [])
        if overlays:
            flat = []
            for t in timeseries:
                flat.extend(t if isinstance(t, list) else [t])
            self._plot_overlays(overlays, list(dict.fromkeys(flat)), df_of_column, output_dir,
                                block_average, grid)

        if self.specs.get('summary-table', False):
            summarized = []
            for t in timeseries:
                summarized.extend(t if isinstance(t, list) else [t])
            for p_ in panels:
                summarized.extend(p_ if isinstance(p_, list) else [p_])
            self._write_summary_table(list(dict.fromkeys(summarized)), df_of_column, output_dir)

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
                    stamp_figure(plt.gcf(), self.build_stamp())
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
    