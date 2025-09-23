# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`MDPlotTask` class for making plots of energy-like quantities from NAMD runs.
This class is a descendant of the :class:`BaseTask <pestifer.core.basetask.BaseTask>` class and is used to extract energy-like data from NAMD log files,
pressure profiles, and XST files, and to generate plots based on this data.
It handles the collection of energy data from multiple NAMD runs, creates CSV files for energy and pressure profile data,
and generates plots for specified traces and profiles.
The plots can include energy traces, pressure profiles, and histograms, with options for units, legends, and grid lines.
It also supports the extraction of data from existing NAMD log files and XST files, allowing for flexible data visualization.

Usage as a task in a build workflow is described in the :ref:`config_ref tasks mdplot` documentation.  This module is also used in standalone form by the :ref:`subs_mdplot` command.

"""
import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from matplotlib import colormaps as cmaps

from .basetask import BaseTask
from .mdtask import MDTask
from ..util.units import g_per_amu, A3_per_cm3
from ..logparsers import NAMDLogParser
from ..util.stringthings import to_latex_math
from ..core.artifacts import *

logger = logging.getLogger(__name__)
logging.getLogger("matplotlib").setLevel(logging.WARNING)

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
    def do(self):
        self.next_basename()
        my_logger(self.specs, logger.debug)
        self.reprocess_logs = self.specs.get('reprocess-logs', False)
        self.explicit_logs = self.specs.get('logs', [])
        self.running_sums = self.specs.get('running_sums', ['cpu_time', 'wall_time'])
        self.dataframes: dict[str, pd.DataFrame] = {}
        self.colormapname = self.specs.get('colormap', 'plasma')
        self.colormap = cmaps.get(self.colormapname, cmaps['plasma'])
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
                            raise FileNotFoundError(f'CSV file {csvs_generated[key]} does not exist.')
                        else:
                            self.csvartifacts.append(artifact)
                    namdlog_objs.append(the_log)
        elif len(self.priortasklist) > 0:
            logger.debug(f'Extracting data from prior tasks: {[pt.index for pt in self.priortasklist]}')
            for pt in self.priortasklist:
                logger.debug(f'Extracting data from prior task {pt.taskname}-{pt.index}')
                artifactfile_collection = list(pt.get_my_artifactfile_collection().filter_by_artifact_type(CSVDataFileArtifact))
                for pt_artifact in artifactfile_collection:
                    self.csvartifacts.append(pt_artifact)
        else:
            raise ValueError('No CSV artifacts found in prior tasks.')

        if len(self.csvartifacts) == 0:
            raise ValueError('No CSV artifacts found.  Cannot extract time series data.')
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
            for csvartifact in csv_artifact_collection:
                logger.debug(f'Collecting data from CSV file {csvartifact.path.name}')
                csvname = csvartifact.path.name
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
                if self.dataframes['xst'].shape[0] != self.dataframes[f'{profiles[0]}profile'].shape[0]:
                    logger.debug(f'xst: {self.dataframes["xst"].shape[0]} rows, {profiles[0]}profile: {self.dataframes[f"{profiles[0]}profile"].shape[0]} rows')
                    # transfer the c_z column from the xst dataframe to the profile dataframe keeping the TS columns aligned
                dp = self.dataframes['xst']
                pp = self.dataframes[f'{profiles[0]}profile']
                # logger.debug(f'pp[TS]: {pp["TS"].values}')
                # logger.debug(f'dp[TS]: {dp[dp["TS"].isin(pp["TS"])]["TS"].values}')
                pp['c_z'] = dp[dp['TS'].isin(pp['TS'])]['c_z'].values
                self.dataframes[f'{profiles[0]}profile'] = pp
        profiles_per_block = self.specs.get('profiles_per_block', 100)
        legend = self.specs.get('legend', False)
        grid = self.specs.get('grid', False)
        if histograms:
            logger.debug(f'Histograms are not yet implemented in the mdplot task.')
        logger.debug(f'Timeseries to plot: {timeseries}')
        for trace in timeseries:
            unitspecs = []
            figsize = self.specs.get('figsize', (9, 6))
            fig, ax = plt.subplots(1, 1, figsize=figsize)
            if type(trace) != list:
                tracelist = [trace]
            else:
                tracelist = trace
            for idx, t_i in enumerate(tracelist):
                units = 1.0
                unitspec = self.specs.get('units', {}).get(t_i, '*')
                if unitspec == '*':
                    units = 1.0
                else:
                    if t_i == 'density':
                        if unitspec in ['g_per_cc', 'g/cc', 'g_per_cm3', 'g/cm3']:
                            units = g_per_amu * A3_per_cm3
                        else:
                            logger.debug(f'Unitspec "{unitspec}" not recognized.')
                            units = 1.0
                unitspecs.append(unitspec)
                key = t_i.upper()
                df = df_of_column.get(key, None)
                if df is None:
                    # try the lowercase key
                    key = t_i.lower()
                    df = df_of_column.get(key, None)
                if df is None:
                    logger.debug(f'No data found for trace {t_i}. Skipping...')
                    continue
                # determine the time step column name
                time_step_column = None
                for tcn in time_step_column_names:
                    if tcn in df.columns:
                        time_step_column = tcn
                        break
                if time_step_column is None:
                    logger.debug(f'No time step column found for trace {t_i}. Skipping...')
                    continue
                color = self.colormap(idx / max(1, len(tracelist)-1))
                if self.colormap_direction == -1 and len(tracelist) > 1:
                    color = self.colormap(1.0 - idx / max(1, len(tracelist)-1))
                label = key
                if label.endswith('_time'):
                    label = label.replace('_time', ' time')
                ax.plot(df[time_step_column], df[key] * units, label=label.title() if '_' not in label else r'$'+label+r'$', color=color)
            ax.set_xlabel('time step')
            axis_labels = self.specs.get('axis-labels', {})
            ax.set_ylabel(', '.join([to_latex_math(axis_labels.get(n, n)) + ' (' + u + ')' for n, u in zip(tracelist, unitspecs)]))
            if legend:
                ax.legend()
            if grid:
                ax.grid(True)
            tracename = '-'.join(tracelist)
            plt.savefig(f'{self.basename}-{tracename}.png', bbox_inches='tight')
            self.register(f'{self.basename}-{tracename}', key=f'{tracename}-timeseries-plot', artifact_type=PNGImageFileArtifact)
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
                    figsize = self.specs.get('figsize', (21, 6))
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
                            xmin = (np.round(min(pprofilex.min()-pprofilex_std.min(), pprofiley.min()-pprofiley_std.min(), pprofilez.min()-pprofilez_std.min()),0)//500-1)*500*units
                            xmax = (np.round(max(pprofilex.max()+pprofilex_std.max(), pprofiley.max()+pprofiley_std.max(), pprofilez.max()+pprofilez_std.max()),0)//500+1)*500*units
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
                    plt.savefig(f'{self.basename}-pressureprofile.png', bbox_inches='tight')
                    self.register(f'{self.basename}-pressureprofile.png', key='pressureprofile-plot', artifact_type=PNGImageFileArtifact)
                    plt.clf()
                else:
                    logger.debug(f'No pressure profile data.  Skipping...')
                    continue
            else:
                logger.debug(f'Profile {profile} not recognized.  Skipping...')
                continue
        self.result = 0
        return self.result
    