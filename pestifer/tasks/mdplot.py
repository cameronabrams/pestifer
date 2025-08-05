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
import pandas as pd

from .basetask import BaseTask
from .md import MDTask
from ..util.units import g_per_amu,A3_per_cm3
from ..util.logparsers import NAMDLog
from ..core.stringthings import to_latex_math
from ..core.artifacts import PNGImageFile, CSVDataFile

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
        is_post_processing = self.specs.get('postprocessing', False)
        running_sums = self.specs.get('running_sums', ['cpu_time', 'wall_time'])
        dataframes = {}
        if not is_post_processing and self.prior:
            # this is an in-build task, so we collect time series data from each prior md task's csv files
            logger.debug(f'Collecting all output data from upstream md tasks...')
            priortaskpointer = self.prior
            priortasklist = []
            logger.debug(f'self={str(self)}; self.prior={str(priortaskpointer)}')
            while priortaskpointer!=None and isinstance(priortaskpointer, MDTask):
                priortasklist.append(priortaskpointer)
                priortaskpointer=priortaskpointer.prior
            priortasklist=priortasklist[::-1]
            for pt in priortasklist:
                artifactfile_collection = pt.get_my_artifactfile_collection()
                csv_artifacts=[a for a in artifactfile_collection if a.key.endswith('-csv')]
                if len(csv_artifacts)==0:
                    raise ValueError(f'No CSV artifacts found in {pt.basename}.  Cannot extract time series data.')
                time_series_titles=[x.key.replace('-csv','') for x in csv_artifacts]
                csv_artifact_dict={x:y for x,y in zip(time_series_titles,csv_artifacts)}
                logger.debug(f'Found {len(time_series_titles)} time series titles in {pt.basename}: {time_series_titles}')
                for tst in time_series_titles:
                    if not tst in dataframes:
                        dataframes[tst]=pd.DataFrame()
                    csvartifact=csv_artifact_dict.get(tst)
                    if csvartifact:
                        csvname=csvartifact.path.name
                        logger.debug(f'Collecting data from CSV file {csvname}')
                        try:
                            newdf=pd.read_csv(csvname,header=0,index_col=None)
                            logger.debug(f'newdf shape: {newdf.shape}')
                            # show the range of the first column
                            if not newdf.empty:
                                logger.debug(f' -> first column range: {newdf.iloc[:,0].min()} - {newdf.iloc[:,0].max()}')
                            if not dataframes[tst].empty and any(col in newdf.columns for col in running_sums):
                                # shift any columns designated as running sums by the final value of the previous dataframe
                                for col in running_sums:
                                    if col in newdf.columns:
                                        newdf[col]=newdf[col]+dataframes[tst][col].iloc[-1]
                            dataframes[tst]=pd.concat([dataframes[tst],newdf],ignore_index=True)
                            logger.debug(f'{tst} dataframe shape: {dataframes[tst].shape}')
                        except:
                            logger.debug(f'For some reason, I could not read this csv file')
                            logger.debug(f'{csvname} does not exist or is not a valid csv file')
        else:
            # this is being run as a standalone, post-processing task
            logger.debug(f'Extracting data from {len(self.specs["logs"])} explicitly named namd logs...')
            self.basename=self.taskname
            logger.debug(self.specs['logs'])
            namdlogs=[]
            for f in self.specs['logs']:
                logger.debug(f'Extracting data from {f}')
                namdlogs.append(NAMDLog.from_file(f))
                nl=namdlogs[-1]
                time_series_titles=list(nl.dataframes.keys())
                logger.debug(f'Found {len(time_series_titles)} time series titles in {nl.basename}: {time_series_titles}')
                for tst in time_series_titles:
                    if not tst in dataframes:
                        logger.debug(f'Creating dataframe for {tst}')
                        dataframes[tst]=pd.DataFrame()
                    newdf=nl.dataframes[tst]
                    if not dataframes[tst].empty and any(col in newdf.columns for col in running_sums):
                        # shift any columns designated as running sums by the final value of the previous dataframe
                        for col in running_sums:
                            if col in newdf.columns:
                                newdf[col]=newdf[col]+dataframes[tst][col].iloc[-1]
                    if not newdf.empty:
                        dataframes[tst]=pd.concat([dataframes[tst],newdf],ignore_index=True)
                        logger.debug(f'{tst} dataframe shape: {dataframes[tst].shape}')

        # save each new dataframe to a csv file
        # TODO: make sure accumulating quantities increment
        for key,df in dataframes.items():
            if not df.empty:
                csvname=f'{self.basename}-{key}.csv'
                df.to_csv(csvname,index=False)
                self.register_current_artifact(CSVDataFile(f'{self.basename}-{key}'),key=f'{key}-csv')

        # build a dictionary of column headings:dataframe pairs
        df_of_column={}
        for key,df in dataframes.items():
            for col in df.columns[1:]: # ignore first column, which is usually 'TS' or 'step'
                if col not in df_of_column:
                    df_of_column[col]=df
                else:
                    logger.debug(f'Column {col} found in multiple dataframes')

        timeseries=self.specs.get('timeseries',[])
        time_step_column_names=self.specs.get('time_step_column_names',['TS','step','steps'])
        histograms=self.specs.get('histograms',[])
        profiles=self.specs.get('profiles',[])
        profiles_per_block=self.specs.get('profiles_per_block',100)
        legend=self.specs.get('legend',False)
        grid=self.specs.get('grid',False)
        if histograms:
            logger.debug(f'Histograms are not yet implemented in the mdplot task.')
        for trace in timeseries:
            unitspecs=[]
            figsize=self.specs.get('figsize',(9,6))
            fig,ax=plt.subplots(1,1,figsize=figsize)
            if type(trace)!=list:
                tracelist=[trace]
            else:
                tracelist=trace
            for t_i in tracelist:
                unitspec=self.specs.get('units',{}).get(t_i,'*')
                if unitspec=='*':
                    units=1.0
                else:
                    if t_i=='density':
                        if unitspec in ['g_per_cc','g/cc','g_per_cm3','g/cm3']:
                            units=g_per_amu*A3_per_cm3
                        else:
                            logger.debug(f'Unitspec "{unitspec}" not recognized.')
                            units=1.0
                unitspecs.append(unitspec)
                key=t_i.upper()
                df=df_of_column.get(key,None)
                if df is None:
                    # try the lowercase key
                    key=t_i.lower()
                    df=df_of_column.get(key,None)
                if df is None:
                    logger.debug(f'No data found for trace {t_i}. Skipping...')
                    continue
                # determine the time step column name
                time_step_column=None
                for tcn in time_step_column_names:
                    if tcn in df.columns:
                        time_step_column=tcn
                        break
                if time_step_column is None:
                    logger.debug(f'No time step column found for trace {t_i}. Skipping...')
                    continue
                ax.plot(df[time_step_column],df[key]*units,label=key.title())
            ax.set_xlabel('time step')
            axis_labels=self.specs.get('axis-labels',{})
            ax.set_ylabel(', '.join([to_latex_math(axis_labels.get(n,n))+' ('+u+')' for n,u in zip(tracelist,unitspecs)]))
            if legend:
                ax.legend()
            if grid:
                ax.grid(True)
            tracename='-'.join(tracelist)
            plt.savefig(f'{self.basename}-{tracename}.png',bbox_inches='tight')
            self.register_current_artifact(PNGImageFile(f'{self.basename}-{tracename}'),key=f'{tracename}-timeseries-plot')
            plt.clf()
        for profile in profiles:
            if profile=='pressureprofile':
                df=dataframes.get('pressureprofile',None)
                if df is None:
                    logger.debug(f'No pressure profile data found. Skipping...')
                    continue
                if not df.empty:
                    unitspec=self.specs.get('units',{}).get('pressure','bar')
                    if unitspec=='*':
                        units=1.0
                    else:
                        if unitspec in ['bar','atm']:
                            units=1.0
                        elif unitspec in ['Pa','pascal']:
                            units=1e-5
                        elif unitspec in ['kPa','kilopascal']:
                            units=1e-2
                        elif unitspec in ['MPa','megapascal']:
                            units=1.0
                        elif unitspec in ['GPa','gigapascal']:
                            units=1e3
                        else:
                            logger.debug(f'Unitspec "{unitspec}" not recognized.')
                            units=1.0
                    figsize=self.specs.get('figsize',(9,6))
                    fig,ax=plt.subplots(1,1,figsize=figsize)
                    ax.set_xlabel(f'Pressure ({unitspec})')
                    ax.set_ylabel('z (Å)')
                    nprofiles=df.shape[0]
                    nblocks=nprofiles//profiles_per_block
                    if nprofiles%profiles_per_block>0:
                        nblocks+=1
                    if nblocks==0:
                        nblocks=1
                    logger.debug(f'Number of profiles: {nprofiles}, profiles per block: {profiles_per_block}, nblocks: {nblocks}')
                    for i in range(nblocks):
                        start=i*profiles_per_block
                        end=start+profiles_per_block
                        if end>nprofiles:
                            end=nprofiles
                        if start<nprofiles:
                            pprofile=df.iloc[start:end,1:]
                            # average all rows of pprofile
                            pressure_profile_mean=pprofile.mean(axis=0)
                            pressure_profile_stdev=pprofile.std(axis=0)
                            # make sure slab_index are integers of the column headings
                            pprofile.columns=pprofile.columns.astype(int)
                            # plot pressure and stdev along x-axis with shaded region and slab index on y axis in reverse numerical order
                            ax.plot(pressure_profile_mean*units,pprofile.columns,label=f'TS {df.iloc[start,0]}-{df.iloc[end-1,0]}')
                            ax.fill_betweenx(pprofile.columns, 
                                             pressure_profile_mean*units-pressure_profile_stdev*units, 
                                             pressure_profile_mean*units+pressure_profile_stdev*units, 
                                             alpha=0.2)
                    # # average all columns of pprofiles
                    # pprofiles['pressure']=pprofiles.mean(axis=1)
                    # pprofiles['stdev']=pprofiles.std(axis=1)
                    # # plot pressure and stdev along x-axis with shaded region and slab index on y axis in reverse numerical order
                    # ax.plot(pprofiles['pressure']*units,pprofiles.index,label='Pressure')
                    # ax.fill_betweenx(pprofiles.index, 
                    #                  pprofiles['pressure']*units-pprofiles['stdev']*units, 
                    #                  pprofiles['pressure']*units+pprofiles['stdev']*units, 
                    #                  alpha=0.2, label='Stdev')
                    # place a vertical line at x=0
                    ax.axvline(x=0, color='k', linestyle='--')
                    if legend:
                        ax.legend()
                    if grid:
                        ax.grid(True)
                    plt.savefig(f'{self.basename}-pressureprofile.png',bbox_inches='tight')
                    self.register_current_artifact(PNGImageFile(f'{self.basename}-pressureprofile.png'),key='pressureprofile-plot')
                    plt.clf()
                else:
                    logger.debug(f'No pressure profile data.  Skipping...')
                    continue
            else:
                logger.debug(f'Profile {profile} not recognized.  Skipping...')
                continue
        self.result=0
        return self.result
    