# Author: Cameron F. Abrams, <cfa22@drexel.edu>
# 
import logging
import matplotlib.pyplot as plt
import pandas as pd
import os

from ..basetask import BaseTask
from ..util.units import g_per_amu,A3_per_cm3
from ..util.logparsers import NAMDLog, NAMDxst
from ..stringthings import to_latex_math

logger=logging.getLogger(__name__)
logging.getLogger("matplotlib").setLevel(logging.WARNING)

class MDPlotTask(BaseTask):
    """ A class for making plots of energy-like quantities from a series of one or more NAMD 
        runs.  Since NAMD runs are always invoked using a log-parser, a csv file is created that 
        contains all energy-like data from the run. """
    yaml_header='mdplot'
    def do(self):
        logger.debug(f'Running {self.__class__.__name__} task with specs {self.specs}')
        self.log_message('initiated')
        self.inherit_state()
        nrg=pd.DataFrame()
        pprofiles=pd.DataFrame()
        xst=None
        if len(self.specs['existing-logs'])==0 and self.prior:
            logger.debug(f'Collecting all output data from upstream md tasks...')
            priortaskpointer=self.prior
            priortasklist=[]
            logger.debug(f'self={str(self)}; self.prior={str(priortaskpointer)}')
            while priortaskpointer!=None and priortaskpointer.yaml_header=='md':
                priortasklist.append(priortaskpointer)
                priortaskpointer=priortaskpointer.prior
            priortasklist=priortasklist[::-1]
            for pt in priortasklist:
                logger.debug(f'Collecting data from {pt.basename}')
                ecsv=pt.statevars.get('energy-csv',None)
                ppcsv=pt.statevars.get('pressureprofile-csv',None)
                if ppcsv:
                    logger.debug(f'Collecting pressure profile data from {ppcsv}')
                    try:
                        newpp=pd.read_csv(ppcsv,header=0,index_col=0)
                        logger.debug(f'newpp shape: {newpp.shape}')
                        pprofiles=pd.concat([pprofiles,newpp],axis=1,ignore_index=False)
                        logger.debug(f'pprofiles dataframe shape: {pprofiles.shape}')
                    except:
                        logger.debug(f'For some reason, I could not read this csv file')
                        logger.debug(f'{ppcsv} does not exist or is not a valid csv file')
                else:
                    logger.debug(f'No pressure profile data found for {pt.basename}')
                if ecsv:
                    logger.debug(f'Collecting data from {ecsv}')
                    try:
                        newcsv=pd.read_csv(ecsv,header=0,index_col=None)
                        logger.debug(f'newcsv shape: {newcsv.shape}')
                        nrg=pd.concat([nrg,newcsv],ignore_index=True)
                        logger.debug(f'nrg dataframe shape: {nrg.shape}')
                    except:
                        logger.debug(f'For some reason, I could not read this csv file')
                        logger.debug(f'{ecsv} does not exist or is not a valid csv file')
                else:
                    logger.debug(f'No energy data found for {pt.basename}')
                xstfile=f'{pt.basename}.xst'
                if os.path.exists(xstfile):
                    logger.debug(f'Collecting data from {xstfile}')
                    try:
                        if not xst:
                            xst=NAMDxst(xstfile)
                        else:
                            xst.add_file(xstfile)
                    except:
                        logger.debug(f'{xstfile} does not exist or is not a valid xst file')
        else:
            logger.debug(f'Extracting data from {len(self.specs["existing-logs"])} explicitly named namd logs...')
            logger.debug(self.specs['existing-logs'])
            for f in self.specs['existing-logs']:
                logger.debug(f'Extracting data from {f}')
                apparent_basename=os.path.splitext(os.path.basename(f))[0]
                l=NAMDLog(basename=apparent_basename)
                l.static(f)
                nrg=pd.concat([nrg,l.energy_df])
                if not l.pressureprofile_df.empty:
                    pprofiles=pd.concat([pprofiles,l.pressureprofile_df],axis=1,ingore_index=False)
            if len(self.specs['existing-xsts'])>0:
                logger.debug(f'Extracting data from {len(self.specs["existing-xsts"])} explicitly named XST files...')
            for f in self.specs['existing-xsts']:
                if not xst:
                    xst=NAMDxst(f)
                else:
                    xst.add_file(f)

        basename=self.specs.get('basename','myplot')
        logger.debug(f'Saving NAMD log data to {basename}-energy.csv.')
        try:
            nrg.to_csv(f'{basename}-energy.csv',header=True,index=False)
        except:
            logger.debug(f'For some reason, I could not write this dataframe to csv')
            logger.debug(nrg.iloc[:3,:].to_string())
        if not pprofiles.empty:
            logger.debug(f'Saving pressure profile data to {basename}-pressureprofile.csv.')
            try:
                pprofiles.to_csv(f'{basename}-pressureprofile',header=True,index=True)
            except:
                logger.debug(f'For some reason, I could not write this dataframe to csv')
                logger.debug(pprofiles.iloc[:3,:].to_string())
        if xst is not None:
            logger.debug(f'Saving cell data to {basename}-cell.csv.')
            try:
                xst.df.to_csv(f'{basename}-cell',header=True,index=False)
            except:
                logger.debug(f'For some reason, I could not write this dataframe to csv')
                logger.debug(xst.df.iloc[:3,:].to_string())
        
        traces=self.specs.get('traces',[])
        histograms=self.specs.get('histograms',[])
        profiles=self.specs.get('profiles',[])
        profiles_per_block=self.specs.get('profiles_per_block',100)
        legend=self.specs.get('legend',False)
        grid=self.specs.get('grid',False)
        if histograms:
            logger.debug(f'Histograms are not yet implemented in the mdplot task.')
        for trace in traces:
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
                if t_i.upper() in nrg:
                    key=t_i.upper()
                    ax.plot(nrg['TS'],nrg[key]*units,label=key.title())
                elif xst is not None and t_i in xst.df:
                    ax.plot(xst.df['step'],xst.df[t_i]*units,label=to_latex_math(t_i))
            ax.set_xlabel('time step')
            ax.set_ylabel(', '.join([to_latex_math(n)+' ('+u+')' for n,u in zip(tracelist,unitspecs)]))
            if legend:
                ax.legend()
            if grid:
                ax.grid(True)
            tracename='-'.join(tracelist)
            plt.savefig(f'{basename}-{tracename}.png',bbox_inches='tight')
            plt.clf()
        for profile in profiles:
            if profile=='pressure':
                if not pprofiles.empty:
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
                    ax.set_ylabel('Slab Index')
                    nprofiles=len(pprofiles.columns)
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
                            pprofile=pprofiles.iloc[:,start:end]
                            # average all columns of pprofile
                            pressure_profile_mean=pprofile.mean(axis=1)
                            pressure_profile_stdev=pprofile.std(axis=1)
                            # plot pressure and stdev along x-axis with shaded region and slab index on y axis in reverse numerical order
                            ax.plot(pressure_profile_mean*units,pprofile.index,label=f'TS {pprofiles.columns[start]}-{pprofiles.columns[end-1]}')
                            ax.fill_betweenx(pprofile.index, 
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
                    plt.savefig(f'{basename}-pressureprofile.png',bbox_inches='tight')
                    plt.clf()
                else:
                    logger.debug(f'No pressure profile data.  Skipping...')
                    continue
            else:
                logger.debug(f'Profile {profile} not recognized.  Skipping...')
                continue
        self.log_message('complete')
        self.result=0
        return super().do()
    