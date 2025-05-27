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
        xst=None
        if len(self.specs['existing-logs'])==0 and self.prior:
            logger.debug(f'Collecting all mdlog.edata tables from upstream md tasks...')
            priortaskpointer=self.prior
            priortasklist=[]
            logger.debug(f'self={str(self)}; self.prior={str(priortaskpointer)}')
            while priortaskpointer!=None and priortaskpointer.yaml_header=='md':
                priortasklist.append(priortaskpointer)
                priortaskpointer=priortaskpointer.prior
            priortasklist=priortasklist[::-1]
            for pt in priortasklist:
                logger.debug(f'Collecting data from {pt.basename}')
                csv=pt.statevars.get('energy-csv',None)
                if csv:
                    logger.debug(f'Collecting data from {csv}')
                    try:
                        newcsv=pd.read_csv(csv,header=0,index_col=None)
                        logger.debug(f'newcsv shape: {newcsv.shape}')
                        nrg=pd.concat([nrg,newcsv],ignore_index=True)
                        logger.debug(f'nrg dataframe shape: {nrg.shape}')
                    except:
                        logger.debug(f'For some reason, I could not read this csv file')
                        logger.debug(f'{csv} does not exist or is not a valid csv file')
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
            if len(self.specs['existing-xsts'])>0:
                logger.debug(f'Extracting data from {len(self.specs["existing-xsts"])} explicitly named XST files...')
            for f in self.specs['existing-xsts']:
                if not xst:
                    xst=NAMDxst(f)
                else:
                    xst.add_file(f)

        savedata=self.specs.get('savedata',None)
        if savedata:
            logger.debug(f'Saving energy-like data to {savedata}.')
            try:
                nrg.to_csv(savedata,header=True,index=False)
            except:
                logger.debug(f'For some reason, I could not write this dataframe to csv')
                logger.debug(nrg.iloc[:3,:].to_string())
            if xst is not None:
                logger.debug(f'Saving cell data to xst-{savedata}.')
                try:
                    xst.df.to_csv(f'xst-{savedata}',header=True,index=False)
                except:
                    logger.debug(f'For some reason, I could not write this dataframe to csv')
                    logger.debug(xst.df.iloc[:3,:].to_string())
        
        traces=self.specs.get('traces',[])
        legend=self.specs.get('legend',False)
        grid=self.specs.get('grid',False)
        basename=self.specs.get('basename','myplot')
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
            # ax.set_ylabel(tracename_label+' ('+','.join([_ for _ in unitspecs if _!='*'])+')')
            if legend:
                ax.legend()
            if grid:
                ax.grid(True)
            tracename='-'.join(tracelist)
            plt.savefig(f'{basename}-{tracename}.png',bbox_inches='tight')
            plt.clf()
        self.log_message('complete')
        self.result=0
        return super().do()
    