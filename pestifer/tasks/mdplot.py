# Author: Cameron F. Abrams, <cfa22@drexel.edu>
# 

import logging
import matplotlib.pyplot as plt
import pandas as pd

from ..basetask import BaseTask
from ..util.units import g_per_amu,A3_per_cm3
from ..util.namdlog import NAMDLog, NAMDxst

logger=logging.getLogger(__name__)
logging.getLogger("matplotlib").setLevel(logging.WARNING)

class MDPlotTask(BaseTask):
    """ A class for making plots of energy-like quantities from a series of one or more NAMD 
        runs """
    yaml_header='mdplot'
    def do(self):
        self.log_message('initiated')
        self.inherit_state()
        datasources=[]
        xstsources=[]
        if len(self.specs['existing-logs'])==0 and self.prior and self.prior.taskname=='md':
            logger.debug(f'Collecting all mdlog.edata tables from upstream md tasks...')
            mdtaskpointer=self.prior
            logger.debug(f'self={str(self)}; self.prior={str(self.prior)}')
            while mdtaskpointer!=None and mdtaskpointer.taskname=='md':
                logger.debug(f'appending log/xst data from prior task {str(mdtaskpointer)}')
                if hasattr(mdtaskpointer,'mdlog'):
                    datasources.append(mdtaskpointer.mdlog.edata)
                else:
                    logger.debug(f'Task {mdtaskpointer.index} has no mdlog.  Is this a bug?')
                if hasattr(mdtaskpointer,'xstlog'):
                    xstsources.append(mdtaskpointer.xstlog.df)
                mdtaskpointer=mdtaskpointer.prior
        else:
            logger.debug(f'Extracting data from {len(self.specs["existing-logs"])} explicitly named namd logs...')
            save_titles=[]
            for f in self.specs['existing-logs'][::-1]:
                if not save_titles:
                    l=NAMDLog(f)
                    l.energy()
                    save_titles=l.etitles
                else:
                    l=NAMDLog(f,inherited_etitles=save_titles)
                    l.energy()
                datasources.append(l.edata)
            if len(self.specs['existing-xsts'])>0:
                logger.debug(f'Extracting data from {len(self.specs["existing-xsts"])} explicitly named XST files...')
            for f in self.specs['existing-xsts'][::-1]:
                l=NAMDxst(f)
                xstsources.append(l.df)

        if len(datasources)==0:
            logger.warning(f'No datasources found for mdplot task.')
            return -1
        logger.debug(f'concatenating energy-like data from {len(datasources)} sequential logs')
        edata=pd.concat(datasources[::-1])
        savedata=self.specs.get('savedata',None)
        xstdata=None
        if len(xstsources)>0: # we instructed the md tasks that only NPT runs write xst files
            xstdata=pd.concat(xstsources[::-1])
        if savedata:
            logger.debug(f'Saving energy-like data to {savedata}.')
            try:
                edata.to_csv(savedata,header=True,index=False)
            except:
                logger.debug(f'For some reason, I could not write this dataframe to csv')
                logger.debug(edata.iloc[:3,:].to_string())
            if len(xstsources)>0:
                logger.debug(f'Saving cell data to xst-{savedata}.')
                try:
                    xstdata.to_csv(f'xst-{savedata}',header=True,index=False)
                except:
                    logger.debug(f'For some reason, I could not write this dataframe to csv')
                    logger.debug(xstdata.iloc[:3,:].to_string())
        
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
                    if unitspec in ['g_per_cc','g/cc','g_per_cm3','g/cm3']:
                        units=g_per_amu*A3_per_cm3
                    else:
                        logger.debug(f'Unitspec "{unitspec}" not recognized.')
                        units=1.0
                unitspecs.append(unitspec)
                if t_i.upper() in edata:
                    key=t_i.upper()
                    ax.plot(edata['TS'],edata[key]*units,label=key.title())
                elif len(xstsources)>0 and t_i in xstdata:
                    ax.plot(xstdata['step'],xstdata[t_i]*units,label=t_i)
            ax.set_xlabel('time step')
            tracename=','.join(tracelist)
            ax.set_ylabel(tracename+' ('+','.join([_ for _ in unitspecs if _!='*'])+')')
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
    