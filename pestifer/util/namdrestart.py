# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
Defines the :class:`NAMDConfig` class for parsing NAMD config files and the :func:`make_namd_restart` function for the NAMD restart task.
"""
import glob
import logging
import os
import shutil
from .logparsers import NAMDLog
from ..core.scripters import Filewriter
from ..core.stringthings import ByteCollector

logger=logging.getLogger(__name__)

class NAMDConfig:
    """
    A class for parsing NAMD configuration files. This class reads a NAMD config file, extracts commands and variable assignments,
    and allows for modification of the configuration by replacing variable assignments and commands.
    
    Parameters
    ----------
    configfile : str
        The path to the NAMD configuration file to be parsed. The file should exist and be readable.
        
    """
    def __init__(self,configfile):
        if not os.path.exists(configfile):
            raise FileNotFoundError(f'{configfile} not found')
        self.filename=configfile
        with open(configfile,'r') as f:
            self.rawlines=f.read().split('\n')
        self.lines=[]
        for n,l in enumerate(self.rawlines):
            if len(l)==0:
                linedict=dict(number=n,ltype='blank')
            elif l[0]=='#':
                if l.count('#')>1:
                    linedict=dict(number=n,ltype='banner',contents=l.replace('#',''))
                else:
                    linedict=dict(
                        number=n,
                        ltype='comment',
                        comment=l[1:].strip())
            else:
                tokens=[_.strip() for _ in l.split()]
                if tokens[0].lower()=='set':
                    linedict=dict(ltype='varassign',
                        number=n,
                        varname=tokens[1],
                        varval=tokens[2])
                else:
                    linedict=dict(
                        ltype='command',
                        number=n,
                        commandname=tokens[0],
                        commandargs=tokens[1:],
                    )
            self.lines.append(linedict)

        self.varsdefined={}
        for line in [x for x in self.lines if x.get('ltype',None)=='varassign']:
            if line['varname'] not in self.varsdefined:
                self.varsdefined[line['varname']]=line['varval']
            else:
                logger.warning(f'{line["varname"]} is set twice in {self.filename}')

        self.varscited=[]
        for command in [x for x in self.lines if x.get('ltype',None)=='command']:
            for arg in command['commandargs']:
                if arg.startswith('$'):
                    if arg[1]==r'{' and arg[-1]==r'}':
                        arg=arg[1:-1]
                    command['varname']=arg[1:]
                    varcited=arg[1:]
                    if not varcited in self.varsdefined:
                        logger.warning(f'TCL referenced variable \'{varcited}\' in {self.filename} does not have an assignment')
                    else:
                        if not varcited in self.varscited:
                            self.varscited.append(varcited)

    def write(self,filename):
        """
        Write the NAMD configuration to a file.
        
        Parameters
        ----------
        filename : str
            The name of the file to which the NAMD configuration will be written. If the file already exists, it will be overwritten.
        """
        W=Filewriter()
        W.filename=filename
        W.banner('pestifer NAMD restart')
        for line in self.lines:
            if line['ltype']=='comment':
                W.comment(line['comment'])
            elif line['ltype']=='varassign':
                W.addline(f'set {line["varname"]} {line["varval"]}')
            elif line['ltype']=='command':
                if len(line["commandargs"])>0:
                    W.addline(f'{line["commandname"]} {" ".join(line["commandargs"])}')
                else:
                    W.addline(f'{line["commandname"]} {line["commandargs"]}')
        W.banner('thanks for using pestifer!')
        W.writefile()

    def backresolve_lines(self):
        """
        Back-resolve variable assignments in the lines of the NAMD configuration.
        """
        for l in self.lines:
            if l['ltype']=='varassign':
                for k,v in self.varsdefined.items():
                    if k==l['varname']:
                        l['varval']=v

    def var_backresolve(self,oldvarname,newvarval):
        """
        Back-resolve a variable assignment in the NAMD configuration.
        """
        if oldvarname in self.varsdefined:
            self.varsdefined[oldvarname]=newvarval
            self.backresolve_lines()
        else:
            logger.warning(f'Cannot back-resolve {oldvarname}')

    def replace_command(self,commandname,args):
        """
        Replace a command in the NAMD configuration with a new command and its arguments.
        
        Parameters
        ----------
        commandname : str
            The name of the command to be replaced. The command name is case-insensitive.
        args : list
            A list of arguments to be passed to the command. The arguments can include variable references (starting with `$`).
            If an argument starts with `$`, it is treated as a variable reference and will be back-resolved to its value.
        """
        commandline=[x for x in self.lines if x.get('ltype',None)=='command' and x.get('commandname','fake').lower()==commandname.lower()][0]
        newargs=[]
        for o,n in zip(commandline['commandargs'],args):
            if o.startswith('$'):
                newargs.append(o)
                self.var_backresolve(o,n)
            else:
                newargs.append(n)
        commandline['commandargs']=newargs

def make_namd_restart(args,**kwargs):
    """
    Create a NAMD restart configuration file based on an existing NAMD log file and configuration file.
    
    Parameters
    ----------
    args : argparse.Namespace
        The command-line arguments parsed by argparse. It should contain the following attributes:

        - ``log``: The path to the NAMD log file.
        - ``config``: The path to the NAMD configuration file.
        - ``new_base``: The new base name for the output files.
        - ``run``: The number of time steps to run. If ``run`` is 0, the script will not run any new time steps.
        - ``slurm``: The path to the SLURM script file, if applicable.
    """
    log=args.log
    config=args.config
    newbasename=args.new_base
    run=args.run
    oldconfig=NAMDConfig(config)
    oldlog=NAMDLog.from_file(log,passfilter=['OUTPUT','RESTART','TCL','TIMESTEP'])
    output_filename=oldlog.metadata.get('output_filename',None)
    if not output_filename:
        logger.error(f'No output filename found in {log}')
    oldconfig.replace_command('outputname',[newbasename])
    last_timestep=oldlog.time_series_data.get('restart',[None])[-1]
    requested_timesteps=oldlog.metadata.get('running_for',None)
    if last_timestep is None:
        raise ValueError(f'No restart time step found in {log}. Please check the log file for errors.')
    if requested_timesteps is None:
        raise ValueError(f'No "TCL: Running for..." metadata found in {log}. Please check the log file for errors.')
    if oldlog.success():
        if run==0:
            logger.warning(f'Run logged in {log} was successful but you did not request any new time steps')
            return
        resstr=''
        oldconfig.replace_command('run',[f'{run}'])
    else:
        remaining_timesteps=requested_timesteps-last_timestep
        oldconfig.replace_command('run',[f'{remaining_timesteps}'])
        resstr='.restart'

    outputfile={}
    for ext,cmd in zip(['coor','vel','xsc'],['bincoordinates','binvelocities','extendedsystem']):
        outputfile[ext]=f'{output_filename}{resstr}.{ext}'
        if os.path.exists(outputfile[ext]):
            oldconfig.replace_command(cmd,[outputfile[ext]])
        else:
            raise FileNotFoundError(f'{outputfile[ext]} not found.')

    oldconfig.replace_command('firsttimestep',[f'{last_timestep}'])

    oldconfig.write(f'{newbasename}.namd')

    if len(args.slurm)>0:
        bc=ByteCollector()
        bc.ingest_file(args.slurm)
        oldscripts=glob.glob(f'%{args.slurm}%-*')
        if len(oldscripts)==0:
            n=1
        else:
            n=max([int(x.split('%-')[-1]) for x in oldscripts])+1
        shutil.copy(args.slurm,f'%{args.slurm}%-{n}')
        bc.reassign('BASENAME',newbasename) 
        bc.write_file(args.slurm)
