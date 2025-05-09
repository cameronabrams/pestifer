# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" Defines the NAMDConfig class for parsing NAMD config files and the make_namd_restart function
"""
import glob
import logging
import os
import shutil
from .util.logparsers import NAMDLog
from .scriptwriters import Filewriter
from .stringthings import ByteCollector

logger=logging.getLogger(__name__)

class NAMDConfig:
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
        for l in self.lines:
            if l['ltype']=='varassign':
                for k,v in self.varsdefined.items():
                    if k==l['varname']:
                        l['varval']=v

    def var_backresolve(self,oldvarname,newvarval):
        if oldvarname in self.varsdefined:
            self.varsdefined[oldvarname]=newvarval
            self.backresolve_lines()
        else:
            logger.warning(f'Cannot back-resolve {oldvarname}')

    def replace_command(self,commandname,args):
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
    log=args.log
    config=args.config
    newbasename=args.new_base
    run=args.run
    oldconfig=NAMDConfig(config)
    oldlog=NAMDLog(log)
    output_filename=oldlog.info.get('OUTPUT FILENAME',None)
    if not output_filename:
        logger.error(f'No output filename found in {log}')
    oldconfig.replace_command('outputname',[newbasename])
    if oldlog.success():
        if run==0:
            logger.warning(f'Run logged in {log} was successful but you did not request any new time steps')
            return
        last_timestep=oldlog.output_timestep
        resstr=''
        oldconfig.replace_command('run',[f'{run}'])
    else:
        last_timestep=oldlog.restart_timestep
        remaining_timesteps=oldlog.requested_timesteps-last_timestep
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
        bc.injest_file(args.slurm)
        oldscripts=glob.glob(f'%{args.slurm}%-*')
        if len(oldscripts)==0:
            n=1
        else:
            n=max([int(x.split('%-')[-1]) for x in oldscripts])+1
        shutil.copy(args.slurm,f'%{args.slurm}%-{n}')
        bc.reassign('BASENAME',newbasename) 
        bc.write_file(args.slurm)
