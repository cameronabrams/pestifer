# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" Defines the PackmolInputWriter class for generating packmol input files 
    Defines the PackmolLog class for parsing packmol log files
"""

import datetime
import logging
import os

from .command import Command
from .util.progress import PackmolProgress
from .scriptwriters import Filewriter
from .stringthings import FileCollector, striplist

logger=logging.getLogger(__name__)

class PackmolInputWriter(Filewriter):
    def __init__(self,config):
        super().__init__(comment_char='#')
        self.indent=4*' '
        self.config=config
        self.progress=self.config.progress
        self.F=FileCollector()
        self.default_ext='.inp'
        self.default_script=f'packmol{self.default_ext}'
        self.scriptname=self.default_script
    
    def newscript(self,basename=None):
        timestampstr=datetime.datetime.today().ctime()
        if basename:
            self.basename=basename
        else:
            self.basename=os.path.splitext(self.default_script)[0]
        self.scriptname=f'{self.basename}{self.default_ext}'
        self.newfile(self.scriptname)
        self.banner(f'{__package__}: {self.basename}{self.default_ext}')
        self.banner(f'Created {timestampstr}')

    def writescript(self):
        self.writefile()

    def runscript(self,*args,**options):
        assert hasattr(self,'scriptname'),f'No scriptname set.'
        self.logname=f'{self.basename}.log'
        logger.debug(f'Log file: {self.logname}')
        cmd=Command(f'{self.config.shell_commands["packmol"]} < {self.scriptname}')
        progress_struct=None
        if self.progress:
            progress_struct=PackmolProgress()
        return cmd.run(ignore_codes=[173],logfile=self.logname,progress=progress_struct)

class PackmolLog:
    def __init__(self,filename):
        self.filename=filename
        if os.path.exists(filename):
            with open(filename,'r') as f:
                lines=f.read().split('\n')
            
            self.poundlines=[]
            for i,l in enumerate(lines):
                ll=l.lstrip()
                if ll.startswith('#'):
                    self.poundlines.append(i)
            self.banners=[]
            nbanners=0
            self.titles={}
            self.tasks={}
            for l,r in zip(self.poundlines[:-1],self.poundlines[1:]):
                if (r-l)<20:
                    candidate=striplist(lines[l+1:r])
                    if candidate!=[]:
                        self.banners.append(candidate)
                        nbanners+=1
                else:
                    dashlines=[]
                    for i,il in enumerate(lines[l+1:r]):
                        if il.lstrip().startswith('-'):
                            dashlines.append(l+1+i)
                    if (dashlines[0]-l)<10:
                        candidate=striplist(lines[l+1:dashlines[0]])
                        if candidate!=[]:
                            self.titles[nbanners]=candidate
                    for tl,tr in zip(dashlines[:-1],dashlines[1:]):
                        candidate=striplist(lines[tl+1:tr])
                        if candidate!=[]:
                            if not nbanners in self.tasks:
                                self.tasks[nbanners]=[]
                            self.tasks[nbanners].append(candidate)

            candidate=[x for x in lines[self.poundlines[-1]:] if x!='' and not x[0] in '-#']
            self.tail=''
            if candidate!=[]:
                self.tail=striplist(candidate)

            self.nbanners=nbanners