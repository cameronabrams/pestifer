# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" Implements progress bar animations for pestifer runs outside of a batch system 
    (for fun)
"""
import progressbar
import logging
from .colors import __plasma__
from .namdlog import getinfo, gettclinfo

logger=logging.getLogger(__name__)
class PestiferProgress:
    def __init__(self,**kwargs): #Elapsed Time: %(elapsed)s\x1b[33mColorful example\x1b[39m
        """ Initialize an instance of PestiferProgress 
        """
        self.name=kwargs.get('name','Elapsed')
        self.color=__plasma__[kwargs.get('colorno',100)]
        self.max_value=kwargs.get('max_value',progressbar.UnknownLength)
        self.timer_format=kwargs.get('timer_format','time: %(elapsed)s')
        if self.max_value==progressbar.UnknownLength:
            self.unmeasured=True
        else:
            self.unmeasured=False
        self.widgets=kwargs.get('widgets',None)
        self.track_stdout=kwargs.get('track_stdout',True)
        self.interval_sec=kwargs.get('interval_sec',1)
        timer_format=f'{self.color}{self.name}{self.color.OFF} {self.timer_format}'
        if not self.widgets:
            if self.max_value==progressbar.UnknownLength:
                self.widgets=[
                    progressbar.Timer(timer_format),' ',progressbar.RotatingMarker()
                ]
            else:
                self.widgets=[
                    progressbar.Timer(timer_format),
                    progressbar.Bar(),' ', progressbar.ETA()
                ]
        self.initialized=False
        self.init_obj=False
        self.meas_obj=False

    def init_f(self,a_string=''):
        self.init_obj=True
        return self.init_obj
    def meas_f(self,a_string=''):
        self.meas_obj=True
        return self.meas_obj
    def comp_f(self):
        return True

    def go(self,a_string=''):
        # if not initialized, attempt to initialize
        if not self.initialized:
            if self.init_f(a_string):
                self.initialized=True
                self.bar=progressbar.ProgressBar(max_value=self.max_value,widgets=self.widgets)
                # print('initialized')
        else:
            if self.meas_f(a_string):
                if not self.unmeasured:
                    progress=self.comp_f()
                    self.bar.update(int(progress*self.max_value))
                else:
                    self.bar.update()

class NAMDProgress(PestiferProgress):
    groupnames=['ETITLE:','ENERGY:','Info:','TCL:']
    infonames=['FIRST TIMESTEP']
    tclnames=['Running for','Minimizing for']
    def __init__(self,**kwargs):
        super().__init__(max_value=200,name='namd',colorno=50,**kwargs)
        self.groups={}
        self.info={}
        self.tcl={}
    def parse_f(self,logstring):
        self.groups={}
        self.info={}
        self.tcl={}
        self.etitles=[]
        loglines=logstring.split('\n')
        for l in loglines:
            if len(l)>0:
                tok=l.split()[0]
                if tok in self.groupnames:
                    if not tok in self.groups:
                        self.groups[tok]=[]
                    self.groups[tok].append(l)
        for etitle in self.groups.get('ETITLE:',[]):
            self.etitles.append(etitle)
        for info in self.groups.get('Info:',[]):
            for il in self.infonames:
                if il in info:
                    self.info[il]=getinfo(il,info)
        for tclline in self.groups.get('TCL:',[]):
            for tl in self.tclnames:
                if tl in tclline:
                    self.tcl[tl]=gettclinfo(tl,tclline)
    def init_f(self,logstring):
        self.init_obj=False
        if not logstring:
            return self.init_obj
        self.parse_f(logstring)
        # if self.info: logger.debug(f'info {self.info}')
        # if self.tcl: logger.debug(f'tcl {self.tcl}')
        first_step=self.info.get('FIRST TIMESTEP',0)
        num_steps=self.tcl.get('Running for',None)
        nummin_steps=self.tcl.get('Minimizing for',None)
        estart=len(self.etitles)>0
        if nummin_steps:
            self.init_obj={'first_step':int(first_step) if first_step else 0}
            self.init_obj['num_steps']=nummin_steps
            return self.init_obj
        elif num_steps:
            if first_step or estart:
                self.init_obj={'first_step':int(first_step)}
                self.init_obj['num_steps']=num_steps
                return self.init_obj
        return self.init_obj
    def meas_f(self,logstring):
        self.meas_obj=False
        if not logstring:
            return self.meas_obj
        self.parse_f(logstring)
        if 'ENERGY:' in self.groups and len(self.groups['ENERGY:'])>0:
            last_line=self.groups['ENERGY:'][-1]
            if len(last_line)<20:
                return False
            tok=last_line.split()
            if len(tok)>1:
                current_time_step=int(tok[1])
                # print(current_time_step)
                self.meas_obj=dict(current_time_step=current_time_step)
                return self.meas_obj
        return self.meas_obj
    def comp_f(self):
        if not self.init_obj or not self.meas_obj:
            return 0
        fac=(self.meas_obj['current_time_step']-self.init_obj['first_step'])/self.init_obj['num_steps']
        assert fac<=1.0,f'error: {self.meas_obj} {self.init_obj}'
        return fac
    
class PackmolProgress(PestiferProgress):
    def __init__(self,**kwargs):
        super().__init__(name='packmol',colorno=150,**kwargs)

class PsfgenProgress(PestiferProgress):
    def __init__(self,**kwargs):
        super().__init__(name='psfgen',colorno=125,**kwargs)

class RingCheckProgress(PestiferProgress):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)


