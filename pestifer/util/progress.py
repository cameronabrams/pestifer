# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" Implements progress bar animations for pestifer runs outside of a batch system 
    (for fun)
"""
import progressbar
import logging
from .colors import __plasma__

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
        self.bar=progressbar.ProgressBar(max_value=self.max_value,widgets=self.widgets)

    def register_update_function(self,func):
        self.comp_f=func
        self.initialized=True

    def go(self):
        if self.initialized:
            if self.unmeasured:
                self.bar.update()
            else:
                progress=self.comp_f()
                self.bar.update(int(progress*self.max_value))
        else:
            self.bar.update()

class NAMDProgress(PestiferProgress):
    def __init__(self,**kwargs):
        super().__init__(max_value=200,name='namd',colorno=50,**kwargs)
    
class PackmolProgress(PestiferProgress):
    def __init__(self,**kwargs):
        super().__init__(max_value=200,name='packmol',colorno=150,**kwargs)

class PsfgenProgress(PestiferProgress):
    def __init__(self,**kwargs):
        super().__init__(name='psfgen',colorno=125,**kwargs)

class RingCheckProgress(PestiferProgress):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)


