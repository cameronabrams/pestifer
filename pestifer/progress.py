# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import progressbar
import logging
logger=logging.getLogger(__name__)
class PestiferProgress:
    def __init__(self,init_f,meas_f,comp_f):
        """ Initialize an instance of PestiferProgress 
        
        Parameters:
        -----------
        init_f: function that analyzes its single string argument and 
        returns either False or an object that describes the iteration
        
        meas_f: function that analyzes its single string argument and
        returns an object that measures the current state
        
        comp_f: function that analyzes the object returned by init_f
        and meas_f to report a number between 0 and 1 measuring
        the progress of whatever process is generating the string
        """
        self.init_f=init_f
        self.meas_f=meas_f
        self.comp_f=comp_f
        self.initialized=False

    def go(self,a_string):
        # if not initialized, attempt to initialize
        if not self.initialized:
            init_obj=self.init_f(a_string)
            logger.debug(f'Progress: init_obj {init_obj}')
            if init_obj:
                self.initialized=True
                self.bar=progressbar.ProgressBar(max_value=50,widgets=[
                    ' [', progressbar.Timer(), '] ',
                          progressbar.Bar(),
                    ' (', progressbar.ETA(), ') ',
                        ])
                # self.bar=progressbar.ProgressBar(max_value=50)
                self.init_obj=init_obj
        else:
            meas_obj=self.meas_f(a_string)
            logger.debug(f'Progress: meas_obj {meas_obj}')
            if meas_obj:
                progress=self.comp_f(self.init_obj,meas_obj)
                logger.debug(f'progress {progress}')
                self.bar.update(int(progress*50))
