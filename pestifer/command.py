# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#
# Class for handling running external commands
#
import atexit
import logging
import os
import signal
import shutil
import subprocess

from glob import glob

logger=logging.getLogger(__name__)

class Command:
    divider_line_length=55
    def __init__(self,command:str,*args,**options):
        self.command=command
        self.args=args
        self.options=options
        self.c=f'{self.command} '+' '.join(args)+' '.join([f'-{k} {v}' for k,v in self.options.items()])
        self.stdout=''
        self.stderr=''

    def run(self,logfile=None,override=(),ignore_codes=[],quiet=True,logparser=None,**kwargs):
        """ runs this Command instance
        
        Parameters:
        -----------

        logfile : str
            name of log file to write process' stdout and stderr. If None (default), stdout and stderr are retained
            only in this instance's stdout and stderr attributes.
        override : tuple
            tuple composed of a "needle" and a "message".  If the needle is found in the stdout or stderr of
            the process, the message is displayed and an error is thrown, halting the program.
        ignore_codes : list
            a list of integer exit codes that are ignored in addition to 0
        quiet : boolean
        logparser: LogParser instance (or descendant) 
            used for progress bar/elapsed time displays and log parsing

        """
        if not quiet:
            logger.debug(f'{self.c}')
        log=None
        if not logfile:
            logger.debug(f'No logfile specified for {self.c}')
        else:
            if os.path.exists(logfile) and not kwargs.get('overwrite_logs',False):
                nlogs=len(glob(f'%{logfile}'))
                shutil.move(logfile,f'%{logfile}-{nlogs+1}%')
                logger.debug(f'Rotating {logfile} to %{logfile}-{nlogs+1}%')
            log=open(logfile,'w')
            logger.debug(f'Opened {logfile} for writing')

        process=subprocess.Popen(self.c,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        global pid
        pid=process.pid
        def kill_child():
            if pid is None:
                pass
            else:
                try:
                    os.kill(pid, signal.SIGTERM)
                except:
                    pass
        atexit.register(kill_child)
        self.stdout=''
        self.stderr=''
        output=''

        while True:
            output=process.stdout.readline()
            self.stdout+=output
            if logparser:
                logparser.update(output)
                logparser.update_progress_bar()
            if log: 
                log.write(output)
                log.flush()
            if output=='' and process.poll() is not None:
                break
        if hasattr(logparser,'progress_bar') and logparser.progress_bar is not None:
            print()
        if logfile:
            logger.debug(f'Log written to {logfile}')
        remaining_stdout,self.stderr=process.communicate()
        if logparser:
            logparser.update(remaining_stdout)
            if hasattr(logparser,'finalize'):
                logparser.finalize()
        if process.returncode!=0 and not process.returncode in ignore_codes:
            logger.error(f'Returncode: {process.returncode}')
            if len(self.stdout)>0:
                logger.error('stdout buffer follows\n'+'*'*self.divider_line_length+'\n'+self.stdout+'\n'+'*'*self.divider_line_length)
            if len(self.stderr)>0:
                logger.error('stderr buffer follows\n'+'*'*self.divider_line_length+'\n'+self.stderr+'\n'+'*'*self.divider_line_length)
            return process.returncode
        if len(override)==2:
            needle,msg=override
            if needle in self.stdout or needle in self.stderr:
                logger.info(f'Returncode: {process.returncode}, but another error was detected:')
                logger.error(msg)
                if len(self.stdout)>0 and needle in self.stdout:
                    logger.error('stdout buffer follows\n'+'*'*self.divider_line_length+'\n'+self.stdout+'\n'+'*'*self.divider_line_length)
                if len(self.stderr)>0 and needle in self.stderr:
                    logger.error('stderr buffer follows\n'+'*'*self.divider_line_length+'\n'+self.stderr+'\n'+'*'*self.divider_line_length)
        return 0