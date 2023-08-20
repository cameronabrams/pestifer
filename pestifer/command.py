import logging
import subprocess
logger=logging.getLogger(__name__)

class Command:
    divider_line_length=55
    def __init__(self,command,*args,**options):
        self.command=command
        self.args=args
        self.options=options
        self.c=f'{self.command} '+' '.join(args)+' '.join([f'-{k} {v}' for k,v in self.options.items()])
        self.stdout=''
        self.stderr=''

    def run(self,override=(),ignore_codes=[],quiet=True):
        if not quiet:
            logger.debug(f'{self.c}')
        process=subprocess.Popen(self.c,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        self.stdout,self.stderr=process.communicate()
        if process.returncode!=0 and not process.returncode in ignore_codes:
            logger.error(f'Returncode: {process.returncode}')
            if len(self.stdout)>0:
                logger.error('stdout buffer follows\n'+'*'*self.divider_line_length+'\n'+self.stdout+'\n'+'*'*self.divider_line_length)
            if len(self.stderr)>0:
                logger.error('stderr buffer follows\n'+'*'*self.divider_line_length+'\n'+self.stderr+'\n'+'*'*self.divider_line_length)
            raise subprocess.SubprocessError(f'Command "{self.c}" failed with returncode {process.returncode}')
        else:
            # logger.info(f'Returncode: {process.returncode}.')
            if len(override)==2:
                needle,msg=override
                if needle in self.stdout or needle in self.stderr:
                    logger.info(f'Returncode: {process.returncode}, but another error was detected:')
                    logger.error(msg)
                    if len(self.stdout)>0:
                        logger.error('stdout buffer follows\n'+'*'*self.divider_line_length+'\n'+self.stdout+'\n'+'*'*self.divider_line_length)
                    if len(self.stderr)>0:
                        logger.error('stderr buffer follows\n'+'*'*self.divider_line_length+'\n'+self.stderr+'\n'+'*'*self.divider_line_length)