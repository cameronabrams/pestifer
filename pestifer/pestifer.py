"""

.. module:: pestifer
   :synopsis: manages command-line interface; sets up logging
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import argparse as ap
import os
import shutil
import logging
logger=logging.getLogger(__name__)

from .stringthings import banner
from .controller import Controller

def _main():
    parser=ap.ArgumentParser()
    parser.add_argument('config',help='input configuration file (yaml)')
    parser.add_argument('-log',help='log file name',default='pestifer.log')
    parser.add_argument('--loglevel',default='info',help='logging level (info)')
    args=parser.parse_args()

    # Set up logging to both a log file and the console
    loglevel=args.loglevel
    loglevel_numeric=getattr(logging, loglevel.upper())
    if os.path.exists(args.log):
        shutil.copyfile(args.log,args.log+'.bak')
    logging.basicConfig(filename=args.log,filemode='w',format='%(asctime)s %(name)s %(message)s',level=loglevel_numeric)
    console=logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter=logging.Formatter('%(levelname)s> %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

    banner(logger.info)
    # Set up the Controller and execute tasks
    logger.info(f'pestifer runtime begins')
    C=Controller(args.config)
    C.do_tasks()    
    logger.info('pestifer runtime ends.')

def cli():
    _main()
