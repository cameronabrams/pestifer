# Author: Cameron F. Abrams <cfa22@drexel.edu>.
"""Automatic generation of PSF files

Pestifer is a python package that facilitates the use
of UIUC's psfgen VMD plug-in for generating PSF files
for use in MD simulations of biomacromolecules using NAMD.

"""
import argparse as ap
import textwrap
import os
import shutil
import yaml
import glob
import logging
logger=logging.getLogger(__name__)
import importlib.metadata
__pestifer_version__ = importlib.metadata.version("pestifer")
from .stringthings import banner, banner_message
from .controller import Controller
from .config import Config

def config_help(args):
    c=Config()
    if not args.no_banner:
        banner(print)
    print(f'Help on user-provided configuration file format')
    directives=args.directives
    c.console_help(*directives)

def config_default(args):
    c=Config()
    if not args.no_banner:
        banner(print)
    print(f'Default directive:')
    directives=args.directives
    specs=c.make_default_specs(*directives)
    print(yaml.dump(specs))

def run(args):
    # Set up logging to both a log file and the console
    loglevel=args.loglevel
    loglevel_numeric=getattr(logging, loglevel.upper())
    if os.path.exists(args.diag):
        shutil.copyfile(args.diag,args.diag+'.bak')
    logging.basicConfig(filename=args.diag,filemode='w',format='%(asctime)s %(name)s %(message)s',level=loglevel_numeric)
    console=logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter=logging.Formatter('%(levelname)s> %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

    if not args.no_banner:
        banner(logger.info)
    # Set up the Controller and execute tasks
    logger.info(f'pestifer runtime begins')
    C=Controller(args.config)
    C.write_complete_config('00-complete.yaml')
    C.do_tasks()    
    logger.info('pestifer runtime ends.')

def list_examples():
    c=Config()
    ex_path=c['Resources']['examples']
    ex_full_paths=glob.glob(ex_path+'/*.yaml')
    exdict={}
    for e in ex_full_paths:
        base=os.path.split(e)[1]
        code=int(base.split('-')[0])
        with open(e,'r') as f:
            d=yaml.safe_load(f)
            desc=d['title']
            exdict[code]=desc
    return exdict

def run_example(args):
    number=args.number
    nstr=f'{number:02d}'
    c=Config()
    ex_path=c['Resources']['examples']
    for ex_yaml in glob.glob(ex_path+'/*.yaml'):
        base=os.path.split(ex_yaml)[1]
        if str(base).startswith(nstr):
            shutil.copy(ex_yaml,'./')
            args.config=base
            run(args)
            break
    else:
        print(f'No example {number} is found.')

def cli():
    commands={
        'config-help':config_help,
        'config-default':config_default,
        'run-example': run_example,
        'run':run
    }
    helps={
        'config-help':'get help on the syntax of input configuration files',
        'config-default':'generate a default input directive',
        'run-example':'build one of the examples provided: '+'; '.join([f'{c}: {d}' for c,d in list_examples().items()]),
        'run':'build a system using instructions in the config file'
    }
    parser=ap.ArgumentParser(description=textwrap.dedent(banner_message),formatter_class=ap.RawDescriptionHelpFormatter)
    subparsers=parser.add_subparsers()
    subparsers.required=True
    command_parsers={}
    for k in commands:
        command_parsers[k]=subparsers.add_parser(k,help=helps[k])
        command_parsers[k].set_defaults(func=commands[k])
        command_parsers[k].add_argument('--no-banner',default=False,action='store_true',help='turn off the banner')
        command_parsers[k].add_argument('--loglevel',type=str,default='debug',help='Log level for messages written to diagnostic log (debug|info)')
        command_parsers[k].add_argument('--diag',type=str,default='pestifer_diagnostics.log',help='diagnostic log file')
    
    command_parsers['run'].add_argument('config',type=str,default=None,help='input configuration file in YAML format')
    command_parsers['run-example'].add_argument('number',type=int,default=None,help='example number')
    command_parsers['config-help'].add_argument('directives',type=str,nargs='*',help='config file directives')
    command_parsers['config-default'].add_argument('directives',type=str,nargs='*',help='config file directives')
    args=parser.parse_args()
    args.func(args)