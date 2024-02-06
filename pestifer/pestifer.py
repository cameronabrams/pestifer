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
import collections
logger=logging.getLogger(__name__)
import importlib.metadata
__pestifer_version__ = importlib.metadata.version("pestifer")
from .stringthings import banner, banner_message
from .controller import Controller
from .config import Config
from .scriptwriters import VMD

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
    c,e=os.path.splitext(os.path.basename(args.config))
    C.write_complete_config(f'{c}-complete.yaml')
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
    return collections.OrderedDict(sorted(exdict.items()))

def run_example(args):
    number=args.number
    nstr=f'{number:02d}'
    c=Config()
    ex_path=c['Resources']['examples']
    for ex_yaml in glob.glob(ex_path+'/*.yaml'):
        base=os.path.split(ex_yaml)[1]
        if str(base).startswith(nstr):
            rebase=base[len(nstr)+1:]
            shutil.copy(ex_yaml,f'./{rebase}')
            args.config=rebase
            run(args)
            break
    else:
        print(f'No example {number} is found.')

# def script(args):
#     scriptname=args.scriptname
#     scriptbase=os.path.basename(scriptname)
#     sp,ss=os.path.splitext(scriptbase)
#     c=Config()
#     if os.path.exists(scriptname): # it's local
#         print(f'Running local script: {scriptname}')
#         vm=VMD(c)
#         vm.basename=sp
#         vm.scriptname=scriptname
#         vm.runscript()
#     else: # assume user is referring to a library script
#         script_dir=c['Resources']['tcl']['scripts']
#         script_path=os.path.join(script_dir,f'{sp}.tcl')
#         if os.path.exists(script_path):
#             print(f'Running library script: {script_path}')
#             vm=VMD(c)
#             vm.newscript(f'local-{sp}')
#             vm.usescript(sp)
#             vm.writescript()
#             vm.runscript()
#         else:
#             print(f'{scriptname} is not found.')

def wheretcl(args):
    c=Config()
    script_dir=c['Resources']['tcl']['scripts']
    proc_dir=c['Resources']['tcl']['proc']
    if args.startup_script_path:
        print(c.vmd_startup_script)
    if args.script_dir:
        print(script_dir)
    if args.proc_dir:
        print(proc_dir)
    if args.verbose:
        print(f'Copy-paste one or more of these commands into a VMD session console:')
        for s_tcl in glob.glob(script_dir+'/*.tcl'):
            print(f'source {s_tcl}')
        for p_tcl in glob.glob(proc_dir+'/*.tcl'):
            print(f'source {p_tcl}')

def cli():
    commands={
        'config-help':config_help,
        'config-default':config_default,
        'run-example': run_example,
        'run':run,
        # 'script':script,
        'wheretcl':wheretcl
    }
    helps={
        'config-help':'get help on the syntax of input configuration files',
        'config-default':'generate a default input directive',
        'run-example':'build one of the examples provided',
        'run':'build a system using instructions in the config file',
        # 'script':'run a single-use special VMD/TcL script',
        'wheretcl':'provides path of TcL scripts for sourcing in interactive VMD'
    }
    descs={
        'config-help':'Use this command to get interactive help on config file directives.',
        'config-default':'This will generate a default config file for you to fill in using a text editor.',
        'run-example':'Build one of the examples:\n'+'\n'.join([f'{c:>3d}: {d}' for c,d in list_examples().items()]),
        'run':'Build a system',
        # 'script':'Run a single-use special VMD/TcL script',
        'wheretcl':'provides path of TcL scripts for sourcing in interactive VMD'
    }
    parser=ap.ArgumentParser(description=textwrap.dedent(banner_message),formatter_class=ap.RawDescriptionHelpFormatter)
    subparsers=parser.add_subparsers()
    subparsers.required=True
    command_parsers={}
    for k in commands:
        command_parsers[k]=subparsers.add_parser(k,description=descs[k],help=helps[k],formatter_class=ap.RawDescriptionHelpFormatter)
        command_parsers[k].set_defaults(func=commands[k])
        command_parsers[k].add_argument('--no-banner',default=False,action='store_true',help='turn off the banner')
        command_parsers[k].add_argument('--loglevel',type=str,default='debug',help='Log level for messages written to diagnostic log (debug|info)')
        command_parsers[k].add_argument('--diag',type=str,default='pestifer_diagnostics.log',help='diagnostic log file')
    
    command_parsers['run'].add_argument('config',type=str,default=None,help='input configuration file in YAML format')
    command_parsers['run-example'].add_argument('number',type=int,default=None,help='example number')
    command_parsers['config-help'].add_argument('directives',type=str,nargs='*',help='config file directives')
    command_parsers['config-default'].add_argument('directives',type=str,nargs='*',help='config file directives')
    # command_parsers['script'].add_argument('scriptname',type=str,default=None,help='name of VMD/TcL script to run')
    command_parsers['wheretcl'].add_argument('--startup-script-path',default=False,action='store_true',help='print full path of VMD startup script used by pestifer')
    command_parsers['wheretcl'].add_argument('--script-dir',default=False,action='store_true',help='print full path of directory of VMD scripts')
    command_parsers['wheretcl'].add_argument('--proc-dir',default=False,action='store_true',help='print full path of directory of VMD/TcL procedure library')
    command_parsers['wheretcl'].add_argument('--verbose',default=False,action='store_true',help='print a long message describing all TcL files available')
    args=parser.parse_args()
    args.func(args)