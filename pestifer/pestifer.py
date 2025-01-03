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
import time
import datetime
import logging
import collections
import importlib.metadata
__pestifer_version__ = importlib.metadata.version("pestifer")
from .stringthings import banner, banner_message, enhanced_banner_message, oxford
from .controller import Controller
from .config import Config, ResourceManager
from .scriptwriters import Psfgen
from .charmmresi import make_RESI_database

logger=logging.getLogger(__name__)

logging.getLogger("pidibble").setLevel(logging.WARNING)
logging.getLogger("ycleptic").setLevel(logging.WARNING)
logging.getLogger("matplotlib").setLevel(logging.WARNING)

def config_help(args):
    c=Config()
    print(f'Help on user-provided configuration file format')
    directives=args.directives
    iprompt='pestifer-help: ' if args.interactive else ''
    c.console_help(directives,interactive_prompt=iprompt,exit=True)

def config_default(args):
    c=Config()
    print(f'Default directive:')
    directives=args.directives
    specs=c.make_default_specs(*directives)
    print(yaml.dump(specs))

def run(args,**kwargs):
    if args.output_dir!='./':
        exec_dir=os.getcwd()
        if not os.path.exists(args.output_dir):
            os.mkdir(args.output_dir)
        shutil.copy(args.config,args.output_dir)
        os.chdir(args.output_dir)
    loglevel_numeric=getattr(logging,args.diagnostic_log_level.upper())
    if args.diagnostic_log_file:
        if os.path.exists(args.diagnostic_log_file):
            shutil.copyfile(args.diagnostic_log_file,args.diagnostic_log_file+'.bak')
        logging.basicConfig(filename=args.diagnostic_log_file,filemode='w',format='%(asctime)s %(name)s %(message)s',level=loglevel_numeric)
    console=logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter=logging.Formatter('%(levelname)s> %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

    # Set up the Controller and execute tasks
    begin_time=time.time()
    logger.info(f'{__package__} runtime begins.')
    allowed_extensions=['.yaml','.yml','.y']
    configname=args.config
    cbase,cext=os.path.splitext(configname)
    if not cext:
        fil=[os.path.exists(f'{cbase}{ext}') for ext in allowed_extensions]
        if any(fil):
            iix=fil.index(True)
            configname=f'{cbase}{allowed_extensions[iix]}'
    C=Controller(configname)
    C.write_complete_config(f'{cbase}-complete.yaml')
    report=C.do_tasks()
    end_time=time.time()
    elapsed_time_s=datetime.timedelta(seconds=(end_time-begin_time))
    logger.info(f'{__package__} runtime ends. Elapsed runtime {time.strftime("%H:%M:%S",time.gmtime(elapsed_time_s.seconds))}.')
    if args.output_dir!='./':
        os.chdir(exec_dir)

def desolvate(args,**kwargs):
    C=Controller(userspecs={'tasks':
                            [
                                {'desolvate':{
                                    'basename':'desolvate',
                                    'keepatselstr':args.keepatselstr,
                                    'psf':args.psf,
                                    'dcd_infiles':args.dcd_infiles,
                                    'dcd_outfile':args.dcd_outfile,
                                    'psf_outfile':args.psf_outfile,
                                    'idx_outfile':args.idx_outfile,
                                    'dcd_stride':args.dcd_stride
                                   }
                                }
                            ]})
    report=C.do_tasks()

def list_examples():
    r=ResourceManager()
    ex_path=r['examples']
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

def fetch_example(args,**kwargs):
    number=args.number
    nstr=f'{number:02d}'
    r=ResourceManager()
    ex_path=r['examples']
    ex_yaml=glob.glob(ex_path+f'/{nstr}*.yaml')
    assert len(ex_yaml)==1,f'No example {number} found'
    ex_yaml=ex_yaml[0]
    with open(ex_yaml,'r') as f:
        ex_dict=yaml.safe_load(f)
    for update in args.config_updates:
        with open(update,'r') as f:
            up_dict=yaml.safe_load(f)
            print(f'ex_dict {ex_dict}')
            print(f'up_dict {up_dict}')
            ex_dict.update(up_dict)
            print(f'ex_dict {ex_dict}')
    
    base=os.path.split(ex_yaml)[1]
    rebase=base[len(nstr)+1:]
    
    with open(rebase,'w') as f:
        yaml.dump(ex_dict,f)
    # shutil.copy(ex_yaml,f'./{rebase}')
    return rebase

def run_example(args,**kwargs):
    args.config=fetch_example(args,**kwargs)
    run(args,**kwargs)

def inittcl(args):
    c=Config()
    assert 'tcl' in c['Resources'],f'No directory \'tcl\' found in resources -- this is a bad {__package__} installation'
    tcl_root=os.path.join(c['Resources']['root'],'tcl')
    macrofile=os.path.join(tcl_root,'macros.tcl')
    if not os.path.exists(macrofile) or args.force:
        logger.info(f'Generating macros.tcl in package directory {tcl_root}')
        p=Psfgen(c)
        p.newfile(macrofile)
        p.addline('# This file generated by pestifer initttcl')
        p.atomselect_macros()
        p.writefile()
        
def wheretcl(args):
    r=ResourceManager()
    assert 'tcl' in r,f'No directory \'tcl\' found in resources -- this is a bad {__package__} installation'
    tcl_root=os.path.join(r['root'],'tcl')
    script_dir=r['tcl']['scripts']
    pkg_dir=r['tcl']['pkg']
    msg=''
    if args.root or args.verbose:
        assert os.path.exists(tcl_root)
        if args.verbose:
            msg='Tcl root: '
        print(f'{msg}{tcl_root}')
    if args.startup_script_path or args.verbose:
        vmd_startup_script=os.path.join(tcl_root,'vmdrc.tcl')
        assert os.path.exists(vmd_startup_script)
        if args.verbose:
            msg='VMD Startup script: '
        print(f'{msg}{vmd_startup_script}')
    if args.script_dir or args.verbose:
        assert os.path.exists(script_dir)
        if args.verbose:
            msg='VMD Script directory: '
        print(f'{msg}{script_dir}')
    if args.pkg_dir or args.verbose:
        assert os.path.exists(pkg_dir)
        if args.verbose:
            msg='TcL package directory: '
        print(f'{msg}{pkg_dir}')

def cli():
    commands={
        'config-help':config_help,
        'config-default':config_default,
        'fetch-example': fetch_example,
        'run-example': run_example,
        'run':run,
        'wheretcl':wheretcl,
        'inittcl':inittcl,
        'make-resi-database':make_RESI_database,
        'desolvate':desolvate
    }
    helps={
        'config-help':'get help on the syntax of input configuration files',
        'config-default':'generate a default input directive',
        'fetch-example':'copy the example\'s YAML config file to the CWD',
        'run-example':'build one of the examples provided',
        'run':'build a system using instructions in the config file',
        'wheretcl':'provides path of TcL scripts for sourcing in interactive VMD',
        'inittcl':'initializes macros from config',
        'make-resi-database':'makes representative psf/pdb files for any CHARMM RESI\'s found in given topology streams',
        'desolvate':'desolvate an existing PSF/DCD'
    }
    descs={
        'config-help':'Use this command to get interactive help on config file directives.',
        'config-default':'This will generate a default config file for you to fill in using a text editor.',
        'fetch-example':'Fetch YAML config for one of the examples:\n'+'\n'.join([f'{c:>3d}: {d}' for c,d in list_examples().items()]),
        'run-example':'Build one of the examples:\n'+'\n'.join([f'{c:>3d}: {d}' for c,d in list_examples().items()]),
        'run':'Build a system',
        'wheretcl':'provides path of TcL scripts for sourcing in interactive VMD',
        'inittcl':'initializes macros from config',
        'make-resi-database':'makes representative psf/pdb files for any CHARMM RESI\'s found in given topology streams',
        'desolvate':'desolvate an existing PSF/DCD'
    }
    parser=ap.ArgumentParser(description=textwrap.dedent(banner_message),formatter_class=ap.RawDescriptionHelpFormatter)
    parser.add_argument('--no-banner',default=False,action='store_true',help='turn off the banner')
    parser.add_argument('--kick-ass-banner',default=False,action='store_true',help=ap.SUPPRESS)
    parser.add_argument('--diagnostic-log-level',type=str,default=None,choices=[None,'info','debug','warning'],help='Log level for messages written to diagnostic log')
    parser.add_argument('--diagnostic-log-file',type=str,default='pestifer_diagnostics.log',help='diagnostic log file')
    subparsers=parser.add_subparsers()
    subparsers.required=False
    command_parsers={}
    for k in commands:
        command_parsers[k]=subparsers.add_parser(k,description=descs[k],help=helps[k],formatter_class=ap.RawDescriptionHelpFormatter)
        command_parsers[k].set_defaults(func=commands[k])
    
    command_parsers['run'].add_argument('config',type=str,default=None,help='input configuration file in YAML format')
    command_parsers['run'].add_argument('--output-dir',type=str,default='./',help='name of output directory relative to CWD (default: %(default)s)')
    command_parsers['run'].add_argument('--diagnostic-log-level',type=str,default='debug',choices=[None,'info','debug','warning'],help='Log level for messages written to diagnostic log')
    command_parsers['fetch-example'].add_argument('number',type=int,default=None,help='example number')
    command_parsers['fetch-example'].add_argument('--config-updates',type=str,nargs='+',default=[],help='yaml files to update example')
    command_parsers['run-example'].add_argument('number',type=int,default=None,help='example number')
    command_parsers['run-example'].add_argument('--output-dir',type=str,default='./',help='name of output directory relative to CWD')
    command_parsers['run-example'].add_argument('--config-updates',type=str,nargs='+',default=[],help='yaml files to update example')
    command_parsers['run-example'].add_argument('--diagnostic-log-level',type=str,default='debug',choices=[None,'info','debug','warning'],help='Log level for messages written to diagnostic log')
    command_parsers['config-help'].add_argument('directives',type=str,nargs='*',help='config file directives')
    command_parsers['config-help'].add_argument('--interactive',default=True,action=ap.BooleanOptionalAction,help='use help in interactive mode')
    command_parsers['config-default'].add_argument('directives',type=str,nargs='*',help='config file directives')
    command_parsers['wheretcl'].add_argument('--startup-script-path',default=False,action='store_true',help='print full path of VMD startup script used by pestifer')
    command_parsers['wheretcl'].add_argument('--verbose',default=False,action='store_true',help='print full paths of all TcL resources of Pestifer')
    command_parsers['wheretcl'].add_argument('--script-dir',default=False,action='store_true',help='print full path of directory of Pestifer\'s VMD scripts')
    command_parsers['wheretcl'].add_argument('--pkg-dir',default=False,action='store_true',help='print full path of directory of Pestifer\'s VMD/TcL package library')
    command_parsers['wheretcl'].add_argument('--root',default=False,action='store_true',help='print full path of Pestifer\'s root TcL directory')
    command_parsers['inittcl'].add_argument('--force',default=False,action='store_true',help='force overwrite of any package-resident tcl files inittcl generates')
    command_parsers['make-resi-database'].add_argument('--streams',type=str,nargs='+',default=['lipid'],help='list of charmmff streams to scan (default: %(default)s)')
    command_parsers['make-resi-database'].add_argument('--force',default=False,action='store_true',help='force overwrite of any existing molecules in the database')
    command_parsers['make-resi-database'].add_argument('--cleanup',default=True,action=ap.BooleanOptionalAction,help='clean up all working files (default: %(default)s)')
    command_parsers['make-resi-database'].add_argument('--resi',nargs='+',type=str,default=[],help='regenerate entry for just these RESIs')
    command_parsers['make-resi-database'].add_argument('--diagnostic-log-level',type=str,default='debug',choices=[None,'info','debug','warning'],help='Log level for messages written to diagnostic log (default: %(default)s)')
    command_parsers['make-resi-database'].add_argument('--lenfac',type=float,default=1.4,help='this factor times topological distance is the cartesian distance to which you want to stretch a molecule (default: %(default)s)')
    command_parsers['make-resi-database'].add_argument('--minimize-steps',type=int,default=500,help='number of minimization steps immediately after each build (default: %(default)s)')
    command_parsers['make-resi-database'].add_argument('--sample-steps',type=int,default=5000,help='number of sample steps (default: %(default)s)')
    command_parsers['make-resi-database'].add_argument('--nsamples',type=int,default=10,help='number of samples (default: %(default)s)')
    command_parsers['make-resi-database'].add_argument('--sample-temperature',type=float,default=300.0,help='number of sample steps (default: %(default)s)')
    command_parsers['make-resi-database'].add_argument('--output-dir',type=str,default='data',help='name of output directory relative to CWD (default: %(default)s)')
    command_parsers['make-resi-database'].add_argument('--fail-dir',type=str,default='fails',help='name of output directory for failed runs relative to CWD (default: %(default)s)')
    command_parsers['make-resi-database'].add_argument('--refic-idx',type=int,default=0,help='index of reference IC to use to build a single molecule (default: %(default)s)')
    command_parsers['desolvate'].add_argument('--psf',type=str,help='name of input PSF file')
    command_parsers['desolvate'].add_argument('--dcd-infiles',type=str,nargs='+',default=[],help='list of input dcd files in chronological order')
    command_parsers['desolvate'].add_argument('--keepatselstr',type=str,default='protein or glycan or lipid',help='VMD atomsel string for atoms you want to keep (default: "%(default)s")')
    command_parsers['desolvate'].add_argument('--psf-outfile',type=str,default='dry.psf',help='name of output PSF file to create (default: %(default)s)')
    command_parsers['desolvate'].add_argument('--dcd-outfile',type=str,default='dry.dcd',help='name of DCD output file to create (default: %(default)s)')
    command_parsers['desolvate'].add_argument('--idx-outfile',type=str,default='dry.idx',help='name of index file for catdcd to create  (default: %(default)s)')
    command_parsers['desolvate'].add_argument('--dcd-stride',type=int,default=1,help='stride in number of frames for catdcd (default: %(default)s)')

    args=parser.parse_args()

    if not args.no_banner:
        banner(print)    
    if args.kick_ass_banner:
        print(enhanced_banner_message)

    if hasattr(args,'func'):
        args.func(args)
    else:
        my_list=oxford(list(commands.keys()))
        print(f'No subcommand found. Expected one of {my_list}')