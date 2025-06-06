# Author: Cameron F. Abrams <cfa22@drexel.edu>.
"""Various utility functions good for pestifer
"""
import glob
import importlib
import inspect
import logging
import os
import sys
import time
# import shutil
import pandas as pd
import numpy as np
from functools import wraps

logger=logging.getLogger(__name__)

_fntiminginfo={}
def countTime(fn):
    @wraps(fn)
    def measure_time(*args,**kwargs):
        keyname=f'{fn.__module__}.{fn.__name__}'
        if hasattr(fn,'__self__'):#.__class__.__name__=='method':
            keyname=f'{fn.__module__}.{type(fn.__self__).__name__}.{fn.__name__}'
        if not keyname in _fntiminginfo:
            _fntiminginfo[keyname]=dict(ncalls=0,totaltime=0.0,avgtimepercall=0.0)
        t1=time.time()
        result=fn(*args, **kwargs)
        t2=time.time()
        _fntiminginfo[keyname]['ncalls']+=1
        _fntiminginfo[keyname]['totaltime']+=t2-t1
        _fntiminginfo[keyname]['avgtimepercall']=_fntiminginfo[keyname]['totaltime']/_fntiminginfo[keyname]['ncalls']
        logger.debug(f'{keyname}: {(t2-t1)*1000:.6f} ms; avg {_fntiminginfo[keyname]["avgtimepercall"]*1000:.6f} ms/call; {_fntiminginfo[keyname]["ncalls"]} calls')
        return result
    return measure_time

def cell_from_xsc(xsc):
    if xsc and os.path.exists(xsc):
        celldf=pd.read_csv(xsc,skiprows=2,header=None,sep=r'\s+',index_col=None)
        col='step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z s_x s_y s_z s_u s_v s_w'.split()[:len(celldf.columns)]
        celldf.columns=col
        if len(celldf.columns)<13:
            return None,None
        avec=np.array(celldf.loc[0,['a_x','a_y','a_z']].to_list())
        bvec=np.array(celldf.loc[0,['b_x','b_y','b_z']].to_list())
        cvec=np.array(celldf.loc[0,['c_x','c_y','c_z']].to_list())
        box=np.array([avec,bvec,cvec])
        orig=np.array(celldf.loc[0,['o_x','o_y','o_z']].to_list())
        return box,orig
    return None,None

def cell_to_xsc(box,orig,xsc):
    with open(xsc,'w') as f:
        f.write('# NAMD extended system configuration output filewritten by pestifer\n')
        f.write('#$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z\n')
        f.write(f'0 {" ".join([f"{_:.6f}" for _ in box.reshape((9,)).tolist()])} {" ".join([f"{_:.6f}" for _ in orig])}\n')

# def is_tool(name):
#     """Checks to see if the object name is an executable"""
#     return shutil.which(name) is not None

def is_periodic(xsc):
    """Checks to see if the contents of the NAMD
    output xsc indicate that the current system is periodic 
    """
    if xsc and os.path.exists(xsc):
        with open(xsc,'r') as f:
            lines=f.read().split('\n')
        specline=lines[1]
        specfields=specline.split()
        reqdfieldlabels='a_x a_y a_z b_x b_y b_z c_x c_y c_z'.split()
        check=all([x in specfields for x in reqdfieldlabels])
        return check
    return False

def special_update(dict1,dict2):
    """Update dict1 with values from dict2 in a "special" way so that
    any list values are appended rather than overwritten
    """
    for k,v in dict2.items():
        ov=dict1.get(k,None)
        if not ov:
            dict1[k]=v
        else:
            if type(v)==list and type(ov)==list:
                for nv in v:
                    if not nv in ov:
                        ov.append(nv)
            elif type(v)==dict and type(ov)==dict:
                ov.update(v)
            else:
                dict1[k]=v # overwrite
    return dict1

def protect_str_arg(arg):
    # returns arg with spaces replaced with underscores
    if arg is None or type(arg) not in [str]:
        return ''
    return arg.replace(' ','_')

def reduce_intlist(L):
    """Generate a "reduced-byte" representation of a list of integers by collapsing 
    runs of adjacent integers into 'i to j' format. Example:

    [1,2,3,4,5,7,8,9,10,12] -> '1 to 5 7 to 10 12'

    Parameters
    ----------
    L: list
        list of integers
    
    Returns
    -------
    str: reduced representation as string

    """
    if not L:
        return ''
    ret=f'{L[0]}'
    if len(L)==2:
        ret+=f' {L[1]}'
        return ret
    inrun=False
    for l,r in zip(L[1:-1],L[2:]):
        adj=(r-l)==1
        if adj and not inrun:
            inrun=True
            ret+=f' to '
        elif not adj and inrun:
            ret+=f'{l} {r}'
            inrun=False
        elif not inrun:
            ret+=f' {l}'
    if inrun:
        ret+=f'{r}'
    return ret

def inspect_package_dir(dirname,key=' ',use_yaml_headers_as_keys=False):
    modules=glob.glob(f'{dirname}/*.py')
    packagename=os.path.split(dirname)[-1]
    for om in modules:
        if '__init__' in om:
            modules.remove(om)
            break
    obj_classes,objlist_classes={},{}
    for om in modules:
        modname=os.path.splitext(os.path.basename(om))[0]
        x,y=inspect_classes(f'pestifer.{packagename}.{modname}',key=key,use_yaml_headers_as_keys=use_yaml_headers_as_keys)
        obj_classes.update(x)
        objlist_classes.update(y)
    return obj_classes,objlist_classes

def inspect_classes(module,key=' ',use_yaml_headers_as_keys=False):
    """Returns the dictionary of names:classes for classes
    defined in the given module
    
    Optionally, if key is given a value, the function returns
    two lists: the first comprises classes that do NOT have
    the key in their names and the second comprises the classes
    that DO.

    Class names are used as keys in the dictionaries that are 
    returned, unless use_yaml_headers_as_keys is True; in that
    case, the string class attribute yaml_header (if the class
    defines one) is used as the key.  I don't know why I did 
    that.

    """
    importlib.import_module(module)
    if key!=' ':
        nonkey_classes={}
        for name,cls in inspect.getmembers(sys.modules[module], lambda x: inspect.isclass(x) and (x.__module__==module) and key not in x.__name__):
            if use_yaml_headers_as_keys and hasattr(cls,'yaml_header'):
                nkey=cls.yaml_header
            else:
                nkey=name
            nonkey_classes[nkey]=cls
        key_classes={}
        for name,cls in inspect.getmembers(sys.modules[module], lambda x: inspect.isclass(x) and (x.__module__==module) and key in x.__name__):
            if use_yaml_headers_as_keys and hasattr(cls,'yaml_header'):
                nkey=cls.yaml_header
            else:
                nkey=name
            key_classes[nkey]=cls
        return nonkey_classes,key_classes
    else:
        classes={}
        for name,cls in inspect.getmembers(sys.modules[module], lambda x: inspect.isclass(x) and (x.__module__==module)):
            if use_yaml_headers_as_keys and hasattr(cls,'yaml_header'):
                nkey=cls.yaml_header
            else:
                nkey=name
            classes[nkey]=cls
        return classes,{}
    
def replace(data,match,repl):
    """Recursive value search-and-replace; data is either list or dictionary; nesting is ok

    Parameters
    ----------
    match: str
        the string we are looking for in any objects in data
    repl: str
        the replacement for the match we are to put into the data

    """
    match_str=r'$('+match+r')'
    if isinstance(data,(dict,list)):
        for k,v in (data.items() if isinstance(data,dict) else enumerate(data)):
            if v==match_str:
                data[k]=repl
            elif type(v)==str and match_str in v:
                data[k]=data[k].replace(match_str,repl)
            replace(v,match,repl)

def flatten(current, key, result):
    """https://stackoverflow.com/questions/24448543/how-to-flatten-a-nested-dictionary
    
    Parameters
    ----------
    current: dict, or not a dict
       the current dictionary or the deepest non-dictionary value
    key:
       running key
    result:
       flattened key and value
    """
    if isinstance(current,dict):
        for k in current:
            new_key="{0}.{1}".format(key,k) if len(key)>0 else k
            flatten(current[k],new_key,result)
    else:
        result[key]=current
    return result

def write_residue_map(the_map,filename,valkeys=['chainID','resseqnum','insertion']):
    flat_map=flatten(the_map,'',{})
    with open(filename,'w') as f:
        for k,v in flat_map.items():
            kl=' '.join(k.split('.'))
            vl=' '.join([str(v.__dict__[x]) for x in valkeys])
            f.write(f'{kl} {vl}\n')

def tarball_walk(tar):
    """Simulates os.walk() for tarball contents. (chatgpt) """
    # This will hold the current directory and all its subdirectories and files
    directories = {}

    # Get all members (files and directories)
    for member in tar.getmembers():
        # Split the member's name into parts, treating '/' as directory separator
        parts = member.name.split('/')
        
        # Reconstruct directory structure by iterating over parts
        for i in range(1, len(parts)):
            dir_path = '/'.join(parts[:i])  # Get the directory path
            if dir_path not in directories:
                directories[dir_path] = []

        # Store the file under its directory
        dir_path = '/'.join(parts[:-1])  # Get the directory of the current file
        if dir_path not in directories:
            directories[dir_path] = []
        directories[dir_path].append(parts[-1])  # Add the file itself

    # Now iterate over directories and their files, similar to os.walk
    for dir_path, files in directories.items():
        yield dir_path, [], files

