import inspect
import sys
import logging
import importlib
import os
from pathlib import Path
logger=logging.getLogger(__name__)
from pestifer import PestiferResources
import shutil

def get_version():
    res_path=Path(PestiferResources.__file__).parent
    logger.debug(f'res_path {res_path}')
    res_parent=Path(res_path).parent
    logger.debug(f'res_parent {res_parent}')
    package_dir=Path(res_parent).parent
    pyproject_toml=os.path.join(package_dir,'pyproject.toml')
    logger.debug(f'pyproject_toml {pyproject_toml}')
    version='UNKNOWN'
    if os.path.exists(pyproject_toml):
        with open(pyproject_toml,'r') as f:
            data=f.read().split('\n')
        for line in data:
            tokens=line.split()
            if tokens[0]=='version':
                version=tokens[2].strip('"')
    return version

def is_tool(name):
    return shutil.which(name) is not None

def is_periodic(cell,xsc):
    if cell and os.path.exists(cell):
        with open(cell,'r') as f:
            lines=f.read().split('\n')
        if len(lines)<4:
            return False
        check=True
        check&=(lines[0].startswith('cellbasisvector'))
        check&=(lines[1].startswith('cellbasisvector'))
        check&=(lines[2].startswith('cellbasisvector'))
        check&=(lines[3].startswith('cellorigin'))
        return check
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

def isidentity(t):
    if t[0][0]==1.0 and t[1][1]==1.0 and t[2][2]==1.0:
        return True
    else:
        return False

def reduce_intlist(L):
    """reduce_intlist generates a "reduced-byte" representation of a list of integers by collapsing runs of adjacent integers into 'i to j' format. Example:

    [1,2,3,4,5,7,8,9,10,12] -> '1 to 5 7 to 10 12'

    :param L: list of integers
    :type L: list
    :return: string of reduced-byte representation
    :rtype: string
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

def inspect_classes(module,key=' ',use_yaml_headers_as_keys=False):
    importlib.import_module(module)
    if key!=' ':
        nonkey_classes={}
        for name,cls in inspect.getmembers(sys.modules[module], lambda x: inspect.isclass(x) and (x.__module__==module) and key not in x.__name__):
            if use_yaml_headers_as_keys:
                nkey=cls.yaml_header
            else:
                nkey=name
            nonkey_classes[nkey]=cls
        key_classes={}
        for name,cls in inspect.getmembers(sys.modules[module], lambda x: inspect.isclass(x) and (x.__module__==module) and key in x.__name__):
            if use_yaml_headers_as_keys:
                nkey=cls.yaml_header
            else:
                nkey=name
            key_classes[nkey]=cls
        return nonkey_classes,key_classes
    else:
        classes={}
        for name,cls in inspect.getmembers(sys.modules[module], lambda x: inspect.isclass(x) and (x.__module__==module)):
            if use_yaml_headers_as_keys:
                nkey=cls.yaml_header
            else:
                nkey=name
            classes[nkey]=cls
        return classes
    
def replace(data,match,repl):
    """Recursive value search-and-replace; data is either list or dictionary; nesting is ok

    :param data: list or dict
    :type data: list or dict
    :param match: replacement string; expected to be found as $(val) in values
    :type match: str
    :param repl: replacement
    :type repl: *
    """
    match_str=r'$('+match+r')'
    if isinstance(data,(dict,list)):
        for k,v in (data.items() if isinstance(data,dict) else enumerate(data)):
            if v==match_str:
                data[k]=repl
            elif type(v)==str and match_str in v:
                data[k]=data[k].replace(match_str,repl)
            replace(v,match,repl)