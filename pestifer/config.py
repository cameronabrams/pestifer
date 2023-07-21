"""

.. module:: config
   :synopsis: Manages all configuration input and parsing
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import os
import logging
logger=logging.getLogger(__name__)
import yaml
from .resources import ResourceManager

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

class Config:
    def __init__(self,userconfigfilename,resman:ResourceManager):
        defaultconfigfilename=os.path.join(resman.resource_paths['config'],'defaults.yaml')
        with open(defaultconfigfilename,'r') as f:
            self.data=yaml.safe_load(f)
        with open(userconfigfilename,'r') as f:
            self.data.update(yaml.safe_load(f))
        vars={'HOME':resman.user_home}
        for v,r in vars.items():
            replace(self.data,v,r)
    def __str__(self):
        return yaml.dump(self.data)        
    