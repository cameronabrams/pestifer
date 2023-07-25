import unittest
import pytest
from pestifer.resources import ResourceManager
from pestifer.config import Config
from pestifer.mods import ModTypesRegistry as mtr
from pestifer.mods import Mutation

def test_moldata():
    r=ResourceManager()
    c=Config('example.yaml',r)
    d=c.data['BuildSteps']
    result_mods=[]
    for step in d:
        if 'mods' in step:
            for k,v in step['mods'].items():
                T=mtr.modtype(k)
                if T:
                    for vv in v:
                        if type(vv)==str: # shortcode
                            result_mods.append(T.from_shortcode(vv))
                        else:
                            result_mods.append(T(vv))
    assert(len(result_mods)==2)
    assert(all([type(x)==Mutation for x in result_mods]))
            

