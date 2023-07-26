import unittest
import pytest
from pestifer.resources import ResourceManager
from pestifer.config import Config
from pestifer.mods import ModTypesRegistry as mtr
from pestifer.mods import Mutation
from pestifer.rcsb import PDBParser

def test_pdbformat():
    p=PDBParser()
    expected_sections=['record_types', 'delimiters', 'record_formats']
    assert(all([x in p.pdb_format_dict.keys() for x in expected_sections]))

def test_parse_type1():
    p=PDBParser(PDBcode='4zmj')
    p.fetch()
    p.read()
    p.parse()
    assert 'HEADER' in p.parsed
    #       classification: [String,[11,50]]
    #   depDate: [String,[51,59]]
    #   idCode: [String,[63,66]]
    assert 'VIRAL PROTEIN'==p.parsed['HEADER'].classification
    assert '04-MAY-15'==p.parsed['HEADER'].depDate
    assert '4ZMJ'==p.parsed['HEADER'].idCode
    assert type(p.parsed['COMPND'].compound)==list
    assert p.parsed['COMPND'].compound[0]=='MOL_ID: 1;'
    assert len(p.parsed['COMPND'].compound)==10

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
            

