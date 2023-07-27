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

def test_parse():
    p=PDBParser(PDBcode='4zmj')
    p.fetch()
    p.read()
    p.parse()
    assert 'HEADER' in p.parsed
    assert 'VIRAL PROTEIN'==p.parsed['HEADER'].classification
    assert '04-MAY-15'==p.parsed['HEADER'].depDate
    assert '4ZMJ'==p.parsed['HEADER'].idCode
    assert type(p.parsed['COMPND'].compound)==list
    assert len(p.parsed['COMPND'].compound)==11
    assert p.parsed['COMPND'].compound[0] =='MOL_ID: 1'
    assert p.parsed['COMPND'].compound[1] =='MOLECULE: ENVELOPE GLYCOPROTEIN GP160'
    assert p.parsed['COMPND'].compound[2]=='A_FAKE_TOKEN: A_FAKE_VALUE'
    assert p.parsed['COMPND'].compound[3]=='CHAIN: G'
    assert p.parsed['TITLE'].title=='CRYSTAL STRUCTURE OF LIGAND-FREE BG505 SOSIP.664 HIV-1 ENV TRIMER THIS IS A FAKE EXTRA LINE THIS IS ANOTHER FAKE EXTRA LINE'
    assert type(p.parsed['COMPND'].tokens)==dict
    # print(p.parsed['COMPND'].tokens)
    assert 'MOL_ID.1' in p.parsed['COMPND'].tokens['compound'].keys()
    assert 'MOL_ID.2' in p.parsed['COMPND'].tokens['compound'].keys()
    assert len(p.parsed['COMPND'].tokens['compound'])==2
    assert type(p.parsed['COMPND'].tokens['compound']['MOL_ID.1'])==dict
    assert type(p.parsed['COMPND'].tokens['compound']['MOL_ID.2'])==dict
    assert len(p.parsed['COMPND'].tokens['compound']['MOL_ID.1'])==5
    assert p.parsed['COMPND'].tokens['compound']['MOL_ID.1']['MOLECULE']=='ENVELOPE GLYCOPROTEIN GP160'
    assert p.parsed['COMPND'].tokens['compound']['MOL_ID.1']['MUTATION']=='YES'
    assert type(p.parsed['SOURCE'].tokens)==dict
    assert 'MOL_ID.1' in p.parsed['SOURCE'].tokens['srcName'].keys()
    assert 'MOL_ID.2' in p.parsed['SOURCE'].tokens['srcName'].keys()
    assert len(p.parsed['SOURCE'].tokens['srcName'])==2
    assert p.parsed['SOURCE'].tokens['srcName']['MOL_ID.1']['ORGANISM_SCIENTIFIC']=='HUMAN IMMUNODEFICIENCY VIRUS 1'
    assert len(p.parsed['ATOM'])==4518
    assert len(p.parsed['ANISOU'])==4518
    assert len(p.parsed['HETATM'])==338
    assert len(p.parsed['HET'])==25
    assert len(p.parsed['LINK'])==25
    assert not 'MODRES' in p.parsed
    assert len(p.parsed['REVDAT'])==6
    print(p.parsed['REVDAT'][0].records)
    assert p.parsed['REVDAT'][0].records==['COMPND', 'REMARK', 'HETNAM', 'LINK', 'SITE', 'ATOM']
    assert len(p.parsed['SEQADV'])==10
    assert p.parsed['SEQADV'][0].resName=='ASN'
    assert p.parsed['SEQADV'][0].chainID=='G'
    assert p.parsed['SEQADV'][0].seqNum==332
    assert p.parsed['SEQADV'][0].iCode==''
    assert p.parsed['SEQADV'][0].database=='UNP'
    assert p.parsed['SEQADV'][0].dbAccession=='Q2N0S6'
    assert p.parsed['SEQADV'][0].dbRes=='THR'
    assert p.parsed['SEQADV'][0].dbSeq==330
    assert p.parsed['SEQADV'][0].conflict=='ENGINEERED MUTATION'
    assert len(p.parsed['SSBOND'])==11
    assert p.parsed['SSBOND'][2].chainID1=='G'
    assert p.parsed['SSBOND'][2].seqNum1==126
    assert p.parsed['SSBOND'][2].chainID2=='G'
    assert p.parsed['SSBOND'][2].seqNum2==196







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
            

