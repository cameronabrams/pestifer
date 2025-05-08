import pytest
import os
import shutil
from pestifer.bilayer import Bilayer,specstrings_builddict
from pestifer.tasks.bilayertask import BilayerEmbedTask
from pestifer.config import Config
from pestifer.charmmtop import CharmmResiDatabase
from pestifer.scriptwriters import Psfgen,VMD,NAMD,Filewriter
from pestifer.util.util import protect_str_arg
from unittest.mock import patch

def test_bilayer_init_empty():
    with patch("pestifer.bilayer.logger.debug") as mock_logger:
        test_bilayer=Bilayer(composition_dict={})
        mock_logger.assert_called_once_with('Empty bilayer')
        assert type(test_bilayer) == Bilayer
        assert test_bilayer.statevars == {}

def test_bilayer_init_nonempty():
    C=Config()
    RM=C.RM        
    pdb_collection=RM.pdb_collection
    resi_database=CharmmResiDatabase()
    resi_database.add_stream('lipid')
    resi_database.add_topology('toppar_all36_moreions.str',streamnameoverride='water_ions')
    test_bilayer=Bilayer(composition_dict={
        'upper_leaflet': [{'name':'POPC','frac':1.0,'conf':0}],
        'lower_leaflet': [{'name':'POPE','frac':1.0,'conf':0}]},
        pdb_collection=pdb_collection,
        resi_database=resi_database)
    assert all([x in test_bilayer.lipid_names for x in ['POPE', 'POPC']])
    assert all([x in ['POPE', 'POPC'] for x in test_bilayer.lipid_names])

def test_bilayer_memgen_to_composition_simple_symm():
    cdict=specstrings_builddict(lipid_specstring='POPC',lipid_ratio_specstring='1',lipid_conformers_specstring='1')
    assert len(cdict['upper_leaflet']) == 1
    assert len(cdict['lower_leaflet']) == 1
    assert len(cdict['upper_chamber']) == 1
    assert len(cdict['lower_chamber']) == 1
    assert cdict['upper_leaflet'][0]['name'] == 'POPC'
    assert cdict['upper_leaflet'][0]['frac'] == 1.0
    assert cdict['upper_leaflet'][0]['conf'] == 1
    assert cdict['lower_leaflet'][0]['name'] == 'POPC'
    assert cdict['lower_leaflet'][0]['frac'] == 1.0
    assert cdict['lower_leaflet'][0]['conf'] == 1
    assert 'patn' not in cdict['upper_leaflet'][0]
    assert 'patn' not in cdict['lower_leaflet'][0]
    assert 'charge' not in cdict['upper_leaflet'][0]
    assert 'charge' not in cdict['lower_leaflet'][0]
    assert 'MW' not in cdict['upper_leaflet'][0]
    assert 'MW' not in cdict['lower_leaflet'][0]
    assert len(cdict['upper_chamber']) == 1
    assert len(cdict['lower_chamber']) == 1
    assert cdict['upper_chamber'][0]['name'] == 'TIP3'
    assert cdict['upper_chamber'][0]['frac'] == 1.0
    assert cdict['lower_chamber'][0]['name'] == 'TIP3'
    assert cdict['upper_chamber'][0]['frac'] == 1.0
    assert 'conf' not in cdict['upper_chamber'][0]
    assert 'conf' not in cdict['lower_chamber'][0]
    assert 'patn' not in cdict['upper_chamber'][0]
    assert 'patn' not in cdict['lower_chamber'][0]
    assert 'charge' not in cdict['upper_chamber'][0]
    assert 'charge' not in cdict['lower_chamber'][0]
    assert 'MW' not in cdict['upper_chamber'][0]
    assert 'MW' not in cdict['lower_chamber'][0]

def test_bilayer_memgen_to_composition_simple_asymm():
    cdict=specstrings_builddict(lipid_specstring='POPC//POPE',lipid_ratio_specstring='1',lipid_conformers_specstring='1')
    assert len(cdict['upper_leaflet']) == 1
    assert len(cdict['lower_leaflet']) == 1
    assert len(cdict['upper_chamber']) == 1
    assert len(cdict['lower_chamber']) == 1
    assert cdict['upper_leaflet'][0]['name'] == 'POPC'
    assert cdict['upper_leaflet'][0]['frac'] == 1.0
    assert cdict['upper_leaflet'][0]['conf'] == 1
    assert cdict['lower_leaflet'][0]['name'] == 'POPE'
    assert cdict['lower_leaflet'][0]['frac'] == 1.0
    assert cdict['lower_leaflet'][0]['conf'] == 1
    assert 'patn' not in cdict['upper_leaflet'][0]
    assert 'patn' not in cdict['lower_leaflet'][0]
    assert 'charge' not in cdict['upper_leaflet'][0]
    assert 'charge' not in cdict['lower_leaflet'][0]
    assert 'MW' not in cdict['upper_leaflet'][0]
    assert 'MW' not in cdict['lower_leaflet'][0]
    assert len(cdict['upper_chamber']) == 1
    assert len(cdict['lower_chamber']) == 1
    assert cdict['upper_chamber'][0]['name'] == 'TIP3'
    assert cdict['upper_chamber'][0]['frac'] == 1.0
    assert cdict['lower_chamber'][0]['name'] == 'TIP3'
    assert cdict['upper_chamber'][0]['frac'] == 1.0
    assert 'conf' not in cdict['upper_chamber'][0]
    assert 'conf' not in cdict['lower_chamber'][0]
    assert 'patn' not in cdict['upper_chamber'][0]
    assert 'patn' not in cdict['lower_chamber'][0]
    assert 'charge' not in cdict['upper_chamber'][0]
    assert 'charge' not in cdict['lower_chamber'][0]
    assert 'MW ' not in cdict['upper_chamber'][0]
    assert 'MW' not in cdict['lower_chamber'][0]

def test_bilayer_memgen_to_composition_complex_symm():
    cdict=specstrings_builddict(lipid_specstring='POPC:POPE',lipid_ratio_specstring='0.25:0.75',lipid_conformers_specstring='1')
    assert len(cdict['upper_leaflet']) == 2
    assert len(cdict['lower_leaflet']) == 2
    assert cdict['upper_leaflet'][0]['name'] == 'POPC'
    assert cdict['upper_leaflet'][0]['frac'] == 0.25
    assert cdict['upper_leaflet'][0]['conf'] == 1
    assert cdict['upper_leaflet'][1]['name'] == 'POPE'
    assert cdict['upper_leaflet'][1]['frac'] == 0.75
    assert cdict['upper_leaflet'][1]['conf'] == 1
    assert cdict['lower_leaflet'][0]['name'] == 'POPC'
    assert cdict['lower_leaflet'][0]['frac'] == 0.25
    assert cdict['lower_leaflet'][0]['conf'] == 1
    assert cdict['lower_leaflet'][1]['name'] == 'POPE'
    assert cdict['lower_leaflet'][1]['frac'] == 0.75
    assert cdict['lower_leaflet'][1]['conf'] == 1
    assert 'patn' not in cdict['upper_leaflet'][0]
    assert 'patn' not in cdict['upper_leaflet'][1]
    assert 'patn' not in cdict['lower_leaflet'][0]
    assert 'patn' not in cdict['lower_leaflet'][1]
    assert 'charge' not in cdict['upper_leaflet'][0]
    assert 'charge' not in cdict['upper_leaflet'][1]
    assert 'charge' not in cdict['lower_leaflet'][0]
    assert 'charge' not in cdict['lower_leaflet'][1]
    assert 'MW' not in cdict['upper_leaflet'][0]
    assert 'MW' not in cdict['upper_leaflet'][1]
    assert 'MW' not in cdict['lower_leaflet'][0]
    assert 'MW' not in cdict['lower_leaflet'][1]
    assert len(cdict['upper_chamber']) == 1
    assert len(cdict['lower_chamber']) == 1
    assert cdict['upper_chamber'][0]['name'] == 'TIP3'
    assert cdict['upper_chamber'][0]['frac'] == 1.0
    assert cdict['lower_chamber'][0]['name'] == 'TIP3'
    assert cdict['upper_chamber'][0]['frac'] == 1.0
    assert 'conf' not in cdict['upper_chamber'][0]
    assert 'conf' not in cdict['lower_chamber'][0]
    assert 'patn' not in cdict['upper_chamber'][0]
    assert 'patn' not in cdict['lower_chamber'][0]
    assert 'charge' not in cdict['upper_chamber'][0]
    assert 'charge' not in cdict['lower_chamber'][0]
    assert 'MW' not in cdict['upper_chamber'][0]
    assert 'MW' not in cdict['lower_chamber'][0]

def test_bilayer_memgen_to_composition_complex_asymm():
    cdict=specstrings_builddict(lipid_specstring='POPC:POPE//POPC:CHL1',lipid_ratio_specstring='0.25:0.75//0.50:0.50',lipid_conformers_specstring='1')
    assert len(cdict['upper_leaflet']) == 2
    assert len(cdict['lower_leaflet']) == 2
    assert cdict['upper_leaflet'][0]['name'] == 'POPC'
    assert cdict['upper_leaflet'][0]['frac'] == 0.25
    assert cdict['upper_leaflet'][0]['conf'] == 1
    assert cdict['upper_leaflet'][1]['name'] == 'POPE'
    assert cdict['upper_leaflet'][1]['frac'] == 0.75
    assert cdict['upper_leaflet'][1]['conf'] == 1
    assert cdict['lower_leaflet'][0]['name'] == 'POPC'
    assert cdict['lower_leaflet'][0]['frac'] == 0.50
    assert cdict['lower_leaflet'][0]['conf'] == 1
    assert cdict['lower_leaflet'][1]['name'] == 'CHL1'
    assert cdict['lower_leaflet'][1]['frac'] == 0.50
    assert cdict['lower_leaflet'][1]['conf'] == 1
    assert 'patn' not in cdict['upper_leaflet'][0]
    assert 'patn' not in cdict['upper_leaflet'][1]
    assert 'patn' not in cdict['lower_leaflet'][0]
    assert 'patn' not in cdict['lower_leaflet'][1]
    assert 'charge' not in cdict['upper_leaflet'][0]
    assert 'charge' not in cdict['upper_leaflet'][1]
    assert 'charge' not in cdict['lower_leaflet'][0]
    assert 'charge' not in cdict['lower_leaflet'][1]
    assert 'MW' not in cdict['upper_leaflet'][0]
    assert 'MW' not in cdict['upper_leaflet'][1]
    assert 'MW' not in cdict['lower_leaflet'][0]
    assert 'MW' not in cdict['lower_leaflet'][1]
    assert len(cdict['upper_chamber']) == 1
    assert len(cdict['lower_chamber']) == 1
    assert cdict['upper_chamber'][0]['name'] == 'TIP3'
    assert cdict['upper_chamber'][0]['frac'] == 1.0
    assert cdict['lower_chamber'][0]['name'] == 'TIP3'
    assert cdict['upper_chamber'][0]['frac'] == 1.0
    assert 'conf' not in cdict['upper_chamber'][0]
    assert 'conf' not in cdict['lower_chamber'][0]
    assert 'patn' not in cdict['upper_chamber'][0]
    assert 'patn' not in cdict['lower_chamber'][0]
    assert 'charge' not in cdict['upper_chamber'][0]
    assert 'charge' not in cdict['lower_chamber'][0]
    assert 'MW' not in cdict['upper_chamber'][0]
    assert 'MW' not in cdict['lower_chamber'][0]


def test_bilayer_init_memgen_style():
    C=Config()
    RM=C.RM        
    pdb_collection=RM.pdb_collection
    resi_database=CharmmResiDatabase()
    resi_database.add_stream('lipid')
    resi_database.add_topology('toppar_all36_moreions.str',streamnameoverride='water_ions')
    cdict=specstrings_builddict(lipid_specstring='POPC',lipid_ratio_specstring='1.0',lipid_conformers_specstring='1')
    test_bilayer=Bilayer(composition_dict=cdict,
        pdb_collection=pdb_collection,resi_database=resi_database)
    assert test_bilayer.lipid_names == ['POPC']
    assert test_bilayer.species_names == ['POPC','TIP3']
    cdict=specstrings_builddict(lipid_specstring='POPC//POPE',lipid_ratio_specstring='1',lipid_conformers_specstring='1')
    test_bilayer=Bilayer(composition_dict=cdict,pdb_collection=pdb_collection,resi_database=resi_database)

    assert all([x in test_bilayer.lipid_names for x in ['POPE', 'POPC']])
    assert all([x in ['POPE', 'POPC'] for x in test_bilayer.lipid_names])
    assert all([x in test_bilayer.species_names for x in ['POPE', 'POPC', 'TIP3']])
    assert all([x in ['POPE', 'POPC', 'TIP3'] for x in test_bilayer.species_names])
    assert test_bilayer.slices['upper_leaflet']['composition'][0]['name'] == 'POPC'
    assert test_bilayer.slices['lower_leaflet']['composition'][0]['name'] == 'POPE'
    assert test_bilayer.slices['upper_leaflet']['composition'][0]['frac'] == 1.0
    assert test_bilayer.slices['lower_leaflet']['composition'][0]['frac'] == 1.0
    assert test_bilayer.slices['upper_leaflet']['composition'][0]['conf'] == 1
    assert test_bilayer.slices['lower_leaflet']['composition'][0]['conf'] == 1
    assert test_bilayer.slices['upper_leaflet']['composition'][0]['patn'] == 100
    assert test_bilayer.slices['lower_leaflet']['composition'][0]['patn'] == 100
    assert test_bilayer.slices['upper_leaflet']['composition'][0]['charge'] == 0.0
    assert test_bilayer.slices['lower_leaflet']['composition'][0]['charge'] == 0.0
    assert test_bilayer.slices['upper_chamber']['composition'][0]['name'] == 'TIP3'
    assert test_bilayer.slices['upper_chamber']['composition'][0]['frac'] == 1.0
    assert test_bilayer.slices['upper_chamber']['composition'][0]['charge'] == 0.0
    assert test_bilayer.slices['upper_chamber']['composition'][0]['patn'] == 3200
    assert test_bilayer.slices['upper_chamber']['composition'][0]['MW'] == 18.0154
    assert test_bilayer.slices['lower_chamber']['composition'][0]['name'] == 'TIP3'
    assert test_bilayer.slices['lower_chamber']['composition'][0]['frac'] == 1.0
    assert test_bilayer.slices['lower_chamber']['composition'][0]['charge'] == 0.0
    assert test_bilayer.slices['lower_chamber']['composition'][0]['patn'] == 3200
    assert test_bilayer.slices['lower_chamber']['composition'][0]['MW'] == 18.0154
    assert test_bilayer.asymmetric == True

    cdict=specstrings_builddict(lipid_specstring='POPC:CHL1//POPE:CHL1',lipid_ratio_specstring='0.5:0.5',lipid_conformers_specstring='1:1')
    test_bilayer=Bilayer(composition_dict=cdict,pdb_collection=pdb_collection,resi_database=resi_database)
    assert all([x in test_bilayer.lipid_names for x in ['POPE', 'POPC','CHL1']])
    assert all([x in ['POPE', 'POPC','CHL1'] for x in test_bilayer.lipid_names])
    assert all([x in test_bilayer.species_names for x in ['POPE', 'POPC','CHL1','TIP3']])
    assert all([x in ['POPE', 'POPC','CHL1','TIP3'] for x in test_bilayer.species_names])
    assert test_bilayer.slices['upper_leaflet']['composition'][0]['name'] == 'POPC'
    assert test_bilayer.slices['upper_leaflet']['composition'][1]['name'] == 'CHL1'
    assert test_bilayer.slices['lower_leaflet']['composition'][0]['name'] == 'POPE'
    assert test_bilayer.slices['lower_leaflet']['composition'][1]['name'] == 'CHL1'
    assert test_bilayer.slices['upper_leaflet']['composition'][0]['frac'] == 0.5
    assert test_bilayer.slices['lower_leaflet']['composition'][0]['frac'] == 0.5
    assert test_bilayer.slices['upper_leaflet']['composition'][0]['conf'] == 1
    assert test_bilayer.slices['lower_leaflet']['composition'][0]['conf'] == 1
    assert test_bilayer.slices['upper_leaflet']['composition'][1]['frac'] == 0.5
    assert test_bilayer.slices['lower_leaflet']['composition'][1]['frac'] == 0.5
    assert test_bilayer.slices['upper_leaflet']['composition'][1]['conf'] == 1
    assert test_bilayer.slices['lower_leaflet']['composition'][1]['conf'] == 1
    assert test_bilayer.slices['upper_leaflet']['composition'][0]['patn'] == 50
    assert test_bilayer.slices['lower_leaflet']['composition'][0]['patn'] == 50
    assert test_bilayer.slices['upper_leaflet']['composition'][1]['patn'] == 50
    assert test_bilayer.slices['lower_leaflet']['composition'][1]['patn'] == 50
    assert test_bilayer.slices['upper_leaflet']['composition'][0]['charge'] == 0.0
    assert test_bilayer.slices['lower_leaflet']['composition'][0]['charge'] == 0.0
    assert test_bilayer.slices['upper_leaflet']['composition'][1]['charge'] == 0.0
    assert test_bilayer.slices['lower_leaflet']['composition'][1]['charge'] == 0.0
    assert test_bilayer.slices['upper_chamber']['composition'][0]['name'] == 'TIP3'
    assert test_bilayer.slices['upper_chamber']['composition'][0]['frac'] == 1.0
    assert test_bilayer.slices['upper_chamber']['composition'][0]['charge'] == 0.0
    assert test_bilayer.slices['upper_chamber']['composition'][0]['patn'] == 3200
    assert test_bilayer.slices['upper_chamber']['composition'][0]['MW'] == 18.0154
    assert test_bilayer.slices['lower_chamber']['composition'][0]['name'] == 'TIP3'
    assert test_bilayer.slices['lower_chamber']['composition'][0]['frac'] == 1.0
    assert test_bilayer.slices['lower_chamber']['composition'][0]['charge'] == 0.0
    assert test_bilayer.slices['lower_chamber']['composition'][0]['patn'] == 3200
    assert test_bilayer.slices['lower_chamber']['composition'][0]['MW'] == 18.0154
    assert test_bilayer.asymmetric == True

    cdict=specstrings_builddict(lipid_specstring='POPS:CHL1//POPE:CHL1',lipid_ratio_specstring='0.75:0.25',lipid_conformers_specstring='3:4//7:2')
    test_bilayer=Bilayer(composition_dict=cdict,pdb_collection=pdb_collection,resi_database=resi_database)
    assert all([x in test_bilayer.lipid_names for x in ['POPS', 'POPE','CHL1']])
    assert all([x in ['POPS', 'POPE','CHL1'] for x in test_bilayer.lipid_names])
    assert all([x in test_bilayer.species_names for x in ['POPS', 'POPE','CHL1','TIP3']])
    assert all([x in ['POPS', 'POPE','CHL1','TIP3','POT'] for x in test_bilayer.species_names])
    assert test_bilayer.slices['upper_leaflet']['composition'][0]['name'] == 'POPS'
    assert test_bilayer.slices['upper_leaflet']['composition'][1]['name'] == 'CHL1'
    assert test_bilayer.slices['lower_leaflet']['composition'][0]['name'] == 'POPE'
    assert test_bilayer.slices['lower_leaflet']['composition'][1]['name'] == 'CHL1'
    assert test_bilayer.slices['upper_leaflet']['composition'][0]['frac'] == 0.75
    assert test_bilayer.slices['lower_leaflet']['composition'][0]['frac'] == 0.75
    assert test_bilayer.slices['upper_leaflet']['composition'][0]['conf'] == 3
    assert test_bilayer.slices['lower_leaflet']['composition'][0]['conf'] == 7
    assert test_bilayer.slices['upper_leaflet']['composition'][1]['frac'] == 0.25
    assert test_bilayer.slices['lower_leaflet']['composition'][1]['frac'] == 0.25
    assert test_bilayer.slices['upper_leaflet']['composition'][1]['conf'] == 4
    assert test_bilayer.slices['lower_leaflet']['composition'][1]['conf'] == 2
    assert test_bilayer.slices['upper_leaflet']['composition'][0]['patn'] == 75
    assert test_bilayer.slices['lower_leaflet']['composition'][0]['patn'] == 75
    assert test_bilayer.slices['upper_leaflet']['composition'][1]['patn'] == 25
    assert test_bilayer.slices['lower_leaflet']['composition'][1]['patn'] == 25
    assert test_bilayer.slices['upper_leaflet']['composition'][0]['charge'] == -1.0
    assert test_bilayer.slices['lower_leaflet']['composition'][0]['charge'] == 0.0
    assert test_bilayer.slices['upper_leaflet']['composition'][1]['charge'] == 0.0
    assert test_bilayer.slices['lower_leaflet']['composition'][1]['charge'] == 0.0
    assert test_bilayer.slices['upper_chamber']['composition'][0]['name'] == 'TIP3'
    assert test_bilayer.slices['upper_chamber']['composition'][0]['frac'] == 1.0
    assert test_bilayer.slices['upper_chamber']['composition'][0]['charge'] == 0.0
    assert test_bilayer.slices['upper_chamber']['composition'][0]['patn'] == 3200
    assert test_bilayer.slices['upper_chamber']['composition'][0]['MW'] == 18.0154
    assert test_bilayer.slices['lower_chamber']['composition'][0]['name'] == 'TIP3'
    assert test_bilayer.slices['lower_chamber']['composition'][0]['frac'] == 1.0
    assert test_bilayer.slices['lower_chamber']['composition'][0]['charge'] == 0.0
    assert test_bilayer.slices['lower_chamber']['composition'][0]['patn'] == 3200
    assert test_bilayer.slices['lower_chamber']['composition'][0]['MW'] == 18.0154
    assert test_bilayer.slices['upper_chamber']['composition'][1]['name'] == 'POT'
    assert test_bilayer.slices['upper_chamber']['composition'][1]['charge'] == 1.0
    assert test_bilayer.slices['upper_chamber']['composition'][1]['patn'] == 38
    assert test_bilayer.slices['upper_chamber']['composition'][1]['MW'] == 39.0983
    assert test_bilayer.slices['lower_chamber']['composition'][1]['name'] == 'POT'
    assert test_bilayer.slices['lower_chamber']['composition'][1]['charge'] == 1.0
    assert test_bilayer.slices['lower_chamber']['composition'][1]['patn'] == 37
    assert test_bilayer.slices['lower_chamber']['composition'][1]['MW'] == 39.0983
    assert test_bilayer.asymmetric == True

def test_bilayer_build_patch():
    if os.path.exists('__test_bilayer_build_patch'):
        shutil.rmtree('__test_bilayer_build_patch')
    os.mkdir('__test_bilayer_build_patch')
    os.chdir('__test_bilayer_build_patch')
    RDB=CharmmResiDatabase()
    RDB.add_stream('lipid')
    RDB.add_topology('toppar_all36_moreions.str',streamnameoverride='water_ions')
    C=Config()
    RM=C.RM
    pdb_collection=RM.pdb_collection
    cdict=specstrings_builddict(lipid_specstring='POPC:CHL1//POPE:CHL1',lipid_ratio_specstring='0.75:0.25//0.33:0.67',lipid_conformers_specstring='3:4//7:2')
    test_bilayer=Bilayer(composition_dict=cdict,pdb_collection=pdb_collection,resi_database=RDB)
    test_bilayer.build_patch()
    assert test_bilayer.patch_ll_corner[0]==pytest.approx(0.0, rel=1e-2)
    assert test_bilayer.patch_ll_corner[1]==pytest.approx(0.0, rel=1e-2)
    assert test_bilayer.patch_ll_corner[2]==pytest.approx(-43.31, rel=1e-2)
    os.chdir('..')
    
def test_bilayer_task_init_symmetric():
    if os.path.exists('__test_bilayer_task_symmetric'):
        shutil.rmtree('__test_bilayer_task_symmetric')
    os.mkdir('__test_bilayer_task_symmetric')
    os.chdir('__test_bilayer_task_symmetric')
    C=Config()
    writers={
            'psfgen': Psfgen(C),
            'vmd':    VMD(C),
            'namd':   NAMD(C),
            'data':   Filewriter()
        }
    idict={'bilayer':{
            'npatch':[2,2],
            'composition':{
                'upper_leaflet': [{'name':'POPC','frac':1.0,'conf':0}],
                'lower_leaflet': [{'name':'POPC','frac':1.0,'conf':0}]
                },
            'relaxation_protocols':{
                'patch':[
                    {'md':{'ensemble':'minimize','nsteps':1000}},
                    {'md':{'ensemble':'NVT','nsteps':1000}},
                    {'md':{'ensemble':'NPT','nsteps':1000}},
                    {'md':{'ensemble':'NPT','nsteps':2000}},
                    {'md':{'ensemble':'NPT','nsteps':4000}},
                    {'md':{'ensemble':'NPT','nsteps':8000}}
                ],
                'bilayer':[
                    {'md':{'ensemble':'minimize'}},
                    {'md':{'ensemble':'NVT','nsteps':1000}},
                    {'md':{'ensemble':'NPT','nsteps':10000}}
                ]}}}
    BET = BilayerEmbedTask(idict,'test_bilayer_task',C,writers,None)
    assert BET.taskname == 'test_bilayer_task'
    result=BET.do()
    os.chdir('..')
    assert result==0

def test_bilayer_task_init_asymmetric():
    if os.path.exists('__test_bilayer_task_asymmetric'):
        shutil.rmtree('__test_bilayer_task_asymmetric')
    os.mkdir('__test_bilayer_task_asymmetric')
    os.chdir('__test_bilayer_task_asymmetric')
    C=Config()
    writers={
            'psfgen': Psfgen(C),
            'vmd':    VMD(C),
            'namd':   NAMD(C),
            'data':   Filewriter()
        }
    idict={'bilayer':{
            'npatch':[2,2],
            'composition':{
                'upper_leaflet': [{'name':'POPC','frac':1.0,'conf':0}],
                'lower_leaflet': [{'name':'POPE','frac':1.0,'conf':0}]
                },
            'relaxation_protocols':{
                'patch':[
                    {'md':{'ensemble':'minimize','nsteps':1000}},
                    {'md':{'ensemble':'NVT','nsteps':1000}},
                    {'md':{'ensemble':'NPT','nsteps':1000}},
                    {'md':{'ensemble':'NPT','nsteps':2000}},
                    {'md':{'ensemble':'NPT','nsteps':4000}},
                    {'md':{'ensemble':'NPT','nsteps':8000}}
                ],
                'bilayer':[
                    {'md':{'ensemble':'minimize','nsteps':1000}},
                    {'md':{'ensemble':'NVT','nsteps':1000}},
                    {'md':{'ensemble':'NPT','nsteps':1000}},
                    {'md':{'ensemble':'NPT','nsteps':10000}}
                ]}}}
    BET = BilayerEmbedTask(idict,'test_bilayer_task',C,writers,None)
    assert BET.taskname == 'test_bilayer_task'
    result=BET.do()
    os.chdir('..')
    assert result==0

def test_bilayer_task_init_asymmetric_multicomponent():
    if os.path.exists('__test_bilayer_task_asymmetric_multicomponent'):
        shutil.rmtree('__test_bilayer_task_asymmetric_multicomponent')
    os.mkdir('__test_bilayer_task_asymmetric_multicomponent')
    os.chdir('__test_bilayer_task_asymmetric_multicomponent')
    C=Config()
    writers={
            'psfgen': Psfgen(C),
            'vmd':    VMD(C),
            'namd':   NAMD(C),
            'data':   Filewriter()
        }
    idict={'bilayer':{
            'npatch':[2,2],
            'composition':{
                'upper_leaflet': [
                    {'name':'POPC','frac':0.5,'conf':0},
                    {'name':'CHL1','frac':0.5,'conf':0}],
                'lower_leaflet': [
                    {'name':'PSM','frac':0.5,'conf':0},
                    {'name':'CHL1','frac':0.5,'conf':0}]},
            'relaxation_protocols':{
                'patch':[
                    {'md':{'ensemble':'minimize','nsteps':1000}},
                    {'md':{'ensemble':'NVT','nsteps':1000}},
                    {'md':{'ensemble':'NPT','nsteps':1000}},
                    {'md':{'ensemble':'NPT','nsteps':2000}},
                    {'md':{'ensemble':'NPT','nsteps':4000}},
                    {'md':{'ensemble':'NPT','nsteps':8000}}
                ],
                'bilayer':[
                    {'md':{'ensemble':'minimize','nsteps':1000}},
                    {'md':{'ensemble':'NVT','nsteps':1000}},
                    {'md':{'ensemble':'NPT','nsteps':1000}},
                    {'md':{'ensemble':'NPT','nsteps':10000}}
                ]}}}
    BET = BilayerEmbedTask(idict,'test_bilayer_task',C,writers,None)
    assert BET.taskname == 'test_bilayer_task'
    result=BET.do()
    os.chdir('..')
    assert result==0

def test_bilayer_task_embed():
    if os.path.exists('__test_bilayer_task_embed'):
        shutil.rmtree('__test_bilayer_task_embed')
    os.mkdir('__test_bilayer_task_embed')
    os.chdir('__test_bilayer_task_embed')
    basename='test_bilayer_embed'
    psf='5e8w-proteinonly.psf'
    pdb='5e8w-proteinonly.pdb'
    bilayer_psf='equilibrate.psf'
    bilayer_pdb='equilibrate.pdb'
    bilayer_xsc='equilibrate.xsc'
    input_data_dir='../../fixtures/embed_inputs'
    for ftype in [psf,pdb,bilayer_psf,bilayer_pdb,bilayer_xsc]:
        shutil.copy(os.path.join(input_data_dir,ftype),'.')
    C=Config()
    pg=Psfgen(C)
    pg.newscript(basename)
    pg.usescript('bilayer_embed')
    pg.writescript(basename,guesscoord=False,regenerate=True,force_exit=True,writepsf=False,writepdb=False)
    result=pg.runscript(psf=psf,
                        pdb=pdb,
                        bilayer_psf=bilayer_psf,
                        bilayer_pdb=bilayer_pdb,
                        bilayer_xsc=bilayer_xsc,
                        z_head_group=protect_str_arg("protein and resid 667"),
                        z_tail_group=protect_str_arg("protein and resid 710"),
                        z_ref_group=protect_str_arg("protein and resid 696"),
                        z_value=0.0,
                        o=basename)

    os.chdir('..')
    assert result==0