import pytest
from pestifer.bilayer import Bilayer
from pestifer.tasks.bilayertask import BilayerEmbedTask
from pestifer.config import Config
from pestifer.charmmtop import CharmmResiDatabase
from pestifer.scriptwriters import Psfgen,VMD,NAMD,Filewriter
from unittest.mock import patch

def test_bilayer_init_empty():
    with patch("pestifer.bilayer.logger.debug") as mock_logger:
        test_bilayer=Bilayer(composition_dict={})
        mock_logger.assert_called_once_with('Empty bilayer')

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

def test_bilayer_init_memgen_style():
    C=Config()
    RM=C.RM        
    pdb_collection=RM.pdb_collection
    resi_database=CharmmResiDatabase()
    resi_database.add_stream('lipid')
    resi_database.add_topology('toppar_all36_moreions.str',streamnameoverride='water_ions')
    test_bilayer=Bilayer(lipid_specstring='POPC',lipid_ratio_specstring='1',lipid_conformers_specstring='1',pdb_collection=pdb_collection,resi_database=resi_database)
    assert test_bilayer.lipid_names == ['POPC']
    assert test_bilayer.species_names == ['POPC','TIP3']

    test_bilayer=Bilayer(lipid_specstring='POPC//POPE',lipid_ratio_specstring='1',lipid_conformers_specstring='1',pdb_collection=pdb_collection,resi_database=resi_database)

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
    test_bilayer=Bilayer(lipid_specstring='POPC:CHL1//POPE:CHL1',lipid_ratio_specstring='0.5:0.5',lipid_conformers_specstring='1:1',pdb_collection=pdb_collection,resi_database=resi_database)
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

    test_bilayer=Bilayer(lipid_specstring='POPS:CHL1//POPE:CHL1',lipid_ratio_specstring='0.75:0.25',lipid_conformers_specstring='3:4//7:2',pdb_collection=pdb_collection,resi_database=resi_database)
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
    RDB=CharmmResiDatabase()
    RDB.add_stream('lipid')
    RDB.add_topology('toppar_all36_moreions.str',streamnameoverride='water_ions')
    C=Config()
    RM=C.RM
    pdb_collection=RM.pdb_collection
    test_bilayer=Bilayer(lipid_specstring='POPC:CHL1//POPE:CHL1',lipid_ratio_specstring='0.75:0.25//0.33:0.67',lipid_conformers_specstring='3:4//7:2',pdb_collection=pdb_collection,resi_database=RDB)
    test_bilayer.build_patch()
    assert test_bilayer.patch_ll_corner[0]==pytest.approx(0.0, rel=1e-2)
    assert test_bilayer.patch_ll_corner[1]==pytest.approx(0.0, rel=1e-2)
    assert test_bilayer.patch_ll_corner[2]==pytest.approx(-43.31, rel=1e-2)

def test_bilayer_task_init():
    C=Config()
    writers={
            'psfgen': Psfgen(C),
            'vmd':    VMD(C),
            'namd':   NAMD(C),
            'data':   Filewriter()
        }
    idict={'composition':{'upper_leaflet': [{'name':'POPC','frac':1.0,'conf':0}],
        'lower_leaflet': [{'name':'POPE','frac':1.0,'conf':0}]}}
    BET = BilayerEmbedTask(idict,'test_bilayer_task',C,writers,None)
    assert BET.taskname == 'test_bilayer_task'
    result=BET.do()
    assert result==0
