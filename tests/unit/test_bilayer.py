import pytest
from pestifer.tasks.bilayer import Bilayer,BilayerEmbedTask
from pestifer.config import Config
from pestifer.charmmtop import CharmmResiDatabase
from pestifer.scriptwriters import Psfgen,VMD,NAMD,Filewriter
from unittest.mock import patch

def test_bilayer_init_empty():
    with patch("pestifer.tasks.bilayer.logger.debug") as mock_logger:
        test_bilayer=Bilayer(composition_dict={})
        mock_logger.assert_called_once_with('Empty bilayer')

def test_bilayer_init_nonempty():
    test_bilayer=Bilayer(composition_dict={
        'upper_leaflet': [{'name':'POPC','frac':1.0,'conf':0}],                                            'lower_leaflet': [{'name':'POPE','frac':1.0,'conf':0}]})
    assert all([x in test_bilayer.lipid_names for x in ['POPE', 'POPC']])
    assert all([x in ['POPE', 'POPC'] for x in test_bilayer.lipid_names])

def test_bilayer_init_missing_upper_leaflet():
    with pytest.raises(AssertionError,match=f'Composition spec missing \'upper_leaflet\' directive'):
        test_bilayer=Bilayer(composition_dict={
        'lower_leaflet': [{'name':'POPE','frac':1.0,'conf':0}]})

def test_bilayer_init_missing_lower_leaflet():
    with pytest.raises(AssertionError,match=f'Composition spec missing \'lower_leaflet\' directive'):
        test_bilayer=Bilayer(composition_dict={
        'upper_leaflet': [{'name':'POPC','frac':1.0,'conf':0}]})
    
def test_bilayer_init_nolist_upper_leaflet():
    with pytest.raises(AssertionError,match=f'upper_leaflet is not a list'):
        test_bilayer=Bilayer(composition_dict={
        'upper_leaflet': {'name':'POPC','frac':1.0,'conf':0},
        'lower_leaflet': [{'name':'POPE','frac':1.0,'conf':0}]})

def test_bilayer_init_nolist_lower_leaflet():
    with pytest.raises(AssertionError,match=f'lower_leaflet is not a list'):
        test_bilayer=Bilayer(composition_dict={
        'upper_leaflet': [{'name':'POPC','frac':1.0,'conf':0}],
        'lower_leaflet': {'name':'POPE','frac':1.0,'conf':0}})

def test_bilayer_init_memgen_style():
    test_bilayer=Bilayer(lipid_specstring='POPC',ratio_specstring='1',conformers_specstring='1')
    assert test_bilayer.lipid_names == ['POPC']
    test_bilayer=Bilayer(lipid_specstring='POPC//POPE',ratio_specstring='1',conformers_specstring='1')
    assert all([x in test_bilayer.lipid_names for x in ['POPE', 'POPC']])
    assert all([x in ['POPE', 'POPC'] for x in test_bilayer.lipid_names])
    assert test_bilayer.composition['upper_leaflet'][0]['name'] == 'POPC'
    assert test_bilayer.composition['lower_leaflet'][0]['name'] == 'POPE'
    assert test_bilayer.composition['upper_leaflet'][0]['frac'] == 1.0
    assert test_bilayer.composition['lower_leaflet'][0]['frac'] == 1.0
    assert test_bilayer.composition['upper_leaflet'][0]['conf'] == 1
    assert test_bilayer.composition['lower_leaflet'][0]['conf'] == 1
    test_bilayer=Bilayer(lipid_specstring='POPC:CHL1//POPE:CHL1',ratio_specstring='0.5:0.5',conformers_specstring='1:1')
    assert all([x in test_bilayer.lipid_names for x in ['POPE', 'POPC','CHL1']])
    assert test_bilayer.composition['upper_leaflet'][0]['name'] == 'POPC'
    assert test_bilayer.composition['upper_leaflet'][1]['name'] == 'CHL1'
    assert test_bilayer.composition['lower_leaflet'][0]['name'] == 'POPE'
    assert test_bilayer.composition['lower_leaflet'][1]['name'] == 'CHL1'
    assert test_bilayer.composition['upper_leaflet'][0]['frac'] == 0.5
    assert test_bilayer.composition['lower_leaflet'][0]['frac'] == 0.5
    assert test_bilayer.composition['upper_leaflet'][0]['conf'] == 1
    assert test_bilayer.composition['lower_leaflet'][0]['conf'] == 1
    assert test_bilayer.composition['upper_leaflet'][1]['frac'] == 0.5
    assert test_bilayer.composition['lower_leaflet'][1]['frac'] == 0.5
    assert test_bilayer.composition['upper_leaflet'][1]['conf'] == 1
    assert test_bilayer.composition['lower_leaflet'][1]['conf'] == 1
    test_bilayer=Bilayer(lipid_specstring='POPC:CHL1//POPE:CHL1',ratio_specstring='0.75:0.25',conformers_specstring='3:4//7:2')
    assert all([x in test_bilayer.lipid_names for x in ['POPE', 'POPC','CHL1']])
    assert test_bilayer.composition['upper_leaflet'][0]['name'] == 'POPC'
    assert test_bilayer.composition['upper_leaflet'][1]['name'] == 'CHL1'
    assert test_bilayer.composition['lower_leaflet'][0]['name'] == 'POPE'
    assert test_bilayer.composition['lower_leaflet'][1]['name'] == 'CHL1'
    assert test_bilayer.composition['upper_leaflet'][0]['frac'] == 0.75
    assert test_bilayer.composition['lower_leaflet'][0]['frac'] == 0.75
    assert test_bilayer.composition['upper_leaflet'][0]['conf'] == 3
    assert test_bilayer.composition['lower_leaflet'][0]['conf'] == 7
    assert test_bilayer.composition['upper_leaflet'][1]['frac'] == 0.25
    assert test_bilayer.composition['lower_leaflet'][1]['frac'] == 0.25
    assert test_bilayer.composition['upper_leaflet'][1]['conf'] == 4
    assert test_bilayer.composition['lower_leaflet'][1]['conf'] == 2
    test_bilayer=Bilayer(lipid_specstring='POPC:CHL1//POPE:CHL1',ratio_specstring='0.75:0.25//0.33:0.67',conformers_specstring='3:4//7:2')
    assert all([x in test_bilayer.lipid_names for x in ['POPE', 'POPC','CHL1']])
    assert test_bilayer.composition['upper_leaflet'][0]['name'] == 'POPC'
    assert test_bilayer.composition['upper_leaflet'][1]['name'] == 'CHL1'
    assert test_bilayer.composition['lower_leaflet'][0]['name'] == 'POPE'
    assert test_bilayer.composition['lower_leaflet'][1]['name'] == 'CHL1'
    assert test_bilayer.composition['upper_leaflet'][0]['frac'] == 0.75
    assert test_bilayer.composition['lower_leaflet'][0]['frac'] == 0.33
    assert test_bilayer.composition['upper_leaflet'][0]['conf'] == 3
    assert test_bilayer.composition['lower_leaflet'][0]['conf'] == 7
    assert test_bilayer.composition['upper_leaflet'][1]['frac'] == 0.25
    assert test_bilayer.composition['lower_leaflet'][1]['frac'] == 0.67
    assert test_bilayer.composition['upper_leaflet'][1]['conf'] == 4
    assert test_bilayer.composition['lower_leaflet'][1]['conf'] == 2

def test_bilayer_set_slice_bounds():
    test_bilayer=Bilayer(lipid_specstring='POPC:CHL1//POPE:CHL1',ratio_specstring='0.75:0.25//0.33:0.67',conformers_specstring='3:4//7:2')
    test_bilayer.set_slice_bounds()
    assert test_bilayer.LC['z-lo']==-27.0
    assert test_bilayer.LC['z-hi']==-22.0
    assert test_bilayer.LL['z-lo']==-22.0
    assert test_bilayer.LL['z-hi']==0.0
    assert test_bilayer.UL['z-lo']==0.0
    assert test_bilayer.UL['z-hi']==22.0
    assert test_bilayer.UC['z-lo']==22.0
    assert test_bilayer.UC['z-hi']==27.0

def test_bilayer_volumizer():
    test_bilayer=Bilayer(lipid_specstring='POPC:CHL1//POPE:CHL1',ratio_specstring='0.75:0.25//0.33:0.67',conformers_specstring='3:4//7:2')
    test_bilayer.set_slice_bounds()
    test_bilayer.volumizer(1000.00)
    assert test_bilayer.LC['INIT-VOLUME'] == pytest.approx(5000.00, rel=1e-2)
    assert test_bilayer.LL['INIT-VOLUME'] == pytest.approx(22000.00, rel=1e-2)
    assert test_bilayer.UL['INIT-VOLUME'] == pytest.approx(22000.00, rel=1e-2)
    assert test_bilayer.UC['INIT-VOLUME'] == pytest.approx(5000.00, rel=1e-2)
    assert test_bilayer.LC['INIT-NWATEREQUIV'] == 167
    assert test_bilayer.LL['INIT-NWATEREQUIV'] == 736
    assert test_bilayer.UL['INIT-NWATEREQUIV'] == 736
    assert test_bilayer.UC['INIT-NWATEREQUIV'] == 167

def test_bilayer_build_patch():
    RDB=CharmmResiDatabase()
    RDB.add_stream('lipid')
    RDB.add_topology('toppar_all36_moreions.str',streamnameoverride='water_ions')
    C=Config()
    RM=C.RM
    pdb_collection=RM.pdb_collection
    test_bilayer=Bilayer(lipid_specstring='POPC:CHL1//POPE:CHL1',ratio_specstring='0.75:0.25//0.33:0.67',conformers_specstring='3:4//7:2',pdb_collection=pdb_collection)
    test_bilayer.set_slice_bounds()
    test_bilayer.volumizer(1000.00)
    test_bilayer.build_patch()
    assert test_bilayer.patch_ll_corner[0]==pytest.approx(0.0, rel=1e-2)
    assert test_bilayer.patch_ll_corner[1]==pytest.approx(0.0, rel=1e-2)
    assert test_bilayer.patch_ll_corner[2]==pytest.approx(-27.0, rel=1e-2)

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
