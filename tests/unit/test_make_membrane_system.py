import os
import shutil
import pytest
from pestifer.tasks.make_membrane_system import MakeMembraneSystemTask
from pestifer.config import Config
from pestifer.controller import Controller
from pestifer.scriptwriters import Psfgen,VMD,NAMD,Filewriter
from pestifer.util.util import protect_str_arg

@pytest.mark.slow
def test_membrane_symmetric_popc():
    if os.path.exists('__test_make_membrane_system_task_symmetric'):
        shutil.rmtree('__test_make_membrane_system_task_symmetric')
    os.mkdir('__test_make_membrane_system_task_symmetric')
    os.chdir('__test_make_membrane_system_task_symmetric')
    C=Config()
    writers={
            'psfgen': Psfgen(C),
            'vmd':    VMD(C),
            'namd':   NAMD(C),
            'data':   Filewriter()
        }
    config_specs={'bilayer':{
            'SAPL': 50,
            'npatch':[3,3],
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
    controller_specs={'controller_index':0,'taskname':'test_make_membrane_system_task','config':C,'writers':writers,'prior':None}
    BET = MakeMembraneSystemTask(config_specs,controller_specs)
    assert BET.taskname == 'test_make_membrane_system_task'
    result=BET.do()
    os.chdir('..')
    assert result==0

@pytest.mark.slow
def test_membrane_symmetric_c6dhpc():
    if os.path.exists('__test_make_membrane_system_task_symmetric'):
        shutil.rmtree('__test_make_membrane_system_task_symmetric')
    os.mkdir('__test_make_membrane_system_task_symmetric')
    os.chdir('__test_make_membrane_system_task_symmetric')
    C=Config()
    writers={
            'psfgen': Psfgen(C),
            'vmd':    VMD(C),
            'namd':   NAMD(C),
            'data':   Filewriter()
        }
    config_specs={'bilayer':{
            'SAPL': 50,
            'npatch':[3,3],
            'composition':{
                'upper_leaflet': [{'name':'C6DHPC','frac':1.0,'conf':0}],
                'lower_leaflet': [{'name':'C6DHPC','frac':1.0,'conf':0}]
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
    controller_specs={'controller_index':0,'taskname':'test_make_membrane_system_task','config':C,'writers':writers,'prior':None}
    BET = MakeMembraneSystemTask(config_specs,controller_specs)
    assert BET.taskname == 'test_make_membrane_system_task'
    result=BET.do()
    os.chdir('..')
    assert result==0

@pytest.mark.slow
def test_membrane_asymmetric_pure_leaflets():
    if os.path.exists('__test_make_membrane_system_task_asymmetric_pure_leaflets'):
        shutil.rmtree('__test_make_membrane_system_task_asymmetric_pure_leaflets')
    os.mkdir('__test_make_membrane_system_task_asymmetric_pure_leaflets')
    os.chdir('__test_make_membrane_system_task_asymmetric_pure_leaflets')
    C=Config()
    writers={
            'psfgen': Psfgen(C),
            'vmd':    VMD(C),
            'namd':   NAMD(C),
            'data':   Filewriter()
        }
    config_specs={'bilayer':{
            'SAPL': 50,
            'npatch':[2,2],
            'composition':{
                'upper_leaflet': [{'name':'POPC','frac':1.0,'conf':0}],
                'lower_leaflet': [{'name':'POPE','frac':1.0,'conf':0}]
                },
            'relaxation_protocols':{
                'patch':[
                    {'md':{'ensemble':'minimize','nsteps':1000}},
                    {'md':{'ensemble':'NVT','nsteps':1000}},
                    {'md':{'ensemble':'NPT','nsteps':1000,'pressure':10.0}},
                    {'md':{'ensemble':'NPT','nsteps':2000,'pressure':10.0}},
                    {'md':{'ensemble':'NPT','nsteps':4000,'pressure':10.0}},
                    {'md':{'ensemble':'NPT','nsteps':8000,'pressure':1.0}},
                    {'md':{'ensemble':'NPT','nsteps':16000,'pressure':1.0}},
                    {'md':{'ensemble':'NPT','nsteps':21000,'pressure':1.0}}
                ],
                'bilayer':[
                    {'md':{'ensemble':'minimize','nsteps':2000}},
                    {'md':{'ensemble':'NVT','nsteps':1000}},
                    {'md':{'ensemble':'NPT','nsteps':1000}},
                    {'md':{'ensemble':'NPT','nsteps':10000}}
                ]}}}
    controller_specs={'controller_index':0,'taskname':'test_make_membrane_system_task','config':C,'writers':writers,'prior':None}
    BET = MakeMembraneSystemTask(config_specs,controller_specs)
    assert BET.taskname == 'test_make_membrane_system_task'
    result=BET.do()
    os.chdir('..')
    assert result==0

@pytest.mark.slow
def test_membrane_asymmetric_multicomponent():
    if os.path.exists('__test_make_membrane_system_task_asymmetric_multicomponent'):
        shutil.rmtree('__test_make_membrane_system_task_asymmetric_multicomponent')
    os.mkdir('__test_make_membrane_system_task_asymmetric_multicomponent')
    os.chdir('__test_make_membrane_system_task_asymmetric_multicomponent')
    C=Config()
    writers={
            'psfgen': Psfgen(C),
            'vmd':    VMD(C),
            'namd':   NAMD(C),
            'data':   Filewriter()
        }
    config_specs={'bilayer':{
            'SAPL': 50,
            'npatch':[3,3],
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
                    {'md':{'ensemble':'NPT','nsteps':1000,'pressure':10.0}},
                    {'md':{'ensemble':'NPT','nsteps':2000,'pressure':10.0}},
                    {'md':{'ensemble':'NPT','nsteps':4000,'pressure':10.0}},
                    {'md':{'ensemble':'NPT','nsteps':8000,'pressure':10.0}},
                    {'md':{'ensemble':'NPT','nsteps':16000}},
                    {'md':{'ensemble':'NPT','nsteps':32000}},
                ],
                'bilayer':[
                    {'md':{'ensemble':'minimize','nsteps':1000}},
                    {'md':{'ensemble':'NVT','nsteps':1000}},
                    {'md':{'ensemble':'NPT','nsteps':1000,'pressure':10.0}},
                    {'md':{'ensemble':'NPT','nsteps':10000,'pressure':10.0}},
                    {'md':{'ensemble':'NPT','nsteps':20000}},
                    {'md':{'ensemble':'NPT','nsteps':40000}}
                ]}}}
    controller_specs={'controller_index':0,'taskname':'test_make_membrane_system_task','config':C,'writers':writers,'prior':None}
    BET = MakeMembraneSystemTask(config_specs,controller_specs)
    assert BET.taskname == 'test_make_membrane_system_task'
    result=BET.do()
    os.chdir('..')
    assert result==0

def test_membrane_embed():
    if os.path.exists('__test_make_membrane_system_task_embed'):
        shutil.rmtree('__test_make_membrane_system_task_embed')
    os.mkdir('__test_make_membrane_system_task_embed')
    os.chdir('__test_make_membrane_system_task_embed')
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
                        z_dist=10.0,
                        o=basename)

    os.chdir('..')
    assert result==0

@pytest.mark.slow
def test_membrane_md_prebuilt():
    if os.path.exists('__test_make_membrane_system_with_md_prebuilt'):
        shutil.rmtree('__test_make_membrane_system_with_md_prebuilt')
    os.mkdir('__test_make_membrane_system_with_md_prebuilt')
    os.chdir('__test_make_membrane_system_with_md_prebuilt')
    psf='5e8w-proteinonly.psf'
    pdb='5e8w-proteinonly.pdb'
    bilayer_psf='equilibrate.psf'
    bilayer_pdb='equilibrate.pdb'
    bilayer_xsc='equilibrate.xsc'
    yaml_file='test.yaml'
    input_data_dir='../../fixtures/embed_inputs'
    for ftype in [psf,pdb,bilayer_psf,bilayer_pdb,bilayer_xsc,yaml_file]:
        shutil.copy(os.path.join(input_data_dir,ftype),'.')
    config=Config(yaml_file)
    C=Controller(config)
    C.do_tasks()
    os.chdir('..')

@pytest.mark.slow
def test_membrane_md_build():
    if os.path.exists('__test_make_membrane_system_with_md_build'):
        shutil.rmtree('__test_make_membrane_system_with_md_build')
    os.mkdir('__test_make_membrane_system_with_md_build')
    os.chdir('__test_make_membrane_system_with_md_build')
    psf='5e8w-proteinonly.psf'
    pdb='5e8w-proteinonly.pdb'
    yaml_file='test2.yaml'
    input_data_dir='../../fixtures/embed_inputs'
    for ftype in [psf,pdb,yaml_file]:
        shutil.copy(os.path.join(input_data_dir,ftype),'.')
    config=Config(yaml_file)
    C=Controller(config)
    C.do_tasks()
    os.chdir('..')

def test_membrane_quilt():
    if os.path.exists('__test_make_membrane_system_task_quilt'):
        shutil.rmtree('__test_make_membrane_system_task_quilt')
    os.mkdir('__test_make_membrane_system_task_quilt')
    os.chdir('__test_make_membrane_system_task_quilt')
    datadir='../../fixtures/quilt_inputs'
    basename='patch'
    for ftype in ['.coor','.psf','.pdb','.xsc']:
        shutil.copy(os.path.join(datadir,basename+ftype),'.')
    psfA=psfB=basename+'.psf'
    pdbA=pdbB=basename+'.pdb'
    xscA=xscB=basename+'.xsc'
    npatchx=npatchy=3
    C=Config()
    pg=Psfgen(C)
    pg.newscript(basename)
    pg.usescript('bilayer_quilt')
    pg.writescript(basename,guesscoord=False,regenerate=True,force_exit=True,writepsf=False,writepdb=False)
    result=pg.runscript(nx=npatchx,ny=npatchy,psfA=psfA,pdbA=pdbA,
                                    psfB=psfB,pdbB=pdbB,xscA=xscA,xscB=xscB,o='quilt_test_sym_3x3')

    os.chdir('..')
    assert result==0

@pytest.mark.slow
def test_membrane_embed_no_orient():
    if os.path.exists('__test_make_membrane_system_task_embed_no_orient'):
        shutil.rmtree('__test_make_membrane_system_task_embed_no_orient')
    os.mkdir('__test_make_membrane_system_task_embed_no_orient')
    os.chdir('__test_make_membrane_system_task_embed_no_orient')
    basename='test_bilayer_embed_no_orient'
    psf='wt.psf'
    pdb='wt-flip.pdb'
    bilayer_psf='big_membrane.psf'
    bilayer_pdb='big_membrane.pdb'
    bilayer_xsc='big_membrane.xsc'
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
                        no_orient=True,
                        z_head_group=protect_str_arg(None),
                        z_tail_group=protect_str_arg(None),
                        z_ref_group=protect_str_arg("protein and resid 696"),
                        z_value=0.0,
                        z_dist=10.0,
                        o=basename)

    os.chdir('..')
    assert result==0