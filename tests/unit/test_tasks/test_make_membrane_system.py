import pytest
import os
import shutil
import unittest

from pestifer.core.artifacts import DataArtifact, StateArtifacts
from pestifer.core.config import Config
from pestifer.core.controller import Controller

from pestifer.scripters import PsfgenScripter

from pestifer.scripters.vmdscripter import VMDScripter
from pestifer.tasks.make_membrane_system import MakeMembraneSystemTask

from pestifer.tasks.validate import ValidateTask
from pestifer.util.util import protect_str_arg

class TestMakeMembraneSystem(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.controller = Controller().configure(Config().configure_new(), terminate=False)
        cls.scripters = cls.controller.config.scripters # shortcut
        cls.common_patch_relaxation_protocols = [
            {'md': {'ensemble': 'minimize', 'nsteps': 1000}},
            {'md': {'ensemble': 'NVT', 'nsteps': 1000}},
            {'md': {'ensemble': 'NPT', 'nsteps': 1000}},
            # {'md': {'ensemble': 'NPT', 'nsteps': 2000}},
            # {'md': {'ensemble': 'NPT', 'nsteps': 4000}}
        ]
        cls.common_quilt_relaxation_protocols = [
            {'md': {'ensemble': 'minimize', 'nsteps': 1000}},
            {'md': {'ensemble': 'NVT', 'nsteps': 1000}},
            {'md': {'ensemble': 'NPT', 'nsteps': 1000}}
        ]

    @pytest.mark.slow
    def test_makemembranesystem_popc(self):
        test_dir = '__test_makemembranesystem_popc'
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        os.mkdir(test_dir)
        os.chdir(test_dir)

        task_list = [
            {'make_membrane_system': {
                'bilayer': {
                    'SAPL': 75,
                    'npatch': [1, 1],
                    'composition': {
                        'upper_leaflet': [{'name': 'POPC', 'frac': 1.0, 'conf': 0}],
                        'lower_leaflet': [{'name': 'POPC', 'frac': 1.0, 'conf': 0}]
                    },
                    'relaxation_protocols':{
                        'patch': self.common_patch_relaxation_protocols,
                        'quilt': self.common_quilt_relaxation_protocols
                    }
                }
            }},
            {'validate': {
                'tests': [
                    {'residue_test': {
                        'name': 'Count of lipid residues',
                        'selection': 'resname POPC',
                        'measure': 'residue_count',
                        'value': 200
                    }}
                ]
            }}]
        self.controller.reconfigure_tasks(task_list)
        self.assertIsInstance(self.controller.tasks[0], MakeMembraneSystemTask)
        self.assertIsInstance(self.controller.tasks[1], ValidateTask)
        self.assertEqual(len(self.controller.tasks), 2)
        result = self.controller.do_tasks()
        task = self.controller.tasks[0]
        validation = self.controller.tasks[1]
        state: StateArtifacts = task.get_current_artifact('state')
        self.assertTrue(state is not None)
        self.assertTrue(state.pdb.exists())
        self.assertTrue(state.psf.exists())
        self.assertTrue(state.coor.exists())
        self.assertTrue(state.xsc.exists())

        validation_results = validation.get_current_artifact('validation_results')
        self.assertIsInstance(validation_results, DataArtifact)
        self.assertIsInstance(validation_results.data, dict)
        self.assertEqual(validation_results.data['npass'], 1)
        self.assertEqual(validation_results.data['nfail'], 0)

        self.assertIsInstance(result, dict)
        os.chdir('..')
        assert result[0]['result'] == 0

    @pytest.mark.slow
    def test_makemembranesystem_popc_pope(self):
        test_dir='__test_makemembranesystem_popc_pope'
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        os.mkdir(test_dir)
        os.chdir(test_dir)
        task_list = [
            {'make_membrane_system': {
                'bilayer': {
                    'SAPL': 75,
                    'patch_nlipids': {
                        'upper': 49,
                        'lower': 49
                    },
                    'npatch': [1, 1],
                    'composition': {
                        'upper_leaflet': [{'name': 'POPC', 'frac': 1.0, 'conf': 0}],
                        'lower_leaflet': [{'name': 'POPE', 'frac': 1.0, 'conf': 0}]
                    },
                    'relaxation_protocols':{
                        'patch': self.common_patch_relaxation_protocols,
                        'quilt': self.common_quilt_relaxation_protocols
                    }
                }
            }},
            {'validate': {
                'tests': [
                    {
                        'residue_test': {
                            'name': 'Count of POPC residues',
                            'selection': 'resname POPC',
                            'measure': 'residue_count',
                            'relation': '<=',
                            'value': 49
                        }
                    },
                    {
                        'residue_test': {
                            'name': 'Count of POPE residues',
                            'selection': 'resname POPE',
                            'measure': 'residue_count',
                            'relation': '<=',
                            'value': 49
                        }
                    }
                ]
            }}]
        self.controller.reconfigure_tasks(task_list)
        self.assertIsInstance(self.controller.tasks[0], MakeMembraneSystemTask)
        self.assertIsInstance(self.controller.tasks[1], ValidateTask)
        self.assertEqual(len(self.controller.tasks), 2)
        result = self.controller.do_tasks()
        task = self.controller.tasks[0]
        validation = self.controller.tasks[1]
        state: StateArtifacts = task.get_current_artifact('state')
        self.assertTrue(state is not None)
        self.assertTrue(state.pdb.exists())
        self.assertTrue(state.psf.exists())
        self.assertTrue(state.coor.exists())
        self.assertTrue(state.xsc.exists())

        validation_results = validation.get_current_artifact('validation_results')
        self.assertIsInstance(validation_results, DataArtifact)
        self.assertIsInstance(validation_results.data, dict)
        self.assertEqual(validation_results.data['npass'], 2)
        self.assertEqual(validation_results.data['nfail'], 0)

        self.assertIsInstance(result, dict)
        os.chdir('..')
        assert result[0]['result'] == 0

    @pytest.mark.slow
    def test_makemembranesystem_popc_chl1_psm_chl1(self):
        test_dir = '__test_makemembranesystem_popc_chl1_psm_chl1'
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        os.mkdir(test_dir)
        os.chdir(test_dir)

        task_list = [
            {'make_membrane_system': {
                'bilayer': {
                    'SAPL': 75,
                    'patch_nlipids': {
                        'upper': 64,
                        'lower': 64
                    },
                    'npatch': [1, 1],
                'composition':{
                    'upper_leaflet': [
                        {'name':'POPC','frac':0.5,'conf':0},
                        {'name':'CHL1','frac':0.5,'conf':0}],
                    'lower_leaflet': [
                        {'name':'PSM','frac':0.5,'conf':0},
                        {'name':'CHL1','frac':0.5,'conf':0}]},
                    'relaxation_protocols':{
                        'patch': self.common_patch_relaxation_protocols,
                        'quilt': self.common_quilt_relaxation_protocols
                    }
                }
            }},
            {'validate': {
                'tests': [
                    {
                        'residue_test': {
                            'name': 'Count of POPC residues',
                            'selection': 'resname POPC',
                            'measure': 'residue_count',
                            'relation': '<=',
                            'value': 32
                        }
                    },
                    {
                        'residue_test': {
                            'name': 'Count of PSM residues',
                            'selection': 'resname PSM',
                            'measure': 'residue_count',
                            'relation': '<=',
                            'value': 32
                        }
                    },
                    {
                        'residue_test': {
                            'name': 'Count of CHL1 residues',
                            'selection': 'resname CHL1',
                            'measure': 'residue_count',
                            'relation': '<=',
                            'value': 64
                        }
                    },
                ]
            }}]
        self.controller.reconfigure_tasks(task_list)
        self.assertIsInstance(self.controller.tasks[0], MakeMembraneSystemTask)
        self.assertIsInstance(self.controller.tasks[1], ValidateTask)
        self.assertEqual(len(self.controller.tasks), 2)
        result = self.controller.do_tasks()
        task = self.controller.tasks[0]
        validation = self.controller.tasks[1]
        state: StateArtifacts = task.get_current_artifact('state')
        self.assertTrue(state is not None)
        self.assertTrue(state.pdb.exists())
        self.assertTrue(state.psf.exists())
        self.assertTrue(state.coor.exists())
        self.assertTrue(state.xsc.exists())

        validation_results = validation.get_current_artifact('validation_results')
        self.assertIsInstance(validation_results, DataArtifact)
        self.assertIsInstance(validation_results.data, dict)
        self.assertEqual(validation_results.data['npass'], 3)
        self.assertEqual(validation_results.data['nfail'], 0)

        self.assertIsInstance(result, dict)
        os.chdir('..')
        assert result[0]['result'] == 0

    def test_makemembranesystem_embed_with_orient(self):
        test_dir = '__test_makemembranesystem_embed_with_orient'
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        os.mkdir(test_dir)
        os.chdir(test_dir)
        # psf = '5e8w-proteinonly.psf'
        # pdb = '5e8w-proteinonly.pdb'
        psf = 'van3.psf'
        pdb = 'van3.pdb'
        bilayer_psf = 'equilibrate.psf'
        bilayer_pdb = 'equilibrate.pdb'
        bilayer_xsc = 'equilibrate.xsc'
        input_data_dir = '../../fixtures/embed_inputs'
        for ftype in [psf, pdb, bilayer_psf, bilayer_pdb, bilayer_xsc]:
            shutil.copy(os.path.join(input_data_dir, ftype), '.')

        vm: VMDScripter = self.scripters['vmd']
        vm.newscript('orient')
        vm.usescript('bilayer_orient')
        vm.writescript('orient')
        result = vm.runscript(psf=psf,
                              pdb=pdb,
                              z_head_group=protect_str_arg("protein and resid 126"),
                              z_tail_group=protect_str_arg("protein and resid 158"),
                              o='oriented')
        opdb = 'oriented.pdb'
        pg: PsfgenScripter = self.scripters['psfgen']
        pg.newscript('embedded')
        pg.usescript('bilayer_embed')
        pg.writescript('embedded', guesscoord=False, regenerate=True, force_exit=True, writepsf=False, writepdb=False)
        result = pg.runscript(psf=psf,
                              pdb=opdb,
                              bilayer_psf=bilayer_psf,
                              bilayer_pdb=bilayer_pdb,
                              bilayer_xsc=bilayer_xsc,
                              z_ref_group=protect_str_arg("protein and resid 142"),
                              z_value=0.0,
                              z_dist=10.0,
                              o='embedded')
        os.chdir('..')
        assert result == 0

    @pytest.mark.slow
    def test_makemembranesystem_quilt_from_patch(self):
        test_dir = '__test_makemembranesystem_quilt_from_patch'
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        os.mkdir(test_dir)
        os.chdir(test_dir)
        datadir = '../../fixtures/quilt_inputs'
        basename = 'patch'
        for ftype in ['.coor','.psf','.pdb','.xsc']:
            shutil.copy(os.path.join(datadir,basename+ftype),'.')
        psfA = psfB = basename+'.psf'
        pdbA = pdbB = basename+'.pdb'
        xscA = xscB = basename+'.xsc'
        npatchx = npatchy = 3
        pg: PsfgenScripter = self.scripters['psfgen']
        pg.newscript(basename)
        pg.usescript('bilayer_quilt')
        pg.writescript(basename, guesscoord=False, regenerate=True, 
                       force_exit=True, writepsf=False, writepdb=False)
        result = pg.runscript(nx=npatchx, ny=npatchy, psfA=psfA, pdbA=pdbA,
                              psfB=psfB, pdbB=pdbB, xscA=xscA, xscB=xscB,
                              o='quilt_test_sym_3x3')
        os.chdir('..')
        assert result == 0

    @pytest.mark.slow
    def test_makemembranesystem_5e8w_psm_chl1_pope_chl1(self):
        test_dir='__test_makemembranesystem_5e8w_psm_chl1_pope_chl1'
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        os.mkdir(test_dir)
        os.chdir(test_dir)
        psf='5e8w-proteinonly.psf'
        pdb='5e8w-proteinonly.pdb'
        task_list = [
            {'continuation': {'psf': '5e8w-proteinonly.psf', 'pdb': '5e8w-proteinonly.pdb'}},
            {'make_membrane_system': {
                'bilayer': {
                    'SAPL': 70.0,
                    'composition': {
                        'upper_leaflet': [
                            {'name': 'PSM', 'frac': 0.5},
                            {'name': 'CHL1', 'frac': 0.5}
                        ],
                        'lower_leaflet': [
                            {'name': 'POPE', 'frac': 0.5},
                            {'name': 'CHL1', 'frac': 0.5}
                        ]
                    },
                    'relaxation_protocols': {
                        'patch': self.common_patch_relaxation_protocols,
                        'quilt': self.common_quilt_relaxation_protocols
                    }
                },
                'embed': {
                    'z_head_group': 'protein and resid 667',
                    'z_tail_group': 'protein and resid 710',
                    'z_ref_group': {
                        'text': 'protein and resid 696',
                        'z_value': 0.0,
                        'z_dist': 10.0
                    }
                }
            }},
            {'md': {
                'ensemble': 'minimize',
                'nsteps': 1000
            }},
            {'validate' : {
                'tests': [
                    {
                        'residue_test': {
                            'name': 'Count of PSM residues',
                            'selection': 'resname PSM',
                            'measure': 'residue_count',
                            'relation': '<=',
                            'value': 200
                        }
                    },
                    {
                        'residue_test': {
                            'name': 'Count of POPE residues',
                            'selection': 'resname POPE',
                            'measure': 'residue_count',
                            'relation': '<=',
                            'value': 200
                        }
                    },
                    {
                        'residue_test': {
                            'name': 'Count of CHL1 residues',
                            'selection': 'resname CHL1',
                            'measure': 'residue_count',
                            'relation': '<=',
                            'value': 400
                        }
                    },
                    {
                        'residue_test' : {
                            'name': 'Test of protein presence',
                            'selection': 'protein',
                            'measure': 'residue_count',
                            'relation': '>=',
                            'value': 1
                        }
                    }
                ]
            }
            }
        ]
        input_data_dir = '../../fixtures/embed_inputs'
        for ftype in [psf, pdb]:
            shutil.copy(os.path.join(input_data_dir,ftype),'.')
        self.controller.reconfigure_tasks(task_list)
        result = self.controller.do_tasks()
        os.chdir('..')
        assert result[0]['result'] == 0

