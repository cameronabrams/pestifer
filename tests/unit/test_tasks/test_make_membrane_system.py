import pytest
import os
import shutil
import unittest

import numpy as np
import pandas as pd

from pestifer.core.artifacts import DataArtifact, StateArtifacts
from pestifer.core.config import Config
from pestifer.core.controller import Controller

from pestifer.scripters import PsfgenScripter

from pestifer.scripters.vmd import VMDScripter
from pestifer.tasks.make_membrane_system import (
    MakeMembraneSystemTask, _AREA_CONVERGENCE_TOL, _per_leaflet_tension)

from pestifer.tasks.validate import ValidateTask
from pestifer.util.util import protect_str_arg

class _StubMDPlot:
    """Minimal stand-in for the trailing mdplot task: just carries a ``dataframes`` dict."""
    def __init__(self, xst=None):
        self.dataframes = {} if xst is None else {'xst': xst}


def _xst_from_area(area):
    """Build an xst-style dataframe whose a_x*b_y equals the given area trajectory."""
    area = np.asarray(area, dtype=float)
    return pd.DataFrame({'a_x': area, 'b_y': np.ones_like(area)})


class TestAreaConvergence(unittest.TestCase):
    """Unit tests for the calibration-cell area-convergence guard (no NAMD required)."""

    # _area_convergence does not touch self, so an unbound call with None is fine
    _check = staticmethod(lambda task, name: MakeMembraneSystemTask._area_convergence(None, task, name))

    def test_plateaued_area_is_converged(self):
        # a flat tail (already equilibrated) -> drift well under tolerance
        area = np.concatenate([np.linspace(8000, 6000, 60), np.full(60, 6000.0)])
        drift = self._check(_StubMDPlot(_xst_from_area(area)), 'flat')
        self.assertIsNotNone(drift)
        self.assertLess(abs(drift), _AREA_CONVERGENCE_TOL)

    def test_still_shrinking_area_flagged(self):
        # area decreasing monotonically right to the end -> drift exceeds tolerance, negative
        area = np.linspace(7000, 6000, 120)
        drift = self._check(_StubMDPlot(_xst_from_area(area)), 'shrinking')
        self.assertIsNotNone(drift)
        self.assertGreater(abs(drift), _AREA_CONVERGENCE_TOL)
        self.assertLess(drift, 0.0)

    def test_missing_trajectory_returns_none(self):
        self.assertIsNone(self._check(_StubMDPlot(None), 'nodata'))
        self.assertIsNone(self._check(None, 'notask'))

    def test_too_short_trajectory_returns_none(self):
        # only 10 frames -> tail of 5 < 8, cannot judge
        drift = self._check(_StubMDPlot(_xst_from_area(np.linspace(7000, 6000, 10))), 'short')
        self.assertIsNone(drift)


def _pp_df_from_dp(dp, nframes=4, c_z_unused=None):
    """Build a pressureprofile-style dataframe whose Pzz-(Pxx+Pyy)/2 equals dp per slab.

    Pxx=Pyy=0, Pzz=dp, repeated over nframes dynamics frames plus a junk TS=0 frame
    (which the parser must discard)."""
    nslabs = len(dp)
    cols = {'TS': [0] + [1000 * (k + 1) for k in range(nframes)]}
    for i in range(nslabs):
        cols[f'x_{i}'] = [0.0] * (nframes + 1)
        cols[f'y_{i}'] = [0.0] * (nframes + 1)
        # TS=0 frame gets a wildly wrong value to prove it is dropped
        cols[f'z_{i}'] = [1e6] + [float(dp[i])] * nframes
    return pd.DataFrame(cols)


class TestPerLeafletTension(unittest.TestCase):
    """Unit tests for the differential-stress (per-leaflet tension) parser (no NAMD)."""

    def test_symmetric_profile_zero_dgamma(self):
        # profile symmetric about the midplane -> upper and lower tensions equal
        nslabs = 40
        idx = np.arange(nslabs)
        dp = (idx - (nslabs - 1) / 2.0) ** 2          # symmetric about the central boundary
        res = _per_leaflet_tension(_pp_df_from_dp(dp), c_z=80.0)
        self.assertIsNotNone(res)
        self.assertEqual(res['nslabs'], nslabs)
        self.assertAlmostEqual(res['upper'], res['lower'], places=6)
        self.assertAlmostEqual(res['dgamma'], 0.0, places=6)

    def test_asymmetric_profile_signed_dgamma(self):
        # make the upper half more tensile than the lower -> nonzero, positive dgamma
        nslabs = 40
        dp = np.ones(nslabs)
        dp[nslabs // 2:] *= 3.0
        res = _per_leaflet_tension(_pp_df_from_dp(dp), c_z=80.0)
        self.assertIsNotNone(res)
        self.assertGreater(res['dgamma'], 0.0)
        # quantitative: gamma = 0.01 * dp * dz, dz = 80/40 = 2
        dz = 80.0 / nslabs
        self.assertAlmostEqual(res['lower'], 0.01 * (nslabs // 2) * 1.0 * dz, places=6)
        self.assertAlmostEqual(res['upper'], 0.01 * (nslabs // 2) * 3.0 * dz, places=6)

    def test_drops_pre_dynamics_frame(self):
        # the TS=0 frame holds 1e6; if it were not dropped the tensions would explode
        dp = np.ones(20)
        res = _per_leaflet_tension(_pp_df_from_dp(dp, nframes=6), c_z=60.0)
        self.assertIsNotNone(res)
        self.assertLess(abs(res['total']), 10.0)      # sane magnitude, not ~1e6

    def test_degenerate_inputs_return_none(self):
        self.assertIsNone(_per_leaflet_tension(None, c_z=80.0))
        self.assertIsNone(_per_leaflet_tension(pd.DataFrame(), c_z=80.0))
        self.assertIsNone(_per_leaflet_tension(pd.DataFrame({'TS': [1000]}), c_z=80.0))


class TestAutostageProtocol(unittest.TestCase):
    """Unit tests for the relaxation auto-staging guard against patch-grid overflow."""

    _stage = staticmethod(MakeMembraneSystemTask._autostage_protocol)

    def test_long_npt_split_into_ramp(self):
        out = self._stage([{'md': {'ensemble': 'NPT', 'nsteps': 2000}}], first=100, cap=2000)
        steps = [s['md']['nsteps'] for s in out]
        self.assertEqual(sum(steps), 2000)          # total preserved
        self.assertEqual(steps[0], 100)             # first run is short
        self.assertEqual(steps[0], min(steps))      # ...and the shortest (ramps up from there)
        self.assertGreater(len(steps), 1)           # actually split
        self.assertTrue(all(s <= 2000 for s in steps))

    def test_short_stage_and_minimize_passthrough(self):
        proto = [{'md': {'ensemble': 'minimize', 'nsteps': 500}},
                 {'md': {'ensemble': 'NVT', 'nsteps': 1000}},
                 {'md': {'ensemble': 'NPT', 'nsteps': 100}}]
        self.assertEqual(self._stage(proto, first=100), proto)   # nothing exceeds first / is barostatted-long

    def test_preserves_other_parameters(self):
        proto = [{'md': {'ensemble': 'NPT', 'nsteps': 400,
                         'other_parameters': {'useflexiblecell': True}}}]
        out = self._stage(proto, first=100, cap=2000)
        self.assertEqual(sum(s['md']['nsteps'] for s in out), 400)
        self.assertTrue(all(s['md']['other_parameters'] == {'useflexiblecell': True} for s in out))

    def test_empty_or_nonlist_passthrough(self):
        self.assertEqual(self._stage({}), {})
        self.assertEqual(self._stage([]), [])
        self.assertIsNone(self._stage(None))


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

    @pytest.mark.slow
    def test_makemembranesystem_5e8w_grid_asymmetric_embed(self):
        # asymmetric (PSM/CHL1 over POPE/CHL1) bilayer embedded around the 5e8w fragment,
        # built with the packmol-free grid packer: calibrate two symmetric patches, size the
        # box to the protein, grid at stress-free per-leaflet counts, stage the relaxation,
        # then embed. Exercises the asymmetric + embedding grid path end to end.
        test_dir = '__test_makemembranesystem_5e8w_grid_asymmetric_embed'
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        os.mkdir(test_dir)
        os.chdir(test_dir)
        psf = '5e8w-proteinonly.psf'
        pdb = '5e8w-proteinonly.pdb'
        task_list = [
            {'continuation': {'psf': psf, 'pdb': pdb}},
            {'make_membrane_system': {
                'bilayer': {
                    'packer': 'grid',
                    'SAPL': 70.0,
                    'patch_nlipids': {'upper': 36, 'lower': 36},
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
            {'validate': {
                'tests': [
                    {'residue_test': {
                        'name': 'Count of PSM residues', 'selection': 'resname PSM',
                        'measure': 'residue_count', 'relation': '<=', 'value': 200}},
                    {'residue_test': {
                        'name': 'Count of POPE residues', 'selection': 'resname POPE',
                        'measure': 'residue_count', 'relation': '<=', 'value': 200}},
                    {'residue_test': {
                        'name': 'Count of CHL1 residues', 'selection': 'resname CHL1',
                        'measure': 'residue_count', 'relation': '<=', 'value': 400}},
                    {'residue_test': {
                        'name': 'Test of protein presence', 'selection': 'protein',
                        'measure': 'residue_count', 'relation': '>=', 'value': 1}}
                ]
            }}
        ]
        input_data_dir = '../../fixtures/embed_inputs'
        for ftype in [psf, pdb]:
            shutil.copy(os.path.join(input_data_dir, ftype), '.')
        self.controller.reconfigure_tasks(task_list)
        result = self.controller.do_tasks()
        task = self.controller.tasks[1]
        validation = self.controller.tasks[2]
        state: StateArtifacts = task.get_current_artifact('state')
        self.assertTrue(state is not None)
        self.assertTrue(state.pdb.exists())
        self.assertTrue(state.psf.exists())
        validation_results = validation.get_current_artifact('validation_results')
        self.assertEqual(validation_results.data['nfail'], 0)
        os.chdir('..')
        assert result[0]['result'] == 0

