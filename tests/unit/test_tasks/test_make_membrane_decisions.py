# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Unit tests for the decisions ``make_membrane_system`` makes around the toolchain calls.

A companion to ``test_make_membrane_system.py``, which covers the physics helpers
(``_per_leaflet_tension``, the area-convergence checks, ``_npatch``, ``_grid_half_mid_zgap``,
``_orientation_align``, the phase warnings) and runs the real end-to-end builds behind
``--runslow``.  This file takes the routing and setup logic between those two: what the task
decides *before* it places a lipid and *after* the relaxation returns.

That is where a mistake is silent.  A grid-placement failure announces itself; misreading the
schema's default ``embed`` block just builds the wrong kind of system, and a wrong stress-free
count ratio builds in exactly the differential stress the calibration exists to remove.

Nothing here needs VMD or NAMD: ``Bilayer`` and every method that shells out are stubbed, and what
is asserted is the arguments they were handed.
"""

import os
import tempfile
import unittest
from unittest import mock

import pandas as pd

from pestifer.tasks.make_membrane_system import (
    MakeMembraneSystemTask, _DEFAULT_AREA_MODULUS,
)

LOGGER = 'pestifer.tasks.make_membrane_system'


def _task(**attrs):
    t = MakeMembraneSystemTask.__new__(MakeMembraneSystemTask)
    t.taskname = 'make_membrane_system'
    t.specs = {}
    t.bilayer_specs = {}
    t.embed_specs = {}
    for k, v in attrs.items():
        setattr(t, k, v)
    return t


def _profile(nslabs=4, nframes=8, pxx=0.0, pyy=0.0, pzz=0.0, first_ts_zero=True):
    """A NAMD pressureprofile dataframe with constant per-slab components."""
    def per_slab(v):
        return list(v) if hasattr(v, '__len__') else [v] * nslabs
    xs, ys, zs = per_slab(pxx), per_slab(pyy), per_slab(pzz)
    ts = list(range(0 if first_ts_zero else 100,
                    100 * nframes + (0 if first_ts_zero else 100), 100))
    data = {'TS': ts}
    for i in range(nslabs):
        data[f'x_{i}'] = [xs[i]] * nframes
        data[f'y_{i}'] = [ys[i]] * nframes
        data[f'z_{i}'] = [zs[i]] * nframes
    return pd.DataFrame(data)


class TestEmbedRequested(unittest.TestCase):
    """The trap this guards is documented in the method and is worth a test of its own: ycleptic
    materializes an ``embed`` block with default scalars even when no embedding was asked for, so
    ``bool(specs.get('embed'))`` is *always* true on the real schema-expanded pipeline.  Getting
    this wrong does not fail loudly -- it declares a protein requirement for a bare-membrane build.
    """

    @staticmethod
    def _schema_default_embed():
        # what ycleptic materializes for a config that requested no embedding
        return {'xydist': 10.0, 'zdist': 10.0, 'margin': 2.0, 'orient': {}, 'z_ref_group': {}}

    def test_the_schema_default_block_is_not_an_embed_request(self):
        self.assertFalse(MakeMembraneSystemTask._embed_requested(
            {'embed': self._schema_default_embed()}))

    def test_no_embed_block_at_all(self):
        self.assertFalse(MakeMembraneSystemTask._embed_requested({}))

    def test_a_non_mapping_embed_is_not_a_request(self):
        self.assertFalse(MakeMembraneSystemTask._embed_requested({'embed': 'yes please'}))

    def test_a_z_axis_selection_is_a_request(self):
        for key in ('z_head_group', 'z_tail_group'):
            with self.subTest(key=key):
                specs = {'embed': {**self._schema_default_embed(), key: 'protein and resid 1'}}
                self.assertTrue(MakeMembraneSystemTask._embed_requested(specs))

    def test_an_explicit_no_orient_opt_in_is_a_request(self):
        specs = {'embed': {**self._schema_default_embed(), 'no_orient': True}}
        self.assertTrue(MakeMembraneSystemTask._embed_requested(specs))

    def test_an_explicit_orient_vector_pair_is_a_request(self):
        for key in ('source', 'target'):
            with self.subTest(key=key):
                specs = {'embed': {**self._schema_default_embed(), 'orient': {key: [0, 0, 1]}}}
                self.assertTrue(MakeMembraneSystemTask._embed_requested(specs))

    def test_a_z_ref_group_anchor_is_a_request_only_when_it_selects_something(self):
        base = self._schema_default_embed()
        self.assertTrue(MakeMembraneSystemTask._embed_requested(
            {'embed': {**base, 'z_ref_group': {'text': 'resid 42'}}}))
        self.assertFalse(MakeMembraneSystemTask._embed_requested(
            {'embed': {**base, 'z_ref_group': {'text': ''}}}))


class TestPipelineContract(unittest.TestCase):
    """What the task promises the pipeline checker."""

    def test_a_bare_membrane_build_requires_nothing(self):
        contract = MakeMembraneSystemTask.pipeline_contract({})
        self.assertFalse(contract.requires)

    def test_an_embedding_build_requires_a_protein_state(self):
        from pestifer.tasks.pipeline_contract import STATE
        contract = MakeMembraneSystemTask.pipeline_contract(
            {'embed': {'z_head_group': 'protein and resid 1'}})
        self.assertIn(STATE, contract.requires)

    def test_the_task_provides_a_solvated_state(self):
        from pestifer.tasks.pipeline_contract import SOLVATED, STATE
        contract = MakeMembraneSystemTask.pipeline_contract({})
        self.assertIn(STATE, contract.provides)
        self.assertIn(SOLVATED, contract.provides)

    def test_an_incoming_protein_is_never_discarded(self):
        """At runtime the task embeds whatever state is present and only builds from scratch when
        none is, so declaring discards_state would warn about a system it does not throw away."""
        self.assertFalse(MakeMembraneSystemTask.pipeline_contract({}).discards_state)

    def test_a_prior_solvation_is_flagged_as_suspect(self):
        from pestifer.tasks.pipeline_contract import SOLVATED
        self.assertIn(SOLVATED, MakeMembraneSystemTask.pipeline_contract({}).warn_if_present)


class TestGuardPiercedLipids(unittest.TestCase):
    """Grid placement leaves no pierced rings, but the opening condensing minimize can pull a tail
    through a sterol ring -- a topological defect that survives further minimization and then
    RATTLE-fails the first dynamics step.  The guard goes in just before dynamics begins."""

    MINIMIZE = {'md': {'ensemble': 'minimize', 'minimize': 1000}}
    NPT = {'md': {'ensemble': 'NPT', 'nsteps': 1000}}
    NVT = {'md': {'ensemble': 'NVT', 'nsteps': 500}}

    def test_a_ring_check_is_inserted_before_the_first_dynamics_stage(self):
        got = _task()._guard_pierced_lipids([self.MINIMIZE, self.NPT])
        self.assertEqual(list(got[0]), ['md'])
        self.assertIn('ring_check', got[1])

    def test_the_ring_check_targets_lipids_and_deletes_the_pierced_one(self):
        got = _task()._guard_pierced_lipids([self.MINIMIZE, self.NPT])
        spec = got[1]['ring_check']
        self.assertEqual(spec['segtypes'], ['lipid'])
        self.assertEqual(spec['delete'], 'piercee')

    def test_a_minimize_follows_the_ring_check(self):
        """Deleting the pierced lipid frees the threading one into a strained pose."""
        got = _task()._guard_pierced_lipids([self.MINIMIZE, self.NPT])
        self.assertEqual(got[2]['md']['ensemble'], 'minimize')

    def test_the_original_stages_are_all_preserved_in_order(self):
        protocol = [self.MINIMIZE, self.NPT, self.NVT]
        got = _task()._guard_pierced_lipids(protocol)
        # the guard adds a ring_check *and* a follow-up minimize, so the original stages are a
        # subsequence of the result rather than everything that is not a ring_check
        it = iter(got)
        self.assertTrue(all(stage in it for stage in protocol),
                        f'original stages not preserved in order: {got}')
        self.assertEqual(len(got), len(protocol) + 2)

    def test_the_guard_is_inserted_only_once(self):
        got = _task()._guard_pierced_lipids([self.MINIMIZE, self.NPT, self.NVT])
        self.assertEqual(sum('ring_check' in s for s in got), 1)

    def test_an_all_minimize_protocol_is_left_alone(self):
        protocol = [self.MINIMIZE, self.MINIMIZE]
        self.assertEqual(_task()._guard_pierced_lipids(protocol), protocol)

    def test_a_protocol_that_opens_with_dynamics_is_still_guarded(self):
        got = _task()._guard_pierced_lipids([self.NPT])
        self.assertIn('ring_check', got[0])

    def test_a_non_list_protocol_passes_through(self):
        for protocol in ({}, None, 'quilt'):
            with self.subTest(protocol=protocol):
                self.assertEqual(_task()._guard_pierced_lipids(protocol), protocol)


class _InitializeCase(unittest.TestCase):
    """Shared machinery for the two initialize suites.

    Named without a ``Test`` prefix so pytest does not collect it, and so the asymmetric suite
    below can reuse the helpers without re-running the symmetric tests under a second name.
    """

    def _run(self, bilayer_specs, asymmetric=False):
        t = _task(bilayer_specs=bilayer_specs)
        t.charmmff_content = mock.Mock()
        t.register = mock.Mock()
        made = []

        def fake_bilayer(composition_dict=None, *a, **kw):
            b = mock.Mock()
            b.register_species_pdbs = []
            b.asymmetric = asymmetric
            # the caller mutates the dict in place between constructions, so snapshot it
            import copy as _copy
            made.append(_copy.deepcopy(composition_dict if composition_dict is not None
                                       else kw.get('composition_dict')))
            return b

        with mock.patch('pestifer.tasks.make_membrane_system.Bilayer', side_effect=fake_bilayer):
            t.initialize()
        return t, made

    @staticmethod
    def _explicit(upper, lower, **extra):
        return {'composition': {'upper_leaflet': upper, 'lower_leaflet': lower, **extra}}


class TestInitialize(_InitializeCase):
    """What ``initialize`` decides before any lipid is placed.

    ``Bilayer`` is stubbed out: constructing a real one needs the force field and can trigger
    conformer generation.  What is under test is the composition dictionary handed to it -- how the
    leaflets are specified, where the solvent chambers come from, which phase each species is
    tagged with.
    """

    def test_an_explicit_composition_is_used_as_given(self):
        _t, made = self._run(self._explicit([{'name': 'POPC', 'frac': 1.0}],
                                            [{'name': 'POPE', 'frac': 1.0}]))
        self.assertEqual([d['name'] for d in made[0]['upper_leaflet']], ['POPC'])
        self.assertEqual([d['name'] for d in made[0]['lower_leaflet']], ['POPE'])

    def test_solvent_chambers_are_filled_in_when_the_composition_names_only_leaflets(self):
        """An explicit ``composition:`` block specifies leaflets only, but the asymmetric branch
        reads and swaps the chambers, so they have to exist by then."""
        _t, made = self._run(self._explicit([{'name': 'POPC', 'frac': 1.0}],
                                            [{'name': 'POPC', 'frac': 1.0}]))
        self.assertIn('upper_chamber', made[0])
        self.assertIn('lower_chamber', made[0])

    def test_an_empty_composition_falls_back_to_the_memgen_specstrings(self):
        _t, made = self._run({'composition': {'upper_leaflet': [], 'lower_leaflet': []},
                              'lipids': 'POPC//POPC', 'mole_fractions': '1.0//1.0'})
        self.assertTrue(made[0]['upper_leaflet'])
        self.assertEqual(made[0]['upper_leaflet'][0]['name'], 'POPC')

    def test_each_species_is_tagged_with_its_leaflet_phase(self):
        """The packer draws a phase-specific conformer ensemble, so the tag has to ride along on
        every species rather than being consulted once."""
        _t, made = self._run(self._explicit(
            [{'name': 'PSM', 'frac': 0.6}, {'name': 'CHL1', 'frac': 0.4}],
            [{'name': 'POPC', 'frac': 1.0}],
            upper_leaflet_phase='Lo', lower_leaflet_phase='Ld'))
        self.assertEqual({d['phase'] for d in made[0]['upper_leaflet']}, {'Lo'})
        self.assertEqual({d['phase'] for d in made[0]['lower_leaflet']}, {'Ld'})

    def test_the_default_phase_is_the_fluid_one(self):
        _t, made = self._run(self._explicit([{'name': 'POPC', 'frac': 1.0}],
                                            [{'name': 'POPC', 'frac': 1.0}]))
        self.assertEqual({d['phase'] for d in made[0]['upper_leaflet']}, {'Ld'})

    def test_an_explicit_species_phase_is_not_overwritten(self):
        _t, made = self._run(self._explicit(
            [{'name': 'POPC', 'frac': 1.0, 'phase': 'Lo'}], [{'name': 'POPC', 'frac': 1.0}],
            upper_leaflet_phase='Ld'))
        self.assertEqual(made[0]['upper_leaflet'][0]['phase'], 'Lo')

    def test_a_mismatched_phase_is_warned_about_during_setup(self):
        with self.assertLogs(LOGGER, level='WARNING') as cm:
            self._run(self._explicit([{'name': 'CHL1', 'frac': 0.5},
                                      {'name': 'POPC', 'frac': 0.5}],
                                     [{'name': 'POPC', 'frac': 1.0}],
                                     upper_leaflet_phase='Ld'))
        self.assertIn('upper_leaflet', ''.join(cm.output))

    def test_a_symmetric_build_makes_exactly_one_bilayer(self):
        t, made = self._run(self._explicit([{'name': 'POPC', 'frac': 1.0}],
                                           [{'name': 'POPC', 'frac': 1.0}]))
        self.assertEqual(len(made), 1)
        self.assertIsNotNone(t.patch)


class TestAsymmetricSymmetrization(_InitializeCase):
    """An asymmetric request is built as two *symmetric* patches -- one mirroring each leaflet --
    which are later combined.  The swapping is done by mutating one dict in place, so the risk is
    that the second patch inherits state from the first.
    """

    UPPER = [{'name': 'PSM', 'frac': 1.0}]
    LOWER = [{'name': 'POPC', 'frac': 1.0}]

    def _asym(self):
        return self._run(self._explicit(list(self.UPPER), list(self.LOWER),
                                        upper_leaflet_phase='Lo', lower_leaflet_phase='Ld'),
                         asymmetric=True)

    def test_two_further_patches_are_built(self):
        t, made = self._asym()
        self.assertEqual(len(made), 3, 'the probe patch plus one per leaflet')
        self.assertIsNone(t.patch, 'the asymmetric patch itself is not used directly')
        self.assertIsNotNone(t.patchA)
        self.assertIsNotNone(t.patchB)

    def test_the_first_patch_mirrors_the_upper_leaflet(self):
        _t, made = self._asym()
        a = made[1]
        self.assertEqual([d['name'] for d in a['upper_leaflet']], ['PSM'])
        self.assertEqual([d['name'] for d in a['lower_leaflet']], ['PSM'])

    def test_the_second_patch_mirrors_the_lower_leaflet(self):
        _t, made = self._asym()
        b = made[2]
        self.assertEqual([d['name'] for d in b['upper_leaflet']], ['POPC'])
        self.assertEqual([d['name'] for d in b['lower_leaflet']], ['POPC'])

    def test_each_patch_carries_its_own_leaflet_phase(self):
        """patchA is the upper leaflet's world and patchB the lower's, so the ordered/fluid
        distinction has to travel with them or both get packed from the same ensemble."""
        _t, made = self._asym()
        self.assertEqual({d['phase'] for d in made[1]['upper_leaflet']}, {'Lo'})
        self.assertEqual({d['phase'] for d in made[2]['upper_leaflet']}, {'Ld'})

    def test_the_stashed_grid_kwargs_are_not_disturbed_by_the_swaps(self):
        """The grid packer rebuilds a Bilayer from these at full box size; the swaps mutate the
        shared composition in place, which is why it is deep-copied when stashed."""
        t, _made = self._asym()
        stashed = t._bilayer_kwargs['composition_dict']
        self.assertEqual([d['name'] for d in stashed['upper_leaflet']], ['PSM'])
        self.assertEqual([d['name'] for d in stashed['lower_leaflet']], ['POPC'])


class TestProteinBoxDims(unittest.TestCase):
    """The gridded membrane must span the oriented protein plus a margin on each side."""

    PDB = ('ATOM      1  N   ALA A   1      -5.000  -3.000   0.000  1.00  0.00\n'
           'ATOM      2  CA  ALA A   1       5.000   3.000  10.000  1.00  0.00\n'
           'HETATM    3  O   HOH A   2       0.000   0.000   5.000  1.00  0.00\n'
           'TER\nEND\n')

    def _dims(self, margin=10.0, pdb=None):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'p.pdb')
            with open(path, 'w') as f:
                f.write(pdb if pdb is not None else self.PDB)
            state = mock.Mock()
            state.pdb.name = path
            t = _task(embed_specs={'xydist': margin})
            t.get_current_artifact = lambda k: state
            return t, t._protein_box_dims()

    def test_the_box_is_the_extent_plus_a_margin_on_each_side(self):
        _t, (Lx, Ly) = self._dims(margin=10.0)
        self.assertAlmostEqual(Lx, 10.0 + 20.0)      # x extent -5..5
        self.assertAlmostEqual(Ly, 6.0 + 20.0)       # y extent -3..3

    def test_the_bare_footprint_is_remembered_for_the_later_span_check(self):
        """embed_protein verifies the condensed box still spans the protein; without this it
        cannot tell an over-condensed quilt from a correctly sized one."""
        t, _dims = self._dims()
        self.assertEqual(t._protein_xy, (10.0, 6.0))

    def test_hetatm_records_are_included(self):
        # a bound ligand or a structural water still occupies space in the membrane
        pdb = self.PDB + 'HETATM    4  O   HOH A   3      50.000   0.000   0.000  1.00  0.00\n'
        _t, (Lx, _Ly) = self._dims(margin=0.0, pdb=pdb)
        self.assertAlmostEqual(Lx, 55.0)

    def test_coordinates_are_read_by_column_not_by_split(self):
        """PDB is a fixed-column format: adjacent coordinates can run together with no space,
        and a whitespace split would silently mis-parse them."""
        pdb = 'ATOM      1  N   ALA A   1    -100.000-100.000   0.000  1.00  0.00\n'
        _t, (Lx, Ly) = self._dims(margin=0.0, pdb=pdb + self.PDB)
        self.assertAlmostEqual(Lx, 105.0)
        self.assertAlmostEqual(Ly, 103.0)


class TestDo(unittest.TestCase):
    """Which branch runs, and in what order.

    Every method it dispatches to needs VMD or NAMD, so they are all stubbed; what is under test
    is the routing, which is decided by two independent facts -- whether a protein state arrived,
    and whether the requested bilayer was symmetric.
    """

    STEPS = ('orient_protein', 'build_grid_membrane', 'build_patch',
             'build_grid_membrane_asymmetric', 'embed_protein')

    def _run(self, state=True, patch=True, prebuilt=False):
        t = _task(using_prebuilt_bilayer=prebuilt, patch=(mock.Mock() if patch else None))
        called = []
        for name in self.STEPS:
            setattr(t, name, (lambda n=name: called.append(n)))
        artifact = None
        if state:
            artifact = mock.Mock()
            artifact.psf = mock.Mock()
            artifact.pdb = mock.Mock()
        t.get_current_artifact = lambda k: artifact
        t.pipeline = mock.Mock()
        rc = t.do()
        return t, called, rc

    def test_a_bare_membrane_build_neither_orients_nor_embeds(self):
        t, called, rc = self._run(state=False)
        self.assertEqual(called, ['build_grid_membrane'])
        self.assertEqual(rc, 0)
        self.assertFalse(t.embedding)

    def test_a_bare_build_promotes_the_membrane_to_the_pipeline_state(self):
        """Nothing else will: with no protein to embed, the quilt *is* the system."""
        t, _called, _rc = self._run(state=False)
        t.pipeline.rekey.assert_called_once_with('quilt_state', 'state')

    def test_an_embedding_build_orients_first_and_embeds_last(self):
        _t, called, _rc = self._run(state=True)
        self.assertEqual(called, ['orient_protein', 'build_grid_membrane', 'embed_protein'])

    def test_an_embedding_build_does_not_rekey_the_membrane(self):
        t, _called, _rc = self._run(state=True)
        t.pipeline.rekey.assert_not_called()

    def test_an_asymmetric_request_calibrates_patches_before_gridding(self):
        """``initialize`` leaves ``patch`` as None for an asymmetric build; the two symmetric
        patches are measured first so the grid can be laid at stress-free counts."""
        _t, called, _rc = self._run(state=False, patch=False)
        self.assertEqual(called, ['build_patch', 'build_grid_membrane_asymmetric'])

    def test_a_prebuilt_bilayer_skips_construction_entirely(self):
        _t, called, _rc = self._run(state=False, prebuilt=True)
        self.assertEqual(called, [])

    def test_a_partial_state_is_not_treated_as_a_protein(self):
        """A state carrying only coordinates cannot be embedded -- both a psf and a pdb are
        needed -- so it must fall through to a bare-membrane build rather than half-embedding."""
        for missing in ('psf', 'pdb'):
            with self.subTest(missing=missing):
                t = _task(using_prebuilt_bilayer=False, patch=mock.Mock())
                for name in self.STEPS:
                    setattr(t, name, mock.Mock())
                artifact = mock.Mock()
                setattr(artifact, missing, None)
                t.get_current_artifact = lambda k: artifact
                t.pipeline = mock.Mock()
                t.do()
                self.assertFalse(t.embedding)


class TestBuildGridMembraneSizing(unittest.TestCase):
    """How many lipids per leaflet the grid is laid with, and at what aspect ratio."""

    def _sizes(self, embedding, bilayer_specs, protein_xy=(100.0, 50.0)):
        t = _task(bilayer_specs=bilayer_specs, embedding=embedding)
        t._protein_box_dims = lambda: protein_xy
        captured = {}
        t._grid_membrane = lambda counts, sapl, aspect: captured.update(
            counts=counts, sapl=sapl, aspect=aspect)
        t.build_grid_membrane()
        return captured

    def test_a_bare_build_is_sized_from_the_patch_tiling(self):
        got = self._sizes(False, {'patch_nlipids': {'upper': 100, 'lower': 100},
                                  'npatch': [2, 3], 'SAPL': 70.0})
        self.assertEqual(got['counts'], {'upper': 600, 'lower': 600})
        self.assertEqual(got['sapl'], 70.0)

    def test_a_bare_build_uses_the_configured_aspect_ratio(self):
        got = self._sizes(False, {'xy_aspect_ratio': 1.5})
        self.assertEqual(got['aspect'], 1.5)

    def test_an_embedding_build_is_sized_from_the_protein_footprint(self):
        # 100 x 50 box at 75 A^2 per lipid -> 5000/75 = 66.7 -> 67 lipids per leaflet
        got = self._sizes(True, {'SAPL': 75.0}, protein_xy=(100.0, 50.0))
        self.assertEqual(got['counts'], {'upper': 67, 'lower': 67})

    def test_an_embedding_build_matches_the_protein_box_aspect(self):
        """A square membrane around an elongated protein would waste lipids and box."""
        got = self._sizes(True, {'SAPL': 75.0}, protein_xy=(100.0, 50.0))
        self.assertAlmostEqual(got['aspect'], 0.5)

    def test_the_grid_is_symmetric_in_both_leaflets(self):
        got = self._sizes(True, {'SAPL': 60.0})
        self.assertEqual(got['counts']['upper'], got['counts']['lower'])


class TestAsymmetricStressFreeCounts(unittest.TestCase):
    """The calibration that makes an asymmetric membrane stress-free.

    Each symmetric calibration patch is tensionless by symmetry, so its relaxed area divided by
    its lipid count is that leaflet's preferred area-per-lipid.  At a common box area both
    leaflets are tensionless when each sits at its own preferred APL, which fixes the count ratio
    ``n_upper/n_lower = apl_lower/apl_upper``.  Getting this ratio wrong builds in exactly the
    differential stress the procedure exists to remove, and nothing downstream would catch it.
    """

    def _build(self, apl_upper=60.0, apl_lower=75.0, n_cal=100, embedding=False,
               protein_xy=(120.0, 60.0), drift=None, diag=None, **bilayer_specs):
        specs = {'patch_nlipids': {'upper': n_cal, 'lower': n_cal}, 'npatch': [1, 1]}
        specs.update(bilayer_specs)
        t = _task(bilayer_specs=specs, embedding=embedding,
                  specs={'diagnose_differential_stress': diag or {}})
        t.patchA = mock.Mock(area=apl_upper * n_cal, area_drift=drift)
        t.patchB = mock.Mock(area=apl_lower * n_cal, area_drift=drift)
        t._protein_box_dims = lambda: protein_xy
        t.quilt = mock.Mock()
        t.diagnose_differential_stress = mock.Mock()
        captured = {}
        t._grid_membrane = lambda counts, sapl, aspect: captured.update(
            counts=counts, sapl=sapl, aspect=aspect)
        t.build_grid_membrane_asymmetric()
        return t, captured

    def test_the_leaflet_with_the_smaller_apl_gets_more_lipids(self):
        _t, got = self._build(apl_upper=60.0, apl_lower=75.0)
        self.assertGreater(got['counts']['upper'], got['counts']['lower'])

    def test_the_count_ratio_is_the_inverse_apl_ratio(self):
        _t, got = self._build(apl_upper=60.0, apl_lower=75.0, n_cal=100)
        # requested upper count is n_cal x npatch = 100; n_lower = 100 * 60/75 = 80
        self.assertEqual(got['counts'], {'upper': 100, 'lower': 80})

    def test_equal_preferred_areas_give_equal_counts(self):
        """A composition that happens to be balanced must not acquire an asymmetry."""
        _t, got = self._build(apl_upper=65.0, apl_lower=65.0)
        self.assertEqual(got['counts']['upper'], got['counts']['lower'])

    def test_the_requested_upper_count_honours_the_patch_tiling(self):
        _t, got = self._build(apl_upper=60.0, apl_lower=60.0, n_cal=100, npatch=[2, 2])
        self.assertEqual(got['counts']['upper'], 400)

    def test_an_embedding_build_sizes_both_leaflets_to_the_protein_box(self):
        # 120 x 60 = 7200 A^2; upper 7200/60 = 120, lower 7200/75 = 96
        _t, got = self._build(apl_upper=60.0, apl_lower=75.0, embedding=True,
                              protein_xy=(120.0, 60.0))
        self.assertEqual(got['counts'], {'upper': 120, 'lower': 96})
        self.assertAlmostEqual(got['aspect'], 0.5)

    def test_the_grid_is_laid_at_the_calibrated_area_not_loose(self):
        """Gridding ~15% loose left the quilt at ~0.64 g/cc and let water seep into the lateral
        gaps during the long condensation; the default is now a token clearance."""
        _t, got = self._build(apl_upper=60.0, apl_lower=75.0, n_cal=100)
        # target area = n_upper * apl_upper = 6000; box_SAPL = 1.05 * 6000 / max(100, 80)
        self.assertAlmostEqual(got['sapl'], 1.05 * 6000.0 / 100)

    def test_the_placement_clearance_is_configurable(self):
        _t, got = self._build(apl_upper=60.0, apl_lower=60.0, quilt_grid_slack=1.20)
        self.assertAlmostEqual(got['sapl'], 1.20 * 6000.0 / 100)

    def test_an_unequilibrated_calibration_cell_is_flagged(self):
        """The classic symptom is two leaflets of different composition reporting near-equal
        APLs because one under-condensed -- which looks fine and is not."""
        with self.assertLogs(LOGGER, level='WARNING') as cm:
            self._build(drift=0.25)
        out = ''.join(cm.output)
        self.assertIn('had not equilibrated', out)
        self.assertIn('upper', out)

    def test_a_settled_calibration_cell_is_silent(self):
        with self.assertNoLogs(LOGGER, level='WARNING'):
            self._build(drift=0.0)

    def test_an_unmeasured_drift_is_not_reported_as_a_problem(self):
        with self.assertNoLogs(LOGGER, level='WARNING'):
            self._build(drift=None)

    def test_the_stress_diagnostic_runs_only_when_asked(self):
        t, _got = self._build()
        t.diagnose_differential_stress.assert_not_called()
        t, _got = self._build(diag={'enabled': True})
        t.diagnose_differential_stress.assert_called_once()


class TestDifferentialStressDiagnostic(unittest.TestCase):
    """Diagnostic only: the membrane is never rebuilt, so what matters is that the advice it
    prints is arithmetically right."""

    def _run(self, pp_df, n_upper=100, n_lower=100, c_z=100.0, processor='cpu', **diag):
        t = _task(provisions={'processor-type': processor})
        t.equilibrate_bilayer = mock.Mock()
        mdplot = mock.Mock()
        mdplot.dataframes = {'pressureprofile': pp_df} if pp_df is not None else {}
        t.subcontroller = mock.Mock(tasks=[mdplot])
        membrane = mock.Mock()
        membrane.box = [[0, 0, 0], [0, 0, 0], [0, 0, c_z]]
        t.diagnose_differential_stress(membrane, n_upper, n_lower, diag)
        return t, membrane

    def test_the_diagnostic_is_skipped_on_a_gpu_build(self):
        """CUDA NAMD does not support pressureProfile; running it anyway would waste the pass."""
        with self.assertLogs(LOGGER, level='WARNING') as cm:
            t, _m = self._run(_profile(), processor='gpu')
        self.assertIn('pressureProfile', ''.join(cm.output))
        t.equilibrate_bilayer.assert_not_called()

    def test_a_missing_profile_is_reported_rather_than_silently_skipped(self):
        with self.assertLogs(LOGGER, level='WARNING') as cm:
            self._run(None)
        self.assertIn('no usable pressure profile', ''.join(cm.output))

    def test_the_measured_difference_is_recorded_on_the_membrane(self):
        _t, membrane = self._run(_profile(nslabs=4, pzz=[0, 0, 20, 20]))
        self.assertNotEqual(membrane.dgamma, 0.0)

    def test_balanced_leaflets_are_reported_as_balanced(self):
        with self.assertLogs(LOGGER, level='INFO') as cm:
            self._run(_profile(pxx=100.0, pyy=100.0, pzz=100.0))
        self.assertIn('look balanced', ''.join(cm.output))

    def test_an_imbalance_names_a_direction_and_a_count(self):
        """Direction is easy to get backwards, so it is worth stating the reasoning.

        gamma_leaflet ~ K_A*(A/A0 - 1), so the leaflet under *more* tension is the stretched one
        -- it has too few lipids for the area it is being held at.  Relieving it means moving
        lipids *into* it.  Here the upper leaflet carries the excess tension, so the advice is to
        move lipids from the lower to the upper.
        """
        with self.assertLogs(LOGGER, level='INFO') as cm:
            self._run(_profile(nslabs=4, pzz=[0, 0, 400, 400]), n_upper=100, n_lower=100)
        out = ''.join(cm.output)
        self.assertIn('gamma_upper          = +200.000', out)
        self.assertIn('move ~', out)
        self.assertIn('from the lower to the upper', out)

    def test_the_suggested_count_follows_the_area_modulus(self):
        """dm = -Dgamma / [K_A * (1/n_up + 1/n_low)] -- a softer membrane needs more lipids
        moved for the same tension difference."""
        def suggested(K_A):
            with self.assertLogs(LOGGER, level='INFO') as cm:
                self._run(_profile(nslabs=4, pzz=[0, 0, 400, 400]), area_modulus=K_A)
            line = next(l for l in cm.output if 'move ~' in l)
            return float(line.split('move ~')[1].split(' ')[0])
        stiff = suggested(_DEFAULT_AREA_MODULUS * 2)
        soft = suggested(_DEFAULT_AREA_MODULUS / 2)
        self.assertGreater(soft, stiff)
        self.assertAlmostEqual(soft / stiff, 4.0, places=1)


if __name__ == '__main__':
    unittest.main()


if __name__ == '__main__':
    unittest.main()
