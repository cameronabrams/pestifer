# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Unit tests for lipid-conformer collection building.

The expensive part of this module -- ``do_psfgen``, which builds and samples an actual lipid --
needs psfgen and NAMD and is exercised by the generator's own acceptance runs.  What is covered
here is everything around it: the phase-order targets, the bisection that tunes a lipid toward an
ordered ensemble, and the directory bookkeeping that decides where a built (or failed) conformer
set ends up.

That bookkeeping is worth testing precisely because it is unglamorous.  A conformer set filed in
the wrong place is not a crash; it is a cache miss that silently rebuilds a lipid on every run, or
worse, a failed build left where a successful one is expected.

``auto_cylinder_apl`` is covered in ``test_athermal_mc.py`` and is not repeated here.
"""

import os
import tempfile
import unittest
from unittest import mock

import numpy as np

from pestifer.charmmff.make_pdb_collection import (
    _PHASE_ORDER_TARGET, _bisect_to_order, all_chains_saturated, do_cleanup, do_resi,
    make_pdb_collection, phase_order_target,
)


class _Sandbox(unittest.TestCase):
    def setUp(self):
        self._cwd = os.getcwd()
        self._tmp = tempfile.TemporaryDirectory()
        os.chdir(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)
        self.addCleanup(os.chdir, self._cwd)


# --- phase targets ---------------------------------------------------------------------------------

class TestPhaseOrderTarget(unittest.TestCase):
    """Which phases ask for an ordering bias, and which are already what they claim to be."""

    def test_no_phase_means_no_tuning(self):
        self.assertIsNone(phase_order_target(None))

    def test_the_fluid_phase_needs_no_bias(self):
        """The pure athermal MC at trans bias 0 already *is* the Ld ensemble -- its order is the
        disordered floor, so tuning toward it would be tuning toward where it already is."""
        self.assertIsNone(phase_order_target('Ld'))

    def test_the_ordered_phase_has_a_numeric_target(self):
        self.assertEqual(phase_order_target('Lo'), _PHASE_ORDER_TARGET['Lo'])
        self.assertGreater(phase_order_target('Lo'), 0.0)

    def test_an_unknown_phase_is_refused_and_the_options_listed(self):
        with self.assertRaises(ValueError) as cm:
            phase_order_target('Lc')
        msg = str(cm.exception)
        self.assertIn('Lc', msg)
        self.assertIn('Lo', msg, 'the error should say what was expected')


class TestAllChainsSaturated(unittest.TestCase):
    """A cheap pre-filter, not a verdict.

    The docstring is explicit that this cannot decide whether a lipid will order: the
    ``unsaturations`` count does not distinguish cis from trans, so a sphingomyelin -- whose
    sphingosine base carries a *trans* double bond that does not kink -- reads as unsaturated here
    yet still orders under the trans bias.  These tests pin the cheap answer, including that it is
    the cheap one.
    """

    def test_a_fully_saturated_lipid(self):
        self.assertTrue(all_chains_saturated([{'unsaturations': 0}, {'unsaturations': 0}]))

    def test_one_unsaturated_chain_is_enough_to_fail(self):
        self.assertFalse(all_chains_saturated([{'unsaturations': 0}, {'unsaturations': 1}]))

    def test_a_missing_count_is_treated_as_saturated(self):
        self.assertTrue(all_chains_saturated([{}, {}]))

    def test_a_lipid_with_no_annotated_chains_is_not_saturated(self):
        """A sterol or detergent annotates empty; calling that 'all saturated' would let it
        through a filter it has no business passing."""
        self.assertFalse(all_chains_saturated([]))
        self.assertFalse(all_chains_saturated(None))

    def test_a_trans_unsaturation_still_reads_as_unsaturated_here(self):
        # the documented false negative: SM's trans double bond does not kink, but this cannot see
        # that, which is why the authoritative verdict is the tuned order rather than this filter
        self.assertFalse(all_chains_saturated([{'unsaturations': 1}, {'unsaturations': 0}]))


# --- order bisection -------------------------------------------------------------------------------

class TestBisectToOrder(unittest.TestCase):
    """Tuning a scalar knob until a lipid's chain order hits a phase target.

    Each ``order_of`` evaluation is a whole MC ensemble, so the iteration budget is small and the
    routine has to be right rather than lucky.  It infers the monotonic direction from the
    endpoints so it can serve either parametrization -- a trans-bias strength (order rises with it)
    or a cylinder looseness (order falls) -- and reports whether the target was bracketed at all.
    """

    @staticmethod
    def _linear(lo_val, hi_val, bounds=(0.0, 1.0)):
        a, b = bounds
        return lambda x: lo_val + (hi_val - lo_val) * (x - a) / (b - a)

    def test_a_bracketed_target_is_found_to_tolerance(self):
        knob, order, bracketed = _bisect_to_order(self._linear(0.0, 0.6), 0.30, (0.0, 1.0))
        self.assertTrue(bracketed)
        self.assertAlmostEqual(order, 0.30, delta=0.02)
        self.assertAlmostEqual(knob, 0.5, delta=0.05)

    def test_a_decreasing_response_works_too(self):
        """The cylinder-looseness parametrization runs the other way; inferring the direction from
        the endpoints is what lets one routine serve both."""
        knob, order, bracketed = _bisect_to_order(self._linear(0.6, 0.0), 0.30, (0.0, 1.0))
        self.assertTrue(bracketed)
        self.assertAlmostEqual(order, 0.30, delta=0.02)

    def test_a_target_above_the_reachable_range_is_reported_unbracketed(self):
        """A cis-unsaturated lipid cannot be ordered into Lo at any bias; saying so is the point,
        because the caller records it as Ld instead of pretending it succeeded."""
        knob, order, bracketed = _bisect_to_order(self._linear(0.0, 0.2), 0.30, (0.0, 1.0))
        self.assertFalse(bracketed)
        self.assertAlmostEqual(order, 0.2, delta=1e-9, msg='the closest reachable order')
        self.assertEqual(knob, 1.0, 'the bound that got closest')

    def test_a_target_below_the_reachable_range_is_reported_unbracketed(self):
        knob, order, bracketed = _bisect_to_order(self._linear(0.4, 0.8), 0.30, (0.0, 1.0))
        self.assertFalse(bracketed)
        self.assertEqual(knob, 0.0)

    def test_a_target_exactly_at_a_bound_counts_as_bracketed(self):
        _knob, order, bracketed = _bisect_to_order(self._linear(0.3, 0.8), 0.30, (0.0, 1.0))
        self.assertTrue(bracketed)
        self.assertAlmostEqual(order, 0.30)

    def test_an_early_hit_stops_evaluating(self):
        """Each evaluation is a full MC ensemble, so landing inside tolerance must end it."""
        calls = []

        def order_of(x):
            calls.append(x)
            return self._linear(0.0, 0.6)(x)

        _bisect_to_order(order_of, 0.30, (0.0, 1.0), tol=0.02)
        self.assertLessEqual(len(calls), 4, f'too many ensembles evaluated: {calls}')

    def test_the_iteration_budget_is_respected(self):
        calls = []

        def order_of(x):
            calls.append(x)
            return self._linear(0.0, 0.6)(x)

        _bisect_to_order(order_of, 0.2999, (0.0, 1.0), tol=1e-9, max_iters=3)
        self.assertEqual(len(calls), 2 + 3, 'two endpoints plus max_iters probes')

    def test_the_best_seen_value_is_returned_even_without_convergence(self):
        knob, order, bracketed = _bisect_to_order(self._linear(0.0, 0.6), 0.30, (0.0, 1.0),
                                                  tol=1e-12, max_iters=2)
        self.assertTrue(bracketed)
        self.assertAlmostEqual(order, 0.30, delta=0.1)
        self.assertGreater(knob, 0.0)

    def test_a_flat_response_does_not_hang(self):
        knob, order, bracketed = _bisect_to_order(lambda x: 0.1, 0.30, (0.0, 1.0))
        self.assertFalse(bracketed)
        self.assertEqual(order, 0.1)


# --- MC conformer sampling ---------------------------------------------------------------------------

class TestSampleAndWriteMcConformers(_Sandbox):
    """The athermal-MC sampler: melt a built lipid's tails into a fluid conformer set.

    The MC itself is covered in ``test_athermal_mc.py``; what is covered here is the wiring around
    it -- how the ordering target selects between a single probe and a tuned bisection, what is
    refused, and how the conformers are written out.

    The write-out deserves particular attention.  It substitutes coordinates into the template PDB
    by fixed columns, so everything outside columns 31-54 must survive byte-for-byte.  This
    repository has already been bitten once by a fixed-column PDB write that shifted fields and
    silently scrambled glycan residues downstream.
    """

    PDB = ('ATOM      1  N   POPC    1      10.000  11.000  12.000  1.00  0.00      MEMB N\n'
           'ATOM      2  C12 POPC    1      13.000  14.000  15.000  1.00  0.00      MEMB C\n'
           'ATOM      3  C13 POPC    1      16.000  17.000  18.000  1.00  0.00      MEMB C\n')

    def _atom(self, serial, name, atomtype='CTL2', wt=12.0):
        a = mock.Mock()
        a.serial, a.atomname, a.atomtype, a.atomicwt = serial, name, atomtype, wt
        return a

    def _sample(self, *, target_order=None, nsamples=3, order_curve=None,
                rotatable=(1,), missing_types=False, pdb=None):
        """Drive the sampler with the MC machinery stubbed out."""
        from pestifer.charmmff import make_pdb_collection as M

        with open('in.pdb', 'w') as f:
            f.write(pdb if pdb is not None else self.PDB)

        psf = mock.Mock()
        psf.atoms.data = [self._atom(1, 'N'), self._atom(2, 'C12'), self._atom(3, 'C13')]
        bond = mock.Mock(serial1=1, serial2=2)
        psf.bonds = [bond]

        params = mock.Mock()
        params.nonbonded = {} if missing_types else {'CTL2': mock.Mock(Rmin_half=2.0)}

        mol = mock.Mock()
        mol.rotatable = list(rotatable)
        mol.confined.sum.return_value = 3

        curve = order_curve if order_curve is not None else (lambda bias: 0.01 * bias)

        def fake_run_mc(m, nsamples=1, **kw):
            bias = kw.get('torsion_bias', 0.0)
            # one distinct coordinate set per sample, tagged by the bias so the writer can be
            # checked against the conformers actually selected
            return [np.full((3, 3), bias + i) for i in range(nsamples)]

        patches = {
            'pestifer.charmmff.make_pdb_collection.PSFContents': mock.Mock(return_value=psf),
            'pestifer.charmmff.charmmffprm.CharmmParamFile': mock.Mock(return_value=params),
            'pestifer.charmmff.athermal_mc.build_lipid_mc': mock.Mock(return_value=mol),
            'pestifer.charmmff.athermal_mc.run_mc': mock.Mock(side_effect=fake_run_mc),
            'pestifer.charmmff.athermal_mc.ensemble_chain_order':
                mock.Mock(side_effect=lambda samples, *a, **k: curve(float(samples[0][0][0]))),
            'pestifer.charmmff.athermal_mc.cylinder_radius_for_apl': mock.Mock(return_value=5.0),
        }
        ctxs = [mock.patch(k, v) for k, v in patches.items()]
        for c in ctxs:
            c.start()
            self.addCleanup(c.stop)

        DB = mock.Mock()
        DB.copy_charmmfile_local.side_effect = RuntimeError('no parameter files in this test')
        return M._sample_and_write_mc_conformers(
            'POPC', 'in.psf', 'in.pdb', heads=['N'], tails=['C13'], par_basenames=['par.prm'],
            DB=DB, nsamples=nsamples, digits=2, cylinder_apl=60.0, cylinder_inflation=1.9,
            n_equil=10, n_decorr=5, seed=1, max_angle=0.5, radius_scale=0.8,
            target_order=target_order)

    def test_missing_vdw_parameters_are_refused_by_name(self):
        """Running without Rmin/2 would silently disable the hard-sphere check."""
        with self.assertRaises(RuntimeError) as cm:
            self._sample(missing_types=True)
        self.assertIn('CTL2', str(cm.exception))

    def test_without_a_target_the_sampler_runs_once_unbiased(self):
        got = self._sample(target_order=None)
        self.assertEqual(got['torsion_bias'], 0.0)

    def test_with_a_target_the_bias_is_tuned(self):
        got = self._sample(target_order=0.20, order_curve=lambda b: 0.01 * b)
        self.assertGreater(got['torsion_bias'], 0.0)
        self.assertAlmostEqual(got['chain_order'], 0.20, delta=0.02)

    def test_an_unreachable_target_warns_and_uses_the_closest_bound(self):
        """A cis-unsaturated lipid cannot be ordered into Lo at any bias."""
        with self.assertLogs('pestifer.charmmff.make_pdb_collection', level='WARNING') as cm:
            got = self._sample(target_order=0.9, order_curve=lambda b: 0.001 * b)
        self.assertIn('not reachable', ''.join(cm.output))
        self.assertLess(got['chain_order'], 0.9)

    def test_a_lipid_with_no_rotatable_torsions_is_flagged(self):
        """Its 'conformers' are copies of one structure, which is worth saying out loud."""
        with self.assertLogs('pestifer.charmmff.make_pdb_collection', level='WARNING') as cm:
            self._sample(rotatable=())
        self.assertIn('no rotatable tail torsions', ''.join(cm.output))

    def test_one_pdb_is_written_per_conformer(self):
        got = self._sample(nsamples=4)
        written = sorted(f for f in os.listdir('.') if f.startswith('POPC-'))
        self.assertEqual(written, ['POPC-00.pdb', 'POPC-01.pdb', 'POPC-02.pdb', 'POPC-03.pdb'])
        self.assertEqual(got['nconf'], 4)

    def test_the_conformer_filenames_are_zero_padded_to_the_requested_width(self):
        self._sample(nsamples=2)
        self.assertIn('POPC-00.pdb', os.listdir('.'))

    def test_every_atom_record_is_written(self):
        self._sample(nsamples=1)
        lines = [l for l in open('POPC-00.pdb') if l.startswith('ATOM')]
        self.assertEqual(len(lines), 3)

    def test_the_file_is_terminated(self):
        self._sample(nsamples=1)
        self.assertEqual(open('POPC-00.pdb').read().splitlines()[-1], 'END')

    def test_everything_outside_the_coordinate_columns_survives_byte_for_byte(self):
        """The failure this guards against is silent: a shifted column changes a residue name or
        a segment ID and nothing downstream notices until the structure is wrong."""
        self._sample(nsamples=1)
        original = self.PDB.splitlines()
        written = [l for l in open('POPC-00.pdb').read().splitlines() if l.startswith('ATOM')]
        for before, after in zip(original, written):
            self.assertEqual(before[:30], after[:30], 'the record prefix moved')
            self.assertEqual(before[54:], after[54:], 'the trailing fields moved')

    def test_the_coordinates_are_replaced_with_the_sampled_ones(self):
        self._sample(nsamples=1, target_order=None)
        line = next(l for l in open('POPC-00.pdb') if l.startswith('ATOM'))
        # the stub returns an all-zeros first sample at bias 0
        self.assertEqual(line[30:54], f'{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}')

    def test_run_together_coordinates_are_read_by_column(self):
        """PDB is fixed-column: adjacent negative coordinates run together with no separating
        space, so a whitespace split would mis-parse them and a naive re-write would shift every
        field after them."""
        pdb = ('ATOM      1  N   POPC    1    -100.000-100.000-100.000  1.00  0.00      MEMB N\n'
               'ATOM      2  C12 POPC    1      13.000  14.000  15.000  1.00  0.00      MEMB C\n'
               'ATOM      3  C13 POPC    1      16.000  17.000  18.000  1.00  0.00      MEMB C\n')
        self._sample(nsamples=1, pdb=pdb)
        written = [l for l in open('POPC-00.pdb').read().splitlines() if l.startswith('ATOM')]
        self.assertEqual(written[0][54:].rstrip(), '  1.00  0.00      MEMB N')
        for before, after in zip(pdb.splitlines(), written):
            self.assertEqual(before[:30], after[:30])


class TestDoPsfgenInputGuard(unittest.TestCase):
    """``do_psfgen`` builds and samples an actual lipid, so it is exercised by the generator's own
    acceptance runs rather than here.  Its input guard fires before any of that and is worth a
    test: asking for more conformers than there are sampling steps cannot be satisfied, and
    finding that out after the build and minimization would waste the expensive part."""

    def test_more_conformers_than_steps_is_refused_up_front(self):
        from pestifer.charmmff.make_pdb_collection import do_psfgen
        DB = mock.Mock()
        with self.assertRaises(ValueError) as cm:
            do_psfgen('POPC', DB, nsamples=100, sample_steps=10)
        self.assertIn('100', str(cm.exception))
        self.assertIn('10', str(cm.exception))
        DB.get_resi.assert_not_called()


# --- directory bookkeeping -------------------------------------------------------------------------

class TestDoCleanup(_Sandbox):
    """What survives in a finished conformer directory: the build script, the metadata, the PSF
    and the conformer PDBs.  Everything else is scratch."""

    def _populate(self, resname='POPC', extra=()):
        os.mkdir('d')
        keep = ['init.tcl', 'info.yaml', f'{resname}-init.psf',
                f'{resname}-01.pdb', f'{resname}-02.pdb']
        for f in list(keep) + list(extra):
            open(os.path.join('d', f), 'w').close()
        return keep

    def test_the_scratch_files_are_removed(self):
        self._populate(extra=['scratch.coor', 'scratch.vel', 'run.log', 'tmp.namd'])
        do_cleanup('POPC', 'd')
        self.assertEqual(sorted(os.listdir('d')),
                         ['POPC-01.pdb', 'POPC-02.pdb', 'POPC-init.psf', 'info.yaml', 'init.tcl'])

    def test_every_conformer_pdb_is_kept(self):
        self._populate(extra=[f'POPC-{i:02d}.pdb' for i in range(3, 11)])
        do_cleanup('POPC', 'd')
        self.assertEqual(len([f for f in os.listdir('d') if f.endswith('.pdb')]), 10)

    def test_the_working_directory_is_restored(self):
        self._populate()
        before = os.getcwd()
        do_cleanup('POPC', 'd')
        self.assertEqual(os.getcwd(), before)

    def test_a_directory_with_nothing_to_clean_is_left_alone(self):
        keep = self._populate()
        do_cleanup('POPC', 'd')
        self.assertEqual(sorted(os.listdir('d')), sorted(keep))


class TestDoResi(_Sandbox):
    """Where a conformer set ends up, and whether it is rebuilt at all.

    A set filed in the wrong place is not a crash: it is a cache miss that rebuilds the lipid on
    every run, or a failed build sitting where a good one is expected.
    """

    def _run(self, result=0, force=False, prebuilt=False, cleanup=False):
        os.makedirs('data', exist_ok=True)
        os.makedirs('fails', exist_ok=True)
        if prebuilt:
            os.makedirs(os.path.join('data', 'POPC'))
            open(os.path.join('data', 'POPC', 'info.yaml'), 'w').close()

        def fake_psfgen(resi, DB, **kw):
            # do_psfgen runs inside the scratch dir do_resi created
            open('init.tcl', 'w').close()
            open('info.yaml', 'w').close()
            open(f'{resi}-init.psf', 'w').close()
            open(f'{resi}-01.pdb', 'w').close()
            open('scratch.log', 'w').close()
            return result

        with mock.patch('pestifer.charmmff.make_pdb_collection.do_psfgen',
                        side_effect=fake_psfgen) as psfgen:
            do_resi('POPC', mock.Mock(), outdir='data', faildir='fails',
                    force=force, cleanup=cleanup)
        return psfgen

    def test_a_successful_build_is_filed_under_the_output_directory(self):
        self._run(result=0)
        self.assertTrue(os.path.isfile(os.path.join('data', 'POPC', 'info.yaml')))
        self.assertFalse(os.path.exists('tmp'), 'the scratch directory should be gone')

    def test_a_failed_build_is_kept_for_inspection(self):
        """Filing a failure under the output directory would look like a cache hit forever."""
        self._run(result=1)
        self.assertTrue(os.path.isdir(os.path.join('fails', 'POPC')))
        self.assertFalse(os.path.exists(os.path.join('data', 'POPC')))

    def test_an_unknown_residue_is_reported_rather_than_filed_as_a_failure(self):
        with self.assertLogs('pestifer.charmmff.make_pdb_collection', level='WARNING') as cm:
            self._run(result=-2)
        self.assertIn('not found', ''.join(cm.output))
        self.assertFalse(os.path.exists(os.path.join('fails', 'POPC')))

    def test_an_existing_build_is_not_repeated(self):
        psfgen = self._run(prebuilt=True)
        psfgen.assert_not_called()

    def test_force_rebuilds_over_an_existing_set(self):
        psfgen = self._run(prebuilt=True, force=True)
        psfgen.assert_called_once()
        self.assertTrue(os.path.isfile(os.path.join('data', 'POPC', 'info.yaml')))

    def test_cleanup_strips_the_scratch_files_from_the_filed_set(self):
        self._run(result=0, cleanup=True)
        self.assertNotIn('scratch.log', os.listdir(os.path.join('data', 'POPC')))

    def test_without_cleanup_the_scratch_files_are_kept(self):
        self._run(result=0, cleanup=False)
        self.assertIn('scratch.log', os.listdir(os.path.join('data', 'POPC')))

    def test_the_working_directory_is_restored(self):
        before = os.getcwd()
        self._run(result=0)
        self.assertEqual(os.getcwd(), before)


# --- the collection driver -------------------------------------------------------------------------

class TestMakePdbCollection(_Sandbox):
    """The CLI entry point: which residues get built, and where."""

    @staticmethod
    def _args(**kw):
        base = dict(streamID=None, substreamID=None, resname='', topfile=None,
                    output_dir='', fail_dir='fails', force=False, cleanup=True,
                    lenfac=1.2, minimize_steps=500, sample_steps=5000, nsamples=10,
                    sample_temperature=300, refic_idx=0, force_constant=1.0,
                    take_ic_from=None, charmmff_release='')
        base.update(kw)
        return mock.Mock(**base)

    def _run(self, args, resnames=('POPC', 'POPE'), contains=True):
        cc = mock.MagicMock()
        cc.get_resnames_of_streamID.return_value = list(resnames)
        cc.__contains__.return_value = contains
        with mock.patch('pestifer.charmmff.make_pdb_collection.ResourceManager',
                        return_value=mock.Mock(charmmff_content=cc)), \
             mock.patch('pestifer.charmmff.make_pdb_collection.do_resi') as do_resi_mock:
            make_pdb_collection(args)
        return do_resi_mock, cc

    def test_a_named_residue_builds_only_that_one(self):
        do_resi_mock, _cc = self._run(self._args(resname='POPC'))
        do_resi_mock.assert_called_once()
        self.assertEqual(do_resi_mock.call_args.args[0], 'POPC')

    def test_a_stream_builds_every_residue_in_it(self):
        do_resi_mock, _cc = self._run(self._args(streamID='lipid'))
        self.assertEqual([c.args[0] for c in do_resi_mock.call_args_list], ['POPC', 'POPE'])

    def test_the_output_directory_defaults_to_the_stream_name(self):
        do_resi_mock, _cc = self._run(self._args(streamID='lipid'))
        self.assertEqual(do_resi_mock.call_args.kwargs['outdir'], 'lipid')
        self.assertTrue(os.path.isdir('lipid'))

    def test_the_output_directory_defaults_to_data_without_a_stream(self):
        do_resi_mock, _cc = self._run(self._args(resname='POPC'))
        self.assertEqual(do_resi_mock.call_args.kwargs['outdir'], 'data')

    def test_an_explicit_output_directory_wins(self):
        do_resi_mock, _cc = self._run(self._args(resname='POPC', output_dir='mine'))
        self.assertEqual(do_resi_mock.call_args.kwargs['outdir'], 'mine')
        self.assertTrue(os.path.isdir('mine'))

    def test_an_empty_failure_directory_is_removed(self):
        """Leaving it behind reads as 'something failed' on every successful run."""
        self._run(self._args(resname='POPC'))
        self.assertFalse(os.path.exists('fails'))

    def test_a_non_empty_failure_directory_is_kept_and_announced(self):
        cc = mock.MagicMock()
        cc.__contains__.return_value = True
        with mock.patch('pestifer.charmmff.make_pdb_collection.ResourceManager',
                        return_value=mock.Mock(charmmff_content=cc)), \
             mock.patch('pestifer.charmmff.make_pdb_collection.do_resi',
                        side_effect=lambda *a, **k: open(os.path.join('fails', 'x'), 'w').close()):
            with self.assertLogs('pestifer.charmmff.make_pdb_collection',
                                 level='WARNING') as cm:
                make_pdb_collection(self._args(resname='POPC'))
        self.assertTrue(os.path.isdir('fails'))
        self.assertIn('Failures', ''.join(cm.output))

    def test_a_stale_scratch_directory_is_cleared_first(self):
        os.mkdir('tmp')
        open(os.path.join('tmp', 'left-over'), 'w').close()
        self._run(self._args(resname='POPC'))
        self.assertFalse(os.path.exists('tmp'))

    def test_a_missing_topology_file_is_refused(self):
        with self.assertRaises(FileNotFoundError):
            self._run(self._args(resname='POPC', topfile='no_such.rtf'))

    def test_a_supplied_topology_file_is_registered(self):
        open('extra.rtf', 'w').close()
        _do_resi, cc = self._run(self._args(resname='POPC', topfile='extra.rtf'))
        cc.add_topology.assert_called_once_with('extra.rtf')

    def test_an_unknown_residue_exits_rather_than_building_nothing(self):
        with self.assertRaises(SystemExit):
            self._run(self._args(resname='NOPE'), contains=False)

    def test_the_sampler_choice_is_passed_through(self):
        do_resi_mock, _cc = self._run(self._args(resname='POPC', sampler='mc',
                                                 cylinder_apl=None, cylinder_inflation=1.9))
        self.assertEqual(do_resi_mock.call_args.kwargs['sampler'], 'mc')
        self.assertIsNone(do_resi_mock.call_args.kwargs['cylinder_apl'],
                          'None means auto-size from the chain count')


if __name__ == '__main__':
    unittest.main()
