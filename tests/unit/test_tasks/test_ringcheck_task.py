# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Unit tests for the ``ring_check`` task's decision logic.

``test_ringcheck.py`` covers the detector -- whether a bond actually threads a ring -- against real
coordinates.  This file covers what the *task* does once it has a list of piercings: which ones it
tries to rotate free, which it deletes, which it refuses to proceed past, and in what order it
tries the rotation strategies.

Those decisions are consequential and none of them is visible in the output when they go wrong.
Deleting the wrong residue silently changes the system's composition; classifying a fixable
piercing as fatal aborts a build that had a way forward; and classifying a fatal one as fixable
hands the next stage a threaded structure that RATTLE-fails minutes later with an unrelated-looking
error.

No VMD or NAMD: the detector and the rotation trials are stubbed, and what is asserted is the
routing and the delete set.
"""

import os
import tempfile
import unittest
from unittest import mock

import numpy as np

from pestifer.core.errors import PestiferBuildError
from pestifer.tasks.ringcheck import RingCheckTask


def _side(segtype, segname='A', resid=42, resname=None, **extra):
    d = {'segtype': segtype, 'segname': segname, 'resid': resid}
    if resname is not None:
        d['resname'] = resname
    d.update(extra)
    return d


def _piercing(piercee, piercer):
    return {'piercee': piercee, 'piercer': piercer}


def _task(**specs):
    t = RingCheckTask.__new__(RingCheckTask)
    t.taskname = 'ring_check'
    t.basename = 'test'
    t.specs = dict(specs)
    return t


# --- reporting -------------------------------------------------------------------------------------

class TestFormatting(unittest.TestCase):
    """Every piercing message names both partners; a report that says only "a ring was pierced"
    is not actionable."""

    def test_a_side_is_named_by_residue_and_position(self):
        self.assertEqual(RingCheckTask._fmt(_side('protein', 'A', 42, 'TYR')), 'TYR A-42')

    def test_a_missing_resname_does_not_break_the_message(self):
        self.assertEqual(RingCheckTask._fmt(_side('lipid', 'M', 7)), '? M-7')

    def test_the_report_names_both_partners_and_the_piercer_segtype(self):
        t = _task()
        lines = []
        t._report([_piercing(_side('protein', 'A', 42, 'TYR'),
                             _side('glycan', 'G1', 3, 'BGLC'))], lines.append)
        self.assertEqual(len(lines), 1)
        self.assertIn('TYR A-42', lines[0])
        self.assertIn('BGLC G1-3', lines[0])
        self.assertIn('glycan', lines[0])

    def test_every_piercing_is_reported(self):
        t = _task()
        lines = []
        t._report([_piercing(_side('lipid'), _side('lipid'))] * 3, lines.append)
        self.assertEqual(len(lines), 3)


# --- what can be rotated free ------------------------------------------------------------------------

class TestRotatable(unittest.TestCase):
    """Which piercings have a rotation worth trying.

    Getting this wrong in either direction is costly: a false negative aborts a build that could
    have been fixed, and a false positive spends the rotation search on something that cannot move
    and then reports failure anyway.
    """

    def test_an_aromatic_side_chain_ring_can_be_rotated(self):
        for resname in ('PHE', 'TYR', 'TRP', 'HIS', 'HSD', 'HSE', 'HSP'):
            with self.subTest(resname=resname):
                p = _piercing(_side('protein', resname=resname), _side('lipid'))
                self.assertTrue(RingCheckTask._rotatable(p))

    def test_a_non_aromatic_protein_ring_speared_by_a_lipid_is_not_rotatable(self):
        """Proline's ring is rigid and a lipid bond carries no glycan pendant to swing."""
        p = _piercing(_side('protein', resname='PRO'), _side('lipid'))
        self.assertFalse(RingCheckTask._rotatable(p))

    def test_a_rigid_ring_speared_by_a_glycan_bond_can_be_rotated(self):
        """The ring cannot move, but the glycan pendant carrying the offending bond can."""
        p = _piercing(_side('protein', resname='PRO'),
                      _side('glycan', bond_serials=[10, 11]))
        self.assertTrue(RingCheckTask._rotatable(p))

    def test_a_glycan_piercer_without_a_bond_offers_nothing_to_rotate(self):
        p = _piercing(_side('protein', resname='PRO'), _side('glycan'))
        self.assertFalse(RingCheckTask._rotatable(p))

    def test_a_glycan_ring_speared_by_a_protein_side_chain_is_rotatable_either_way(self):
        p = _piercing(_side('glycan', ring_serials=[1, 2, 3]), _side('protein'))
        self.assertTrue(RingCheckTask._rotatable(p))
        p = _piercing(_side('glycan'), _side('protein', bond_serials=[4, 5]))
        self.assertTrue(RingCheckTask._rotatable(p))

    def test_a_lipid_ring_is_not_rotated(self):
        """A sterol ring speared by a tail is resolved by deletion, not rotation."""
        p = _piercing(_side('lipid', resname='CHL1'), _side('lipid'))
        self.assertFalse(RingCheckTask._rotatable(p))


# --- the task's decision tree --------------------------------------------------------------------------

class _DoCase(unittest.TestCase):
    """Drives ``do()`` with the detector and the rotation search stubbed."""

    def _run(self, piercings, *, after_rotation=None, rotation_fails=None, **specs):
        specs.setdefault('segtypes', ['lipid', 'glycan', 'protein'])
        t = _task(**specs)
        state = mock.Mock()
        state.psf.name, state.pdb.name = 'in.psf', 'in.pdb'
        state.xsc = None
        t.get_current_artifact = lambda k: state
        t.next_basename = mock.Mock()
        t.register = mock.Mock()
        t._cleanup_declash_scratch = mock.Mock()
        self.scripter = mock.Mock()
        t.get_scripter = mock.Mock(return_value=self.scripter)

        if rotation_fails is not None:
            t._resolve_by_rotation = mock.Mock(return_value=(None, rotation_fails))
        else:
            t._resolve_by_rotation = mock.Mock(return_value=('fixed.pdb', None))

        # first call returns the initial list; a second (post-rotation) call returns whatever
        # the rotation left behind
        rounds = [list(piercings)]
        if after_rotation is not None:
            rounds.append(list(after_rotation))
        self.checks = []

        def fake_ring_check(psf, pdb, xsc, **kw):
            self.checks.append(kw)
            return rounds.pop(0) if rounds else []

        # assigned before do() so a test that expects an abort can still inspect the task
        self.task = t
        with mock.patch('pestifer.tasks.ringcheck.ring_check', side_effect=fake_ring_check), \
             mock.patch('pestifer.tasks.ringcheck.shutil.copy'):
            self.result = t.do()
        return t

    def _delatom_lines(self):
        return [c.args[0] for c in self.scripter.addline.call_args_list
                if c.args and c.args[0].startswith('delatom')]


class TestNothingToDo(_DoCase):

    def test_a_clean_structure_passes(self):
        self._run([])
        self.assertEqual(self.result, 0)
        self.task.register.assert_not_called()

    def test_no_segtypes_checks_nothing_and_says_so(self):
        """``segtypes`` has no default, so an unset one would otherwise scan nothing in silence."""
        with self.assertLogs('pestifer.tasks.ringcheck', level='WARNING') as cm:
            t = self._run([_piercing(_side('lipid'), _side('lipid'))], segtypes=[])
        self.assertIn('no segtypes specified', ''.join(cm.output))
        self.assertEqual(self.result, 0)
        self.assertEqual(self.checks, [], 'the detector should not have run at all')


class TestRotationPath(_DoCase):

    AROMATIC = _piercing(_side('protein', 'A', 42, 'TYR'), _side('lipid', 'M', 7))

    def test_a_rotatable_piercing_is_sent_to_the_rotation_search(self):
        t = self._run([self.AROMATIC], after_rotation=[])
        t._resolve_by_rotation.assert_called_once()
        self.assertEqual(self.result, 0)

    def test_a_successful_rotation_keeps_the_psf_and_registers_a_new_pdb(self):
        """A rotation moves coordinates only; re-running psfgen would be both wasteful and a
        chance to change the topology."""
        t = self._run([self.AROMATIC], after_rotation=[])
        t.register.assert_called_once()
        registered = t.register.call_args.args[0]
        self.assertIs(registered['psf'], t.get_current_artifact('state').psf)

    def test_the_structure_is_re_checked_after_rotating(self):
        """The rotation is only believed if the detector agrees afterwards."""
        self._run([self.AROMATIC], after_rotation=[])
        self.assertEqual(len(self.checks), 2)

    def test_scratch_from_the_rotation_search_is_cleaned_up_on_success(self):
        t = self._run([self.AROMATIC], after_rotation=[])
        t._cleanup_declash_scratch.assert_called_once()

    def test_an_unresolvable_rotation_aborts_with_both_partners_named(self):
        with self.assertRaises(PestiferBuildError) as cm:
            self._run([self.AROMATIC], rotation_fails=self.AROMATIC)
        msg = str(cm.exception)
        self.assertIn('TYR A-42', msg)
        self.assertIn('? M-7', msg)
        self.assertIn('rotate', msg, 'the error should suggest a fix')

    def test_scratch_is_kept_when_the_rotation_fails(self):
        """The trial PDBs are the evidence for debugging why nothing cleared it."""
        with self.assertRaises(PestiferBuildError):
            self._run([self.AROMATIC], rotation_fails=self.AROMATIC)
        self.task._cleanup_declash_scratch.assert_not_called()


class TestFatalPiercings(_DoCase):
    """A protein or glycan ring that no rotation cleared cannot be resolved by deletion -- you
    cannot drop a residue out of a chain -- so the build stops."""

    def test_an_unrotatable_protein_ring_is_fatal(self):
        p = _piercing(_side('protein', resname='PRO'), _side('lipid'))
        with self.assertRaises(PestiferBuildError) as cm:
            self._run([p])
        self.assertIn('protein/glycan', str(cm.exception))

    def test_an_unrotatable_glycan_ring_is_fatal(self):
        p = _piercing(_side('glycan'), _side('lipid'))
        with self.assertRaises(PestiferBuildError):
            self._run([p])

    def test_a_piercing_that_survives_rotation_is_fatal(self):
        rot = _piercing(_side('protein', resname='TYR'), _side('lipid'))
        still = _piercing(_side('protein', resname='PRO'), _side('lipid'))
        with self.assertRaises(PestiferBuildError):
            self._run([rot], after_rotation=[still])

    def test_the_error_suggests_a_rebuild(self):
        p = _piercing(_side('glycan'), _side('lipid'))
        with self.assertRaises(PestiferBuildError) as cm:
            self._run([p])
        self.assertIn('rebuild', str(cm.exception))


class TestLipidDeletion(_DoCase):
    """A lipid ring speared by a tail is resolved by deleting a lipid, which changes the system's
    composition -- so which residues go is worth pinning exactly."""

    PIERCING = _piercing(_side('lipid', 'MEMB', 10, 'CHL1'),
                         _side('lipid', 'MEMB', 20, 'POPC'))

    def test_by_default_the_pierced_lipid_is_deleted(self):
        self._run([self.PIERCING])
        self.assertEqual(self._delatom_lines(), ['delatom MEMB 10'])

    def test_the_piercer_can_be_deleted_instead(self):
        self._run([self.PIERCING], delete='piercer')
        self.assertEqual(self._delatom_lines(), ['delatom MEMB 20'])

    def test_both_removes_the_pair(self):
        """A tail threaded through a sterol ring contorts *both* lipids, so removing only one can
        leave the other in a pose that RATTLE-fails at the first dynamics step."""
        self._run([self.PIERCING], delete='both')
        self.assertEqual(sorted(self._delatom_lines()),
                         ['delatom MEMB 10', 'delatom MEMB 20'])

    def test_both_never_deletes_a_non_lipid_piercer(self):
        """A glycan or protein bond cannot be dropped without breaking the molecule."""
        p = _piercing(_side('lipid', 'MEMB', 10), _side('protein', 'A', 5))
        self._run([p], delete='both')
        self.assertEqual(self._delatom_lines(), ['delatom MEMB 10'])

    def test_delete_none_refuses_to_proceed(self):
        with self.assertRaises(PestiferBuildError) as cm:
            self._run([self.PIERCING], delete='none')
        self.assertIn('delete=none', str(cm.exception))

    def test_a_residue_named_twice_is_deleted_once(self):
        two = [self.PIERCING, self.PIERCING]
        self._run(two)
        self.assertEqual(self._delatom_lines(), ['delatom MEMB 10'])

    def test_the_deletion_is_reported_at_info(self):
        """A silent composition change is misleading; the rotation path already reports."""
        with self.assertLogs('pestifer.tasks.ringcheck', level='INFO') as cm:
            self._run([self.PIERCING])
        self.assertIn('Deleting', ''.join(cm.output))

    def test_deleting_rewrites_the_topology(self):
        """Unlike a rotation, a deletion changes the PSF, so both files are re-registered."""
        t = self._run([self.PIERCING])
        registered = t.register.call_args_list[0].args[0]
        self.assertIn('psf', registered)
        self.assertIn('pdb', registered)
        self.scripter.runscript.assert_called_once()


class TestDetectorArguments(_DoCase):

    def test_the_configured_cutoff_and_ring_size_reach_the_detector(self):
        self._run([], cutoff=4.0, max_ring_size=6)
        self.assertEqual(self.checks[0]['cutoff'], 4.0)
        self.assertEqual(self.checks[0]['max_ring_size'], 6)

    def test_the_defaults_are_the_documented_ones(self):
        self._run([])
        self.assertEqual(self.checks[0]['cutoff'], 3.5)
        self.assertEqual(self.checks[0]['max_ring_size'], 7)

    def test_the_requested_segtypes_reach_the_detector(self):
        self._run([], segtypes=['lipid'])
        self.assertEqual(self.checks[0]['segtypes'], ['lipid'])


# --- side-chain rotation scoring -----------------------------------------------------------------------

class TestTrySidechain(unittest.TestCase):
    """Choosing among side-chain rotations that clear the ring.

    Any rotation that un-threads the ring is a candidate, but they are not equivalent: swinging a
    tyrosine free of one ring can bury it in something else.  So candidates are scored by the new
    heavy-atom clashes they introduce and the cheapest wins, with a clash-free one ending the
    search immediately -- the trials are the expensive part.
    """

    def setUp(self):
        self._cwd = os.getcwd()
        self._tmp = tempfile.TemporaryDirectory()
        os.chdir(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)
        self.addCleanup(os.chdir, self._cwd)

    def _try(self, *, clears=(), clashes=None, written=None):
        """``clears`` names the candidate files that un-thread the ring; ``clashes`` maps a
        candidate to its new-clash count."""
        clashes = clashes or {}
        t = _task()
        self.written = []

        def fake_write(psf, pdb, seg, resid, chi, degrees, prefix):
            self.written.append(chi)
            for deg in degrees:
                name = f'{prefix}_{deg}.pdb'
                if written is None or name in written:
                    open(name, 'w').close()

        t._write_sidechain_candidates = fake_write

        checker = mock.Mock()
        self.loaded = []

        def load_coords(path):
            self.loaded.append(path)
            return np.zeros((3, 3))

        checker.load_coords.side_effect = load_coords
        # check_coords returns truthy while still pierced
        checker.check_coords.side_effect = lambda coords, box, targets: (
            [] if self.loaded[-1] in clears else ['still pierced'])
        checker.clash_count.side_effect = lambda coords, moved: clashes.get(self.loaded[-1], 0)

        return t._try_sidechain(checker, 'in.psf', 'in.pdb', np.ones((3, 3)), None, None,
                                'A', 42, 'TYR', ('A', 42), 'pfx')

    def test_no_rotation_clearing_the_ring_yields_nothing(self):
        self.assertIsNone(self._try(clears=()))

    def test_chi2_is_tried_before_chi1(self):
        """chi2 swings the ring itself; chi1 moves the whole side chain and disturbs more."""
        self._try(clears=())
        self.assertEqual(self.written, [2, 1])

    def test_a_clearing_rotation_is_returned(self):
        got = self._try(clears={'pfx-chi2_180.pdb'})
        self.assertEqual(got, 'pfx-chi2_180.pdb')

    def test_a_clash_free_candidate_ends_the_search(self):
        self._try(clears={'pfx-chi2_180.pdb'}, clashes={'pfx-chi2_180.pdb': 0})
        self.assertEqual(self.written, [2], 'chi1 should not have been attempted')
        self.assertEqual(len(self.loaded), 1, 'no further candidates should be loaded')

    def test_the_candidate_with_the_fewest_new_clashes_wins(self):
        got = self._try(clears={'pfx-chi2_180.pdb', 'pfx-chi2_120.pdb', 'pfx-chi2_240.pdb'},
                        clashes={'pfx-chi2_180.pdb': 5, 'pfx-chi2_120.pdb': 2,
                                 'pfx-chi2_240.pdb': 9})
        self.assertEqual(got, 'pfx-chi2_120.pdb')

    def test_the_search_continues_into_chi1_when_chi2_only_offers_clashy_poses(self):
        got = self._try(clears={'pfx-chi2_180.pdb', 'pfx-chi1_120.pdb'},
                        clashes={'pfx-chi2_180.pdb': 4, 'pfx-chi1_120.pdb': 0})
        self.assertEqual(self.written, [2, 1])
        self.assertEqual(got, 'pfx-chi1_120.pdb')

    def test_a_candidate_that_was_not_written_is_skipped(self):
        """A rotation the writer could not produce is absent, not empty."""
        got = self._try(clears={'pfx-chi2_120.pdb'}, written={'pfx-chi2_120.pdb'})
        self.assertEqual(got, 'pfx-chi2_120.pdb')
        self.assertEqual(self.loaded, ['pfx-chi2_120.pdb'])

    def test_the_accepted_rotation_is_announced(self):
        with self.assertLogs('pestifer.tasks.ringcheck', level='INFO') as cm:
            self._try(clears={'pfx-chi2_180.pdb'}, clashes={'pfx-chi2_180.pdb': 0})
        out = ''.join(cm.output)
        self.assertIn('TYR A-42', out)
        self.assertIn('chi2', out)


class TestSidechainSerials(unittest.TestCase):
    """Which atoms move when a side chain is rotated: everything that is not backbone."""

    @staticmethod
    def _atom(segname, resid, atomname):
        a = mock.Mock()
        a.segname, a.atomname = segname, atomname
        a.serial = hash((segname, resid, atomname)) % 10000
        a.resid = mock.Mock(resid=resid)
        return a

    def _serials(self, atoms, segname='A', resid=42):
        checker = mock.Mock()
        checker._atoms = atoms
        return RingCheckTask._sidechain_serials(checker, segname, resid)

    def test_backbone_atoms_are_excluded(self):
        atoms = [self._atom('A', 42, n) for n in ('N', 'CA', 'C', 'O', 'CB', 'CG')]
        got = self._serials(atoms)
        self.assertEqual(len(got), 2, 'only CB and CG are side chain')

    def test_other_residues_are_excluded(self):
        atoms = [self._atom('A', 42, 'CB'), self._atom('A', 43, 'CB')]
        self.assertEqual(len(self._serials(atoms)), 1)

    def test_other_segments_are_excluded(self):
        atoms = [self._atom('A', 42, 'CB'), self._atom('B', 42, 'CB')]
        self.assertEqual(len(self._serials(atoms)), 1)

    def test_the_resid_comparison_survives_a_type_mismatch(self):
        """Resids arrive as ints from one source and strings from another; comparing them
        directly would quietly select nothing."""
        atoms = [self._atom('A', '42', 'CB')]
        self.assertEqual(len(self._serials(atoms, resid=42)), 1)


class TestGlycanPendantDelegation(unittest.TestCase):
    """The common case -- a rigid ring speared by a glycan bond -- hands the branch carrying that
    bond straight to the shared engine, the same one used at glycan graft time."""

    def test_the_piercing_bond_branch_is_handed_to_the_engine(self):
        t = _task()
        piercer = _side('glycan', bond_serials=[10, 11])
        with mock.patch('pestifer.tasks.ringcheck.try_glycan_pendant',
                        return_value='cand.pdb') as engine:
            got = t._try_glycan_pendant(mock.Mock(), 'in.pdb', None, None, piercer,
                                        ('A', 42), 'PRO A-42', 'pfx')
        self.assertEqual(got, 'cand.pdb')
        self.assertIs(engine.call_args.args[4], piercer)
        self.assertIn('pfx-glycan.pdb', engine.call_args.args)


class TestCleanupScratchIsBestEffort(unittest.TestCase):
    """The accepted pose is already registered, so a file that cannot be removed is untidy rather
    than a problem -- it must not take the build down after the piercing was successfully fixed."""

    def test_an_unremovable_file_does_not_raise(self):
        t = _task()
        t.taskname = 'ring_check'
        with mock.patch('pestifer.tasks.ringcheck.glob.glob',
                        return_value=['ring_check-declash-A-42-chi1_120.pdb']), \
             mock.patch('pestifer.tasks.ringcheck.os.remove',
                        side_effect=OSError('permission denied')):
            t._cleanup_declash_scratch()          # must not raise


class TestPiercedGlycanAndPiercerSidechain(unittest.TestCase):
    """The two rotations that move something other than the obvious partner.

    Both delegate to the shared glycan-pendant engine, which moves whichever graph branch the
    serials it is handed belong to.  What each of these decides is *which* serials to hand it.
    """

    def _pierced_glycan(self, piercee):
        t = _task()
        with mock.patch('pestifer.tasks.ringcheck.try_glycan_pendant',
                        return_value='cand.pdb') as engine:
            got = t._try_pierced_glycan(mock.Mock(), 'in.pdb', None, None, piercee,
                                        ('G1', 3), 'pfx')
        return got, engine

    def test_a_ring_needs_two_atoms_to_identify_its_branch(self):
        got, engine = self._pierced_glycan(_side('glycan', ring_serials=[7]))
        self.assertIsNone(got)
        engine.assert_not_called()

    def test_no_ring_serials_means_nothing_to_swing(self):
        got, engine = self._pierced_glycan(_side('glycan'))
        self.assertIsNone(got)
        engine.assert_not_called()

    def test_two_ring_atoms_stand_in_for_a_bonded_pair(self):
        """A ring is 2-connected, so no bridge separates its atoms: any two of them identify the
        branch carrying the ring exactly as a bonded pair does."""
        got, engine = self._pierced_glycan(_side('glycan', ring_serials=[7, 8, 9, 10, 11, 12]))
        self.assertEqual(got, 'cand.pdb')
        self.assertEqual(engine.call_args.args[4]['bond_serials'], [7, 8])

    def test_the_mover_kind_is_reported_to_the_engine(self):
        _got, engine = self._pierced_glycan(_side('glycan', ring_serials=[1, 2]))
        self.assertEqual(engine.call_args.kwargs['mover_kind'], 'pierced glycan')

    def _piercer_sidechain(self, piercer, movers=(101, 102)):
        t = _task()
        t._sidechain_serials = mock.Mock(return_value=list(movers))
        with mock.patch('pestifer.tasks.ringcheck.try_glycan_pendant',
                        return_value='cand.pdb') as engine:
            got = t._try_piercer_sidechain(mock.Mock(), 'in.pdb', None, None, piercer,
                                           ('G1', 3), 'label', 'pfx')
        return got, engine

    def test_a_piercer_with_no_bond_offers_no_hinge(self):
        got, engine = self._piercer_sidechain(_side('protein'))
        self.assertIsNone(got)
        engine.assert_not_called()

    def test_a_residue_with_no_side_chain_cannot_be_rotated(self):
        """Glycine: nothing but backbone, so there is no branch to swing."""
        got, engine = self._piercer_sidechain(_side('protein', bond_serials=[1, 2]), movers=())
        self.assertIsNone(got)
        engine.assert_not_called()

    def test_the_movers_are_restricted_to_the_side_chain(self):
        """Restricting the eligible hinges to bonds whose moving branch stays inside the side
        chain leaves exactly that residue's chi bonds -- no per-residue chi table needed."""
        _got, engine = self._piercer_sidechain(_side('protein', bond_serials=[1, 2]))
        self.assertEqual(engine.call_args.kwargs['mover_serials'], [101, 102])
        self.assertEqual(engine.call_args.kwargs['mover_kind'], 'side chain')


class TestWriteSidechainCandidates(unittest.TestCase):
    """One PDB per candidate rotation, each measured from the *input* pose."""

    def _write(self, degrees=(120, 240, 180)):
        t = _task()
        cm = mock.Mock()
        cm.coords = np.zeros((3, 3))
        cm.residue_of_segname.return_value = 'RES'
        applied, written, coords_at_apply = [], [], []

        def apply_scrot(chi, rn, deg):
            applied.append((chi, rn, deg))
            coords_at_apply.append(cm.coords.copy())
            cm.coords = cm.coords + deg          # a visible, cumulative-if-wrong change

        cm.apply_scrot.side_effect = apply_scrot
        cm.write_pdb.side_effect = written.append
        with mock.patch('pestifer.molecule.coordmanip.CoordManipulator', return_value=cm):
            t._write_sidechain_candidates('in.psf', 'in.pdb', 'A', 42, 2, list(degrees), 'pfx')
        return applied, written, coords_at_apply

    def test_one_file_per_requested_rotation(self):
        _a, written, _c = self._write(degrees=(120, 240, 180))
        self.assertEqual(written, ['pfx_120.pdb', 'pfx_240.pdb', 'pfx_180.pdb'])

    def test_every_rotation_is_applied_about_the_requested_chi(self):
        applied, _w, _c = self._write(degrees=(120, 240))
        self.assertEqual([(chi, deg) for chi, _rn, deg in applied], [(2, 120), (2, 240)])

    def test_each_rotation_starts_from_the_original_pose(self):
        """If the coordinates were not reset between candidates the rotations would compound,
        so 'rotate by 240' would silently mean 'rotate by 360'."""
        _a, _w, coords_at_apply = self._write(degrees=(120, 240, 180))
        for c in coords_at_apply:
            self.assertTrue((c == 0).all(), 'a candidate was built on top of the previous one')


# --- rotation strategy order ---------------------------------------------------------------------------

class TestResolveByRotationStrategy(unittest.TestCase):
    """Which rotation is tried, and in which order.

    The ordering is a deliberate choice worth pinning: for a glycan ring speared by a protein side
    chain, the *glycan* is moved first because it is the model-built partner, where the protein
    side chain generally is not.  Reversing that would disturb experimentally-determined
    coordinates in preference to modeled ones.
    """

    def _resolve(self, piercing, wins=None):
        t = _task()
        state = mock.Mock()
        state.psf.name, state.pdb.name = 'in.psf', 'in.pdb'
        state.xsc = None
        tried = []

        def strategy(name):
            def f(*a, **kw):
                tried.append(name)
                return 'cand.pdb' if wins == name else None
            return f

        for name in ('_try_sidechain', '_try_glycan_pendant', '_try_pierced_glycan',
                     '_try_piercer_sidechain'):
            setattr(t, name, strategy(name))

        with mock.patch('pestifer.tasks.ringcheck.RingChecker'):
            fixed, failed = t._resolve_by_rotation(state, [piercing], 3.5, ['protein'], 7)
        return tried, fixed, failed

    def test_an_aromatic_ring_rotates_its_own_side_chain(self):
        p = _piercing(_side('protein', resname='TYR'), _side('lipid'))
        tried, fixed, failed = self._resolve(p, wins='_try_sidechain')
        self.assertEqual(tried, ['_try_sidechain'])
        self.assertEqual(fixed, 'cand.pdb')
        self.assertIsNone(failed)

    def test_a_rigid_ring_speared_by_a_glycan_rotates_the_pendant(self):
        p = _piercing(_side('protein', resname='PRO'), _side('glycan', bond_serials=[1, 2]))
        tried, _f, _x = self._resolve(p, wins='_try_glycan_pendant')
        self.assertEqual(tried, ['_try_glycan_pendant'])

    def test_an_aromatic_ring_falls_back_to_the_glycan_pendant(self):
        p = _piercing(_side('protein', resname='TYR'), _side('glycan', bond_serials=[1, 2]))
        tried, _f, _x = self._resolve(p, wins='_try_glycan_pendant')
        self.assertEqual(tried, ['_try_sidechain', '_try_glycan_pendant'])

    def test_a_pierced_glycan_moves_the_glycan_before_the_protein(self):
        p = _piercing(_side('glycan'), _side('protein'))
        tried, _f, _x = self._resolve(p, wins='_try_piercer_sidechain')
        self.assertEqual(tried, ['_try_pierced_glycan', '_try_piercer_sidechain'],
                         'the model-built glycan should be disturbed before the protein')

    def test_the_protein_is_left_alone_when_the_glycan_clears_it(self):
        p = _piercing(_side('glycan'), _side('protein'))
        tried, _f, _x = self._resolve(p, wins='_try_pierced_glycan')
        self.assertEqual(tried, ['_try_pierced_glycan'])

    def test_the_first_piercing_that_cannot_be_cleared_is_returned(self):
        p = _piercing(_side('protein', resname='TYR'), _side('lipid'))
        _tried, fixed, failed = self._resolve(p, wins=None)
        self.assertIsNone(fixed)
        self.assertIs(failed, p)


if __name__ == '__main__':
    unittest.main()
