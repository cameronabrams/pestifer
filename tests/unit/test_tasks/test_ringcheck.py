import unittest
from pestifer.psfutil.psfring import ring_check, RingChecker
from pestifer.tasks.ringcheck import RingCheckTask
from pestifer.util.util import cell_from_xsc
import os
import tempfile
import types
import logging
import pytest
import numpy as np
logger=logging.getLogger(__name__)


class TestDeclashScratchCleanup(unittest.TestCase):
    def test_cleanup_removes_only_declash_scratch(self):
        # the rotation-resolution path leaves per-candidate trial PDBs and VMD gen
        # scripts/logs (all prefixed '<taskname>-declash-'); cleanup must remove exactly
        # those and nothing else (not the registered ring_check output, not other files)
        d = tempfile.mkdtemp()
        cwd = os.getcwd()
        os.chdir(d)
        try:
            scratch = [
                'ring_check-declash-G-330-chi1_120.pdb',
                'ring_check-declash-G-330-chi2-gen.tcl',
                'ring_check-declash-G-330-chi2-gen.log',
                'ring_check-declash-I-417-glycan.pdb',
            ]
            keep = ['00-04-00_ring_check.pdb', 'unrelated.dat']
            for n in scratch + keep:
                open(n, 'w').close()
            # call the method with a stand-in self (it only reads self.taskname)
            RingCheckTask._cleanup_declash_scratch(types.SimpleNamespace(taskname='ring_check'))
            self.assertEqual(sorted(os.listdir('.')), sorted(keep))
        finally:
            os.chdir(cwd)

class TestRingCheck(unittest.TestCase):
    def test_ring_check_coords_1(self):
        dir='1'
        pdb=os.path.join(dir,'S2.pdb')
        psf=os.path.join(dir,'S2.psf')
        xsc=os.path.join(dir,'test.xsc')
        result=ring_check(psf,pdb,xsc,cutoff=3.5)
        self.assertEqual(len(result),1)
        p=result[0]
        self.assertEqual(p['piercee']['resid'],67)
        self.assertEqual(p['piercee']['segname'],'S2')
        self.assertEqual(p['piercee']['segtype'],'lipid')
        self.assertEqual(p['piercer']['resid'],645)
        self.assertEqual(p['piercer']['segname'],'S2')
        self.assertEqual(p['piercer']['segtype'],'lipid')

    def test_cutoff_invariance(self):
        # the search cutoff only prunes candidates before the 3.5 A pierced_by gate, so any
        # cutoff >= 3.5 must find exactly the same piercings; the schema default is 4.0
        dir='1'
        pdb=os.path.join(dir,'S2.pdb')
        psf=os.path.join(dir,'S2.psf')
        xsc=os.path.join(dir,'test.xsc')
        base=ring_check(psf,pdb,xsc,cutoff=3.5)
        self.assertEqual(len(base),1)
        key=lambda r:sorted(repr(p) for p in r)
        for cut in (4.0,10.0):
            self.assertEqual(key(ring_check(psf,pdb,xsc,cutoff=cut)),key(base),
                             f'cutoff={cut} changed the detected piercings')

    @pytest.mark.slow
    def test_ring_check_coords_2(self):
        dir='2'
        pdb=os.path.join(dir,'S2.pdb')
        psf=os.path.join(dir,'S2.psf')
        xsc=os.path.join(dir,'test.xsc')
        result=ring_check(psf,pdb,xsc,cutoff=3.5)
        self.assertEqual(len(result),2)
        p=result[1]
        self.assertEqual(p['piercee']['resid'],45)
        self.assertEqual(p['piercee']['segname'],'S2')
        self.assertEqual(p['piercer']['resid'],595)
        self.assertEqual(p['piercer']['segname'],'S2')
        p=result[0]
        self.assertEqual(p['piercee']['resid'],67)
        self.assertEqual(p['piercee']['segname'],'S2')
        self.assertEqual(p['piercer']['resid'],645)
        self.assertEqual(p['piercer']['segname'],'S2')

    @pytest.mark.slow
    def test_ring_check_coords_3(self):
        dir='3'
        pdb=os.path.join(dir,'S2.pdb')
        psf=os.path.join(dir,'S2.psf')
        xsc=os.path.join(dir,'test.xsc')
        result=ring_check(psf,pdb,xsc,cutoff=3.5)
        self.assertEqual(len(result),3)
        p=result[0]
        self.assertEqual(p['piercee']['resid'],67)
        self.assertEqual(p['piercee']['segname'],'S2')
        self.assertEqual(p['piercer']['resid'],645)
        self.assertEqual(p['piercer']['segname'],'S2')
        p=result[1]
        self.assertEqual(p['piercee']['resid'],45)
        self.assertEqual(p['piercee']['segname'],'S2')
        self.assertEqual(p['piercer']['resid'],595)
        self.assertEqual(p['piercer']['segname'],'S2')
        p=result[2]
        self.assertEqual(p['piercee']['resid'],243)
        self.assertEqual(p['piercee']['segname'],'S1')
        self.assertEqual(p['piercer']['resid'],529)
        self.assertEqual(p['piercer']['segname'],'S1')

    def test_ring_check_coords_5(self):
        # synthetic: BGLC pyranose ring (glycan) in XY plane, CHL1 bond (lipid) along Z piercing through center
        dir='5'
        pdb=os.path.join(dir,'test.pdb')
        psf=os.path.join(dir,'test.psf')
        xsc=os.path.join(dir,'test.xsc')
        result=ring_check(psf,pdb,xsc,cutoff=3.5,segtypes=['lipid','glycan'])
        self.assertEqual(len(result),1)
        p=result[0]
        self.assertEqual(p['piercee']['segname'],'AG01')
        self.assertEqual(p['piercee']['resid'],1)
        self.assertEqual(p['piercee']['segtype'],'glycan')
        self.assertEqual(p['piercer']['segname'],'LIP1')
        self.assertEqual(p['piercer']['resid'],2)
        self.assertEqual(p['piercer']['segtype'],'lipid')

    def test_ring_check_vacuum_full_scan(self):
        # whole-system scan with no XSC (vacuum): the box is taken from the coordinate
        # extents.  Guards the else-branch of check() where ll/ur must both be set (a
        # regression: 'ur' was left unbound, raising UnboundLocalError on the first
        # non-periodic ring_check of a build).
        dir='5'
        pdb=os.path.join(dir,'test.pdb')
        psf=os.path.join(dir,'test.psf')
        result=ring_check(psf,pdb,None,cutoff=3.5,segtypes=['lipid','glycan'])
        self.assertEqual(len(result),1)
        self.assertEqual(result[0]['piercee']['segname'],'AG01')
        self.assertEqual(result[0]['piercer']['segname'],'LIP1')

    def test_ring_check_cell_less_xsc(self):
        # an xsc that EXISTS and parses but carries no cell -- a vacuum run writes an
        # origin-only `step o_x o_y o_z`, 4 columns.  check() used to guard on the xsc *path*,
        # so this file entered the periodic branch with box=None and died in np.diagonal(None)
        # ("diag requires an array of at least two dimensions").  It must scan as vacuum
        # instead, giving the same answer as xsc=None.
        dir='5'
        pdb=os.path.join(dir,'test.pdb')
        psf=os.path.join(dir,'test.psf')
        with tempfile.NamedTemporaryFile('w',suffix='.xsc',delete=False) as f:
            f.write('# NAMD extended system configuration output file\n')
            f.write('#$LABELS step o_x o_y o_z\n')
            f.write('1000 0 0 0\n')
            cell_less_xsc=f.name
        try:
            self.assertEqual(cell_from_xsc(cell_less_xsc),(None,None))
            result=ring_check(psf,pdb,cell_less_xsc,cutoff=3.5,segtypes=['lipid','glycan'])
            vacuum=ring_check(psf,pdb,None,cutoff=3.5,segtypes=['lipid','glycan'])
            self.assertEqual(len(result),1)
            self.assertEqual(result,vacuum)
            # the targeted re-check path takes box=None through _check_fast as well
            c=RingChecker(psf,cutoff=3.5,segtypes=['lipid','glycan'])
            self.assertEqual(len(c.check(pdb,xsc=cell_less_xsc,only_piercees=[('AG01',1)])),1)
        finally:
            os.unlink(cell_less_xsc)

    def test_ring_check_resname_and_only_piercees(self):
        # new fields/filters used by the protein side-chain auto-rotate path
        dir='5'
        pdb=os.path.join(dir,'test.pdb')
        psf=os.path.join(dir,'test.psf')
        xsc=os.path.join(dir,'test.xsc')
        result=ring_check(psf,pdb,xsc,cutoff=3.5,segtypes=['lipid','glycan'])
        self.assertEqual(len(result),1)
        # resname is now reported and actually resolved (needed to name the fix, e.g.
        # 'HIS G-330'); it must not fall back to '?' (a resid int/str key mismatch bug)
        self.assertNotEqual(result[0]['piercee'].get('resname'), '?')
        self.assertNotEqual(result[0]['piercer'].get('resname'), '?')
        self.assertTrue(result[0]['piercee']['resname'])
        # only_piercees restricts which rings are checked (targeted re-check after a rotation)
        hit=ring_check(psf,pdb,xsc,cutoff=3.5,segtypes=['lipid','glycan'],only_piercees=[('AG01',1)])
        self.assertEqual(len(hit),1)
        miss=ring_check(psf,pdb,xsc,cutoff=3.5,segtypes=['lipid','glycan'],only_piercees=[('ZZZZ',999)])
        self.assertEqual(len(miss),0)

    def test_ringchecker_reuse_is_consistent(self):
        # the auto-rotate path builds one RingChecker and calls check() per candidate;
        # repeated checks (including a targeted only_piercees check) must not corrupt state
        dir='5'
        pdb=os.path.join(dir,'test.pdb')
        psf=os.path.join(dir,'test.psf')
        xsc=os.path.join(dir,'test.xsc')
        c=RingChecker(psf,cutoff=3.5,segtypes=['lipid','glycan'])
        self.assertEqual(len(c.check(pdb,xsc=xsc)),1)
        self.assertEqual(len(c.check(pdb,xsc=xsc)),1)                       # reuse -> same
        self.assertEqual(len(c.check(pdb,xsc=xsc,only_piercees=[('AG01',1)])),1)
        self.assertEqual(len(c.check(pdb,xsc=xsc,only_piercees=[('ZZZZ',9)])),0)
        self.assertEqual(len(c.check(pdb,xsc=xsc)),1)                       # filter didn't persist

    @pytest.mark.slow
    def test_ring_check_coords_4(self):
        # checks when molecules are in different periodic images
        dir='4'
        pdb=os.path.join(dir,'S2.pdb')
        psf=os.path.join(dir,'S2.psf')
        xsc=os.path.join(dir,'test.xsc')
        result=ring_check(psf,pdb,xsc,cutoff=3.5)
        self.assertEqual(len(result),1)
        p=result[0]
        self.assertEqual(p['piercee']['resid'],907)
        self.assertEqual(p['piercee']['segname'],'S1')
        self.assertEqual(p['piercer']['resid'],354)
        self.assertEqual(p['piercer']['segname'],'S2')

class TestScanSeesAPiercerHangingOutsideTheBox(unittest.TestCase):
    """NAMD's ``wrapAll`` keeps each molecule whole, so a lipid straddling a box face leaves
    bonds hanging outside the box -- by up to a molecule's length.  The whole-system scan
    filed each bond in a link cell by its midpoint, CLIPPED to the grid rather than wrapped
    into it, so an overhanging bond sat in the cell at the face it overhung and no ring beside
    its periodic image ever tested it.

    Found 2026-11-22 on example 17 rep-02's patchB: SOPE L1-191 ``C217-C218`` at y 71-72 A in
    a 65.8 A box threaded CHL1 L1-101's C ring dead center from the raw placement on, through
    both resolving passes, invisible to each; a minimize then re-wrapped the lipid, and the
    verify-only pass -- which may not move anything -- saw it and stopped the build (twice, on
    two commits, from the same seeded placement).  An independent ray-cast found it at every
    stage.  The resolver's targeted re-check (``_check_fast``) was already periodic; only the
    scan was not.

    Here the fixture's one piercer is moved WHOLE by one box length, exactly what wrapAll does
    to a molecule whose center crosses a face.  The physics is unchanged, so the answer must be.
    """

    DIR = '1'
    PIERCER = ('POPE', 645)

    def _shifted(self, axis, sign):
        box, _ = cell_from_xsc(os.path.join(self.DIR, 'test.xsc'))
        L = float(np.diagonal(box)[axis])
        out = []
        with open(os.path.join(self.DIR, 'S2.pdb')) as f:
            for ln in f:
                if ln.startswith(('ATOM', 'HETATM')) and ln[17:21].strip() == self.PIERCER[0] \
                        and int(ln[22:26]) == self.PIERCER[1]:
                    xyz = [float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]
                    xyz[axis] += sign * L
                    ln = ln[:30] + ''.join(f'{v:8.3f}' for v in xyz) + ln[54:]
                out.append(ln)
        fd, path = tempfile.mkstemp(suffix='.pdb', dir='.')
        with os.fdopen(fd, 'w') as f:
            f.writelines(out)
        self.addCleanup(os.unlink, path)
        return path

    def test_the_fixture_has_one_piercing_by_that_molecule(self):
        r = ring_check(os.path.join(self.DIR, 'S2.psf'), os.path.join(self.DIR, 'S2.pdb'),
                       os.path.join(self.DIR, 'test.xsc'), cutoff=4.0, segtypes=['lipid', 'glycan'])
        self.assertEqual([(p['piercer']['resname'], p['piercer']['resid']) for p in r], [self.PIERCER])

    def test_a_piercer_one_box_length_away_is_still_found(self):
        for axis in range(3):
            for sign in (+1, -1):
                with self.subTest(axis='xyz'[axis], sign=sign):
                    r = ring_check(os.path.join(self.DIR, 'S2.psf'), self._shifted(axis, sign),
                                   os.path.join(self.DIR, 'test.xsc'), cutoff=4.0,
                                   segtypes=['lipid', 'glycan'])
                    self.assertEqual([(p['piercer']['resname'], p['piercer']['resid']) for p in r],
                                     [self.PIERCER])


class TestVerifyOnlyPass(unittest.TestCase):
    """``resolve: false`` -- the terminator of the rotate-then-minimize loop.

    A rotation that clears a piercing is a rigid displacement, and it is chosen for *fewest*
    new heavy-atom clashes, not *few*: example 17's patchB accepted poses leaving 5, 2 and 3
    of them on 2026-10-19 (and 24 in an earlier pass), then went straight into NVT and died
    with six RATTLE failures whose atoms were all in or beside a rotated branch.  So the last
    check before dynamics must not be allowed to rotate: by then the loop has already rotated
    and relaxed twice, and anything left is something it cannot fix.
    """

    SPECS = dict(cutoff=4.0, segtypes=['lipid'], delete='none', max_ring_size=7)

    def _stub(self, **overrides):
        """A stand-in ``self`` carrying only what ``do()`` touches, same trick as
        TestDeclashScratchCleanup above."""
        d = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_ringcheck', '1')
        named = lambda f: types.SimpleNamespace(name=os.path.join(d, f))
        state = types.SimpleNamespace(psf=named('S2.psf'), pdb=named('S2.pdb'),
                                      xsc=named('test.xsc'))
        task = types.SimpleNamespace(specs={**self.SPECS, **overrides}, outcome={},
                                     result=None, rotated=False)
        task.get_current_artifact = lambda key: state if key == 'state' else None
        task._fmt = RingCheckTask._fmt
        task._report = lambda p, lvl: RingCheckTask._report(task, p, lvl)
        task.record_outcome = lambda **f: task.outcome.update(f)
        task._rotatable = RingCheckTask._rotatable
        # if do() ever reaches the rotation path, this records it instead of moving atoms
        def _spy(*a, **k):
            task.rotated = True
            return (None, [])          # "nothing moved", so do() breaks out of the loop
        task._resolve_by_rotation = _spy
        return task

    def test_the_fixture_really_is_pierced(self):
        """Positive control: without it, a verify pass that raised for the wrong reason --
        or one that found nothing -- would look the same."""
        task = self._stub()
        d = os.path.dirname(task.get_current_artifact('state').psf.name)
        found = ring_check(os.path.join(d, 'S2.psf'), os.path.join(d, 'S2.pdb'),
                           os.path.join(d, 'test.xsc'), cutoff=4.0, segtypes=['lipid'])
        self.assertEqual(len(found), 1)

    def test_verify_only_raises_and_does_not_rotate(self):
        from pestifer.core.errors import PestiferBuildError
        task = self._stub(resolve=False)
        with self.assertRaises(PestiferBuildError) as cm:
            RingCheckTask.do(task)
        self.assertIn('survived the rotate-and-minimize loop', str(cm.exception))
        self.assertFalse(task.rotated, 'a verify-only pass moved atoms')
        self.assertEqual(task.outcome['rotations'], 0)
        self.assertEqual(task.outcome['piercings_initial'], 1)
        self.assertIs(task.outcome['resolve'], False)

    def test_resolving_pass_does_try_to_rotate(self):
        """The negative control for the test above: with ``resolve`` left at its default the
        same fixture must reach the rotation path, or ``test_verify_only_raises_and_does_not
        _rotate`` would pass against a task that never rotates at all."""
        from pestifer.core.errors import PestiferBuildError
        task = self._stub()
        # PestiferBuildError specifically, not Exception: a bare `assertRaises(Exception)`
        # here passed against an AttributeError from an incomplete stub, which is this
        # repository's "a check that passes is not evidence the check ran" in one line.
        with self.assertRaises(PestiferBuildError) as cm:
            RingCheckTask.do(task)      # delete='none' and nothing rotated -> still fatal
        self.assertIn('delete is none', str(cm.exception))
        self.assertTrue(task.rotated, 'a resolving pass never attempted a rotation')
