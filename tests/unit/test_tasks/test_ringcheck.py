import unittest
from pestifer.psfutil.psfring import ring_check, RingChecker
from pestifer.tasks.ringcheck import RingCheckTask
import os
import tempfile
import types
import logging
import pytest
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