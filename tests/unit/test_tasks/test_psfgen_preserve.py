# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""P2.0 tests for psfgen readpsf-preserve mode: a psfgen after a STATE-provider (continuation)
with no fetched SOURCE edits the incoming topology in place instead of rebuilding from segments."""
import unittest
import os
from pathlib import Path

from pestifer.core.controller import Controller
from pestifer.core.config import Config
from pestifer.core.errors import PestiferBuildError
from pestifer.tasks.psfgen import PsfgenTask
from pestifer.tasks.continuation import ContinuationTask
from pestifer.psfutil.psfcontents import PSFContents
from pestifer.core.artifacts import StateArtifacts


class TestPsfgenPreserveMode(unittest.TestCase):

    def setUp(self):
        self.controller = Controller().configure(Config().configure_new(), terminate=False)
        # runs in the test_psfgen_preserve/ subdir (see conftest change_test_dir); fixtures live one
        # level up under test_tasks/fixtures/
        input_dir = Path('../fixtures/continuation_inputs')
        for ext in ('psf', 'pdb', 'xsc', 'coor', 'vel'):
            dest = Path(f'my_6pti.{ext}')
            if dest.exists() or dest.is_symlink():
                dest.unlink()
            os.symlink((input_dir / f'my_6pti.{ext}').resolve(), dest)

    def _cont(self):
        return dict(psf='my_6pti.psf', pdb='my_6pti.pdb', xsc='my_6pti.xsc',
                    coor='my_6pti.coor', vel='my_6pti.vel', verify_parameters=False)

    def test_preserve_reproduces_the_system(self):
        # continuation (STATE, no SOURCE) -> psfgen must run in preserve mode and carry the topology
        # through unchanged (readpsf pass-through): same atom count in and out.
        self.controller.reconfigure_tasks([{'continuation': self._cont()}, {'psfgen': {}}])
        self.assertIsInstance(self.controller.tasks[0], ContinuationTask)
        self.assertIsInstance(self.controller.tasks[1], PsfgenTask)
        self.controller.do_tasks()
        pg = self.controller.tasks[1]
        self.assertEqual(pg.result, 0)
        state: StateArtifacts = pg.get_current_artifact('state')
        self.assertTrue(state.psf.exists())
        n_in = len(PSFContents('my_6pti.psf').atoms)
        n_out = len(PSFContents(state.psf.name).atoms)
        self.assertEqual(n_in, n_out)
        # the produced state is a psfgen output, distinct from the incoming file
        self.assertNotEqual(state.psf.name, 'my_6pti.psf')

    def test_preserve_applies_patch(self):
        # P2.1: a `patches` mod is applied to the incoming topology (readpsf + patch + regenerate).
        # ASPP protonates ASP A:3 (12 atoms) -> adds one proton (13 atoms); the total grows by one.
        def _asp3_atoms(psf):
            return [a for a in psf.atoms if a.segname == 'A' and str(a.resid.resid) == '3'
                    and a.resname.startswith('ASP')]
        n_before = len(PSFContents('my_6pti.psf').atoms)
        self.assertEqual(len(_asp3_atoms(PSFContents('my_6pti.psf'))), 12)
        self.controller.reconfigure_tasks(
            [{'continuation': self._cont()}, {'psfgen': {'mods': {'patches': ['ASPP:A:3']}}}])
        self.controller.do_tasks()
        pg = self.controller.tasks[1]
        self.assertEqual(pg.result, 0)
        out = PSFContents(pg.get_current_artifact('state').psf.name)
        self.assertEqual(len(out.atoms), n_before + 1)
        self.assertEqual(len(_asp3_atoms(out)), 13)

    def test_preserve_applies_ssbond(self):
        # P2.2: an `ssbonds` mod routes to a `patch DISU` on the loaded project. BPTI's cysteines are
        # all already disulfide-bonded (no free thiol to form a NEW bond), so this verifies routing +
        # acceptance + a clean run: the mod is accepted (not rejected), the DISU patch is emitted, and
        # psfgen completes. (The atom-changing mechanism itself is covered by test_preserve_applies_patch.)
        self.controller.reconfigure_tasks(
            [{'continuation': self._cont()}, {'psfgen': {'mods': {'ssbonds': ['A_5-A_55']}}}])
        self.controller.do_tasks()
        pg = self.controller.tasks[1]
        self.assertEqual(pg.result, 0)
        import glob
        tcl = glob.glob('*psfgen-build.tcl')
        self.assertTrue(tcl, 'psfgen build script not found')
        self.assertIn('patch DISU A:5 A:55', Path(tcl[0]).read_text())

    def test_preserve_applies_crotation(self):
        # P2.3: a coord torsion rotation (irotations) applies on the just-written state via the
        # lazily-built base molecule. Rotating CHI1 of ASP A:3 pivots the sidechain about CA-CB:
        # CB (the rotation axis base) stays put; CG/OD1/OD2 swing out. Topology is unchanged.
        import numpy as np
        from pestifer.molecule.coordmanip import CoordManipulator

        def sidechain(psf, pdb):
            names = ['CB', 'CG', 'OD1', 'OD2']
            cm = CoordManipulator(psf, pdb)
            atoms = PSFContents(psf).atoms
            idx = {a.atomname: i for i, a in enumerate(atoms)
                   if a.segname == 'A' and str(a.resid.resid) == '3' and a.atomname in names}
            return {n: cm.coords[idx[n]] for n in names}

        before = sidechain('my_6pti.psf', 'my_6pti.pdb')
        self.controller.reconfigure_tasks(
            [{'continuation': self._cont()}, {'psfgen': {'mods': {'irotations': ['CHI1,A,3,60.0']}}}])
        self.controller.do_tasks()
        pg = self.controller.tasks[1]
        self.assertEqual(pg.result, 0)
        out_psf = pg.get_current_artifact('state').psf.name
        self.assertEqual(len(PSFContents(out_psf).atoms), len(PSFContents('my_6pti.psf').atoms))
        after = sidechain(out_psf, pg.get_current_artifact('state').pdb.name)
        # CB is on the CHI1 axis -> essentially fixed; the distal carboxylate swings out.
        self.assertLess(float(np.linalg.norm(after['CB'] - before['CB'])), 0.05)
        self.assertGreater(float(np.linalg.norm(after['OD2'] - before['OD2'])), 1.0)

    def test_preserve_applies_link(self):
        # P2.4: a `links` mod's patch is resolved from the residues' geometry (assign_residues ->
        # set_patchname) using the lazily-built molecule, then emitted. On the committed glycoprotein
        # fixture, the ASN A:61 -- glycan V:1304 attachment resolves to NGLA and is emitted as
        # `patch NGLA A:61 V:1304`. (The fixture's glycan is already linked, so this verifies
        # resolution + emission; the atom-changing mechanism is covered by the patches test.)
        gdir = Path('../fixtures/cleave_inputs')
        for ext in ('psf', 'pdb'):
            dest = Path(f'in.{ext}')
            if dest.exists() or dest.is_symlink():
                dest.unlink()
            os.symlink((gdir / f'in.{ext}').resolve(), dest)
        self.controller.reconfigure_tasks(
            [{'continuation': dict(psf='in.psf', pdb='in.pdb', verify_parameters=False)},
             {'psfgen': {'mods': {'links': ['A_61_ND2-V_1304_C1']}}}])
        self.controller.do_tasks()
        pg = self.controller.tasks[1]
        self.assertEqual(pg.result, 0)
        import glob
        tcl = glob.glob('*psfgen-build.tcl')
        self.assertTrue(tcl, 'psfgen build script not found')
        self.assertIn('patch NGLA A:61 V:1304', Path(tcl[0]).read_text())

    def _setup_glycoprotein(self):
        gdir = Path('../fixtures/cleave_inputs')
        for ext in ('psf', 'pdb'):
            dest = Path(f'in.{ext}')
            if dest.exists() or dest.is_symlink():
                dest.unlink()
            os.symlink((gdir / f'in.{ext}').resolve(), dest)
        import shutil
        shutil.copy((Path('../fixtures/graft_inputs/mannose_donor.pdb')).resolve(), 'mannose_donor.pdb')

    def test_preserve_applies_additive_graft(self):
        # P2.5: an additive graft extends a glycan at a *terminal* receiver (no downstream). The donor's
        # source-index sugar (mannose_donor C:1) aligns onto terminal mannose V:1314 and its donor
        # residue (C:2) is carried onto it in a fresh segment, patched across. Result: the glycan grows.
        self._setup_glycoprotein()
        self.controller.reconfigure_tasks(
            [{'continuation': dict(psf='in.psf', pdb='in.pdb', verify_parameters=False)},
             {'psfgen': {'mods': {'grafts': ['V_1314:mannose_donor,C_1-2']}}}])
        self.controller.do_tasks()
        pg = self.controller.tasks[1]
        self.assertEqual(pg.result, 0)
        out = PSFContents(pg.get_current_artifact('state').psf.name)
        base = PSFContents('in.psf')
        grafted = [a for a in out.atoms if a.segname.startswith('GRF')]
        self.assertGreater(len(grafted), 0)                       # a fresh graft segment was added
        self.assertEqual(len(out.atoms), len(base.atoms) + len(grafted))

    def test_preserve_rejects_replace_graft(self):
        # A graft onto an *internal* glycan residue (V:1301, which has downstream sugars) would
        # *replace* that glycan -- removing residues, i.e. re-segmenting -- which readpsf-preserve
        # cannot do. It must hard-error as P3 rather than silently produce a clashing double glycan.
        self._setup_glycoprotein()
        self.controller.reconfigure_tasks(
            [{'continuation': dict(psf='in.psf', pdb='in.pdb', verify_parameters=False)},
             {'psfgen': {'mods': {'grafts': ['V_1301:mannose_donor,C_1-2']}}}])
        with self.assertRaises(PestiferBuildError) as ctx:
            self.controller.do_tasks()
        msg = str(ctx.exception)
        self.assertIn('downstream', msg)
        self.assertIn('P3', msg)

    def test_preserve_rejects_unsupported_mod(self):
        # A mod the preserve path does not yet support (e.g. mutations, which re-segment -- P3) must
        # hard-error rather than be silently ignored, naming it so the user knows what was dropped.
        self.controller.reconfigure_tasks(
            [{'continuation': self._cont()}, {'psfgen': {'mods': {'mutations': ['A:ALA,3,GLY']}}}])
        with self.assertRaises(PestiferBuildError) as ctx:
            self.controller.do_tasks()
        self.assertIn('not yet supported', str(ctx.exception))
        self.assertIn('mutations', str(ctx.exception))


if __name__ == '__main__':
    unittest.main()
