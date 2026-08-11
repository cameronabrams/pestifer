# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""P3 tests for re-segmenting edits on an incoming PSF: a mutation/deletion/insertion routes psfgen
to build-mode re-segmentation (reconstruct a Molecule from the incoming PSF+coords, apply the edit,
rebuild), instead of the readpsf-preserve path."""
import unittest
import os

import pytest
from pathlib import Path

from pestifer.core.controller import Controller
from pestifer.core.config import Config
from pestifer.core.errors import PestiferBuildError
from pestifer.tasks.psfgen import PsfgenTask
from pestifer.psfutil.psfcontents import PSFContents
from pestifer.core.artifacts import StateArtifacts

# drives external VMD/psfgen; skipped when those binaries are not on PATH
pytestmark = pytest.mark.needs_tools


class TestPsfgenResegment(unittest.TestCase):

    def setUp(self):
        self.controller = Controller().configure(Config().configure_new(), terminate=False)
        input_dir = Path('../fixtures/continuation_inputs')
        for ext in ('psf', 'pdb', 'xsc', 'coor', 'vel'):
            dest = Path(f'my_6pti.{ext}')
            if dest.exists() or dest.is_symlink():
                dest.unlink()
            os.symlink((input_dir / f'my_6pti.{ext}').resolve(), dest)

    def _cont(self, psf='my_6pti.psf'):
        return dict(psf=psf, pdb='my_6pti.pdb', xsc='my_6pti.xsc', coor='my_6pti.coor',
                    vel='my_6pti.vel', verify_parameters=False)

    def test_mutation_routes_to_rebuild(self):
        # A mutation on an incoming PSF re-segments (build-mode): the residue changes, everything else
        # is preserved, the disulfides are re-derived from the PSF, and the box (xsc) carries forward.
        base = PSFContents('my_6pti.psf')
        self.controller.reconfigure_tasks(
            [{'continuation': self._cont()}, {'psfgen': {'mods': {'mutations': ['A:ASN,24,ALA']}}}])
        self.assertIsInstance(self.controller.tasks[1], PsfgenTask)
        self.controller.do_tasks()
        pg = self.controller.tasks[1]
        self.assertEqual(pg.result, 0)
        state: StateArtifacts = pg.get_current_artifact('state')
        out = PSFContents(state.psf.name)
        self.assertEqual({a.resname for a in out.atoms if a.segname == 'A' and str(a.resid.resid) == '24'},
                         {'ALA'})
        self.assertLess(len(out.atoms), len(base.atoms))                       # ASN->ALA drops sidechain
        # disulfides re-derived: all 6 CYS remain bonded (no free HG1 thiol hydrogen)
        self.assertEqual(sum(1 for a in out.atoms if a.resname == 'CYS' and a.atomname == 'HG1'), 0)
        self.assertTrue(getattr(state, 'xsc', None) and state.xsc.name == 'my_6pti.xsc')  # box carried

    def test_deletion_routes_to_rebuild(self):
        # Deleting the C-terminal two residues re-segments chain A two residues shorter.
        def n_A(psf):
            return len({(a.segname, str(a.resid.resid)) for a in psf.atoms if a.segname == 'A'})
        before = n_A(PSFContents('my_6pti.psf'))
        self.controller.reconfigure_tasks(
            [{'continuation': self._cont()}, {'psfgen': {'mods': {'deletions': ['A:56-57']}}}])
        self.controller.do_tasks()
        pg = self.controller.tasks[1]
        self.assertEqual(pg.result, 0)
        out = PSFContents(pg.get_current_artifact('state').psf.name)
        self.assertEqual(n_A(out), before - 2)
        self.assertFalse(any(a.segname == 'A' and str(a.resid.resid) in ('56', '57') for a in out.atoms))

    def test_rebuildability_guard_names_unbuildable_residue(self):
        # A residue psfgen has no RESI for cannot be rebuilt; a re-segmenting edit must hard-error up
        # front, naming it, rather than silently drop it. Rename one residue to a bogus name in the PSF.
        lines = Path('my_6pti.psf').read_text().splitlines()
        seen, patched, out = False, False, []
        for line in lines:
            if '!NATOM' in line:
                seen = True
                out.append(line)
                continue
            toks = line.split()
            if seen and not patched and len(toks) == 9 and toks[0].isdigit() and toks[3] == 'ASN':
                line = line.replace(' ASN ', ' ZZZ ', 1)
                patched = True
            out.append(line)
        self.assertTrue(patched)
        bad = Path('my_6pti_bad.psf')
        bad.write_text('\n'.join(out) + '\n')
        self.controller.reconfigure_tasks(
            [{'continuation': self._cont(psf=str(bad))},
             {'psfgen': {'mods': {'mutations': ['A:ALA,10,GLY']}}}])
        with self.assertRaises(PestiferBuildError) as ctx:
            self.controller.do_tasks()
        self.assertIn('ZZZ', str(ctx.exception))
        self.assertIn('RESI', str(ctx.exception))
        bad.unlink()


if __name__ == '__main__':
    unittest.main()
