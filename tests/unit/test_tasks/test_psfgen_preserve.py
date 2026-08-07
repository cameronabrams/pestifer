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

    def test_preserve_rejects_mods_for_now(self):
        # Applying mods to a pre-built topology is P2.1+; the preserve path must hard-error rather
        # than silently ignore them.
        self.controller.reconfigure_tasks(
            [{'continuation': self._cont()}, {'psfgen': {'mods': {'ssbonds': ['A_5 A_55']}}}])
        with self.assertRaises(PestiferBuildError) as ctx:
            self.controller.do_tasks()
        self.assertIn('not yet supported', str(ctx.exception))


if __name__ == '__main__':
    unittest.main()
