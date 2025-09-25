import unittest
import os
from pathlib import Path
from pestifer.core.controller import Controller
from pestifer.core.config import Config 
from pestifer.tasks.continuation import ContinuationTask
from pestifer.tasks.terminate import TerminateTask
from pestifer.molecule.molecule import Molecule
from pestifer.core.artifacts import StateArtifacts, YAMLFileArtifact
class TestContinuationTask(unittest.TestCase):

    def setUp(self):
        self.controller=Controller().configure(Config().configure_new(), terminate=False)    
        input_dir = Path('../fixtures/continuation_inputs')
        # copy inputs to working directory
        psf = input_dir / 'my_6pti.psf'
        pdb = input_dir / 'my_6pti.pdb'
        xsc = input_dir / 'my_6pti.xsc'
        coor = input_dir / 'my_6pti.coor'
        vel = input_dir / 'my_6pti.vel'
        # copy to cwd
        dest_psf = Path('my_6pti.psf')
        dest_pdb = Path('my_6pti.pdb')
        dest_xsc = Path('my_6pti.xsc')
        dest_coor = Path('my_6pti.coor')
        dest_vel = Path('my_6pti.vel')
        dest_chainmap = Path('chainmap.yaml')
        if dest_chainmap.exists():
            dest_chainmap.unlink()
        dest_tarball = Path('my_system.tar.gz')
        if dest_tarball.exists():
            dest_tarball.unlink()
        if dest_psf.exists():
            dest_psf.unlink()
        if dest_pdb.exists():
            dest_pdb.unlink()
        if dest_xsc.exists():
            dest_xsc.unlink()
        if dest_coor.exists():
            dest_coor.unlink()
        if dest_vel.exists():
            dest_vel.unlink()
        os.symlink(psf.resolve(), dest_psf)
        os.symlink(pdb.resolve(), dest_pdb)
        os.symlink(xsc.resolve(), dest_xsc)
        os.symlink(coor.resolve(), dest_coor)
        os.symlink(vel.resolve(), dest_vel)

    def test_continuation_task(self):
        task_list = [{'continuation': dict(pdb='my_6pti.pdb', psf='my_6pti.psf', xsc='my_6pti.xsc', coor='my_6pti.coor',vel='my_6pti.vel')}]
        self.controller.reconfigure_tasks(task_list)
        self.assertEqual(len(self.controller.tasks), 1)
        self.assertIsInstance(self.controller.tasks[0], ContinuationTask)
        self.controller.do_tasks()
        task = self.controller.tasks[0]
        state: StateArtifacts = task.get_current_artifact('state')
        psf = state.psf.name
        self.assertEqual(psf, 'my_6pti.psf')
        base_molecule: Molecule = task.get_current_artifact_data('base_molecule')
        self.assertIsInstance(base_molecule, Molecule)

    def test_termination_task(self):
        task_list = [
            {'continuation':
              dict(pdb='my_6pti.pdb', psf='my_6pti.psf', xsc='my_6pti.xsc', coor='my_6pti.coor', vel='my_6pti.vel')}, 
            {'terminate':dict(cleanup=False,chainmapfile='chainmap.yaml')}]
        self.controller.reconfigure_tasks(task_list)
        self.assertEqual(len(self.controller.tasks), 2)
        self.assertIsInstance(self.controller.tasks[0], ContinuationTask)
        self.assertIsInstance(self.controller.tasks[1], TerminateTask)
        self.controller.do_tasks()
        task = self.controller.tasks[0]
        state = self.controller.pipeline.get_state_artifact(produced_by=task)
        # state: StateArtifacts = task.get_current_artifact('state')
        psf = state.psf
        self.assertEqual(psf.name, 'my_6pti.psf')
        pdb = state.pdb
        self.assertEqual(pdb.name, 'my_6pti.pdb')
        vel = state.vel
        self.assertEqual(vel.name, 'my_6pti.vel')
        xsc = state.xsc
        self.assertEqual(xsc.name, 'my_6pti.xsc')
        coor = state.coor
        self.assertEqual(coor.name, 'my_6pti.coor')
        self.assertTrue(psf.exists())
        self.assertTrue(pdb.exists())
        self.assertTrue(vel.exists())
        self.assertTrue(xsc.exists())
        self.assertTrue(coor.exists())
        psf.path.unlink()
        pdb.path.unlink()
        vel.path.unlink()
        xsc.path.unlink()
        coor.path.unlink()
        base_molecule: Molecule = task.get_current_artifact_data('base_molecule')
        self.assertIsInstance(base_molecule, Molecule)
        chainmapfile: YAMLFileArtifact = task.get_current_artifact('chainmapfile')
        self.assertIsInstance(chainmapfile, YAMLFileArtifact)
        self.assertTrue(Path('my_system.tar.gz').exists())
        self.assertTrue(chainmapfile.exists())
        Path('my_system.tar.gz').unlink()
        chainmapfile.path.unlink()

