import unittest

from pestifer.core.controller import Controller
from pestifer.core.config import Config 
from pestifer.tasks.continuation import ContinuationTask
from pestifer.tasks.terminate import TerminateTask
from pestifer.molecule.molecule import Molecule
from pestifer.core.artifacts import StateArtifacts, PSFFileArtifact, PDBFileArtifact, NAMDCoorFileArtifact, NAMDXscFileArtifact, NAMDVelFileArtifact, YAMLFileArtifact
class TestContinuationTask(unittest.TestCase):

    def setUp(self):
        self.controller=Controller().configure(Config().configure_new())    

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
        task_list = [{'continuation': dict(pdb='my_6pti.pdb', psf='my_6pti.psf', xsc='my_6pti.xsc', coor='my_6pti.coor', vel='my_6pti.vel')}, {'terminate':dict(cleanup=False)}]
        self.controller.reconfigure_tasks(task_list)
        self.assertEqual(len(self.controller.tasks), 2)
        self.assertIsInstance(self.controller.tasks[0], ContinuationTask)
        self.assertIsInstance(self.controller.tasks[1], TerminateTask)
        self.controller.do_tasks()
        task = self.controller.tasks[0]
        state: StateArtifacts = task.get_current_artifact('state')
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
        base_molecule: Molecule = task.get_current_artifact_data('base_molecule')
        self.assertIsInstance(base_molecule, Molecule)
        chainmapfile: YAMLFileArtifact = task.get_current_artifact('chainmapfile')
        self.assertIsInstance(chainmapfile, YAMLFileArtifact)
        self.assertTrue(chainmapfile.exists())
