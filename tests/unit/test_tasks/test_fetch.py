from pestifer.tasks.fetch import FetchTask
from pestifer.core.pipeline import PipelineContext
from pestifer.core.artifacts import Artifact, PDBFile, CIFFile
import unittest
import os

class TestFetchTask(unittest.TestCase):

    def test_fetch_task_create_empty(self):
        pipeline = PipelineContext()
        task = FetchTask()
        self.assertIsInstance(task, FetchTask)
        self.assertFalse(task.is_provisioned)
        execute_attempt_result = task.execute()
        self.assertEqual(execute_attempt_result, -1, "Executing an unprovisioned task should return -1.")
        self.assertFalse(hasattr(task, 'pipeline'))
        task.provision(packet=dict(pipeline=pipeline))
        self.assertTrue(hasattr(task, 'pipeline'))
        with self.assertRaises((ValueError, AssertionError)):
            task.execute()

    def test_fetch_task_do_fetch_pdb(self):
        pipeline = PipelineContext()
        task = FetchTask(specs={'source': 'pdb', 'sourceID': '1A2B'})
        task.provision(packet=dict(pipeline=pipeline))
        self.assertIsInstance(task, FetchTask)
        task.execute()
        self.assertIn('base_coordinates', pipeline.head)
        artifact = pipeline.get_artifact('base_coordinates')
        self.assertIsInstance(artifact, Artifact)
        self.assertIsInstance(artifact, PDBFile)
        self.assertEqual(artifact.key, 'base_coordinates')
        self.assertEqual(artifact.produced_by, task)
        self.assertTrue(artifact.exists())
        os.remove(artifact.path)  # Clean up the fetched file
        self.assertFalse(os.path.exists(artifact.path), "Fetched PDB file should be removed after test.")

    def test_fetch_task_do_fetch_mmcif(self):
        pipeline = PipelineContext()
        task = FetchTask(specs={'source': 'pdb', 'sourceID': '4tvp','source_format': 'cif'})
        task.provision(packet=dict(pipeline=pipeline))
        self.assertIsInstance(task, FetchTask)
        task.execute()
        self.assertIn('base_coordinates', pipeline.head)
        artifact = pipeline.get_artifact('base_coordinates')
        self.assertIsInstance(artifact, Artifact)
        self.assertIsInstance(artifact, CIFFile)
        self.assertEqual(artifact.key, 'base_coordinates')
        self.assertEqual(artifact.produced_by, task)
        self.assertTrue(artifact.exists())
        os.remove(artifact.path)  # Clean up the fetched file
        self.assertFalse(os.path.exists(artifact.path), "Fetched CIF file should be removed after test.")