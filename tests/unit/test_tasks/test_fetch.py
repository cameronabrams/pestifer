from pestifer.tasks.fetch import FetchTask
from pestifer.core.pipeline import PipelineContext
from pestifer.core.artifacts import Artifact, PDBFile, CIFFile
import unittest
import os

class TestFetchTask(unittest.TestCase):

    def test_fetch_task_create_empty(self):
        pipeline = PipelineContext()
        task = FetchTask(pipeline)
        self.assertIsInstance(task, FetchTask)

    def test_fetch_task_do_fetch_pdb(self):
        pipeline = PipelineContext()
        task = FetchTask(pipeline, config_specs={'source': 'pdb', 'sourceID': '1A2B'})
        self.assertIsInstance(task, FetchTask)
        task.do()
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
        task = FetchTask(pipeline, config_specs={'source': 'pdb', 'sourceID': '4tvp','source_format': 'cif'})
        self.assertIsInstance(task, FetchTask)
        task.do()
        self.assertIn('base_coordinates', pipeline.head)
        artifact = pipeline.get_artifact('base_coordinates')
        self.assertIsInstance(artifact, Artifact)
        self.assertIsInstance(artifact, CIFFile)
        self.assertEqual(artifact.key, 'base_coordinates')
        self.assertEqual(artifact.produced_by, task)
        self.assertTrue(artifact.exists())
        os.remove(artifact.path)  # Clean up the fetched file
        self.assertFalse(os.path.exists(artifact.path), "Fetched CIF file should be removed after test.")