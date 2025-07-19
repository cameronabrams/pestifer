from pestifer.tasks.fetch import FetchTask
from pestifer.core.pipeline import PipelineContext, Artifact, PDBFile, CIFFile
import unittest
import os

class TestFetchTask(unittest.TestCase):

    def test_initiate_task_create_empty(self):
        pipeline = PipelineContext()
        task = FetchTask(pipeline)
        self.assertIsInstance(task, FetchTask)

    def test_initiate_task_do_fetch_pdb(self):
        pipeline = PipelineContext()
        task = FetchTask(pipeline, config_specs={'source': 'pdb', 'sourceID': '1A2B'})
        self.assertIsInstance(task, FetchTask)
        task.do()
        self.assertIn('base_coordinates', pipeline.artifacts)
        artifact = pipeline.artifacts['base_coordinates']
        self.assertIsInstance(artifact, Artifact)
        self.assertEqual(artifact.key, 'base_coordinates')
        self.assertEqual(artifact.value_type, PDBFile)
        self.assertEqual(artifact.produced_by, task)
        self.assertTrue(os.path.isfile(artifact.value.path))
        # os.remove(artifact.value.path)  # Clean up the fetched file

    def test_initiate_task_do_fetch_mmcif(self):
        pipeline = PipelineContext()
        task = FetchTask(pipeline, config_specs={'source': 'pdb', 'sourceID': '4tvp','source_format': 'cif'})
        self.assertIsInstance(task, FetchTask)
        task.do()
        self.assertIn('base_coordinates', pipeline.artifacts)
        artifact = pipeline.artifacts['base_coordinates']
        self.assertIsInstance(artifact, Artifact)
        self.assertEqual(artifact.key, 'base_coordinates')
        self.assertEqual(artifact.value_type, CIFFile)
        self.assertEqual(artifact.produced_by, task)
        self.assertTrue(os.path.isfile(artifact.value.path))
        # os.remove(artifact.value.path)  # Clean up the fetched file