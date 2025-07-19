from pestifer.core.pipeline import *
import unittest

class TestPipeline(unittest.TestCase):
    def test_pipeline_initialization(self):
        pipeline = PipelineContext()
        self.assertIsInstance(pipeline, PipelineContext)

    def test_register_artifact(self):
        pipeline = PipelineContext()
        pipeline.register('test_key', 'test_value', str, 'test_task')
        returned_artifact_value = pipeline.get('test_key')
        self.assertIsNotNone(returned_artifact_value)
        self.assertEqual(returned_artifact_value, 'test_value')
        added_artifact = pipeline.artifacts['test_key']
        self.assertIsInstance(added_artifact, Artifact)
        self.assertEqual(added_artifact.key, 'test_key')
        self.assertEqual(added_artifact.value, 'test_value')
        self.assertEqual(added_artifact.value_type, str)
        self.assertEqual(added_artifact.produced_by, 'test_task')
