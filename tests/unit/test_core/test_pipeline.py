from pestifer.core.pipeline import *
from pestifer.core.artifacts import *
import unittest

class TestPipeline(unittest.TestCase):
    def test_pipeline_initialization(self):
        pipeline = PipelineContext()
        self.assertIsInstance(pipeline, PipelineContext)

    def test_register_artifact(self):
        pipeline = PipelineContext()
        pipeline.register(PDBFile('1gc1').stamp('I am the maker'))
        added_artifact = pipeline.get_artifact('pdb')
        self.assertIsInstance(added_artifact, Artifact)
        self.assertIsInstance(added_artifact, PDBFile)
        self.assertEqual(added_artifact.key, 'pdb')
        self.assertEqual(added_artifact.value, '1gc1')
        self.assertEqual(added_artifact.path, Path('1gc1.pdb'))
        self.assertEqual(added_artifact.produced_by, 'I am the maker')
        returned_artifact_value = pipeline.get_artifact_value('pdb')
        self.assertIsNotNone(returned_artifact_value)
        self.assertEqual(returned_artifact_value, '1gc1')

    def test_pipeline_history(self):
        pipeline = PipelineContext()
        artifact1 = PDBFile('1gc1').stamp('First maker')
        artifact2 = PDBFile('2gc2').stamp('Second maker')
        pipeline.register(artifact1)
        pipeline.register(artifact2)
        
        self.assertIn('pdb', pipeline.head)
        self.assertEqual(len(pipeline.history), 1)
        self.assertEqual(pipeline.history[0], artifact1)

    def test_get_artifact_series(self):
        pipeline = PipelineContext()
        artifact1 = PDBFile('1gc1').stamp('First maker')
        artifact2 = PDBFile('2gc2').stamp('Second maker')
        pipeline.register(artifact1, key='pdb_series')
        pipeline.register(artifact2, key='pdb_series')
        self.assertEqual(len(pipeline.head), 1)
        self.assertEqual(len(pipeline.history), 1)
        self.assertEqual(pipeline.get_artifact('pdb_series'), artifact2)
        self.assertEqual(pipeline.history[0], artifact1)
        self.assertEqual(pipeline.history[0].key, 'pdb_series')
        series = pipeline.get_artifact_series('pdb_series')
        self.assertIsInstance(series, ArtifactList)
        self.assertIsInstance(series, ArtifactFileList)
        # self.assertEqual(len(series), 2)
        self.assertIn(artifact1, series)
        self.assertIn(artifact2, series)

    def test_stash_current_artifact(self):
        pipeline = PipelineContext()
        artifact = PDBFile('1gc1').stamp('Stash maker')
        pipeline.register(artifact)
        self.assertIn('pdb', pipeline.head)
        pipeline.stash('pdb')
        self.assertNotIn('pdb', pipeline.head)
        self.assertEqual(len(pipeline.history), 1)
        self.assertEqual(pipeline.history[0], artifact)