import unittest
from pestifer.tasks.basetask import BaseTask, VMDTask
from pestifer.core.pipeline import PipelineContext
from pestifer.core.artifacts import TXTFileArtifact

class TestBaseTask(unittest.TestCase):
    
    def test_basetask_is_abstract(self):
        with self.assertRaises(TypeError):
            BaseTask()
    
    def test_basetask_subclass(self):
        class ConcreteTask(BaseTask):
            def do(self):
                pass
        
        task = ConcreteTask()
        self.assertIsInstance(task, BaseTask)
        self.assertTrue(hasattr(task, 'do'))

    def test_basetask_subclass_with_pipeline(self):
        class ConcreteTaskWithPipeline(BaseTask):            
            def do(self):
                pass
        pipeline_context = PipelineContext()
        task = ConcreteTaskWithPipeline(specs={})
        task.provision(packet=dict(pipeline=pipeline_context))
        
        self.assertTrue(task.is_provisioned)
        self.assertIsInstance(task, BaseTask)
        self.assertTrue(hasattr(task, 'do'))
        self.assertEqual(task.pipeline, pipeline_context)
        task.register('afile', key='test_artifact', artifact_type=TXTFileArtifact)
        self.assertIn('test_artifact', task.pipeline.head)

class TestVMDTask(unittest.TestCase):
    
    def test_vmdtask_is_abstract(self):
        with self.assertRaises(TypeError):
            VMDTask()
    
    def test_vmdtask_subclass(self):
        class ConcreteVMDTask(VMDTask):
            def do(self):
                pass