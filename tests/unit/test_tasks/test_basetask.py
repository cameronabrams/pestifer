import unittest
from pathlib import Path
from pestifer.tasks.basetask import BaseTask, VMDTask, restore_chain_ids_from_reference
from pestifer.core.pipeline import PipelineContext
from pestifer.core.artifacts import TXTFileArtifact


def _atom(serial, name, resid, chain, segid, x=0.0):
    # A minimally valid PDB ATOM line: chainID at col 22 (idx 21), resSeq cols 23-26,
    # segid cols 73-76. Matches the columns restore_chain_ids_from_reference reads.
    return (f'ATOM  {serial:>5d} {name:<4s} ALA {chain}{resid:>4d}    '
            f'{x:>8.3f}{0.0:>8.3f}{0.0:>8.3f}  1.00  0.00      {segid:<4s}\n')


class TestRestoreChainIds(unittest.TestCase):
    """restore_chain_ids_from_reference: chain column is restored by (segid, resid)."""

    def _write(self, path, lines):
        with open(path, 'w') as f:
            f.writelines(lines)

    def test_restores_chain_from_reference(self):
        import tempfile, os
        with tempfile.TemporaryDirectory() as d:
            ref = os.path.join(d, 'ref.pdb')
            tgt = os.path.join(d, 'tgt.pdb')
            # Reference: two segments deliberately given DIFFERENT chains (A and B).
            self._write(ref, [_atom(1, 'CA', 1, 'A', 'PROA'),
                              _atom(2, 'CA', 1, 'B', 'PROB')])
            # Target as catdcd would emit it: chain collapsed to the segid leading char.
            self._write(tgt, [_atom(1, 'CA', 1, 'P', 'PROA'),
                              _atom(2, 'CA', 1, 'P', 'PROB')])
            restore_chain_ids_from_reference(tgt, ref)
            chains = [ln[21] for ln in Path(tgt).read_text().splitlines()]
            self.assertEqual(chains, ['A', 'B'])

    def test_unmatched_atoms_keep_catdcd_chain(self):
        import tempfile, os
        with tempfile.TemporaryDirectory() as d:
            ref = os.path.join(d, 'ref.pdb')
            tgt = os.path.join(d, 'tgt.pdb')
            self._write(ref, [_atom(1, 'CA', 1, 'A', 'PROA')])
            # Target has an extra residue absent from the reference; it must be left alone.
            self._write(tgt, [_atom(1, 'CA', 1, 'P', 'PROA'),
                              _atom(2, 'CA', 9, 'W', 'WATX')])
            restore_chain_ids_from_reference(tgt, ref)
            chains = [ln[21] for ln in Path(tgt).read_text().splitlines()]
            self.assertEqual(chains, ['A', 'W'])

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