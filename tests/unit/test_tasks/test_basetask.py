import unittest
import tempfile
from pathlib import Path
from pestifer.tasks.basetask import BaseTask, VMDTask, restore_pdb_chains
from pestifer.core.pipeline import PipelineContext
from pestifer.core.artifacts import TXTFileArtifact


def _atom(serial, name, resname, chain, resseq, x, segid, icode=' ', altloc=' '):
    # PDB ATOM record with segid in cols 73-76.  resName occupies cols 18-21
    # (left-justified, so 4-char CHARMM resnames like TIP3 fit and chain stays at col 22).
    return (f"ATOM  {serial:>5d} {name:>4s}{altloc}{resname:<4s}{chain}{resseq:>4d}{icode}"
            f"   {x:8.3f}{0.0:8.3f}{0.0:8.3f}{0.0:6.2f}{0.0:6.2f}      {segid:<4s}\n")


class TestRestorePdbChains(unittest.TestCase):
    """restore_pdb_chains: reinstate the authoritative chain a PSF round-trip loses."""

    def _write(self, lines):
        f = tempfile.NamedTemporaryFile('w', suffix='.pdb', delete=False)
        f.writelines(lines); f.close()
        return f.name

    def _chains(self, path):
        return [l[21] for l in Path(path).read_text().splitlines()
                if l.startswith(('ATOM', 'HETATM'))]

    def test_restores_chain_when_segid_lead_differs(self):
        # reference: chain 'A' but segid 'PROA' (segid[0]='P' != chain)
        ref = self._write([_atom(1, 'N', 'ALA', 'A', 1, 1.0, 'PROA'),
                           _atom(2, 'CA', 'ALA', 'A', 1, 2.0, 'PROA')])
        # target: same atoms but chain re-derived to 'P' (segid leading char)
        tgt = self._write([_atom(1, 'N', 'ALA', 'P', 1, 9.0, 'PROA'),
                           _atom(2, 'CA', 'ALA', 'P', 1, 8.0, 'PROA')])
        n = restore_pdb_chains(tgt, ref)
        self.assertEqual(n, 2)
        self.assertEqual(self._chains(tgt), ['A', 'A'])

    def test_noop_when_chains_already_match(self):
        ref = self._write([_atom(1, 'N', 'ALA', 'A', 1, 1.0, 'A')])
        tgt = self._write([_atom(1, 'N', 'ALA', 'A', 1, 9.0, 'A')])
        self.assertEqual(restore_pdb_chains(tgt, ref), 0)
        self.assertEqual(self._chains(tgt), ['A'])

    def test_atoms_absent_from_reference_keep_their_chain(self):
        # target has an extra water not in the reference -> its chain is untouched
        ref = self._write([_atom(1, 'N', 'ALA', 'A', 1, 1.0, 'PROA')])
        tgt = self._write([_atom(1, 'N', 'ALA', 'P', 1, 9.0, 'PROA'),
                           _atom(2, 'OH2', 'TIP3', 'W', 1, 5.0, 'WT1')])
        restore_pdb_chains(tgt, ref)
        self.assertEqual(self._chains(tgt), ['A', 'W'])   # solute fixed, water kept

    def test_matches_by_segid_so_same_resid_atomname_across_chains_dont_collide(self):
        # two segments share resid+atomname; segid disambiguates
        ref = self._write([_atom(1, 'CA', 'ALA', 'A', 1, 1.0, 'PROA'),
                           _atom(2, 'CA', 'ALA', 'B', 1, 2.0, 'PROB')])
        tgt = self._write([_atom(1, 'CA', 'ALA', 'P', 1, 9.0, 'PROA'),
                           _atom(2, 'CA', 'ALA', 'P', 1, 8.0, 'PROB')])
        restore_pdb_chains(tgt, ref)
        self.assertEqual(self._chains(tgt), ['A', 'B'])

    def test_preserves_coordinates_and_only_touches_chain_column(self):
        ref = self._write([_atom(1, 'N', 'ALA', 'A', 1, 1.0, 'PROA')])
        tgt_lines = [_atom(1, 'N', 'ALA', 'P', 1, 9.0, 'PROA')]
        tgt = self._write(tgt_lines)
        restore_pdb_chains(tgt, ref)
        new = Path(tgt).read_text().splitlines()[0]
        old = tgt_lines[0].rstrip('\n')
        # every column except 22 (index 21) is unchanged
        self.assertEqual(new[:21], old[:21])
        self.assertEqual(new[22:], old[22:])
        self.assertEqual(new[21], 'A')

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