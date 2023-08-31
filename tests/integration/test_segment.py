import unittest
from pestifer.molecule import Molecule
from pestifer.config import Config
from pestifer.chainids import ChainIDManager
class TestSegment(unittest.TestCase):
    def test_segment(self):
        chainIDmanager=ChainIDManager()
        m=Molecule(source='4zmj',config=Config(),chainIDmanager=chainIDmanager).activate_biological_assembly(1)
        au=m.asymmetric_unit
        protein_segments=[x for x in au.Segments if x.segtype=='PROTEIN']
        self.assertTrue(len(protein_segments),2)
        expected_set_of_chainIDs={'G','B','A','D','C','E','F'}
        self.assertEqual(set(au.Segments.segnames),expected_set_of_chainIDs)
        e=au.Segments.get(segname='E')
        chk1=all([x.chainID=='E' for x in e.residues])
        self.assertTrue(chk1)
        chk2=[]
        for x in e.residues:
            chk2.append(all([y.chainID=='E' for y in x.atoms]))
        self.assertTrue(all(chk2))
 