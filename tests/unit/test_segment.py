import unittest
from pestifer.molecule import Molecule
from pestifer.config import ConfigSetup
from pestifer.segment import SegmentList
import numpy as np
class TestSegment(unittest.TestCase):

    def test_segment(self):
        self.config=ConfigSetup('')
        m=Molecule.from_rcsb(pdb_code='4zmj').activate_biological_assembly(1)
        au=m.asymmetric_unit
        protein_segments=[x for x in au.Segments if x.segtype=='PROTEIN']
        with open('psfgentest.tcl','w') as f:
            for s in protein_segments:
                f.write(s.protein_stanza(np.identity(3),{}))
        self.assertTrue(len(protein_segments),2)
        self.assertEqual(au.chainIDs,{'G','B','A','D','C'})


 