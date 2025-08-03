import unittest
from pestifer.molecule.bioassemb import BioAssemb, Transform, TransformList, BioAssembList
from mmcif.api.PdbxContainers import DataContainer
from pestifer.molecule.asymmetricunit import AsymmetricUnit
from pidibble.pdbparse import PDBRecord, PDBRecordDict
import numpy as np

class TestBioAssemb(unittest.TestCase):
    def test_bioassemb_init(self):
        ba = BioAssemb()
        self.assertIsInstance(ba.transforms, TransformList)
        self.assertEqual(ba.name, 'Assembly1')
        self.assertEqual(ba.parent_molecule, None)

    def test_bioassemb_init_with_transforms(self):
        transforms = TransformList([Transform(index=1, tmat=np.eye(4), applies_chainIDs=['A'])])
        self.assertEqual(len(transforms), 1)
        ba = BioAssemb(transforms=transforms, name='TestAssembly')
        self.assertEqual(ba.name, 'TestAssembly')
        self.assertEqual(len(ba.transforms), 1)
        self.assertEqual(ba.transforms[0].index, 1)
