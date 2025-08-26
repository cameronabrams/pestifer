
import unittest
from pestifer.molecule.molecule import Molecule
from pestifer.molecule.bioassemb import BioAssemb, Transform, TransformList, BioAssembList
from mmcif.api.PdbxContainers import DataContainer
from pestifer.molecule.asymmetricunit import AsymmetricUnit
import numpy as np

class TestBioAssemb(unittest.TestCase):

    def test_bioassemb_init(self):
        ba = BioAssemb()
        self.assertIsInstance(ba.transforms, TransformList)
        self.assertEqual(ba.name, 'Assembly0')
        self.assertEqual(ba.parent_molecule, None)

    def test_bioassemb_init_with_transforms(self):
        transforms = TransformList([Transform(index=1, tmat=np.eye(4), applies_chainIDs=['A'])])
        self.assertEqual(len(transforms), 1)
        ba = BioAssemb(transforms=transforms, name='TestAssembly')
        self.assertEqual(ba.name, 'TestAssembly')
        self.assertEqual(len(ba.transforms), 1)
        self.assertEqual(ba.transforms[0].index, 1)

    def test_bioassemb_set_parent_molecule(self):
        m = Molecule()
        ba = BioAssemb()
        ba.set_parent_molecule(m)
        self.assertEqual(ba.parent_molecule, m)
