import unittest
from pestifer.molecule.residue import Residue, EmptyResidue
from pestifer.molecule.atom import AtomList, Atom
from pestifer.core.config import Config
from pestifer.core.labels import Labels
from pestifer.util.cifutil import CIFdict, CIFload
from pestifer.objs.insertion import InsertionList, Insertion
from io import StringIO
from pidibble.pdbparse import PDBParser
import yaml

class TestEmptyResidue(unittest.TestCase):

    def test_empty_residue(self):
        e = EmptyResidue(
            resname='ALA',
            resseqnum=1,
            insertion='',
            chainID='A',
            resolved=False,
            segtype='UNSET'
        )
        self.assertEqual(e.resname, 'ALA')
        self.assertEqual(e.resseqnum, 1)
        self.assertEqual(e.insertion, '')
        self.assertEqual(e.chainID, 'A')
        self.assertFalse(e.resolved)
        self.assertEqual(e.segtype, 'UNSET')
        self.assertFalse(hasattr(e, 'atoms'))

    def test_empty_residue_from_shortcode(self):
        shortcode = "1:A:ALA1:"
        r = EmptyResidue.new(shortcode)
        self.assertEqual(r.model, 1)
        self.assertEqual(r.chainID, 'A')
        self.assertEqual(r.resname, 'ALA')
        self.assertEqual(r.resseqnum, 1)
        self.assertEqual(r.insertion, '')
        self.assertFalse(r.resolved)
        self.assertEqual(r.segtype, 'UNSET')
        self.assertEqual(r.auth_asym_id, None)
        self.assertEqual(r.auth_comp_id, None)
        self.assertEqual(r.auth_seq_id, None)

class TestResidue(unittest.TestCase):

    def test_residue_creation(self):
        r = Residue(
            resname='ALA',
            resseqnum=1,
            insertion='',
            chainID='A',
            resolved=True,
            segtype='UNSET',
            atoms=AtomList([])  # Empty list for atoms
        )
        self.assertIsInstance(r, Residue)
        self.assertEqual(r.resname, 'ALA')
        self.assertEqual(r.resseqnum, 1)
        self.assertEqual(r.insertion, '')
        self.assertEqual(r.chainID, 'A')
        self.assertTrue(r.resolved)
        self.assertEqual(r.segtype, 'UNSET')
        self.assertIsInstance(r.atoms, AtomList)  # Atoms should be an AtomList