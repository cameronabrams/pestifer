# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import unittest
from pathlib import Path
from pidibble.pdbparse import PDBParser

from pestifer.molecule.atom import AtomList, Atom
from pestifer.molecule.residue import Residue, ResiduePlaceholder, ResidueList, ResiduePlaceholderList
from pestifer.objs.deletion import Deletion, DeletionList
from pestifer.objs.insertion import Insertion, InsertionList
from pestifer.objs.link import Link, LinkList
from pestifer.objs.resid import ResID
from pestifer.objs.substitution import Substitution, SubstitutionList
from pestifer.util.cifutil import CIFload

class TestResiduePlaceholder(unittest.TestCase):

    def test_empty_residue(self):
        e = ResiduePlaceholder(
            resname='ALA',
            resid=ResID(1),
            chainID='A',
            resolved=False,
            segtype='UNSET',
            segname='A'
        )
        self.assertEqual(e.resname, 'ALA')
        self.assertEqual(e.resid.resseqnum, 1)
        self.assertEqual(e.resid.insertion, None)
        self.assertEqual(e.chainID, 'A')
        self.assertFalse(e.resolved)
        self.assertEqual(e.segtype, 'UNSET')
        self.assertFalse(hasattr(e, 'atoms'))

    def test_empty_residue_from_shortcode(self):
        shortcode = "1:A:ALA1:"
        r = ResiduePlaceholder(shortcode)
        self.assertEqual(r.model, 1)
        self.assertEqual(r.chainID, 'A')
        self.assertEqual(r.resname, 'ALA')
        self.assertEqual(r.resid, ResID(1))
        self.assertFalse(r.resolved)
        self.assertEqual(r.segtype, 'UNSET')
        self.assertEqual(r.auth_asym_id, None)
        self.assertEqual(r.auth_comp_id, None)
        self.assertEqual(r.auth_seq_id, None)

class TestResiduePlaceholderList(unittest.TestCase):

    inputs_dir = Path(__file__).parents[2] / "inputs"


    def test_empty_residue_list(self):
        rl = ResidueList()
        self.assertIsInstance(rl, ResidueList)
        self.assertEqual(len(rl), 0)

    def test_empty_residue_list_from_pdb(self):
        p = PDBParser(filepath=self.inputs_dir / '4zmj.pdb').parse().parsed
        empty_residue_list = ResiduePlaceholderList.from_pdb(p)
        self.assertGreater(len(empty_residue_list), 0)

    def test_empty_residue_list_from_cif(self):
        p = CIFload(Path(self.inputs_dir / '4zmj.cif'))
        empty_residue_list = ResiduePlaceholderList.from_cif(p)
        self.assertGreater(len(empty_residue_list), 0)

class TestResidue(unittest.TestCase):

    def test_residue_creation(self):
        r = Residue(
            resname='ALA',
            resid=ResID(1),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        self.assertIsInstance(r, Residue)
        self.assertEqual(r.resname, 'ALA')
        self.assertEqual(r.resid, ResID(1))
        self.assertEqual(r.chainID, 'A')
        self.assertTrue(r.resolved)
        self.assertEqual(r.segtype, 'UNSET')
        self.assertIsInstance(r.atoms, AtomList)  # Atoms should be an AtomList
    
    def test_residue_from_residueplaceholder(self):
        e = ResiduePlaceholder(
            resname='ALA',
            resid=ResID(1),
            chainID='A',
            resolved=False,
            segname='A',
            segtype='UNSET'
        )
        r = Residue(e)
        self.assertIsInstance(r, Residue)
        self.assertEqual(r.resname, 'ALA')
        self.assertEqual(r.resid, ResID(1))
        self.assertEqual(r.chainID, 'A')
        self.assertFalse(r.resolved)
        self.assertEqual(r.segtype, 'UNSET')
        self.assertIsInstance(r.atoms, AtomList)

    def test_residue_create_and_linking(self):
        r1 = Residue(
            resname='ALA',
            resid=ResID(1),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        r2 = Residue(
            resname='GLY',
            resid=ResID(2),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        r1.down.append(r2)
        r2.up.append(r1)
        self.assertEqual(r1.down[0], r2)
        self.assertEqual(r2.up[0], r1)
        self.assertIsInstance(r1.down, ResidueList)
        self.assertIsInstance(r2.up, ResidueList)

        link1 = Link('A_1_C-G_2_N')
        r1.uplink.append(link1)
        r2.downlink.append(link1)
        self.assertEqual(r1.uplink[0], link1)
        self.assertEqual(r2.downlink[0], link1)
        self.assertIsInstance(r1.uplink, LinkList)
        self.assertIsInstance(r2.downlink, LinkList)

    def test_residue_create_and_link_to(self):
        r1 = Residue(
            resname='ALA',
            resid=ResID(1),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        r2 = Residue(
            resname='GLY',
            resid=ResID(2),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        link1 = Link('A_1_C-G_2_N')
        r1.link_to(r2,link1)
        self.assertIn(r2, r1.down)
        self.assertIn(r1, r2.up)
        self.assertIn(link1, r1.downlink)
        self.assertIn(link1, r2.uplink)
        r1.unlink(r2,link1)
        self.assertNotIn(r2, r1.down)
        self.assertNotIn(r1, r2.up)
        self.assertNotIn(link1, r1.downlink)
        self.assertNotIn(link1, r2.uplink)

    def test_residue_pair_relations(self):
        r1 = Residue(
            resname='ALA',
            resid=ResID(1),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        r2 = Residue(
            resname='GLY',
            resid=ResID(2),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        self.assertTrue(r1<r2)
        r3 = Residue(
            resname='GLY',
            resid=ResID('2A'),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        r4 = Residue(
            resname='ALA',
            resid=ResID(1),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        self.assertTrue(r2<r3)
        self.assertFalse(r2.same_resid(r3))
        self.assertEqual(r1,r4)

    def test_residue_add_atom(self):
        r = Residue(
            resname='ALA',
            resid=ResID(1),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        atom = Atom(
            serial=1,
            name='C',
            altloc=' ',
            resname='ALA',
            chainID='A',
            resid=ResID(1),
            x=1.0,
            y=2.0,
            z=3.0,
            occ=1.0,
            beta=0.0,
            elem='C',
            charge=' ',
            segname='A'
        )
        self.assertTrue(r.add_atom(atom))
        self.assertIn(atom, r.atoms)
        self.assertEqual(len(r.atoms), 1)
        another_atom = Atom(
            serial=2,
            name='O',
            altloc=' ',
            resname='GLY',
            chainID='A',
            resid=ResID(1),
            x=1.5,
            y=2.5,
            z=3.5,
            occ=1.0,
            beta=0.0,
            elem='O',
            charge=' ',
            segname='A'
        )
        self.assertFalse(r.add_atom(another_atom))

        r.set_chainID('B')
        self.assertEqual(r.chainID, 'B')
        self.assertEqual(r.atoms[0].chainID, 'B')

    def test_residue_get_down_group(self):
        r1 = Residue(
            resname='ALA',
            resid=ResID(1),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        r2 = Residue(
            resname='GLY',
            resid=ResID(2),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        link1=  Link('A_1_C-G_2_N')
        r1.link_to(r2,link1)
        downstream_residues, downstream_links = r1.get_down_group()
        self.assertIn(r2, downstream_residues)
        self.assertIn(link1, downstream_links)
        self.assertIsInstance(downstream_residues, ResidueList)
        self.assertIsInstance(downstream_links, LinkList)

class TestResidueList(unittest.TestCase):
    def test_residue_list_creation(self):
        rl = ResidueList()
        self.assertIsInstance(rl, ResidueList)
        self.assertEqual(len(rl), 0)

    def test_residue_list_add(self):
        r1 = Residue(
            resname='ALA',
            resid=ResID(1),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        rl = ResidueList()
        rl.append(r1)
        self.assertIn(r1, rl)
        self.assertEqual(len(rl), 1)

    def test_residue_list_from_atomlist(self):
        atom1 = Atom(
            serial=1,
            name='C',
            altloc=' ',
            resname='ALA',
            chainID='A',
            resid=ResID(1),
            x=1.0,
            y=2.0,
            z=3.0,
            occ=1.0,
            beta=0.0,
            elem='C',
            charge=' ',
            segname='A'
        )
        atom2 = Atom(
            serial=2,
            name='O',
            altloc=' ',
            resname='GLY',
            chainID='A',
            resid=ResID(2),
            x=1.5,
            y=2.5,
            z=3.5,
            occ=1.0,
            beta=0.0,
            elem='O',
            charge=' ',
            segname='A'
        )
        atom_list = AtomList([atom1, atom2])
        residue_list = ResidueList.from_residuegrouped_atomlist(atom_list)
        self.assertEqual(len(residue_list), 2)
        self.assertIsInstance(residue_list[0], Residue)
        self.assertIsInstance(residue_list[1], Residue)

    def test_residue_list_get_sublist(self):
        r1 = Residue(
            resname='ALA',
            resid=ResID(1),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        r2 = Residue(
            resname='GLY',
            resid=ResID(2),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        r3 = Residue(
            resname='SER',
            resid=ResID(3),
            chainID='B',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        rl = ResidueList([r1, r2, r3])
        sublist = rl.get(lambda x: x.chainID == 'A')
        self.assertEqual(len(sublist), 2)
        self.assertIn(r1, sublist)
        self.assertIn(r2, sublist)
        self.assertNotIn(r3, sublist)

    def test_residue_list_get_by_resname(self):
        r1 = Residue(
            resname='ALA',
            resid=ResID(1),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        r2 = Residue(
            resname='GLY',
            resid=ResID(2),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        r3 = Residue(
            resname='ALA',
            resid=ResID(3),
            chainID='B',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        rl = ResidueList([r1, r2, r3])
        alist = rl.get(lambda x: x.resname == 'ALA')
        self.assertEqual(len(alist), 2)
        self.assertIn(r1, alist)
        self.assertIn(r3, alist)
        self.assertNotIn(r2, alist)

    def test_residue_list_get_by_resseqnum_insertion(self):
        r1 = Residue(
            resname='ALA',
            resid=ResID(1),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        r2 = Residue(
            resname='GLY',
            resid=ResID(1, 'A'),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        r3 = Residue(
            resname='SER',
            resid=ResID(3, 'A'),
            chainID='B',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        rl = ResidueList([r1, r2, r3])
        singleton = rl.get(lambda x: x.resid == ResID(1, 'A'))
        self.assertIsInstance(singleton, Residue)
        self.assertEqual(singleton, r2)

    def test_residue_list_apply_segtypes(self):
        r1 = Residue(
            resname='ALA',
            resid=ResID(1),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        r2 = Residue(
            resname='GLY',
            resid=ResID(2),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        rl = ResidueList([r1, r2])
        rl.apply_segtypes()
        self.assertEqual(r1.segtype, 'protein')
        self.assertEqual(r2.segtype, 'protein')

    def test_residue_list_deletions(self):
        r1 = Residue(
            resname='ALA',
            resid=ResID(1),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        r2 = Residue(
            resname='GLY',
            resid=ResID(2),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        r3 = Residue(
            resname='SER',
            resid=ResID(3),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        r4 = Residue(
            resname='THR',
            resid=ResID(4),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        rl = ResidueList([r1, r2, r3, r4])
        deletion1 = Deletion(
            chainID="A",
            resid1=ResID(2),
            resid2=ResID(3)
        )
        dlist=DeletionList([deletion1])
        excised=rl.deletion(dlist)
        self.assertEqual(len(excised), 2)
        self.assertIn(r2, excised)
        self.assertIn(r3, excised)
        self.assertNotIn(r1, excised)
        self.assertNotIn(r4, excised)
        self.assertIn(r1, rl)
        self.assertIn(r4, rl)

    def test_residue_list_substitutions(self):
        r1 = Residue(
            resname='ALA',
            resid=ResID(1),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        r2 = Residue(
            resname='GLY',
            resid=ResID(2),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        r3 = Residue(
            resname='SER',
            resid=ResID(3),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        r4 = Residue(
            resname='THR',
            resid=ResID(4, 'A'),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        rl = ResidueList([r1, r2, r3, r4])
        substitution1 = Substitution(
            chainID="A",
            resid1=ResID(2),
            resid2=ResID(3),
            subseq="WAVE"
        )
        slist=SubstitutionList([substitution1])
        chain: ResidueList = rl.get(lambda x: x.chainID == substitution1.chainID)
        self.assertTrue(len(chain) > 0)
        seqadvs, substituted=rl.substitutions(slist)
        self.assertEqual(len(substituted), 0) # subseq is greater in length that the number of residues substituted
        self.assertEqual(len(seqadvs), 2)
        self.assertEqual(len(rl), 4)

    def test_residue_list_insertions(self):
        r1 = Residue(
            resname='ALA',
            resid=ResID(1),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        r2 = Residue(
            resname='GLY',
            resid=ResID(2),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        r3 = Residue(
            resname='SER',
            resid=ResID(3),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        r4 = Residue(
            resname='THR',
            resid=ResID(4),
            chainID='A',
            resolved=True,
            segtype='UNSET',
            segname='A',
            atoms=AtomList([])  # Empty list for atoms
        )
        rl = ResidueList([r1, r2, r3, r4])
        insertion1 = Insertion(
            chainID="A",
            resid=ResID(3),
            sequence="GGG",
            integer_increment=False
        )
        ilist=InsertionList([insertion1])
        self.assertEqual(len(rl), 4)
        rl.apply_insertions(ilist)
        self.assertEqual(len(rl), 7)