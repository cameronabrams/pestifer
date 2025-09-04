import unittest
from collections import UserDict
from pestifer.objs.graft import Graft, GraftList
from pestifer.objs.resid import ResID
from pestifer.objs.link import Link, LinkList
from pestifer.molecule.residue import Residue, ResidueList
from pestifer.molecule.atom import AtomList
from pestifer.molecule.molecule import Molecule
from pestifer.molecule.segment import Segment, SegmentList
from pestifer.core.objmanager import ObjManager
class TestGraft(unittest.TestCase):
    def test_graft_creation(self):
        graft = Graft(
            chainID="A",
            target_root=ResID(20),
            source_pdbid="4ijk",
            source_chainID="A",
            source_root=ResID(1),
        )
        self.assertIsInstance(graft, Graft)
        self.assertEqual(graft.chainID, "A")
        self.assertEqual(graft.target_root.resid, 20)
        self.assertEqual(graft.source_pdbid, "4ijk")
        self.assertEqual(graft.source_chainID, "A")
        self.assertEqual(graft.source_root.resid, 1)
        self.assertEqual(graft.obj_id, 0)
        self.assertEqual(graft._yaml_header, 'grafts')
        self.assertEqual(graft._objcat, 'seq')

    def test_graft_from_shortcode(self):
        """
        The shortcode format is "target:source", where:
            - target is in the format "C_RRR[-SSS]"
              - C is the chainID
              - RRR is the residue number and insertion code of the first residue in the target alignment basis
              - SSS is the optional residue number and insertion code of the last residue in the target alignment basis
            - source is in the format "pdbid,C_RRR[#SSS][-TTT]"
              - pdbid is the basename of the PDB file or PDB ID
              - C is the chainID in the source PDB
              - RRR is the residue number and insertion code of the first residue in the source alignment basis
              - SSS is optional, as above
              - TTT is optional, completing RRR-TTT range for entire graft source
        """
        shortcode = "A_520:4ijk,A_1-10"
        graft = Graft(shortcode)
        self.assertIsInstance(graft, Graft)
        self.assertEqual(graft.chainID, "A")
        self.assertEqual(graft.target_root.resid, 520)
        self.assertEqual(graft.source_pdbid, "4ijk")
        self.assertEqual(graft.source_chainID, "A")
        self.assertEqual(graft.source_root.resid, 1)
        self.assertEqual(graft.source_end.resid, 10)
        self.assertEqual(graft.obj_id, 0)
        # do it again to make sure the counter increments
        graft2 = Graft(shortcode)
        self.assertEqual(graft2.obj_id, 1)

    def test_graft_activation(self):
        mol = Molecule() # empty donor molecule
        self.assertTrue(hasattr(mol, 'asymmetric_unit'))
        au = mol.asymmetric_unit
        self.assertTrue(hasattr(au, 'segments'))
        self.assertIsInstance(au.segments, SegmentList)
        null_seg = au.segments.get(lambda x: x.segname == 'A')
        self.assertIsNone(null_seg)
        om = mol.objmanager
        self.assertTrue(isinstance(om, ObjManager))
        self.assertTrue(isinstance(om, UserDict))
        null_topol = om.get('topol', {})
        self.assertEqual(null_topol, {})
        # make a donor segment
        r = ResidueList([
                Residue(chainID='G',
                        resname='NAG', 
                        resid=ResID(1),
                        atoms=AtomList([]),
                        segtype='glycan',
                        segname='G',
                        resolved=True),
                Residue(chainID='G', 
                        resname='NAG', 
                        resid=ResID(2),
                        atoms=AtomList([]),
                        segtype='glycan',
                        segname='G',
                        resolved=True)
            ])
        s = Segment(
            chainID="G",
            segname='G', 
            residues=r
        )
        # make a mol that has one segment in 
        # its segment list
        mol.asymmetric_unit.segments = SegmentList([s])
        # make a graft that should match the
        # first residue of the segment
        g = Graft(
            chainID="A",
            target_root=ResID(666),
            source_pdbid="4ijk",
            source_chainID="G",
            source_root=ResID(1),
            source_end=ResID(2)
        )
        self.assertIsInstance(g, Graft)
        self.assertEqual(g.source_end, ResID(2))
        g.activate(mol)
        self.assertEqual(g.source_molecule, mol)
        self.assertEqual(g.source_seg, s)
        self.assertEqual(g.source_root, ResID(1))
        self.assertEqual(len(g.index_residues), 1)
        self.assertEqual(g.index_residues[0], r[0])
        self.assertEqual(len(g.donor_residues), 1)
        self.assertEqual(g.donor_residues[0], r[1])

        recv = Molecule()  # receiver molecule
        # make a receiver segment
        r = ResidueList([
                Residue(chainID='A',
                        resname='ASN', 
                        resid=ResID(666),
                        atoms=AtomList([]),
                        segtype='protein',
                        segname='A',
                        resolved=True),
                Residue(chainID='A', 
                        resname='SER', 
                        resid=ResID(667),
                        atoms=AtomList([]),
                        segtype='protein',
                        segname='A',
                        resolved=True)
            ])
        test_res = r.get(lambda x: x.chainID == 'A' and x.resid == ResID(666))
        self.assertEqual(test_res, r[0])
        s = Segment(
            chainID="A",
            segname='A', 
            residues=r
        )
        recv.asymmetric_unit.segments = SegmentList([s])
        g.assign_receiver_residues(r)
        self.assertEqual(len(g.residues), 1)

class TestGraftList(unittest.TestCase):
    def test_graft_list_creation(self):
        graft_list = GraftList()
        self.assertIsInstance(graft_list, GraftList)
        self.assertEqual(len(graft_list), 0)

    def test_graft_list_addition(self):
        graft = Graft(
            chainID="A",
            target_root=ResID(520),
            source_pdbid="4ijk",
            source_chainID="A",
            source_root=ResID(1)
        )
        graft_list = GraftList()
        graft_list.append(graft)
        self.assertEqual(len(graft_list), 1)
        self.assertEqual(graft_list[0], graft)

    def test_graft_list_assignment(self):
        graft_list = GraftList()
        shortcode = "$_520:4ijk,A_1-10"
        for c in 'ABCDEFG':
            graft = Graft(shortcode.replace('$', c))
            self.assertEqual(graft.chainID, c)
            graft_list.append(graft)
        self.assertEqual(len(graft_list), 7)