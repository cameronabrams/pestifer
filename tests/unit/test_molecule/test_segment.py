import unittest
import logging
logger = logging.getLogger(__name__)
from pestifer.molecule.segment import Segment, SegmentList
from pestifer.molecule.residue import ResidueList, Residue, ResiduePlaceholder
from pestifer.objs.resid import ResID
from pestifer.objs.cleavagesite import CleavageSite
from pestifer.molecule.chainidmanager import ChainIDManager
from pestifer.core.objmanager import ObjManager
from pestifer.objs.insertion import Insertion, InsertionList
from pestifer.objs.mutation import Mutation, MutationList

class TestSegment(unittest.TestCase):
    def test_segment_creation_from_dict(self):
        input_dict = {
            'segtype': 'protein',
            'segname': 'A',
            'chainID': 'A',
            'residues': ResidueList([]),
            'subsegments': [],
            'parent_chain': 'A',
            'specs': {}
        }
        segment = Segment(input_dict)
        self.assertIsInstance(segment, Segment)
        self.assertEqual(segment.segtype, 'protein')
        self.assertEqual(segment.segname, 'A')
        self.assertEqual(segment.chainID, 'A')
        self.assertEqual(segment.parent_chain, 'A')
        self.assertEqual(segment.specs, {})

    def test_segment_creation_from_residue(self):
        residue = Residue(ResiduePlaceholder(chainID='A', resid=ResID(100), resname='ALA', segtype='protein', resolved=True))
        segment = Segment(residue)
        self.assertIsInstance(segment, Segment)
        self.assertEqual(segment.segtype, 'protein')
        self.assertEqual(segment.segname, 'A')
        self.assertEqual(segment.chainID, 'A')
        self.assertEqual(len(segment.residues), 1)
        self.assertEqual(segment.residues[0].resname, 'ALA')

    def test_segment_creation_from_residue_list(self):
        residues = ResidueList([
            Residue(ResiduePlaceholder(chainID='A', resid=ResID(100), resname='ALA', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='A', resid=ResID(101), resname='GLY', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='A', resid=ResID(102), resname='SER', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='A', resid=ResID(103), resname='THR', segtype='protein', resolved=True))
        ])
        segment = Segment(residues)
        self.assertIsInstance(segment, Segment)
        self.assertEqual(segment.segtype, 'protein')
        self.assertEqual(segment.segname, 'A')
        self.assertEqual(segment.chainID, 'A')
        self.assertEqual(len(segment.residues), 4)
        self.assertEqual(segment.residues[0].resname, 'ALA')

    def test_segment_cleave(self):
        residues = ResidueList([
            Residue(ResiduePlaceholder(chainID='A', resid=ResID(100), resname='ALA', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='A', resid=ResID(101), resname='GLY', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='A', resid=ResID(102), resname='SER', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='A', resid=ResID(103), resname='THR', segtype='protein', resolved=True))
        ])
        clv = CleavageSite('A:101-102')
        segment = Segment(residues)
        cleaved_segment = segment.cleave(clv,'B')
        self.assertIsInstance(cleaved_segment, Segment)
        self.assertEqual(cleaved_segment.segtype, 'protein')
        self.assertEqual(cleaved_segment.segname, 'B')
        self.assertEqual(cleaved_segment.chainID, 'B')
        self.assertEqual(len(cleaved_segment.residues), 2)
        logger.debug(f"Cleaved segment resid: {cleaved_segment.residues[0].resid}")
        self.assertEqual(cleaved_segment.residues[0].resid, ResID(102))

class TestSegmentList(unittest.TestCase):
    def test_segment_list_creation(self):
        segment_list = SegmentList()
        self.assertIsInstance(segment_list, SegmentList)
        self.assertEqual(len(segment_list), 0)

    def test_segment_list_append(self):
        segment = Segment({
            'segtype': 'protein',
            'segname': 'A',
            'chainID': 'A',
            'residues': ResidueList([]),
            'subsegments': [],
            'parent_chain': 'A',
            'specs': {}
        })
        segment_list = SegmentList()
        segment_list.append(segment)
        self.assertEqual(len(segment_list), 1)
        self.assertEqual(segment_list[0].segname, 'A')

    def test_segment_list_generate(self):
        residue_list = [
            Residue(ResiduePlaceholder(chainID='A', resid=ResID(1), resname='ALA', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='A', resid=ResID(2), resname='GLY', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='A', resid=ResID(3), resname='SER', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='A', resid=ResID(4), resname='THR', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='B', resid=ResID(1), resname='ALA', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='B', resid=ResID(2), resname='GLY', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='B', resid=ResID(3), resname='SER', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='B', resid=ResID(4), resname='THR', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='S', resid=ResID(404), resname='BMA', segtype='glycan', resolved=True))
        ]
        chainIDmanager = ChainIDManager()
        seq_spec={}

        segment_list = SegmentList.generate(seq_spec=seq_spec,residues=ResidueList(residue_list), chainIDmanager=chainIDmanager)
        self.assertIsInstance(segment_list, SegmentList)
        self.assertEqual(len(segment_list), 3)
        self.assertEqual(segment_list[0].segname, 'A')
        self.assertEqual(len(segment_list[0].residues), 4)
        for i in range(4):
            self.assertEqual(segment_list[0].residues[i], residue_list[i])
        self.assertEqual(segment_list[1].segname, 'B')
        self.assertEqual(len(segment_list[1].residues), 4)
        for i in range(4, 8):
            self.assertEqual(segment_list[1].residues[i-4], residue_list[i])
        self.assertEqual(segment_list[2].segname, 'S')
        self.assertEqual(len(segment_list[2].residues), 1)
        for i in range(8, 9):
            self.assertEqual(segment_list[2].residues[i-8], residue_list[i])

        self.assertEqual(segment_list[0],segment_list.get_segment_of_residue(residue_list[0]))
        self.assertEqual(segment_list[1],segment_list.get_segment_of_residue(residue_list[5]))
        self.assertEqual(segment_list[2],segment_list.get_segment_of_residue(residue_list[8]))

        segment_list.remove(segment_list[0])
        self.assertEqual(len(segment_list), 2)
        self.assertEqual(segment_list[0].segname, 'B')
        self.assertEqual(segment_list[1].segname, 'S')

    def test_segment_list_inherit_objs(self):
        residue_list = [
            Residue(ResiduePlaceholder(chainID='A', resid=ResID(1), resname='ALA', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='A', resid=ResID(2), resname='GLY', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='A', resid=ResID(3), resname='SER', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='A', resid=ResID(4), resname='THR', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='B', resid=ResID(1), resname='ALA', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='B', resid=ResID(2), resname='GLY', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='B', resid=ResID(3), resname='SER', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='B', resid=ResID(4), resname='THR', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='S', resid=ResID(404), resname='BMA', segtype='glycan', resolved=True))
        ]
        chainIDmanager = ChainIDManager()
        seq_spec={}

        segment_list = SegmentList.generate(seq_spec=seq_spec,residues=ResidueList(residue_list), chainIDmanager=chainIDmanager)

        O = ObjManager()
        mut = Mutation('A:ALA,1,GLY')
        O.ingest(mut)
        self.assertIn('seq', O)
        self.assertIn('mutations', O['seq'])
        self.assertIsInstance(O['seq']['mutations'], MutationList)

        insertion = Insertion('B,3,GGG')
        O.ingest(insertion)
        self.assertIn('seq', O)
        self.assertIn('insertions', O['seq'])
        self.assertIsInstance(O['seq']['insertions'], InsertionList)

        self.assertTrue(hasattr(segment_list[0], '_obj_manager'))
        self.assertEqual(len(segment_list[0]._obj_manager), 0)
        self.assertTrue(hasattr(segment_list[1], '_obj_manager'))
        self.assertTrue(hasattr(segment_list[2], '_obj_manager'))

        segment_list.inherit_objs(O)
        self.assertEqual(len(segment_list[0]._obj_manager), 1)
        self.assertIn('seq', segment_list[0]._obj_manager)
        self.assertIn('mutations', segment_list[0]._obj_manager['seq'])
        self.assertIsInstance(segment_list[0]._obj_manager['seq']['mutations'], MutationList)
        self.assertEqual(len(segment_list[0]._obj_manager['seq']['mutations']), 1)

        # insertions are not inheritable
        self.assertNotIn('seq', segment_list[1]._obj_manager)

    def test_segment_list_collect_residues(self):
        residue_list = [
            Residue(ResiduePlaceholder(chainID='A', resid=ResID(1), resname='ALA', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='A', resid=ResID(2), resname='GLY', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='A', resid=ResID(3), resname='SER', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='A', resid=ResID(4), resname='THR', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='B', resid=ResID(1), resname='ALA', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='B', resid=ResID(2), resname='GLY', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='B', resid=ResID(3), resname='SER', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='B', resid=ResID(4), resname='THR', segtype='protein', resolved=True)),
            Residue(ResiduePlaceholder(chainID='S', resid=ResID(404), resname='BMA', segtype='glycan', resolved=True))
        ]
        chainIDmanager = ChainIDManager()
        seq_spec={}

        segment_list = SegmentList.generate(seq_spec=seq_spec,residues=ResidueList(residue_list), chainIDmanager=chainIDmanager)

        collected_residues = segment_list.collect_residues()
        self.assertIsInstance(collected_residues, ResidueList)
        self.assertEqual(len(collected_residues), 9)
        self.assertEqual(collected_residues[0].resid, ResID(1))
