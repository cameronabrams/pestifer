import unittest
from pestifer.molecule.segment import Segment, SegmentList
from pestifer.molecule.residue import ResidueList, Residue, EmptyResidue
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
        segment = Segment.new(input_dict)
        self.assertIsInstance(segment, Segment)
        self.assertEqual(segment.segtype, 'protein')
        self.assertEqual(segment.segname, 'A')
        self.assertEqual(segment.chainID, 'A')
        self.assertEqual(segment.parent_chain, 'A')
        self.assertEqual(segment.specs, {})

    def test_segment_creation_from_residue(self):
        residue = Residue.new(EmptyResidue(chainID='A', resseqnum=100, resname='ALA', insertion='', segtype='protein', resolved=True))
        segment = Segment.new(residue)
        self.assertIsInstance(segment, Segment)
        self.assertEqual(segment.segtype, 'protein')
        self.assertEqual(segment.segname, 'A')
        self.assertEqual(segment.chainID, 'A')
        self.assertEqual(len(segment.residues), 1)
        self.assertEqual(segment.residues[0].resname, 'ALA')

    def test_segment_creation_from_residue_list(self):
        residues = ResidueList([
            Residue.new(EmptyResidue(chainID='A', resseqnum=100, resname='ALA', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='A', resseqnum=101, resname='GLY', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='A', resseqnum=102, resname='SER', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='A', resseqnum=103, resname='THR', insertion='', segtype='protein', resolved=True))
        ])
        segment = Segment.new(residues)
        self.assertIsInstance(segment, Segment)
        self.assertEqual(segment.segtype, 'protein')
        self.assertEqual(segment.segname, 'A')
        self.assertEqual(segment.chainID, 'A')
        self.assertEqual(len(segment.residues), 4)
        self.assertEqual(segment.residues[0].resname, 'ALA')

    def test_segment_cleave(self):
        residues = ResidueList([
            Residue.new(EmptyResidue(chainID='A', resseqnum=100, resname='ALA', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='A', resseqnum=101, resname='GLY', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='A', resseqnum=102, resname='SER', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='A', resseqnum=103, resname='THR', insertion='', segtype='protein', resolved=True))
        ])
        clv = CleavageSite.new('A:101-102')
        segment = Segment.new(residues)
        cleaved_segment = segment.cleave(clv,'B')
        self.assertIsInstance(cleaved_segment, Segment)
        self.assertEqual(cleaved_segment.segtype, 'protein')
        self.assertEqual(cleaved_segment.segname, 'B')
        self.assertEqual(cleaved_segment.chainID, 'B')
        self.assertEqual(len(cleaved_segment.residues), 2)
        self.assertEqual(cleaved_segment.residues[0].resseqnum, 102)

class TestSegmentList(unittest.TestCase):
    def test_segment_list_creation(self):
        segment_list = SegmentList()
        self.assertIsInstance(segment_list, SegmentList)
        self.assertEqual(len(segment_list), 0)

    def test_segment_list_append(self):
        segment = Segment.new({
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
            Residue.new(EmptyResidue(chainID='A', resseqnum=1, resname='ALA', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='A', resseqnum=2, resname='GLY', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='A', resseqnum=3, resname='SER', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='A', resseqnum=4, resname='THR', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='B', resseqnum=1, resname='ALA', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='B', resseqnum=2, resname='GLY', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='B', resseqnum=3, resname='SER', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='B', resseqnum=4, resname='THR', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='S', resseqnum=404, resname='BMA', insertion='', segtype='glycan', resolved=True))
        ]
        chainIDmanager = ChainIDManager()
        seq_spec={}

        segment_list = SegmentList.generate(seq_spec=seq_spec,residues=ResidueList(residue_list), chainIDmanager=chainIDmanager)
        self.assertIsInstance(segment_list, SegmentList)
        self.assertEqual(len(segment_list), 3)
        self.assertEqual(segment_list[0].segname, 'A')
        self.assertEqual(len(segment_list[0].residues), 4)
        self.assertEqual(segment_list[1].segname, 'B')
        self.assertEqual(len(segment_list[1].residues), 4)
        self.assertEqual(segment_list[2].segname, 'S')
        self.assertEqual(len(segment_list[2].residues), 1)

        self.assertEqual(segment_list[0],segment_list.get_segment_of_residue(residue_list[0]))
        self.assertEqual(segment_list[1],segment_list.get_segment_of_residue(residue_list[4]))
        self.assertEqual(segment_list[2],segment_list.get_segment_of_residue(residue_list[8]))

        segment_list.remove(segment_list[0])
        self.assertEqual(len(segment_list), 2)
        self.assertEqual(segment_list[0].segname, 'B')
        self.assertEqual(segment_list[1].segname, 'S')

    def test_segment_list_inherit_objs(self):
        residue_list = [
            Residue.new(EmptyResidue(chainID='A', resseqnum=1, resname='ALA', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='A', resseqnum=2, resname='GLY', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='A', resseqnum=3, resname='SER', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='A', resseqnum=4, resname='THR', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='B', resseqnum=1, resname='ALA', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='B', resseqnum=2, resname='GLY', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='B', resseqnum=3, resname='SER', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='B', resseqnum=4, resname='THR', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='S', resseqnum=404, resname='BMA', insertion='', segtype='glycan', resolved=True))
        ]
        chainIDmanager = ChainIDManager()
        seq_spec={}

        segment_list = SegmentList.generate(seq_spec=seq_spec,residues=ResidueList(residue_list), chainIDmanager=chainIDmanager)

        O = ObjManager()
        mut = Mutation.new('A:ALA,1,GLY')
        O.ingest(mut)
        self.assertIn('seq', O)
        self.assertIn('mutations', O['seq'])
        self.assertIsInstance(O['seq']['mutations'], MutationList)

        insertion = Insertion.new('B,3,GGG')
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
            Residue.new(EmptyResidue(chainID='A', resseqnum=1, resname='ALA', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='A', resseqnum=2, resname='GLY', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='A', resseqnum=3, resname='SER', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='A', resseqnum=4, resname='THR', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='B', resseqnum=1, resname='ALA', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='B', resseqnum=2, resname='GLY', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='B', resseqnum=3, resname='SER', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='B', resseqnum=4, resname='THR', insertion='', segtype='protein', resolved=True)),
            Residue.new(EmptyResidue(chainID='S', resseqnum=404, resname='BMA', insertion='', segtype='glycan', resolved=True))
        ]
        chainIDmanager = ChainIDManager()
        seq_spec={}

        segment_list = SegmentList.generate(seq_spec=seq_spec,residues=ResidueList(residue_list), chainIDmanager=chainIDmanager)

        collected_residues = segment_list.collect_residues()
        self.assertIsInstance(collected_residues, ResidueList)
        self.assertEqual(len(collected_residues), 9)
        self.assertEqual(collected_residues[0].resseqnum, 1)
        self.assertEqual(collected_residues[8].resseqnum, 404)
