# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import unittest
from pestifer.core.objmanager import ObjManager
from pestifer.objs.insertion import Insertion, InsertionList
from pestifer.objs.deletion import Deletion, DeletionList
from pestifer.objs.link import Link, LinkList
from pestifer.objs.mutation import Mutation, MutationList
from pestifer.molecule.residue import Residue, ResidueList, ResiduePlaceholder
from pestifer.objs.resid import ResID

class TestObjManager(unittest.TestCase):

    def test_objmanager_initialization(self):
        O=ObjManager()
        self.assertEqual(len(O._obj_classes),16)
        sample_yaml = 'mutations'
        self.assertIn(sample_yaml, O._obj_classes_byYAML)
        Cls=O._obj_classes_byYAML[sample_yaml]
        self.assertEqual(sample_yaml, Cls._yaml_header)
        self.assertEqual(Cls, Mutation)
        LCls=O._objlist_classes_byYAML[sample_yaml]
        self.assertEqual(LCls, MutationList)
        self.assertEqual(LCls._item_type, Mutation)

    def test_objmanager_ingest_obj(self):
        O=ObjManager()
        mut = Mutation('C:THR,154,GLY')
        O.ingest(mut)
        self.assertIn('seq', O)
        self.assertIn('mutations', O['seq'])
        self.assertIsInstance(O['seq']['mutations'], MutationList)

        insertion = Insertion('C,154A,GGG')
        O.ingest(insertion)
        self.assertIn('seq', O)
        self.assertIn('insertions', O['seq'])
        self.assertIsInstance(O['seq']['insertions'], InsertionList)

        link = Link('C_154_ND2-C_1540_C1')
        O.ingest(link)
        self.assertIn('topol', O)
        self.assertIn('links', O['topol'])
        self.assertIsInstance(O['topol']['links'], LinkList)

    def test_objmanager_ingest_dict(self):
        O=ObjManager()
        input_specs = {
            'mutations': MutationList([Mutation('C:THR,154,GLY'), Mutation('C:ALA,155,VAL')]),
            'insertions': InsertionList([Insertion('C,154A,GGG'), Insertion('C,155B,AAA')]),
            'links': ['C_154_ND2-C_1540_C1']
        }
        O.ingest(input_specs)
        self.assertIn('seq', O)
        self.assertIn('mutations', O['seq'])
        self.assertIsInstance(O['seq']['mutations'], MutationList)
        self.assertEqual(len(O['seq']['mutations']), 2)
        self.assertIn('insertions', O['seq'])
        self.assertIsInstance(O['seq']['insertions'], InsertionList)
        self.assertEqual(len(O['seq']['insertions']), 2)
        self.assertIn('topol', O)
        self.assertIn('links', O['topol'])
        self.assertIsInstance(O['topol']['links'], LinkList)

    def test_objmanager_ingest_objlist(self):
        O=ObjManager()
        mut_list = MutationList([Mutation('C:THR,154,GLY'), Mutation('C:ALA,155,VAL')])
        O.ingest(mut_list)
        self.assertIn('seq', O)
        self.assertIn('mutations', O['seq'])
        self.assertIsInstance(O['seq']['mutations'], MutationList)
        self.assertEqual(len(O['seq']['mutations']), 2)

        insertion_list = InsertionList([Insertion('C,154A,GGG'), Insertion('C,155B,AAA')])
        O.ingest(insertion_list)
        self.assertIn('seq', O)
        self.assertIn('insertions', O['seq'])
        self.assertIsInstance(O['seq']['insertions'], InsertionList)
        self.assertEqual(len(O['seq']['insertions']), 2)

    def test_objmanager_filter_copy(self):
        O=ObjManager()
        for C in 'ABCDE':
            for i in range(5):
                mut = Mutation(f'{C}:THR,{150+i},GLY')
                deletion = Deletion(f"{C}:{10+i}-20")
                O.ingest(mut)
                O.ingest(deletion)
        chainE = O.filter_copy(lambda x: x.chainID == 'E')
        self.assertIn('seq', chainE)
        self.assertIn('mutations', chainE['seq'])
        self.assertIsInstance(chainE['seq']['mutations'], MutationList)
        self.assertEqual(len(chainE['seq']['mutations']), 5)
        self.assertIn('deletions', chainE['seq'])
        self.assertIsInstance(chainE['seq']['deletions'], DeletionList)
        self.assertEqual(len(chainE['seq']['deletions']), 5)

        mutations_only = O.filter_copy(lambda x: isinstance(x, Mutation))
        self.assertIn('seq', mutations_only)
        self.assertIn('mutations', mutations_only['seq'])
        self.assertIsInstance(mutations_only['seq']['mutations'], MutationList)
        self.assertEqual(len(mutations_only['seq']['mutations']), 25)
        self.assertNotIn('deletions', mutations_only['seq'])

    def test_objmanager_expel(self):
        O=ObjManager()
        for C in 'ABCDE':
            for i in range(5):
                mut = Mutation(f'{C}:THR,{150+i},GLY')
                deletion = Deletion(f"{C}:{10+i}-20")
                O.ingest(mut)
                O.ingest(deletion)
        counts = O.counts_by_header()
        self.assertEqual(counts['mutations'], 25)
        self.assertEqual(counts['deletions'], 25)
        R=ResidueList()
        for C in 'ABCDE':
            R.append(Residue(ResiduePlaceholder(f'{C}:T151')))
        self.assertEqual(len(R), 5)
        self.assertEqual(R[0].resname, 'THR')
        self.assertEqual(R[0].resid, ResID(151))

        ex_om=O.expel(R)
        self.assertIsInstance(ex_om, ObjManager)
        # all mutations of T151 from all chains should be expelled
        self.assertIn('seq', ex_om)
        self.assertIn('mutations', ex_om['seq'])
        self.assertIsInstance(ex_om['seq']['mutations'], MutationList)
        expelled_counts = ex_om.counts_by_header()
        self.assertEqual(expelled_counts['mutations'], 5)
        post_expel_counts = O.counts_by_header()
        self.assertEqual(post_expel_counts['mutations'], 20)
        self.assertEqual(post_expel_counts['deletions'], 25)
