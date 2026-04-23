# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import unittest

from pestifer.core.objmanager import ObjManager
from pestifer.objs.align import Align, AlignList


class TestAlignObj(unittest.TestCase):

    def test_required_field_only(self):
        a = Align(ref_pdb='ref.pdb')
        self.assertEqual(a.ref_pdb, 'ref.pdb')
        self.assertIsNone(a.ref_psf)
        self.assertEqual(a.mobile_sel, 'all')
        self.assertIsNone(a.ref_sel)
        self.assertEqual(a.apply_to, 'all')

    def test_all_fields(self):
        a = Align(
            ref_pdb='ref.pdb',
            ref_psf='ref.psf',
            mobile_sel='backbone',
            ref_sel='backbone and chain A',
            apply_to='protein',
        )
        self.assertEqual(a.ref_pdb, 'ref.pdb')
        self.assertEqual(a.ref_psf, 'ref.psf')
        self.assertEqual(a.mobile_sel, 'backbone')
        self.assertEqual(a.ref_sel, 'backbone and chain A')
        self.assertEqual(a.apply_to, 'protein')

    def test_yaml_header_and_category(self):
        a = Align(ref_pdb='x.pdb')
        self.assertEqual(a._yaml_header, 'align')
        self.assertEqual(a._objcat, 'coord')

    def test_align_list(self):
        lst = AlignList()
        lst.append(Align(ref_pdb='a.pdb'))
        lst.append(Align(ref_pdb='b.pdb'))
        self.assertEqual(len(lst), 2)

    def test_objmanager_ingest_single(self):
        om = ObjManager()
        om.ingest({'align': [{'ref_pdb': 'ref.pdb', 'mobile_sel': 'backbone'}]})
        aligns = om.get('coord', {}).get('align')
        self.assertIsNotNone(aligns)
        self.assertEqual(len(aligns), 1)
        a = aligns[0]
        self.assertEqual(a.ref_pdb, 'ref.pdb')
        self.assertEqual(a.mobile_sel, 'backbone')
        self.assertIsNone(a.ref_sel)
        self.assertEqual(a.apply_to, 'all')

    def test_objmanager_ingest_multiple(self):
        om = ObjManager()
        om.ingest({'align': [
            {'ref_pdb': 'a.pdb'},
            {'ref_pdb': 'b.pdb', 'ref_psf': 'b.psf', 'mobile_sel': 'name CA'},
        ]})
        aligns = om.get('coord', {}).get('align')
        self.assertEqual(len(aligns), 2)
        self.assertIsNone(aligns[0].ref_psf)
        self.assertEqual(aligns[1].ref_psf, 'b.psf')
        self.assertEqual(aligns[1].mobile_sel, 'name CA')
