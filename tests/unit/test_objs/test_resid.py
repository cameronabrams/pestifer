import unittest
from pestifer.objs.resid import ResID, ResIDList

class TestResID(unittest.TestCase):
    def test_resid_creation(self):
        resid = ResID(resseqnum=123, insertion='A')
        self.assertIsInstance(resid, ResID)
        self.assertEqual(resid.resseqnum, 123)
        self.assertEqual(resid.insertion, 'A')
        self.assertEqual(repr(resid), "ResID(resseqnum=123, insertion='A')")

    def test_resid_from_args(self):
        resid = ResID(123, 'A')
        self.assertIsInstance(resid, ResID)
        self.assertEqual(resid.resseqnum, 123)
        self.assertEqual(resid.insertion, 'A')
        self.assertEqual(repr(resid), "ResID(resseqnum=123, insertion='A')")
        self.assertEqual(resid.resid, '123A')

    def test_resid_from_str(self):
        resid = ResID('123A')
        self.assertEqual(resid.resseqnum, 123)
        self.assertEqual(resid.insertion, 'A')
        self.assertEqual(resid.resid, '123A')

    def test_resid_relational(self):
        resid1 = ResID(123, 'A')
        resid2 = ResID(123, 'B')
        resid3 = ResID(124, None)

        self.assertTrue(resid1 < resid2)
        self.assertTrue(resid1 < resid3)
        self.assertTrue(resid2 < resid3)

    def test_resid_resrange(self):
        res_range = ResIDList('123-125A')
        self.assertEqual(len(res_range), 2)
        self.assertEqual(res_range[0].resid, 123)
        self.assertEqual(res_range[1].resid, '125A')

        res_range = ResIDList('123-125A#125C')
        self.assertEqual(len(res_range), 3)
        self.assertEqual(res_range[0].resid, 123)
        self.assertEqual(res_range[1].resid, '125A')
        self.assertEqual(res_range[2].resid, '125C')

        res_range = ResIDList('123C')
        self.assertEqual(len(res_range), 1)
        self.assertEqual(res_range[0].resid, '123C')