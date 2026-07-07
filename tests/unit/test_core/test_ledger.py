import tempfile
import unittest

from pestifer.core import ledger


class TestLedger(unittest.TestCase):
    def test_append_read_ids_and_fields(self):
        d = tempfile.mkdtemp()
        self.assertEqual(ledger.read(d), [])                     # missing ledger -> empty
        e1 = ledger.append(d, category='charmmff', verb='add-residue',
                           summary='add MYLIG', files=['a/b.str'], author='X <x@x>',
                           branch='modpkg/x', timestamp='2026-07-07T00:00:00+00:00')
        e2 = ledger.append(d, category='example', verb='add', summary='add ex 19',
                           timestamp='2026-07-07T00:01:00+00:00')
        self.assertEqual(e1['id'], 1)
        self.assertEqual(e2['id'], 2)                            # ids increment
        self.assertEqual(e1['schema'], ledger.SCHEMA_VERSION)
        self.assertEqual(e1['reverted'], False)
        self.assertEqual(e1['files'], ['a/b.str'])
        got = ledger.read(d)
        self.assertEqual([e['summary'] for e in got], ['add MYLIG', 'add ex 19'])
        self.assertEqual(ledger.get(d, 2)['verb'], 'add')
        self.assertIsNone(ledger.get(d, 99))

    def test_write_all_marks_reverted(self):
        d = tempfile.mkdtemp()
        ledger.append(d, category='charmmff', verb='add-residue', summary='add MYLIG',
                      timestamp='2026-07-07T00:00:00+00:00')
        entries = ledger.read(d)
        entries[0]['reverted'] = True
        entries[0]['reverted_by'] = 2
        ledger.write_all(d, entries)
        self.assertTrue(ledger.get(d, 1)['reverted'])
        self.assertEqual(ledger.get(d, 1)['reverted_by'], 2)

    def test_format_filters_and_status(self):
        d = tempfile.mkdtemp()
        ledger.append(d, category='charmmff', verb='add-residue', summary='add MYLIG',
                      timestamp='2026-07-07T00:00:00+00:00')
        ledger.append(d, category='example', verb='add', summary='add ex',
                      timestamp='2026-07-07T00:01:00+00:00')
        out = ledger.format_entries(ledger.read(d))
        self.assertIn('#1', out)
        self.assertIn('charmmff/add-residue', out)
        self.assertIn('example/add', out)
        # category filter
        only = ledger.format_entries(ledger.read(d), category='example')
        self.assertNotIn('add-residue', only)
        self.assertIn('example/add', only)
        # empty
        self.assertEqual(ledger.format_entries([]), '(no ledger entries)')


if __name__ == '__main__':
    unittest.main()
