import unittest
from pestifer.stringthings import ByteCollector, FileCollector
class TestStringthings(unittest.TestCase):
    def test_byte_collector(self):
        bc=ByteCollector()
        bc.write('hello')
        a_string=bc.byte_collector
        self.assertEqual(a_string,'hello')
    def test_file_collector(self):
        fc=FileCollector()
        fc.append('delete_me.xyz')
        fc.append('and_me.abc')
        self.assertEqual(len(fc),2)
        fc.flush()
