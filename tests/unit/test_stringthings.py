import unittest
from pestifer.stringthings import ByteCollector
class TestStringthings(unittest.TestCase):
    def test_byte_collector(self):
        bc=ByteCollector()
        bc.write('hello')
        a_string=bc.byte_collector
        self.assertEqual(a_string,'hello')
