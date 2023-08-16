import unittest
from pestifer.stringthings import ByteCollector, FileCollector
class TestStringthings(unittest.TestCase):
    def test_byte_collector_write(self):
        bc=ByteCollector()
        bc.write('hello')
        a_string=bc.byte_collector
        self.assertEqual(a_string,'hello')
    def test_byte_collector_addline(self):
        bc=ByteCollector()
        bc.addline('this is a line')
        self.assertEqual('this is a line\n',str(bc))
    def test_byte_collector_comment(self):
        bc=ByteCollector()
        bc.comment('this is a comment')
        self.assertEqual('# this is a comment\n',str(bc))
        bc.reset()
        bc.comment('this is a multiline comment. What is the deal with potato chips? I mean you\'re crispy, you\'re greasy... just make up your mind, potato chips!')
        self.assertEqual('# this is a multiline comment. What is the deal with potato chips? I mean you\'re\n# crispy, you\'re greasy... just make up your mind, potato chips!\n',str(bc))
    def test_file_collector(self):
        fc=FileCollector()
        fc.append('delete_me.xyz')
        fc.append('and_me.abc')
        self.assertEqual(len(fc),2)
        fc.flush()
