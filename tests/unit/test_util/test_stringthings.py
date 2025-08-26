import unittest
from pestifer.util.stringthings import ByteCollector, FileCollector, to_latex_math, parse_filter_expression

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

    def test_byte_collector_reassign(self):
        bc=ByteCollector()
        bc.ingest_file('reassigntest.sh')
        bc.reassign('MYVAR1',64)
        bc.reassign('EXMYVAR2',63)
        bc.reassign('ntasks',8,style='SLURM')
        bc.reassign('p','debug',style='SLURM')
        bc.write_file('test.sh')
        with open('test.sh','r') as f:
            lines=f.read().split('\n')
        self.assertTrue('MYVAR1=64' in lines)
        self.assertTrue('export EXMYVAR2=63' in lines)
        self.assertTrue('#SBATCH --ntasks=8' in lines)
        self.assertTrue('#SBATCH -p debug' in lines)

    def test_file_collector(self):
        fc=FileCollector()
        fc.append('delete_me.xyz')
        fc.append('and_me.abc')
        self.assertEqual(len(fc),2)
        fc.flush()

    def test_to_latex_math(self):
        s='a_x'
        res=to_latex_math(s)
        self.assertEqual(res,r'a$_{x}$')

    def test_parse_filter_expression(self):
        class MockObject:
            def __init__(self, **kwargs):
                self.__dict__.update(kwargs)

        expr = "segtype == 'protein' and resid > 100"
        filter_func = parse_filter_expression(expr)
        self.assertTrue(filter_func(MockObject(segtype='protein', resid=150)))
        self.assertFalse(filter_func(MockObject(segtype='lipid', resid=50)))

        expr = "(segtype == 'protein' and resid > 100) or (segtype == 'lipid' and resid < 50)"
        filter_func = parse_filter_expression(expr)
        self.assertTrue(filter_func(MockObject(segtype='protein', resid=150)))
        self.assertTrue(filter_func(MockObject(segtype='lipid', resid=25)))

        expr = "resname == 'PO4'"
        filter_func = parse_filter_expression(expr)
        self.assertTrue(filter_func(MockObject(resname='PO4')))
        self.assertFalse(filter_func(MockObject(resname='PO5')))