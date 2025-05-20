import unittest
from pestifer.stringthings import ByteCollector, FileCollector, ri_range, to_latex_math
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
        bc.injest_file('reassigntest.sh')
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
    def test_ri_range(self):
        test_string='123'
        rs=ri_range(test_string)
        self.assertTrue(len(rs)==1)
        r1,i1=rs[0]
        self.assertEqual(r1,123)
        self.assertEqual(i1,'')
        test_string='123A#123B#144C-199'
        rs=ri_range(test_string)
        self.assertTrue(len(rs)==4)
        r1,i1=rs[0]
        r2,i2=rs[1]
        r3,i3=rs[2]
        r4,i4=rs[3]
        self.assertEqual(r1,123)
        self.assertEqual(i1,'A')
        self.assertEqual(r2,123)
        self.assertEqual(i2,'B')
        self.assertEqual(r3,144)
        self.assertEqual(i3,'C')
        self.assertEqual(r4,199)
        self.assertEqual(i4,'')

    def test_to_latex_math(self):
        s='a_x'
        res=to_latex_math(s)
        self.assertEqual(res,r'a$_{x}$')
