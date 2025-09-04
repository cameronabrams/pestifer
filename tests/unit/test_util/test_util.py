from pestifer.util.util import *
import unittest
from pestifer.tasks.psfgen import *
import os
from pestifer.objs.mutation import Mutation
from pestifer.objs.patch import Patch
from pestifer import objs
from pestifer.core.command import Command
from pestifer.util.util import protect_str_arg

import logging
logger = logging.getLogger(__name__)

class TestUtil(unittest.TestCase):
    def test_str_arg(self):
        with open('test_str_arg.tcl','w') as f:
            f.write("""
pestifer_init

proc abc { arg } {
    puts "arg: $arg"
}
catch {set a "hello world"}
for { set i 0 } { $i < [llength $argv] } { incr i } {
   if { [lindex $argv $i] == "-a"} {
      incr i
      set a [lindex $argv $i]
   }
}
set a_arg [deprotect_str_arg $a]
vmdcon -info $a_arg
exit
"""
                    )
        arg='abc 123'
        c=Command(f"vmd -dispdev text -e test_str_arg.tcl -args -a {arg}")
        c.run()
        self.assertFalse('abc 123' in c.stdout)
        c=Command(f"vmd -dispdev text -e test_str_arg.tcl -args -a {protect_str_arg(arg)}")
        c.run()
        logger.debug(f'Output: {c.stdout}')
        self.assertTrue('abc 123' in c.stdout)
        os.remove('test_str_arg.tcl')
       
    def test_reduce_intlist(self):
        l=[1,2,3,4,5,7,8,9,10,12]
        L=reduce_intlist(l)
        self.assertEqual(L,'1 to 5 7 to 10 12')
        l=list(range(100))
        l.remove(23)
        l.remove(73)
        l.remove(74)
        l.remove(75)
        L=reduce_intlist(l)
        self.assertEqual(L,'0 to 22 24 to 72 76 to 99')
        l=[1,2]
        L=reduce_intlist(l)
        self.assertEqual(L,'1 2')

    def test_replace(self):
        starting_dict={
            'SomeIntegers':[1,2,'$(VARIABLE3)'],
            'SomeWords':{'subdata':['$(VARIABLE2)','hello']},
            'YetMoreData':[0.1,'$(VARIABLE4)'],
            'NestyMcNestface':{'level1':{'level2':{'ALIST':'$(VAR6)'}},'level1-str':'$(VAR5)'},
            'SearchReplace': '$(VAR7)!!!'
        }

        expected_dict={
            'SomeIntegers':[1,2,3],
            'SomeWords':{'subdata':['word','hello']},
            'YetMoreData':[0.1,0.99],
            'NestyMcNestface':{'level1':{'level2':{'ALIST':['this','is','a','list']}},'level1-str':'A string'},
            'SearchReplace':'replaced!!!'
        }

        replacements={
            'VARIABLE1':'SomeIntegers',
            'VARIABLE3':3,
            'VARIABLE2':'word',
            'VARIABLE4':0.99,
            'VAR5':'A string',
            'VAR6': ['this','is','a','list'],
            'VAR7': 'replaced'
        }
        for s,r in replacements.items():
            replace(starting_dict,s,r)
        self.assertEqual(starting_dict,expected_dict)
    
    def test_is_periodic(self):
        self.assertFalse(is_periodic('no.xsc'))
        self.assertTrue(is_periodic('yes.xsc'))

    def test_cell_from_xsc(self):
        box,orig=cell_from_xsc('test.xsc')
        self.assertTrue(box.shape==(3,3))
        self.assertTrue(orig.shape==(3,))
        self.assertEqual(box[0][0],191.901209569)
        self.assertEqual(box[0][1],0)
        self.assertEqual(box[0][2],0)
        self.assertEqual(box[1][1],191.901209569)
        self.assertEqual(box[1][0],0)
        self.assertEqual(box[1][2],0)
        self.assertEqual(box[2][2],206.131438627)
        self.assertEqual(box[2][0],0)
        self.assertEqual(box[2][1],0)
        self.assertEqual(orig[0],0)
        self.assertEqual(orig[1],0)
        self.assertEqual(orig[2],0)
