"""

.. module:: test_config
   :synopsis: tests pestifer.Config()
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""

import unittest
import platform
import os
import yaml
from pestifer.config import Config, ResourceManager

class ConfigTest(unittest.TestCase):
    def test_resource_manager(self):
        rm=ResourceManager()
        for r in ['charmmff','config','examples','tcl']:
            self.assertTrue(r in rm)
        for cd in ['custom','toppar']:
            self.assertTrue(cd in rm['charmmff'])
            
    def test_config_nouser(self):
        c=Config()
        self.assertTrue('Resources' in c)
        self.assertEqual(type(c['Resources']),ResourceManager)
        rr=c['Resources']['tcl']['scripts']
        print(rr)
        self.assertEqual(os.path.join(rr,'pestifer-vmd.tcl'),c.vmd_startup_script)
        self.assertEqual(c.segtype_by_resname('ALA'),'protein')
        self.assertEqual(c.segtype_by_resname('MAN'),'glycan')
        self.assertEqual(c.segtype_by_resname('CL'),'ion')
        self.assertEqual(c.charmmify_resname('MAN'),'AMAN')
        self.assertEqual(c.charmmify_resname('ALA'),'ALA')
        self.assertEqual(c.res_123('A'),'ALA')
        self.assertEqual(c.res_321('PHE'),'F')

    def test_config_boolean(self):
        c=Config()
        self.assertTrue(c)

    def test_config_user(self):
        rm=ResourceManager()
        configfile=rm['examples']+'/example1.yaml'
        c=Config(userconfigfile=configfile)
        self.assertTrue('user' in c)
        self.assertTrue('tasks' in c['user'])
        self.assertEqual(len(c['user']['tasks'][0].keys()),1)
        task1=c['user']['tasks'][0]
        name,specs=[(x,y) for x,y in task1.items()][0]
        self.assertEqual(name,'psfgen')
        self.assertTrue('source' in specs)
        self.assertEqual(specs['source']['rcsb'],'6pti')
        self.assertTrue('minimize' in specs)
        self.assertTrue('cleanup' in specs)
        task2=c['user']['tasks'][1]
        name,specs=[(x,y) for x,y in task2.items()][0]
        self.assertEqual(name,'solvate')
        task3=c['user']['tasks'][2]
        name,specs=[(x,y) for x,y in task3.items()][0]
        self.assertEqual(name,'relax')
