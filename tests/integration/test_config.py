"""

.. module:: test_config
   :synopsis: tests pestifer.Config()
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""

import unittest
import glob
import os
import yaml
from pestifer.config import Config, segtype_of_resname, charmm_resname_of_pdb_resname, res_123, res_321
from pestifer.resourcemanager import ResourceManager

class ConfigTest(unittest.TestCase):
            
    def test_config_nouser_globals(self):
        c=Config()
        self.assertEqual(os.path.join(c.tcl_root,'vmdrc.tcl'),c.vmd_startup_script)
        self.assertEqual(segtype_of_resname['ALA'],'protein')
        self.assertEqual(segtype_of_resname['MAN'],'glycan')
        self.assertEqual(segtype_of_resname['CL'],'ion')
        self.assertEqual(charmm_resname_of_pdb_resname.get('MAN'),'AMAN')
        self.assertEqual(charmm_resname_of_pdb_resname.get('ALA'),None)
        self.assertEqual(res_123['A'],'ALA')
        self.assertEqual(res_321['PHE'],'F')

    def test_config_boolean(self):
        c=Config()
        self.assertTrue(c)

    def test_config_user(self):
        RM=ResourceManager()
        configfile=RM.get_example_yaml_by_index(1)
        c=Config(userfile=configfile)
        self.assertTrue('user' in c)
        self.assertTrue('tasks' in c['user'])
        self.assertEqual(len(c['user']['tasks'][0].keys()),1)
        task1=c['user']['tasks'][0]
        name,specs=[(x,y) for x,y in task1.items()][0]
        self.assertEqual(name,'psfgen')
        self.assertTrue('source' in specs)
        self.assertEqual(specs['source']['id'],'6pti')
        self.assertTrue('cleanup' in specs)
        task2=c['user']['tasks'][1]
        name,specs=[(x,y) for x,y in task2.items()][0]
        self.assertEqual(name,'md')
        task3=c['user']['tasks'][2]
        name,specs=[(x,y) for x,y in task3.items()][0]
        self.assertEqual(name,'solvate')

    def test_config_task_source(self):
        RM=ResourceManager()
        configfile=RM.get_example_yaml_by_index(5)
        c=Config(userfile=configfile)
        self.assertTrue('user' in c)
        self.assertTrue('tasks' in c['user'])
        self.assertEqual(len(c['user']['tasks'][0].keys()),1)
        task1=c['user']['tasks'][0]
        name,specs=[(x,y) for x,y in task1.items()][0]
        self.assertEqual(name,'psfgen')
        self.assertTrue('source' in specs)
        self.assertEqual(specs['source']['id'],'4zmj')
        self.assertTrue('cleanup' in specs)
        sourcespecs=specs['source']
        self.assertTrue('transform_reserves' in sourcespecs)
        resm=sourcespecs['transform_reserves']
        self.assertTrue('G' in resm)
        self.assertTrue('B' in resm)
        self.assertTrue(resm['G']==['H','J'])
        self.assertTrue(resm['B']==['C','D'])

    def test_config_help_examples(self):
        RM=ResourceManager()
        configfiles=RM.get_examples_as_list(fullpaths=True)
        for configfile in configfiles:
            d,n=os.path.split(configfile)
            b,e=os.path.splitext(n)
            c=Config(userfile=configfile)
            self.assertTrue('base' in c)
            self.assertTrue('user' in c)
            D=c['base']['directives']
            self.assertEqual(len(D),6)
            tld=[x['name'] for x in D]
            self.assertEqual(tld,['charmmff', 'psfgen', 'namd', 'title', 'paths', 'tasks'])
            T_idx=tld.index('title')
            T=D[T_idx]
            self.assertEqual(T['name'],'title')
            self.assertEqual(T['type'],'str')
            self.assertTrue('text' in T)
            self.assertEqual(T['default'],'Pestifer')
            U=c['user']
            self.assertTrue('title' in U)
            self.assertTrue('paths' in U)
            self.assertTrue('tasks' in U)
            paths=U['paths']
            self.assertEqual(paths['namd3'],'namd3')
            self.assertEqual(paths['charmrun'],'charmrun')
            self.assertEqual(paths['vmd'],'vmd')
            self.assertEqual(paths['packmol'],'packmol')
            self.assertEqual(paths['catdcd'],'catdcd')
            self.assertFalse('charmff' in paths)
            with open(f'{b}-complete.yaml','w') as f:
                yaml.dump(U,f)
            tasks=U['tasks']
            print(tasks)
            self.assertTrue('psfgen' in tasks[0])
            specs=tasks[0]['psfgen']
            self.assertTrue('source' in specs)
            self.assertTrue('mods' in specs)
            self.assertTrue('ssbondsdelete' in specs['mods'])
            source_specs=specs['source']
            self.assertTrue('id' in source_specs)

    def test_config_pdb_depot(self):
        C=Config('test.yaml')
        self.assertEqual(C['user']['paths']['pdb_depot'],'pdb_depot')