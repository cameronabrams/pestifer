"""

.. module:: test_config
   :synopsis: tests pestifer.Config()
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""

import unittest
import glob
import os
import yaml
from pestifer.config import Config, ResourceManager, segtype_of_resname, charmm_resname_of_pdb_resname, res_123, res_321

class ConfigTest(unittest.TestCase):
    def test_resource_manager(self):
        rm=ResourceManager()
        for r in ['charmmff','config','examples','tcl']:
            self.assertTrue(r in rm)
        for cd in ['custom','toppar']:
            self.assertTrue(cd in rm['charmmff'])
            
    def test_config_nouser_globals(self):
        c=Config()
        self.assertTrue('Resources' in c)
        self.assertEqual(type(c['Resources']),ResourceManager)
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
        rm=ResourceManager()
        configfile=rm['examples']+'/01-bpti.yaml'
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
        rm=ResourceManager()
        configfile=os.path.join(rm['examples'],'05-hiv-env-4zmj-package.yaml')
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
        rm=ResourceManager()
        configfiles=glob.glob(rm['examples']+'/*.yaml')
        for configfile in configfiles:
            d,n=os.path.split(configfile)
            b,e=os.path.splitext(n)
            c=Config(userfile=configfile)
            self.assertTrue('base' in c)
            self.assertTrue('user' in c)
            D=c['base']['directives']
            self.assertEqual(len(D),6)
            tld=[x['name'] for x in D]
            self.assertEqual(tld,['charmmff', 'psfgen', 'namd2', 'title', 'paths', 'tasks'])
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
            self.assertEqual(paths['namd2'],'/usr/local/bin/namd2')
            self.assertEqual(paths['charmrun'],'/usr/local/bin/charmrun')
            self.assertEqual(paths['vmd'],'/usr/local/bin/vmd')
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
