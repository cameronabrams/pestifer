# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import unittest
import os
import shutil
import yaml
from pestifer.core.config import Config
from pestifer.core.labels import Labels #segtype_of_resname, charmm_resname_of_pdb_resname, res_123, res_321
from pestifer.core.resourcemanager import ResourceManager

class ConfigTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.RM=ResourceManager()
        cls.c=Config(quiet=True)

    def test_labels(self):
        c=self.c
        self.assertEqual(os.path.join(c.tcl_root,'vmdrc.tcl'),c.vmd_startup_script)
        self.assertEqual(Labels.segtype_of_resname['ALA'],'protein')
        self.assertEqual(Labels.segtype_of_resname['MAN'],'glycan')
        self.assertEqual(Labels.segtype_of_resname['CL'],'ion')
        self.assertEqual(Labels.charmm_resname_of_pdb_resname.get('MAN'),'AMAN')
        self.assertEqual(Labels.charmm_resname_of_pdb_resname.get('ALA'),None)
        self.assertEqual(Labels.res_123['A'],'ALA')
        self.assertEqual(Labels.res_321['PHE'],'F')

    def test_config_boolean(self):
        c=self.c
        self.assertTrue(c)

    def test_config_userdict(self):
        C=self.c
        ud=C['user']
        D=Config(userdict=ud)
        self.assertTrue('user' in D)
        
    def test_config_user(self):
        tmpdir='__test_config_user'
        if os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)
        os.mkdir(tmpdir)
        os.chdir(tmpdir)
        EM=self.RM.example_manager
        configfile=EM.checkout_example(1)
        c=Config(userfile=configfile)
        self.assertTrue('user' in c)
        self.assertTrue('tasks' in c['user'])
        self.assertEqual(len(c['user']['tasks'][0].keys()),1)
        task1=c['user']['tasks'][0]
        name,specs=[(x,y) for x,y in task1.items()][0]
        self.assertEqual(name,'psfgen')
        self.assertTrue('source' in specs)
        self.assertEqual(specs['source']['id'],'6pti')
        task2=c['user']['tasks'][1]
        name,specs=[(x,y) for x,y in task2.items()][0]
        self.assertEqual(name,'md')
        task3=c['user']['tasks'][2]
        name,specs=[(x,y) for x,y in task3.items()][0]
        self.assertEqual(name,'solvate')
        os.chdir('..')

    def test_config_task_source(self):
        tmpdir='__test_config_task_source'
        if os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)
        os.mkdir(tmpdir)
        os.chdir(tmpdir)
        EM=self.RM.example_manager
        configfile=EM.checkout_example(7)
        c=Config(userfile=configfile)
        self.assertTrue('user' in c)
        self.assertTrue('tasks' in c['user'])
        self.assertEqual(len(c['user']['tasks'][0].keys()),1)
        task1=c['user']['tasks'][0]
        name,specs=[(x,y) for x,y in task1.items()][0]
        self.assertEqual(name,'psfgen')
        self.assertTrue('source' in specs)
        self.assertEqual(specs['source']['id'],'4zmj')
        sourcespecs=specs['source']
        self.assertTrue('transform_reserves' in sourcespecs)
        resm=sourcespecs['transform_reserves']
        self.assertTrue('G' in resm)
        self.assertTrue('B' in resm)
        self.assertTrue(resm['G']==['H','J'])
        self.assertTrue(resm['B']==['C','D'])
        os.chdir('..')

    def test_config_help_examples(self):
        EM=self.RM.example_manager
        configfolders=[x.name for x in EM.examples_list]
        configfiles=[f'{x.name}.yaml' for x in EM.examples_list]
        tmpdir='__test_config_help_examples'
        if os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)
        os.mkdir(tmpdir)
        os.chdir(tmpdir)
        for index,(configfolder,configfile) in enumerate(zip(configfolders,configfiles)):
            configfile_path=os.path.join(EM.path,configfolder,configfile)
            self.assertTrue(os.path.exists(configfile_path),f'Config file {configfile_path} does not exist')
            EM.checkout_example(index+1)
            self.assertTrue(os.path.exists(configfile),f'Config file {configfile} does not exist after checkout')
            c=Config(userfile=configfile_path)
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
            with open(f'{configfile}-complete.yaml','w') as f:
                yaml.dump(U,f)
            tasks=U['tasks']
            self.assertTrue('psfgen' in tasks[0])
            specs=tasks[0]['psfgen']
            self.assertTrue('source' in specs)
            self.assertTrue('mods' in specs)
            self.assertTrue('ssbondsdelete' in specs['mods'])
            source_specs=specs['source']
            self.assertTrue('id' in source_specs or 'alphafold' in source_specs)
        os.chdir('..')
