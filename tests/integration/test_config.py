"""

.. module:: test_config
   :synopsis: tests pestifer.Config()
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""

import unittest
import platform
import os
import yaml
from pestifer.config import Config, ResourceManager, dwalk

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

    def test_config_help_example1(self):
        rm=ResourceManager()
        configfile=rm['examples']+'/example1.yaml'
        c=Config(userconfigfile=configfile)
        self.assertTrue('user' in c)
        self.assertTrue('help' in c)
        H=c['help']
        D=H['directives']
        self.assertEqual(len(D),3)
        tld=[x['name'] for x in D]
        self.assertEqual(tld,['title','paths','tasks'])
        T_idx=tld.index('title')
        T=D[T_idx]
        self.assertEqual(T['name'],'title')
        self.assertEqual(T['type'],'str')
        self.assertTrue('text' in T)
        self.assertEqual(T['default'],'Pestifer')
        U=c['user']
        self.assertTrue('title' in U)
        self.assertFalse('paths' in U)
        self.assertTrue('tasks' in U)
        dwalk(H,U)
        utitle=U['title']
        self.assertEqual(utitle,'BPTI')
        paths=U['paths']
        self.assertEqual(paths['namd2'],'/usr/local/bin/namd2')
        self.assertEqual(paths['charmrun'],'/usr/local/bin/charmrun')
        self.assertEqual(paths['vmd'],'/usr/local/bin/vmd')
        self.assertFalse('charmff' in paths)
        with open('example1_complete.yaml','w') as f:
            yaml.dump(U,f)
        tasks=U['tasks']
        print(tasks)
        self.assertTrue('psfgen' in tasks[0])
        specs=tasks[0]['psfgen']
        self.assertTrue('source' in specs)
        source_specs=specs['source']
        self.assertTrue('id' in source_specs)

    def test_config_help_example2(self):
        rm=ResourceManager()
        configfile=rm['examples']+'/example2.yaml'
        c=Config(userconfigfile=configfile)
        self.assertTrue('user' in c)
        self.assertTrue('help' in c)
        H=c['help']
        D=H['directives']
        self.assertEqual(len(D),3)
        tld=[x['name'] for x in D]
        self.assertEqual(tld,['title','paths','tasks'])
        T_idx=tld.index('title')
        T=D[T_idx]
        self.assertEqual(T['name'],'title')
        self.assertEqual(T['type'],'str')
        self.assertTrue('text' in T)
        self.assertEqual(T['default'],'Pestifer')
        U=c['user']
        self.assertTrue('title' in U)
        self.assertFalse('paths' in U)
        self.assertTrue('tasks' in U)
        dwalk(H,U)
        utitle=U['title']
        self.assertEqual(utitle,'BPTI')
        paths=U['paths']
        self.assertEqual(paths['namd2'],'/usr/local/bin/namd2')
        self.assertEqual(paths['charmrun'],'/usr/local/bin/charmrun')
        self.assertEqual(paths['vmd'],'/usr/local/bin/vmd')
        self.assertFalse('charmff' in paths)
        with open('example2_complete.yaml','w') as f:
            yaml.dump(U,f)
        tasks=U['tasks']
        print(tasks)
        self.assertTrue('psfgen' in tasks[0])
        specs=tasks[0]['psfgen']
        self.assertTrue('source' in specs)
        source_specs=specs['source']
        self.assertTrue('id' in source_specs)

    def test_config_help_example3(self):
        rm=ResourceManager()
        configfile=rm['examples']+'/example3.yaml'
        c=Config(userconfigfile=configfile)
        self.assertTrue('user' in c)
        self.assertTrue('help' in c)
        H=c['help']
        D=H['directives']
        self.assertEqual(len(D),3)
        tld=[x['name'] for x in D]
        self.assertEqual(tld,['title','paths','tasks'])
        T_idx=tld.index('title')
        T=D[T_idx]
        self.assertEqual(T['name'],'title')
        self.assertEqual(T['type'],'str')
        self.assertTrue('text' in T)
        self.assertEqual(T['default'],'Pestifer')
        U=c['user']
        self.assertTrue('title' in U)
        self.assertFalse('paths' in U)
        self.assertTrue('tasks' in U)
        dwalk(H,U)
        utitle=U['title']
        self.assertEqual(utitle,'BPTI')
        paths=U['paths']
        self.assertEqual(paths['namd2'],'/usr/local/bin/namd2')
        self.assertEqual(paths['charmrun'],'/usr/local/bin/charmrun')
        self.assertEqual(paths['vmd'],'/usr/local/bin/vmd')
        self.assertFalse('charmff' in paths)
        with open('example3_complete.yaml','w') as f:
            yaml.dump(U,f)
        tasks=U['tasks']
        print(tasks)
        self.assertTrue('psfgen' in tasks[0])
        specs=tasks[0]['psfgen']
        self.assertTrue('source' in specs)
        source_specs=specs['source']
        self.assertTrue('id' in source_specs)

    def test_config_help_example4(self):
        rm=ResourceManager()
        configfile=rm['examples']+'/example4.yaml'
        c=Config(userconfigfile=configfile)
        self.assertTrue('user' in c)
        self.assertTrue('help' in c)
        H=c['help']
        D=H['directives']
        self.assertEqual(len(D),3)
        tld=[x['name'] for x in D]
        self.assertEqual(tld,['title','paths','tasks'])
        T_idx=tld.index('title')
        T=D[T_idx]
        self.assertEqual(T['name'],'title')
        self.assertEqual(T['type'],'str')
        self.assertTrue('text' in T)
        self.assertEqual(T['default'],'Pestifer')
        U=c['user']
        self.assertTrue('title' in U)
        self.assertFalse('paths' in U)
        self.assertTrue('tasks' in U)
        dwalk(H,U)
        utitle=U['title']
        self.assertEqual(utitle,'Example 4 HIV-1 Env Trimer 4zmj')
        paths=U['paths']
        self.assertEqual(paths['namd2'],'/usr/local/bin/namd2')
        self.assertEqual(paths['charmrun'],'/usr/local/bin/charmrun')
        self.assertEqual(paths['vmd'],'/usr/local/bin/vmd')
        self.assertFalse('charmff' in paths)
        with open('example4_complete.yaml','w') as f:
            yaml.dump(U,f)
        tasks=U['tasks']
        print(tasks)
        self.assertTrue('psfgen' in tasks[0])
        specs=tasks[0]['psfgen']
        self.assertTrue('source' in specs)
        source_specs=specs['source']
        self.assertTrue('id' in source_specs)

    def test_config_help_example5(self):
        rm=ResourceManager()
        configfile=rm['examples']+'/example5.yaml'
        c=Config(userconfigfile=configfile)
        self.assertTrue('user' in c)
        self.assertTrue('help' in c)
        H=c['help']
        D=H['directives']
        self.assertEqual(len(D),3)
        tld=[x['name'] for x in D]
        self.assertEqual(tld,['title','paths','tasks'])
        T_idx=tld.index('title')
        T=D[T_idx]
        self.assertEqual(T['name'],'title')
        self.assertEqual(T['type'],'str')
        self.assertTrue('text' in T)
        self.assertEqual(T['default'],'Pestifer')
        U=c['user']
        self.assertTrue('title' in U)
        self.assertFalse('paths' in U)
        self.assertTrue('tasks' in U)
        dwalk(H,U)
        utitle=U['title']
        self.assertEqual(utitle,'HIV-1 Env Trimer 4tvp, Fabs and phosphate ions removed')
        paths=U['paths']
        self.assertEqual(paths['namd2'],'/usr/local/bin/namd2')
        self.assertEqual(paths['charmrun'],'/usr/local/bin/charmrun')
        self.assertEqual(paths['vmd'],'/usr/local/bin/vmd')
        self.assertFalse('charmff' in paths)
        with open('example5_complete.yaml','w') as f:
            yaml.dump(U,f)
        tasks=U['tasks']
        print(tasks)
        self.assertTrue('psfgen' in tasks[0])
        specs=tasks[0]['psfgen']
        self.assertTrue('source' in specs)
        source_specs=specs['source']
        self.assertTrue('id' in source_specs)

    def test_config_help_example6(self):
        rm=ResourceManager()
        configfile=rm['examples']+'/example6.yaml'
        c=Config(userconfigfile=configfile)
        self.assertTrue('user' in c)
        self.assertTrue('help' in c)
        H=c['help']
        D=H['directives']
        self.assertEqual(len(D),3)
        tld=[x['name'] for x in D]
        self.assertEqual(tld,['title','paths','tasks'])
        T_idx=tld.index('title')
        T=D[T_idx]
        self.assertEqual(T['name'],'title')
        self.assertEqual(T['type'],'str')
        self.assertTrue('text' in T)
        self.assertEqual(T['default'],'Pestifer')
        U=c['user']
        self.assertTrue('title' in U)
        self.assertFalse('paths' in U)
        self.assertTrue('tasks' in U)
        dwalk(H,U)
        utitle=U['title']
        self.assertEqual(utitle,'HIV-1 Env Trimer 8fad')
        paths=U['paths']
        self.assertEqual(paths['namd2'],'/usr/local/bin/namd2')
        self.assertEqual(paths['charmrun'],'/usr/local/bin/charmrun')
        self.assertEqual(paths['vmd'],'/usr/local/bin/vmd')
        self.assertFalse('charmff' in paths)
        with open('example6_complete.yaml','w') as f:
            yaml.dump(U,f)
        tasks=U['tasks']
        print(tasks)
        self.assertTrue('psfgen' in tasks[0])
        specs=tasks[0]['psfgen']
        self.assertTrue('source' in specs)
        source_specs=specs['source']
        self.assertTrue('id' in source_specs)

    def test_config_help_example7(self):
        rm=ResourceManager()
        configfile=rm['examples']+'/example7.yaml'
        c=Config(userconfigfile=configfile)
        self.assertTrue('user' in c)
        self.assertTrue('help' in c)
        H=c['help']
        D=H['directives']
        self.assertEqual(len(D),3)
        tld=[x['name'] for x in D]
        self.assertEqual(tld,['title','paths','tasks'])
        T_idx=tld.index('title')
        T=D[T_idx]
        self.assertEqual(T['name'],'title')
        self.assertEqual(T['type'],'str')
        self.assertTrue('text' in T)
        self.assertEqual(T['default'],'Pestifer')
        U=c['user']
        self.assertTrue('title' in U)
        self.assertFalse('paths' in U)
        self.assertTrue('tasks' in U)
        dwalk(H,U)
        utitle=U['title']
        self.assertEqual(utitle,'HIV-1 Env Trimer 8fae')
        paths=U['paths']
        self.assertEqual(paths['namd2'],'/usr/local/bin/namd2')
        self.assertEqual(paths['charmrun'],'/usr/local/bin/charmrun')
        self.assertEqual(paths['vmd'],'/usr/local/bin/vmd')
        self.assertFalse('charmff' in paths)
        with open('example7_complete.yaml','w') as f:
            yaml.dump(U,f)
        tasks=U['tasks']
        print(tasks)
        self.assertTrue('psfgen' in tasks[0])
        specs=tasks[0]['psfgen']
        self.assertTrue('source' in specs)
        source_specs=specs['source']
        self.assertTrue('id' in source_specs)
