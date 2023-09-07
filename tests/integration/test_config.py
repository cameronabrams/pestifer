"""

.. module:: test_config
   :synopsis: tests pestifer.Config()
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""

import unittest
import glob
import os
import yaml
from pestifer.config import Config, ResourceManager, userhelp, segtype_of_resname, charmm_resname_of_pdb_resname, res_123, res_321

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
        rr=c['Resources']['tcl']['scripts']
        self.assertEqual(os.path.join(rr,'pestifer-vmd.tcl'),c.vmd_startup_script)
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

    def test_config_help(self):
        c=Config()
        with open('help.out','w') as f:
            L=c['help']['directives']
            userhelp(L,f.write,end='\n')
            userhelp(L,f.write,'title',end='\n')
            userhelp(L,f.write,'paths',end='\n')
            userhelp(L,f.write,'paths','vmd',end='\n')
            userhelp(L,f.write,'tasks',end='\n')
            userhelp(L,f.write,'tasks','psfgen',end='\n')
            userhelp(L,f.write,'tasks','psfgen','source',end='\n')
            userhelp(L,f.write,'tasks','psfgen','source','id',end='\n')
            userhelp(L,f.write,'tasks','psfgen','source','file_format',end='\n')
            userhelp(L,f.write,'tasks','psfgen','source','sequence',end='\n')
            userhelp(L,f.write,'tasks','psfgen','source','sequence','fix_conflicts',end='\n')
            userhelp(L,f.write,'tasks','psfgen','source','sequence','loops',end='\n')
            userhelp(L,f.write,'tasks','psfgen','source','sequence','loops','min_loop_length',end='\n')
            userhelp(L,f.write,'tasks','ligate',end='\n')
            userhelp(L,f.write,'tasks','ligate','steer',end='\n')
            userhelp(L,f.write,'tasks','solvate',end='\n')
            userhelp(L,f.write,'tasks','relax',end='\n')
            userhelp(L,f.write,'tasks','relax','ensemble',end='\n')
            userhelp(L,f.write,'tasks','terminate',end='\n')
        self.assertTrue(os.path.isfile('help.out'))
            

    def test_config_user(self):
        rm=ResourceManager()
        configfile=rm['examples']+'/01-bpti.yaml'
        c=Config(userconfigfile=configfile)
        self.assertTrue('user' in c)
        self.assertTrue('tasks' in c['user'])
        self.assertEqual(len(c['user']['tasks'][0].keys()),1)
        task1=c['user']['tasks'][0]
        name,specs=[(x,y) for x,y in task1.items()][0]
        self.assertEqual(name,'psfgen')
        self.assertTrue('source' in specs)
        self.assertEqual(specs['source']['id'],'6pti')
        self.assertTrue('minimize' in specs)
        self.assertTrue('cleanup' in specs)
        task2=c['user']['tasks'][1]
        name,specs=[(x,y) for x,y in task2.items()][0]
        self.assertEqual(name,'solvate')
        task3=c['user']['tasks'][2]
        name,specs=[(x,y) for x,y in task3.items()][0]
        self.assertEqual(name,'relax')

    def test_config_help_examples(self):
        rm=ResourceManager()
        configfiles=glob.glob(rm['examples']+'/*.yaml')
        for configfile in configfiles:
            d,n=os.path.split(configfile)
            b,e=os.path.splitext(n)
            c=Config(userconfigfile=configfile)
            self.assertTrue('base' in c)
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
