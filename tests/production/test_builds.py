# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import unittest
from pestifer.config import Config
import os
import shutil
from glob import glob
import subprocess
import yaml
        
def cleanup(excl_exts=['log','png','csv']):
    excl=set()
    for ex in excl_exts:
        excl=excl.union(set(glob(f'*.{ex}')))
    rmv=list(set(glob('*'))-excl)
    for r in rmv:
        os.remove(r)

def do_it(exnumber):
    nstr=f'{exnumber:02d}'
    c=Config()
    ex_path=c['Resources']['examples']
    ex_yaml=glob(ex_path+f'/{nstr}*.yaml')
    assert len(ex_yaml)==1,f'No example {exnumber} found'
    ex_yaml=ex_yaml[0]
    base=os.path.split(ex_yaml)[1]
    rebase=base[len(nstr)+1:]
    d=str(exnumber)
    if os.path.exists(d):
        shutil.rmtree(d)
    os.mkdir(d)
    os.chdir(d)
    shutil.copy(ex_yaml,f'./{rebase}')
    with open(rebase,'r') as f:
        specs=yaml.safe_load(f)
    prod_basename=specs['tasks'][-1]['terminate']['package']['basename']
    prod_tgz=f'{prod_basename}.tgz'
    cmd=f'pestifer run {rebase}'
    res=subprocess.run(cmd,shell=True,executable='/bin/bash',check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    with open('console.log','w') as f:
        f.write(res.stdout)
    assert(os.path.exists(prod_tgz))
        # end_message=res.stderr.split('\n')[-2]
        # assert(end_message=='INFO> pestifer runtime ends.')
        # cleanup()

class TestBuild(unittest.TestCase):
    def test_example_build01(self):
        do_it(1)
    def test_example_build02(self):
        do_it(2)
    def test_example_build03(self):
        do_it(3)
    def test_example_build04(self):
        do_it(4)
    def test_example_build05(self):
        do_it(5)
    def test_example_build06(self):
        do_it(6)
    def test_example_build07(self):
        do_it(7)
    def test_example_build08(self):
        do_it(8)
    def test_example_build09(self):
        do_it(9)
    def test_example_build10(self):
        do_it(10)
    def test_example_build11(self):
        do_it(11)
    def test_example_build12(self):
        do_it(12)
