# Author: Cameron F. Abrams, <cfa22@drexel.edu>
# import unittest
from pestifer.config import Config
from pestifer.resourcemanager import ResourceManager
import os
import shutil
from glob import glob
import subprocess
import yaml
import pytest
        
def cleanup(excl_exts=['log','png','csv']):
    excl=set()
    for ex in excl_exts:
        excl=excl.union(set(glob(f'*.{ex}')))
    rmv=list(set(glob('*'))-excl)
    for r in rmv:
        os.remove(r)

def do_it(exnumber):
    nstr=f'{exnumber:02d}'
    RM=ResourceManager()
    configfile=RM.get_example_yaml_by_index(exnumber)
    configfile_basename=os.path.basename(configfile)
    d=str(exnumber)
    if os.path.exists(d):
        shutil.rmtree(d)
    os.mkdir(d)
    os.chdir(d)
    shutil.copy(configfile,f'./{configfile_basename}')
    with open(configfile_basename,'r') as f:
        specs=yaml.safe_load(f)
    prod_basename=specs['tasks'][-1]['terminate']['package']['basename']
    prod_tgz=f'{prod_basename}.tgz'
    cmd=f'pestifer run {configfile_basename}'
    res=subprocess.run(cmd,shell=True,executable='/bin/bash',check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    with open('console.log','w') as f:
        f.write(res.stdout)
    assert(os.path.exists(prod_tgz))
        # end_message=res.stderr.split('\n')[-2]
        # assert(end_message=='INFO> pestifer runtime ends.')
    cleanup()

@pytest.mark.slow
def test_example_build01():
    do_it(1)
@pytest.mark.slow
def test_example_build02():
    do_it(2)
@pytest.mark.slow
def test_example_build03():
    do_it(3)
@pytest.mark.slow
def test_example_build04():
    do_it(4)
@pytest.mark.slow
def test_example_build05():
    do_it(5)
@pytest.mark.slow
def test_example_build06():
    do_it(6)
@pytest.mark.slow
def test_example_build07():
    do_it(7)
@pytest.mark.slow
def test_example_build08():
    do_it(8)
@pytest.mark.slow
def test_example_build09():
    do_it(9)
@pytest.mark.slow
def test_example_build10():
    do_it(10)
@pytest.mark.slow
def test_example_build11():
    do_it(11)
@pytest.mark.slow
def test_example_build12():
    do_it(12)
@pytest.mark.slow
def test_example_build13():
    do_it(13)
@pytest.mark.slow
def test_example_build14():
    do_it(14)
