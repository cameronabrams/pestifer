# Author: Cameron F. Abrams, <cfa22@drexel.edu>

from pestifer.core.pestifer import run_example
import unittest
import pestifer.resources
from argparse import Namespace
import os
import shutil
import subprocess
import yaml
import pytest

import logging

logger = logging.getLogger(__name__)

def do_it(exnumber):
    nstr=f'{exnumber:02d}'
    if os.path.exists(f'__test_build_{nstr}'):
        shutil.rmtree(f'__test_build_{nstr}')
    os.mkdir(f'__test_build_{nstr}')
    os.chdir(f'__test_build_{nstr}')
    args = Namespace(index=exnumber, config=None, output_dir='./', log_file=f'diagnostics-{nstr}.log', gpu=False, complete_config=False, log_level='debug')
    controller = run_example(args)
    prod_basename = controller.tasks[-1].specs['package']['basename']
    logging.debug(f'prod_basename {prod_basename}')
    prod_tgz = f'{prod_basename}.tar.gz'
    assert(os.path.exists(prod_tgz))
    os.chdir('..')

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
@pytest.mark.slow
def test_example_build15():
    do_it(15)
@pytest.mark.slow
def test_example_build16():
    do_it(16)
@pytest.mark.slow
def test_example_build17():
    do_it(17)
