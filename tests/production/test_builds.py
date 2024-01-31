# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import unittest
from pestifer.pestifer import run_example
import argparse
import os
import shutil
from glob import glob

def do_it(exnumber):
    args=argparse.Namespace
    args.number=exnumber
    args.loglevel='DEBUG'
    args.diag='diagnostics.log'
    args.no_banner=False
    d=str(exnumber)
    if os.path.exists(d):
        shutil.rmtree(d)
    os.mkdir(d)
    os.chdir(d)
    run_example(args)

def cleanup():
    logs=set(glob('*.log'))
    images=set(glob('*.png'))
    data=set(glob('*.csv'))
    nonlogs=list(set(glob('*')) - logs - images - data)
    for n in nonlogs:
        os.remove(n)

class TestBuild(unittest.TestCase):

    def test_example_build01(self):
        do_it(1)
        cmpfile='prod_6pti.tgz'
        self.assertTrue(os.path.exists(cmpfile))
        cleanup()
    def test_example_build02(self):
        do_it(2)
        cmpfile='prod_6pti.tgz'
        self.assertTrue(os.path.exists(cmpfile))
        cleanup()
    def test_example_build03(self):
        do_it(3)
        cmpfile='prod_6pti.tgz'
        self.assertTrue(os.path.exists(cmpfile))
        cleanup()
    def test_example_build04(self):
        do_it(4)
        cmpfile='prod_6pti.tgz'
        self.assertTrue(os.path.exists(cmpfile))
        cleanup()
    def test_example_build05(self):
        do_it(5)
        cmpfile='prod_4zmj.tgz'
        self.assertTrue(os.path.exists(cmpfile))
        cleanup()
    def test_example_build06(self):
        do_it(6)
        cmpfile='prod_4tvp.tgz'
        self.assertTrue(os.path.exists(cmpfile))
        cleanup()
    def test_example_build07(self):
        do_it(7)
        cmpfile='prod_8fad.tgz'
        self.assertTrue(os.path.exists(cmpfile))
        cleanup()
    def test_example_build08(self):
        do_it(8)
        cmpfile='prod_8fae.tgz'
        self.assertTrue(os.path.exists(cmpfile))
        cleanup()
    def test_example_build09(self):
        do_it(9)
        cmpfile='prod_7txd.tgz'
        self.assertTrue(os.path.exists(cmpfile))
        cleanup()
    def test_example_build10(self):
        do_it(10)
        cmpfile='prod_5vn3.tgz'
        self.assertTrue(os.path.exists(cmpfile))
        cleanup()
    def test_example_build05(self):
        do_it(5)
        cmpfile='prod_4zmj.tgz'
        self.assertTrue(os.path.exists(cmpfile))
        cleanup()
    def test_example_build11(self):
        do_it(11)
        cmpfile='prod_2ins.tgz'
        self.assertTrue(os.path.exists(cmpfile))
        cleanup()
    def test_example_build12(self):
        do_it(12)
        cmpfile='prod_4zxb.tgz'
        self.assertTrue(os.path.exists(cmpfile))
        cleanup()
