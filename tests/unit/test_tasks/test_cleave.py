from pestifer.core.controller import Controller
from pestifer.core.config import Config
import os
from pathlib import Path
import unittest

class TestCleaveTask(unittest.TestCase):
    def setUp(self):
        self.controller = Controller().configure(Config().configure_new())
        self.scripters = self.controller.config.scripters # shortcut
        input_dir = Path('../fixtures/cleave_inputs')
        # copy inputs to working directory
        psf = input_dir / 'in.psf'
        pdb = input_dir / 'in.pdb'
        # copy to cwd
        dest_psf = Path('in.psf')
        dest_pdb = Path('in.pdb')
        if dest_psf.exists():
            dest_psf.unlink()
        if dest_pdb.exists():
            dest_pdb.unlink()
        os.symlink(psf.resolve(), dest_psf)
        os.symlink(pdb.resolve(), dest_pdb)
    def test_cleave(self):
       tasklist = [
           {'continuation': {'psf':'in.psf','pdb':'in.pdb'}},
           {'cleave': { 
               'sites': [
                   'A:685-686',
                   'B:685-686',
                   'C:685-686'
               ]
           }},
           {'validate': {
               'tests': [
                   {'connection_test': {
                       'name': 'cleavages',
                       'selection': 'protein and chain A B C and resid 685 686',
                       'connection_type': 'interresidue',
                       'connection_count': 0
                   }},
               ]
           }},
           {'terminate': {'cleanup':True}}
       ]
       self.controller.reconfigure_tasks(tasklist)
       result = self.controller.do_tasks()
       Path('my_system.tar.gz').unlink()
       Path('artifacts.tar.gz').unlink()
