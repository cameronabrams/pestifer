from pestifer.core.controller import Controller
from pestifer.core.config import Config
import os
from pathlib import Path
import unittest
import pytest

pytestmark = pytest.mark.needs_tools

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
       # Cleaving the protein must not change the glycans.  This test once passed while every GlcNAc
       # came out as glucose (126 BGLCNA -> BGLC, 756 acetyl atoms gone): the segment PDBs psfgen
       # reads keep four columns of a residue name, and nothing here looked at what was built.
       self.assertEqual(self._glycan_composition('my_system.psf'), self._glycan_composition('../fixtures/cleave_inputs/in.psf'))
       self.assertEqual(self._glycan_composition('my_system.psf')['BGLCNA'], 126)
       Path('my_system.tar.gz').unlink(missing_ok=True)
       Path('artifacts.tar.gz').unlink(missing_ok=True)

    @staticmethod
    def _glycan_composition(psf):
       import collections
       from pestifer.core.labels import Labels
       lines = open(psf).read().splitlines()
       i = next(k for k, l in enumerate(lines) if '!NATOM' in l)
       n = int(lines[i].split()[0])
       residues = {(t[1], t[2]): t[3] for t in (l.split() for l in lines[i + 1:i + 1 + n])}
       return collections.Counter(r for r in residues.values() if Labels.segtype_of_resname.get(r) == 'glycan')
