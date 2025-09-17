import os
import unittest
from pidibble.pdbparse import PDBParser
from pestifer.core.config import Config
from pathlib import Path

class TestTCL(unittest.TestCase):

    def setUp(self):
        self.config = Config().configure_new()
        self.vmd = self.config.scripters['vmd']
        inputs_dir = Path(__file__).parent.parent.parent / 'inputs'
        pti = inputs_dir / '6pti.pdb'
        dest_pti = Path('6pti.pdb')
        if dest_pti.exists():
            dest_pti.unlink()
        os.symlink(pti.resolve(), dest_pti)
        input_dir = Path('../fixtures/tcl_inputs')
        bad = input_dir / 'bad.pdb'
        dest_bad = Path('bad.pdb')
        if dest_bad.exists():
            dest_bad.unlink()
        os.symlink(bad.resolve(), dest_bad)
        good = input_dir / 'good.pdb'
        dest_good = Path('good.pdb')
        if dest_good.exists():
            dest_good.unlink()
        os.symlink(good.resolve(), dest_good)

    def tearDown(self):
        dest_pti = Path('6pti.pdb')
        if dest_pti.exists():
            dest_pti.unlink()
        dest_bad = Path('bad.pdb')
        if dest_bad.exists():
            dest_bad.unlink()
        dest_good = Path('good.pdb')
        if dest_good.exists():
            dest_good.unlink()
        scripts = Path('.').glob('*.tcl')
        for script in scripts:
            script.unlink()
        logs = Path('.').glob('*.log')
        for log in logs:
            log.unlink()
        cycled_logs = Path('.').glob('%*%')
        for log in cycled_logs:
            log.unlink()

    def test_backup(self):
        # test backup and restor procs in saverestore.tcl
        source='6pti'
        o=PDBParser(source_id=source, source_db='rcsb').parse()
        oatoms=o.parsed['ATOM']
        ox=[a.x for a in oatoms if a.name=='CA']
        self.assertEqual(len(ox),58)
        self.vmd.newscript('testtcl')
        self.vmd.addline(f'mol new {source}.pdb')
        self.vmd.addline(f'set a [atomselect top all]')
        self.vmd.addline(r'set attributes { name x }')
        self.vmd.addline(f'set data [ backup $a $attributes ]')
        self.vmd.addline(f'[atomselect top "name CA"] set x 0')
        self.vmd.addline(f'[atomselect top "name CA"] writepdb bad.pdb')
        self.vmd.addline(f'restore $a $attributes $data')
        self.vmd.addline(f'[atomselect top "name CA"] writepdb good.pdb')
        self.vmd.writescript()
        self.vmd.runscript()
        p=PDBParser(filepath='bad.pdb').parse()
        badatoms=p.parsed['ATOM']
        self.assertEqual(len(badatoms),58)
        badx=[a.x for a in badatoms]
        self.assertTrue(not any(badx))
        g=PDBParser(filepath='good.pdb').parse()
        goodatoms=g.parsed['ATOM']
        goodx=[a.x for a in goodatoms]
        goodchk=all([x==y for x,y in zip(ox,goodx)])
        self.assertTrue(goodchk)

    def test_backup_scriptwriter(self):
        # test backup and restor procs in saverestore.tcl
        source='6pti'
        o=PDBParser(source_id=source, source_db='rcsb').parse()
        oatoms=o.parsed['ATOM']
        ox=[a.x for a in oatoms if a.name=='CA']
        oy=[a.y for a in oatoms if a.name=='CA']
        self.assertEqual(len(ox),58)
        self.vmd.newscript('testtcl')
        self.vmd.addline(f'mol new {source}.pdb')
        self.vmd.addline(f'set a [atomselect top all]')
        self.vmd.backup_selection('a')
        self.vmd.addline(f'[atomselect top "name CA"] set x 0')
        self.vmd.addline(f'[atomselect top "name CA"] set y -1')
        self.vmd.addline(f'[atomselect top "name CA"] set z 1')
        self.vmd.addline(f'[atomselect top "name CA"] writepdb bad.pdb')
        self.vmd.restore_selection('a')
        self.vmd.addline(f'[atomselect top "name CA"] writepdb good.pdb')
        self.vmd.writescript()
        self.vmd.runscript()
        p=PDBParser(filepath='bad.pdb').parse()
        badatoms=p.parsed['ATOM']
        self.assertEqual(len(badatoms),58)
        badx=[a.x for a in badatoms]
        self.assertTrue(not any(badx))
        bady=[a.y for a in badatoms]
        ychk=all([y==-1 for y in bady])
        self.assertTrue(ychk)
        g=PDBParser(filepath='good.pdb').parse()
        goodatoms=g.parsed['ATOM']
        goodx=[a.x for a in goodatoms]
        goodchk=all([x==y for x,y in zip(ox,goodx)])
        self.assertTrue(goodchk)
        goody=[a.y for a in goodatoms]
        goodchk=all([x==y for x,y in zip(oy,goody)])
        self.assertTrue(goodchk)
