import unittest
from pidibble.pdbparse import PDBParser
from pestifer.scriptwriters import *
from pestifer.config import Config

class TestTCL(unittest.TestCase):
    def test_backup(self):
        # test backup and restor procs in saverestore.tcl
        source='6pti'
        config=Config()
        vmd=VMD(config)
        o=PDBParser(PDBcode=source).parse()
        oatoms=o.parsed['ATOM']
        ox=[a.x for a in oatoms if a.name=='CA']
        self.assertEqual(len(ox),58)
        vmd.newscript('testtcl')
        vmd.addline(f'mol new {source}.pdb')
        vmd.addline(f'set a [atomselect top all]')
        vmd.addline(r'set attributes { name x }')
        vmd.addline(f'set data [ backup $a $attributes ]')
        vmd.addline(f'[atomselect top "name CA"] set x 0')
        vmd.addline(f'[atomselect top "name CA"] writepdb bad.pdb')
        vmd.addline(f'restore $a $attributes $data')
        vmd.addline(f'[atomselect top "name CA"] writepdb good.pdb')
        vmd.writescript()
        vmd.runscript()
        p=PDBParser(PDBcode='bad').parse()
        badatoms=p.parsed['ATOM']
        self.assertEqual(len(badatoms),58)
        badx=[a.x for a in badatoms]
        self.assertTrue(not any(badx))
        g=PDBParser(PDBcode='good').parse()
        goodatoms=g.parsed['ATOM']
        goodx=[a.x for a in goodatoms]
        goodchk=all([x==y for x,y in zip(ox,goodx)])
        self.assertTrue(goodchk)

    def test_backup_scriptwriter(self):
        # test backup and restor procs in saverestore.tcl
        source='6pti'
        config=Config()
        vmd=VMD(config)
        o=PDBParser(PDBcode=source).parse()
        oatoms=o.parsed['ATOM']
        ox=[a.x for a in oatoms if a.name=='CA']
        oy=[a.y for a in oatoms if a.name=='CA']
        self.assertEqual(len(ox),58)
        vmd.newscript('testtcl')
        vmd.addline(f'mol new {source}.pdb')
        vmd.addline(f'set a [atomselect top all]')
        vmd.backup_selection('a')
        vmd.addline(f'[atomselect top "name CA"] set x 0')
        vmd.addline(f'[atomselect top "name CA"] set y -1')
        vmd.addline(f'[atomselect top "name CA"] set z 1')
        vmd.addline(f'[atomselect top "name CA"] writepdb bad.pdb')
        vmd.restore_selection('a')
        vmd.addline(f'[atomselect top "name CA"] writepdb good.pdb')
        vmd.writescript()
        vmd.runscript()
        p=PDBParser(PDBcode='bad').parse()
        badatoms=p.parsed['ATOM']
        self.assertEqual(len(badatoms),58)
        badx=[a.x for a in badatoms]
        self.assertTrue(not any(badx))
        bady=[a.y for a in badatoms]
        ychk=all([y==-1 for y in bady])
        self.assertTrue(ychk)
        g=PDBParser(PDBcode='good').parse()
        goodatoms=g.parsed['ATOM']
        goodx=[a.x for a in goodatoms]
        goodchk=all([x==y for x,y in zip(ox,goodx)])
        self.assertTrue(goodchk)
        goody=[a.y for a in goodatoms]
        goodchk=all([x==y for x,y in zip(oy,goody)])
        self.assertTrue(goodchk)
