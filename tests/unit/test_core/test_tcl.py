import unittest
from pidibble.pdbparse import PDBParser
from pestifer.core.config import Config

class TestTCL(unittest.TestCase):

    def setUp(self):
        self.config = Config().configure_new()
        self.vmd = self.config.scripters['vmd']

    def test_backup(self):
        # test backup and restor procs in saverestore.tcl
        source='6pti'
        o=PDBParser(PDBcode=source).parse()
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
        o=PDBParser(PDBcode=source).parse()
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
