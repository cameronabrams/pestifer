import unittest
from pestifer.molecule import Molecule
from pestifer.config import Config
from pestifer.chainids import ChainIDManager

class TestMolecule(unittest.TestCase):

    def test_molecule_au(self):
        c=Config()
        self.assertTrue(hasattr(c,'defs'))
        m=Molecule(config=c,source='1gc1')
        au=m.asymmetric_unit
        self.assertEqual(m.source,'1gc1')
        self.assertEqual(len(au.Atoms),7877)
        self.assertEqual(len(au.SSBonds),14)
        self.assertEqual(len(au.Mutations),2)
        self.assertEqual(len(au.Conflicts),1)
        self.assertEqual(len(au.Missings),28)
        self.assertEqual(len(au.Links),15)
        self.assertEqual(len(au.Residues),1566)
        r0=au.Residues[0]
        self.assertEqual(r0.name,'THR')
        self.assertEqual(r0.segtype,'PROTEIN')
        waterres=[r for r in au.Residues if r.segtype=='WATER']
        self.assertEqual(len(waterres),603)
        glycanres=[r for r in au.Residues if r.segtype=='GLYCAN']
        self.assertEqual(len(glycanres),15)
        g0=glycanres[0]
        self.assertEqual(g0.name,'NAG')

    def test_molecule_links(self):
        m=Molecule(source='4zmj',config=Config())
        au=m.asymmetric_unit
        l=au.Links[0]
        self.assertEqual(l.residue1.segtype,'PROTEIN')
        self.assertEqual(l.residue1.name,'ASN')
        self.assertEqual(l.residue1.chainID,'G')
        self.assertEqual(l.residue1.resseqnum,156)
        self.assertEqual(l.atom1.name,'ND2')
        self.assertEqual(l.atom1.altloc,'')
        self.assertEqual(l.residue2.segtype,'GLYCAN')
        self.assertEqual(l.residue2.name,'NAG')
        self.assertEqual(l.residue2.chainID,'G')
        self.assertEqual(l.residue2.resseqnum,615)
        self.assertEqual(l.atom2.name,'C1')
        self.assertEqual(l.atom2.altloc,'')
        self.assertTrue(l.residue2 in l.residue1.down)
        self.assertTrue(l.residue1 in l.residue2.up)
        self.assertTrue(l in l.residue1.downlink)
        self.assertTrue(l in l.residue2.uplink)
        l=au.Links[-1]
        self.assertEqual(l.residue1.segtype,'GLYCAN')
        self.assertEqual(l.residue1.name,'NAG')
        self.assertEqual(l.residue1.chainID,'D')
        self.assertEqual(l.residue1.resseqnum,1)
        self.assertEqual(l.atom1.name,'O4')
        self.assertEqual(l.atom1.altloc,'')
        self.assertEqual(l.residue2.segtype,'GLYCAN')
        self.assertEqual(l.residue2.name,'NAG')
        self.assertEqual(l.residue2.chainID,'D')
        self.assertEqual(l.residue2.resseqnum,2)
        self.assertEqual(l.atom2.name,'C1')
        self.assertEqual(l.atom2.altloc,'')

    def test_molecule_bioassemb(self):
        m=Molecule(source='4zmj',config=Config())
        self.assertEqual(1,len(m.biological_assemblies))
        m.activate_biological_assembly(1,ChainIDManager())
        ba=m.active_biological_assembly
        self.assertEqual(len(ba.transforms),3)
        cm=ba.transforms[0].chainIDmap
        self.assertEqual(cm['G'],'G')
        self.assertTrue(ba.transforms[0].is_identity())
        cm=ba.transforms[1].chainIDmap
        self.assertEqual(cm['G'],'J')
        cm=ba.transforms[2].chainIDmap
        self.assertEqual(cm['G'],'O')
        m=Molecule(source='4tvp',config=Config(),excludes={'chains':['H','L','E','D'],'resnames':['PO4']})
        self.assertEqual(1,len(m.biological_assemblies))
        m.activate_biological_assembly(1,ChainIDManager())
        ba=m.active_biological_assembly
        self.assertEqual(len(ba.transforms),3)
        cm=ba.transforms[0].chainIDmap
        self.assertEqual(cm['G'],'G')
        self.assertTrue(ba.transforms[0].is_identity())
        cm=ba.transforms[1].chainIDmap
        self.assertEqual(cm['G'],'U')
        cm=ba.transforms[2].chainIDmap
        self.assertEqual(cm['G'],'k')

    def test_molecule_ancestry(self):
        m=Molecule(source='4zmj',reset_counter=True,config=Config())
        au=m.asymmetric_unit
        auao=au.ancestor_obj
        self.assertEqual(auao,m)
        for s in au.Segments:
            print(str(s))
            sao=s.ancestor_obj
            self.assertEqual(sao,m)
            self.assertEqual(sao.molid,0)

    def test_molecule_adjust_serials(self):
        m=Molecule(source='4zmj',config=Config())
        au=m.asymmetric_unit
        self.assertTrue(hasattr(au,'Ters'))
        ters=au.Ters
        self.assertEqual(len(ters),2)
        atom_serials=[x.serial for x in au.Atoms]
        orig_atom_serials=[]
        for a in au.Atoms:
            if '_ORIGINAL_' in a.__dict__:
                orig_atom_serials.append(a._ORIGINAL_['serial'])
            else:
                orig_atom_serials.append(a.serial)
        self.assertEqual(len(atom_serials),4856)
        self.assertEqual(len(orig_atom_serials),4856)
        self.assertFalse(all([x==y for x,y in zip(atom_serials,orig_atom_serials)]))
        genuine_atoms=[au.Atoms.get(serial=x.serial) for x in ters]
        self.assertEqual(len(genuine_atoms),2)
        self.assertEqual(atom_serials[-1],4856)

    def test_molecule_cif(self):
        source='4zmj'
        config=Config()
        config['rcsb_file_format']='mmCIF'
        m=Molecule(source=source,config=config)
        au=m.asymmetric_unit
        m.activate_biological_assembly(1,ChainIDManager(format=config['rcsb_file_format']))
        self.assertEqual(len(au.Residues),659)
        ba=m.active_biological_assembly
        self.assertEqual(len(ba.transforms),3)
        self.assertEqual(len(au.Mutations),4)
        for m in au.Mutations:
            self.assertTrue(m.chainID in ['A','B'])