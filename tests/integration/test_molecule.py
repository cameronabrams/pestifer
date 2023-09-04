import unittest
from pestifer.molecule import Molecule
from pestifer.config import Config, segtype_of_resname
from pestifer.chainids import ChainIDManager
from io import StringIO
import yaml
class TestMolecule(unittest.TestCase):
    def get_source_dict(self,pdbid):
        source=f"""
source:
    biological_assembly: 0
    exclude: {{}}
    file_format: PDB
    id: {pdbid}
    sequence:
    fix_conflicts: true
    fix_engineered_mutations: true
    include_terminal_loops: false
    loops:
        declash:
        maxcycles: 20
        min_loop_length: 4
        sac_res_name: GLY
"""
        f=StringIO(source)
        directive=yaml.safe_load(f)
        return directive
    def test_molecule_au(self):
        c=Config()
        directive=self.get_source_dict('1gc1')
        m=Molecule(source=directive["source"])
        au=m.asymmetric_unit
        self.assertEqual(m.source['id'],'1gc1')
        self.assertEqual(len(au.Atoms),7877)
        self.assertEqual(len(au.SSBonds),14)
        self.assertEqual(len(au.Mutations),2)
        self.assertEqual(len(au.Conflicts),1)
        self.assertEqual(len(au.Missings),28)
        self.assertEqual(len(au.Links),15)
        self.assertEqual(len(au.Residues),1566)
        self.assertEqual(len(au.Segments),13)
        r0=au.Residues[0]
        self.assertEqual(r0.name,'THR')
        self.assertEqual(r0.segtype,'protein')
        for st in c.segtypes:
            res=[r for r in au.Residues if r.segtype==st]
            cbs=0
            for x in au.Segments:
                if x.segtype==st:
                    cbs+=len(x.residues)
            self.assertEqual(len(res),cbs)
        glycanres=[r for r in au.Residues if r.segtype=='glycan']
        g0=glycanres[0]
        self.assertEqual(g0.name,'NAG')

    def test_molecule_links(self):
        c=Config()
        directive=self.get_source_dict('4zmj')
        m=Molecule(source=directive["source"])
        au=m.asymmetric_unit
        l=au.Links[0]
        self.assertEqual(l.residue1.segtype,'protein')
        self.assertEqual(l.residue1.name,'ASN')
        self.assertEqual(l.residue1.chainID,'G')
        self.assertEqual(l.residue1.resseqnum,156)
        self.assertEqual(l.atom1.name,'ND2')
        self.assertEqual(l.atom1.altloc,'')
        self.assertEqual(l.residue2.segtype,'glycan')
        self.assertEqual(l.residue2.name,'NAG')
        self.assertEqual(l.residue2.chainID,'E')
        self.assertEqual(l.residue2.resseqnum,615)
        self.assertEqual(l.atom2.name,'C1')
        self.assertEqual(l.atom2.altloc,'')
        self.assertTrue(l.residue2 in l.residue1.down)
        self.assertTrue(l.residue1 in l.residue2.up)
        self.assertTrue(l in l.residue1.downlink)
        self.assertTrue(l in l.residue2.uplink)
        l=au.Links[-1]
        self.assertEqual(l.residue1.segtype,'glycan')
        self.assertEqual(l.residue1.name,'NAG')
        self.assertEqual(l.residue1.chainID,'D')
        self.assertEqual(l.residue1.resseqnum,1)
        self.assertEqual(l.atom1.name,'O4')
        self.assertEqual(l.atom1.altloc,'')
        self.assertEqual(l.residue2.segtype,'glycan')
        self.assertEqual(l.residue2.name,'NAG')
        self.assertEqual(l.residue2.chainID,'D')
        self.assertEqual(l.residue2.resseqnum,2)
        self.assertEqual(l.atom2.name,'C1')
        self.assertEqual(l.atom2.altloc,'')

    def test_molecule_bioassemb_4zmj(self):
        c=Config()
        cidm=ChainIDManager()
        directive=self.get_source_dict('4zmj')
        directive["source"]["biological_assembly"]=1
        m=Molecule(source=directive["source"],chainIDmanager=cidm)
        self.assertEqual(1,len(m.biological_assemblies))
        m.activate_biological_assembly(directive["source"]["biological_assembly"])
        ba=m.active_biological_assembly
        self.assertEqual(len(ba.transforms),3)
        cm=ba.transforms[0].chainIDmap
        self.assertEqual(cm['G'],'G')
        self.assertTrue(ba.transforms[0].is_identity())
        cm=ba.transforms[1].chainIDmap
        self.assertEqual(cm['G'],'H')
        cm=ba.transforms[2].chainIDmap
        self.assertEqual(cm['G'],'O')

    def test_molecule_bioassemb_4tvp(self):
        c=Config()
        cidm=ChainIDManager()
        directive=self.get_source_dict('4tvp')
        directive["source"]["biological_assembly"]=1
        directive["source"]["exclude"]["chains"]=['H','L','E','D']
        directive["source"]["exclude"]["resnames"]=['SO4']
        m=Molecule(source=directive["source"],chainIDmanager=cidm)
        self.assertEqual(1,len(m.biological_assemblies))
        m.activate_biological_assembly(directive["source"]["biological_assembly"])
        ba=m.active_biological_assembly
        self.assertEqual(len(ba.transforms),3)
        cm=ba.transforms[0].chainIDmap
        self.assertEqual(cm['G'],'G')
        self.assertTrue(ba.transforms[0].is_identity())
        cm=ba.transforms[1].chainIDmap
        self.assertEqual(cm['G'],'U')
        cm=ba.transforms[2].chainIDmap
        self.assertEqual(cm['G'],'o')

    def test_molecule_ancestry(self):
        directive=self.get_source_dict('4zmj')
        directive["source"]["biological_assembly"]=1
        m=Molecule(source=directive["source"],reset_counter=True)
        au=m.asymmetric_unit
        auao=au.ancestor_obj
        self.assertEqual(auao,m)
        for s in au.Segments:
            sao=s.ancestor_obj
            self.assertEqual(sao,m)
            self.assertEqual(sao.molid,0)

    def test_molecule_adjust_serials(self):
        directive=self.get_source_dict('4zmj')
        directive["source"]["biological_assembly"]=1
        m=Molecule(source=directive["source"],reset_counter=True)
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
        directive=self.get_source_dict('4zmj')
        directive["source"]["biological_assembly"]=1
        directive["source"]["file_format"]="mmCIF"
        cidm=ChainIDManager(format=directive["source"]["file_format"])
        m=Molecule(source=directive["source"],chainIDmanager=cidm)
        au=m.asymmetric_unit
        m.activate_biological_assembly(directive["source"]["biological_assembly"])
        self.assertEqual(len(au.Residues),659)
        ba=m.active_biological_assembly
        self.assertEqual(len(ba.transforms),3)
        self.assertEqual(len(au.Mutations),4)
        for m in au.Mutations:
            self.assertTrue(m.chainID in ['A','B'])