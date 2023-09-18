import unittest
from pestifer.molecule import Molecule
from pestifer.config import Config, segtype_of_resname
from pestifer.chainids import ChainIDManager
from pestifer.modcontainer import ModContainer
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
        mods=au.mods
        topomods=mods.topomods
        seqmods=mods.seqmods
        self.assertEqual(m.sourcespecs['id'],'1gc1')
        self.assertEqual(len(au.atoms),7877)
        self.assertEqual(len(topomods.ssbonds),14)
        self.assertEqual(len(seqmods.mutations),3)
        self.assertEqual(len(topomods.links),15)
        self.assertEqual(len(au.residues),1566)
        self.assertEqual(len(au.segments),13)
        r0=au.residues[0]
        self.assertEqual(r0.resname,'THR')
        self.assertEqual(r0.segtype,'protein')
        for st in c.segtypes:
            res=[r for r in au.residues if r.segtype==st]
            cbs=0
            for x in au.segments:
                if x.segtype==st:
                    cbs+=len(x.residues)
            self.assertEqual(len(res),cbs)
        glycanres=[r for r in au.residues if r.segtype=='glycan']
        g0=glycanres[0]
        self.assertEqual(g0.resname,'NAG')

    def test_molecule_links(self):
        c=Config()
        directive=self.get_source_dict('4zmj')
        m=Molecule(source=directive["source"])
        au=m.asymmetric_unit
        l=au.mods.topomods.links[0]
        self.assertEqual(l.residue1.segtype,'protein')
        self.assertEqual(l.residue1.resname,'ASN')
        self.assertEqual(l.residue1.chainID,'G')
        self.assertEqual(l.residue1.resseqnum,156)
        self.assertEqual(l.atom1.name,'ND2')
        self.assertEqual(l.atom1.altloc,'')
        self.assertEqual(l.residue2.segtype,'glycan')
        self.assertEqual(l.residue2.resname,'NAG')
        self.assertEqual(l.residue2.chainID,'E')
        self.assertEqual(l.residue2.resseqnum,615)
        self.assertEqual(l.atom2.name,'C1')
        self.assertEqual(l.atom2.altloc,'')
        self.assertTrue(l.residue2 in l.residue1.down)
        self.assertTrue(l.residue1 in l.residue2.up)
        self.assertTrue(l in l.residue1.downlink)
        self.assertTrue(l in l.residue2.uplink)
        l=au.mods.topomods.links[-1]
        self.assertEqual(l.residue1.segtype,'glycan')
        self.assertEqual(l.residue1.resname,'NAG')
        self.assertEqual(l.residue1.chainID,'C')
        self.assertEqual(l.residue1.resseqnum,1)
        self.assertEqual(l.atom1.name,'O4')
        self.assertEqual(l.atom1.altloc,'')
        self.assertEqual(l.residue2.segtype,'glycan')
        self.assertEqual(l.residue2.resname,'NAG')
        self.assertEqual(l.residue2.chainID,'C')
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
        self.assertEqual(cm['G'],'N')

    def test_molecule_ancestry(self):
        directive=self.get_source_dict('4zmj')
        directive["source"]["biological_assembly"]=1
        m=Molecule(source=directive["source"],usermods=ModContainer(),reset_counter=True)
        au=m.asymmetric_unit
        auao=au.ancestor_obj
        self.assertEqual(auao,m)
        for s in au.segments:
            sao=s.ancestor_obj
            self.assertEqual(sao,m)
            self.assertEqual(sao.molid,0)

    def test_molecule_adjust_serials(self):
        directive=self.get_source_dict('4zmj')
        directive["source"]["biological_assembly"]=1
        m=Molecule(source=directive["source"],reset_counter=True)
        au=m.asymmetric_unit
        print(au.mods.seqmods)
        self.assertTrue(hasattr(au.mods.seqmods,'terminals'))
        ters=au.mods.seqmods.terminals
        self.assertEqual(len(ters),2)
        atom_serials=[x.serial for x in au.atoms]
        orig_atom_serials=[]
        for a in au.atoms:
            if '_ORIGINAL_' in a.__dict__:
                orig_atom_serials.append(a._ORIGINAL_['serial'])
            else:
                orig_atom_serials.append(a.serial)
        self.assertEqual(len(atom_serials),4856)
        self.assertEqual(len(orig_atom_serials),4856)
        self.assertFalse(all([x==y for x,y in zip(atom_serials,orig_atom_serials)]))
        genuine_atoms=[au.atoms.get(serial=x.serial) for x in ters]
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
        self.assertEqual(len(au.residues),659)
        ba=m.active_biological_assembly
        self.assertEqual(len(ba.transforms),3)
        self.assertEqual(len(au.mods.seqmods.mutations),4)
        for m in au.mods.seqmods.mutations:
            self.assertTrue(m.chainID in ['A','B'])