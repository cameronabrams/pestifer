import unittest
from pestifer.molecule import Molecule
from pestifer.config import Config, segtype_of_resname
from pestifer.chainidmanager import ChainIDManager
from pestifer.modmanager import ModManager
from pestifer.mods import *
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
"""
        f=StringIO(source)
        directive=yaml.safe_load(f)
        return directive
    def get_source_dict_alphafold(self,ac):
        source=f"""
source:
    biological_assembly: 0
    exclude: {{}}
    file_format: PDB
    alphafold: {ac}
    sequence:
        fix_conflicts: true
        fix_engineered_mutations: true
        include_terminal_loops: false
"""
        f=StringIO(source)
        directive=yaml.safe_load(f)
        return directive
    def test_molecule_au(self):
        c=Config()
        directive=self.get_source_dict('1gc1')
        m=Molecule(source=directive["source"])
        au=m.asymmetric_unit
        mods=au.modmanager
        topomods=mods.get('topomods',{})
        seqmods=mods.get('seqmods',{})
        ssbonds=topomods.get('ssbonds',[])
        links=topomods.get('links',[])
        mutations=seqmods.get('mutations',[])
        self.assertEqual(m.sourcespecs['id'],'1gc1')
        self.assertEqual(len(au.atoms),7877)
        self.assertEqual(len(ssbonds),14)
        self.assertEqual(len(mutations),3)
        self.assertEqual(len(links),15)
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

    def test_molecule_alphafold(self):
        c=Config()
        ac='O46077'
        directive=self.get_source_dict_alphafold(ac)
        m=Molecule(source=directive["source"])
        au=m.asymmetric_unit
        self.assertEqual(len(au.atoms),3230)
        self.assertEqual(len(au.residues),397)
        mods=au.modmanager
        topomods=mods.get('topomods',{})
        seqmods=mods.get('seqmods',{})
        print(seqmods)
        os.remove(f'{ac}.pdb')
        os.remove(f'{ac}.json')
        self.assertEqual(len(topomods),0)
        self.assertEqual(len(seqmods),1)
        self.assertTrue('terminals' in seqmods)
        self.assertEqual(len(seqmods['terminals']),1)


    def test_molecule_links(self):
        c=Config()
        directive=self.get_source_dict('4zmj')
        m=Molecule(source=directive["source"])
        au=m.asymmetric_unit
        links=au.modmanager['topomods']['links']
        l=links[0]
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
        c=Config()
        directive=self.get_source_dict('4zmj')
        directive["source"]["biological_assembly"]=1
        m=Molecule(source=directive["source"],modmanager=ModManager(),reset_counter=True)
        au=m.asymmetric_unit
        auao=au.ancestor_obj
        self.assertEqual(auao,m)
        for s in au.segments:
            sao=s.ancestor_obj
            self.assertEqual(sao,m)
            self.assertEqual(sao.molid,0)

    def test_molecule_adjust_serials(self):
        c=Config()
        directive=self.get_source_dict('4zmj')
        directive["source"]["biological_assembly"]=1
        m=Molecule(source=directive["source"],reset_counter=True)
        au=m.asymmetric_unit
        ters=au.modmanager.get('seqmods',{}).get('terminals',[])
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
        c=Config()
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
        mutations=au.modmanager.get('seqmods',{}).get('mutations',{})
        self.assertEqual(len(mutations),4)
        for m in mutations:
            self.assertTrue(m.chainID in ['A','B'])

    def test_molecule_existing(self):
        c=Config()
        source={
            'prebuilt':{
                'psf':'existing.psf',
                'pdb':'existing.pdb'
            }
        }
        m=Molecule(source=source).activate_biological_assembly(0)
        au=m.asymmetric_unit
        ssbonds=au.modmanager.get('topomods',{}).get('ssbonds',SSBondList([]))
        self.assertEqual(len(ssbonds),27)
        links=au.modmanager.get('topomods',{}).get('links',LinkList([]))
        self.assertEqual(len(links),129)
        fl=links[0]
        self.assertEqual(fl.residue1.resname,'ASN')
        self.assertEqual(fl.residue1.resseqnum,611)
        self.assertEqual(fl.residue2.resname,'BGLCNA')
        self.assertEqual(fl.residue2.resseqnum,701)
        self.assertEqual(len(au.segments),45)
        self.assertEqual(len(au.segments.filter(segtype='protein')),6)
        self.assertEqual(len(au.segments.filter(segtype='glycan')),39)

