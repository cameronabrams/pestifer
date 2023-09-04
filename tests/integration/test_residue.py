import unittest
from pestifer.residue import Residue
from pestifer.config import Config
from pestifer.molecule import Molecule
from io import StringIO
import yaml
class TestResidue(unittest.TestCase):
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
    def test_from_atom(self):
        c=Config()
        directive=self.get_source_dict('1gc1')
        m=Molecule(source=directive["source"])
        au=m.asymmetric_unit
        residues=[]
        a=au.Atoms[0]
        r=Residue(a)
        self.assertEqual(a.chainID,r.chainID)
        residues.append(r)
        for a in au.Atoms[1:]:
            if not any([x.add_atom(a) for x in residues]):
                r=Residue(a)
                residues.append(r)
        self.assertEqual(len(residues),1538)

    def test_segtypes(self):
        c=Config()
        directive=self.get_source_dict('4tvp')
        m=Molecule(source=directive["source"])
        au=m.asymmetric_unit
        p=au.Residues.filter(segtype='protein')
        pc=[]
        for r in p:
            if not r.chainID in pc:
                pc.append(r.chainID)
        self.assertEqual(pc,['G', 'B', 'L', 'H', 'D', 'E'])
        i=au.Residues.filter(segtype='ion')
        ic=[]
        for r in i:
            if not r.chainID in ic:
                ic.append(r.chainID)
        self.assertEqual(ic,['X', 'Y', 'Z'])
        w=au.Residues.filter(segtype='water')
        wc=[]
        for r in w:
            if not r.chainID in wc:
                wc.append(r.chainID)
        self.assertEqual(wc,['a', 'b', 'c', 'd', 'e', 'f'])
        g=au.Residues.filter(segtype='glycan')
        gc=[]
        for r in g:
            if not r.chainID in gc:
                gc.append(r.chainID)
        self.assertEqual(gc,['A', 'C', 'F', 'I', 'J', 'K', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W'])
        l=au.Residues.filter(segtype='ligand')
        lc=[]
        for r in l:
            if not r.chainID in lc:
                lc.append(r.chainID)
        self.assertEqual(lc,[])
        o=au.Residues.filter(segtype='other')
        oc=[]
        for r in o:
            if not r.chainID in oc:
                oc.append(r.chainID)
        self.assertEqual(oc,[])

    def test_residuelist(self):
        c=Config()
        directive=self.get_source_dict('6m0j')
        m=Molecule(source=directive["source"])
        au=m.asymmetric_unit
        r=au.Residues.get_residue(resseqnum=427,chainID='A')
        self.assertEqual(r.name,'ASP')
        a=au.Residues.get_atom('OD2',resseqnum=427,chainID='A')
        self.assertEqual(a.serial,3302)
