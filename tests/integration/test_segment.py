import unittest
from pestifer.molecule import Molecule
from pestifer.config import Config
from pestifer.chainidmanager import ChainIDManager
from io import StringIO
import yaml
class TestSegment(unittest.TestCase):
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
    def test_segment(self):
        c=Config()
        cidm=ChainIDManager()
        directive=self.get_source_dict('4zmj')
        directive["source"]["biological_assembly"]=1
        m=Molecule(source=directive["source"],chainIDmanager=cidm).activate_biological_assembly(1)
        au=m.asymmetric_unit
        protein_segments=[x for x in au.segments if x.segtype=='protein']
        self.assertTrue(len(protein_segments),2)
        expected_set_of_chainIDs={'G','B','A','C','E','F'}
        self.assertEqual(set(au.segments.segnames),expected_set_of_chainIDs)
        e=au.segments.get(segname='E')
        chk1=all([x.chainID=='E' for x in e.residues])
        self.assertTrue(chk1)
        chk2=[]
        for x in e.residues:
            chk2.append(all([y.chainID=='E' for y in x.atoms]))
        self.assertTrue(all(chk2))
 