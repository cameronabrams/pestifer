import unittest
from pestifer.cifutil import CIFdict, CIFload
from pidibble.pdbparse import PDBParser
from mmcif.io.IoAdapterCore import IoAdapterCore
from mmcif.api.PdbxContainers import DataContainer
from pestifer.residue import Atom, AtomList, ResidueList
from pestifer.config import Config
from pestifer.scriptwriters import VMD
from pestifer.util import reduce_intlist
import os
class TestCIF(unittest.TestCase):

    def test_fetch(self):
        source='8fae'
        PDBParser(PDBcode=source,input_format='mmCIF').fetch()
        self.assertTrue(os.path.isfile(f'{source}.cif'))
        os.remove(f'{source}.cif')

    def test_load(self):
        source='8fae'
        p_struct=CIFload(source)
        os.remove(f'{source}.cif')
        self.assertEqual(type(p_struct),DataContainer)

    def test_dict(self):
        source='8fae'
        p_struct=CIFload(source)
        os.remove(f'{source}.cif')
        recname='pdbx_unobs_or_zero_occ_residues'
        obj=p_struct.getObj(recname)
        self.assertEqual(len(obj),22)
        cd=CIFdict(obj,0,lowercase=True)
        attrlist=['id', 'PDB_model_num', 'polymer_flag', 'occupancy_flag', 'auth_asym_id', 'auth_comp_id', 'auth_seq_id', 'PDB_ins_code', 'label_asym_id', 'label_comp_id', 'label_seq_id']
        for attr in attrlist:
            self.assertTrue(attr.lower() in cd)
        cd=CIFdict(obj,0,lowercase=False)
        for attr in attrlist:
            self.assertTrue(attr in cd)

    def test_atoms(self):
        source='8fae'
        p_struct=CIFload(source)
        obj=p_struct.getObj('atom_site')
        atoms=AtomList([Atom(CIFdict(obj,i)) for i in range(len(obj))])
        self.assertEqual(len(atoms),17693)

    def test_residues(self):
        source='8fae'
        p_struct=CIFload(source)
        obj=p_struct.getObj('atom_site')
        atoms=AtomList([Atom(CIFdict(obj,i)) for i in range(len(obj))])
        self.assertEqual(len(atoms),17693)
        residues=ResidueList(atoms)
        self.assertEqual(len(residues),2082)
        uCIDs=residues.unique_chainIDs()
        print(uCIDs)
        self.assertEqual(len(uCIDs),82)
        nres=0
        for c in uCIDs:
            chain=residues.filter(chainID=c)
            nres+=len(chain)
        self.assertEqual(len(residues),nres)

    def test_vmd(self):
        source='8fae'
        config=Config()
        config.rcsb_file_format='mmCIF'
        vmd=VMD(config)
        vmd.newscript('testcif')
        vmd.addline(f'mol new {source}.cif')
        p_struct=CIFload(source)
        obj=p_struct.getObj('atom_site')
        atoms=AtomList([Atom(CIFdict(obj,i)) for i in range(len(obj))])
        residues=ResidueList(atoms)
        uCIDs=residues.unique_chainIDs()
        nres=0
        for c in uCIDs:
            chain=residues.filter(chainID=c)
            resids=[]
            for x in chain:
                resids.extend([str(y.resseqnum) for y in x.atoms])
            residlist=' '.join(resids)
            serials=chain.atom_serials(as_type=int)
            vmd_red_list=reduce_intlist(serials)
            vmd.addline(f'set a [atomselect top "serial {vmd_red_list}"]')
            vmd.addline(f'set c [lsort -unique [$a get chain]]')
            vmd.addline(f'$a set chain {c}')
            vmd.addline(f'$a set resid [ list {residlist} ]')
            vmd.addline(f'set cn [lsort -unique [$a get chain]]')
            vmd.addline(f'set resids [$a get resid]')
            vmd.addline(f'puts "LOOK $cn {c}"')
            vmd.addline(f'set b [atomselect top "chain {c}"]')
            vmd.addline(f'puts "COUNTS [$b num] {len(serials)}"')
            vmd.addline(f'puts "RESIDS && $resids && {residlist}"')
            nres+=len(chain)
        vmd.endscript()
        vmd.writescript()
        vmd.runscript()
        with open('testcif.log','r') as f:
            output=f.read().split('\n')
        for l in output:
            if l.startswith('LOOK'):
                fields=l.split()
                cvmd=fields[1]
                ccif=fields[2]
                # if len(ccif)==2:
                #     ccif=f'{ccif}1'
                self.assertEqual(ccif,cvmd)
            if l.startswith('COUNTS'):
                fields=l.split()
                cvmd=fields[1]
                ccif=fields[2]
                self.assertEqual(ccif,cvmd)
            if l.startswith('RESIDS'):
                fields=l.split('&&')
                rvmd=fields[1].split()
                rcif=fields[2].split()
                for i,j in zip(rvmd,rcif):
                    self.assertEqual(i,j)
 