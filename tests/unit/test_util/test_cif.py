import unittest
import logging
from pathlib import Path
logger=logging.getLogger(__name__)
from pestifer.util.cifutil import CIFdict, CIFload
from mmcif.api.PdbxContainers import DataContainer
from pestifer.molecule.residue import ResidueList, ResiduePlaceholder, ResiduePlaceholderList
from pestifer.molecule.atom import Atom, AtomList
from pestifer.objs.seqadv import Seqadv, SeqadvList
from pestifer.core.config import Config
from pestifer.molecule.bioassemb import Transform, TransformList, BioAssemb, BioAssembList
from pestifer.util.util import reduce_intlist
from pestifer.objs.resid import ResID
import os
import numpy as np
import glob

class TestCIF(unittest.TestCase):

    def setUp(self):
        input_dir = Path(__file__).parents[2] / "inputs"
        zmj = input_dir / '4zmj.cif'
        fae = input_dir / '8fae.cif'
        # copy to cwd
        dest_zmj = Path('4zmj.cif')
        dest_fae = Path('8fae.cif')
        if dest_zmj.exists():
            dest_zmj.unlink()
        if dest_fae.exists():
            dest_fae.unlink()
        os.symlink(zmj.resolve(), dest_zmj)
        os.symlink(fae.resolve(), dest_fae)

    def tearDown(self):
        dest_zmj = Path('4zmj.cif')
        dest_fae = Path('8fae.cif')
        if dest_zmj.exists():
            dest_zmj.unlink()
        if dest_fae.exists():
            dest_fae.unlink()
        logs = Path('.').glob('*.log')
        for log in logs:
            log.unlink()

    def test_load(self):
        source='8fae'
        p_struct=CIFload(source)
        self.assertEqual(type(p_struct),DataContainer)
        self.assertFalse(type(p_struct)==dict)

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
        residues=ResidueList.from_residuegrouped_atomlist(atoms)
        self.assertEqual(len(residues),2082)
        uCIDs=residues.uniqattrs(['chainID'])['chainID']
        print(uCIDs)
        self.assertEqual(len(uCIDs),82)
        nres=0
        for c in uCIDs:
            chain=residues.filter(lambda x: x.chainID == c)
            nres+=len(chain)
        self.assertEqual(len(residues),nres)

    def test_vmd(self):
        source='8fae'
        config=Config().configure_new()
        # config['rcsb_file_format']='mmCIF'
        vmd=config.get_scripter('vmd')
        vmd.newscript('testcif')
        vmd.addline(f'mol new {source}.cif')
        p_struct=CIFload(source)
        obj=p_struct.getObj('atom_site')
        atoms=AtomList([Atom(CIFdict(obj,i)) for i in range(len(obj))])
        raw_residues=ResidueList.from_residuegrouped_atomlist(atoms)
        raw_residues.apply_segtypes()
        residues=raw_residues.get(lambda x: x.segtype == 'protein')  # 8fae has some glycan residues with resid 0
        uCIDs=residues.uniqattrs(['chainID'])['chainID']
        nres=0
        for c in uCIDs:
            chain=residues.filter(lambda x: x.chainID == c)
            resids=[]
            for x in chain:
                resids.extend([str(y.resid) for y in x.atoms])
            residlist=' '.join(resids)
            serials=chain.atom_serials(as_type=int)
            vmd_red_list=reduce_intlist(serials)
            vmd.addline(f'set a [atomselect top "serial {vmd_red_list}"]')
            vmd.addline(f'set c [lsort -unique [$a get chain]]')
            vmd.addline(f'$a set chain {c}')
            vmd.addline(f'$a set resid [ list {residlist} ]')
            vmd.addline(f'set cn [lsort -unique [$a get chain]]')
            vmd.addline(f'set resids [$a get resid]')
            vmd.addline(f'puts "CHAIN $cn {c}"')
            vmd.addline(f'set b [atomselect top "chain {c}"]')
            vmd.addline(f'puts "COUNTS [$b num] {len(serials)}"')
            vmd.addline(f'puts "RESIDS && $resids && {residlist}"')
            vmd.addline(f'foreach avmd $resids apyt [list {residlist}] ' + r'{')
            vmd.addline(f'puts "   $avmd $apyt"')
            vmd.addline(r'}')
            nres+=len(chain)
        vmd.writescript()
        vmd.runscript()
        with open('testcif.log','r') as f:
            output=f.read().split('\n')
        old_logs=glob.glob('%*%')
        for ol in old_logs:
            os.remove(ol)
        self.assertTrue(len(output)>0)
        for l in output:
            if l.startswith('CHAIN'):
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
                    self.assertEqual(i,str(ResID.split_ri(j)[0]))  # CIF files do not have insertion codes on their native residue sequence numbers
        Path('testcif.tcl').unlink()

    def test_biomolassemb_cif(self):
        source='8fae'
        p_struct=CIFload(source)
        Assemblies=p_struct.getObj('pdbx_struct_assembly')
        gen=p_struct.getObj('pdbx_struct_assembly_gen')
        oper=p_struct.getObj('pdbx_struct_oper_list')
        BAList=BioAssembList()
        for i in range(len(Assemblies)):
            assemb_id=Assemblies.getValue('id',i)
            self.assertEqual(assemb_id,'1')
            this_gen_idx=gen.selectIndices(assemb_id,'assembly_id')[0]
            self.assertEqual(this_gen_idx,0)
            this_opers=gen.getValue('oper_expression',this_gen_idx).split(',')
            this_asyms=gen.getValue('asym_id_list',this_gen_idx).split(',')
            idx=0
            transforms=TransformList()
            for k,opere in enumerate(this_opers):
                oper_idx=oper.selectIndices(opere,'id')[0]
                self.assertEqual(k,oper_idx)
                m=np.identity(3)
                v=np.zeros(3)
                for i in range(3):
                    I=i+1
                    vlabel=f'vector[{I}]'
                    v[i]=float(oper.getValue(vlabel,oper_idx))
                    for j in range(3):
                        J=j+1
                        mlabel=f'matrix[{I}][{J}]'
                        m[i][j]=float(oper.getValue(mlabel,oper_idx))
                T=Transform(m,v,this_asyms,idx)
                transforms.append(T)
                idx+=1
            BA=BioAssemb(transforms)
            BAList.append(BA)
        self.assertEqual(len(BAList),1)

    def test_seqadv_cif(self):
        source='4zmj'
        pr=CIFload(source)
        obj=pr.getObj('atom_site')
        Atoms=AtomList([Atom(CIFdict(obj,i)) for i in range(len(obj))])
        obj=pr.getObj('struct_ref_seq_dif')
        Seqadvs=SeqadvList([Seqadv(CIFdict(obj,i)) for i in range(len(obj))])
        obj=pr.getObj('pdbx_unobs_or_zero_occ_residues')
        EmptyResidues=ResiduePlaceholderList([ResiduePlaceholder(CIFdict(obj,i)) for i in range(len(obj))])
        fromAtoms=ResidueList.from_residuegrouped_atomlist(Atoms)
        fromEmptyResidues=ResidueList.from_ResiduePlaceholderlist(EmptyResidues)
        Residues=fromAtoms+fromEmptyResidues
        Seqadvs.assign_residues(Residues)
        for s in Seqadvs:
            myres=Residues.get(lambda x: x.chainID == s.chainID and x.resid == s.resid)
            logger.debug(f'Seqadv {s.chainID}:{s.resid} typekey {s.typekey} residue {myres.resid if myres else None}')
            if not myres:
                self.assertTrue('engineered' not in s.typekey and 'conflict' not in s.typekey)
            else:
                self.assertTrue('engineered' in s.typekey or 'conflict' not in s.typekey)
