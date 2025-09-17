import unittest
import logging
from pestifer.molecule.molecule import Molecule
from pestifer.core.config import Config
from pestifer.core.labels import Labels #segtype_of_resname

from pestifer.core.pipeline import PipelineContext
from pestifer.tasks.fetch import FetchTask
from pestifer.tasks.continuation import ContinuationTask
from pestifer.objs.resid import ResID
from pestifer.objs.ssbond import SSBondList
from pestifer.objs.link import LinkList
from io import StringIO
import os
import yaml
from pathlib import Path

logger=logging.getLogger(__name__)

class TestMolecule(unittest.TestCase):

    segtypes = list(Labels.segtypes.keys())
    inputs_dir = Path(__file__).parents[2] / "inputs"

    def tearDown(self):
        ac='O46077'
        Path('1gc1.pdb').unlink(missing_ok=True)
        Path('4zmj.pdb').unlink(missing_ok=True)
        for f in Path('.').glob(ac+'*'):
            f.unlink()

    # def tearDown(self):
        
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
    source_id: {ac}
    source_db: alphafold
    sequence:
        fix_conflicts: true
        fix_engineered_mutations: true
        include_terminal_loops: false
"""
        f=StringIO(source)
        directive=yaml.safe_load(f)
        return directive
    
    def test_molecule_au(self):
        pipeline = PipelineContext()
        task = FetchTask(specs={'source': 'rcsb', 'sourceID': '1gc1'})
        task.provision(packet=dict(pipeline=pipeline))
        self.assertIsInstance(task, FetchTask)
        task.execute()
        m=Molecule(source={
                        'source_id': '1gc1',
                        'source_db': 'rcsb',
                        'file_format': 'PDB',
                        'sequence': {
                            'fix_conflicts': True,
                            'fix_engineered_mutations': True,
                            'include_terminal_loops': False
                        }
                    }, molid=0)
        au=m.asymmetric_unit
        objs=au.objmanager
        topol=objs.get('topol',{})
        seq=objs.get('seq',{})
        ssbonds=topol.get('ssbonds',[])
        links=topol.get('links',[])
        mutations=seq.get('mutations',[])
        self.assertEqual(m.sourcespecs['source_id'],'1gc1')
        self.assertEqual(len(au.atoms),7877)
        self.assertEqual(len(ssbonds),14)
        self.assertEqual(len(mutations),3)
        self.assertEqual(len(links),15)
        self.assertEqual(len(au.residues),1566)
        self.assertEqual(len(au.segments),13)
        r0=au.residues[0]
        self.assertEqual(r0.resname,'THR')
        self.assertEqual(r0.segtype,'protein')
        for st in self.segtypes:
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
        ac='O46077'
        pipeline = PipelineContext()
        task = FetchTask(specs={'source': 'alphafold', 'sourceID': ac})
        task.provision(packet=dict(pipeline=pipeline))
        self.assertIsInstance(task, FetchTask)
        task.execute()
        m=Molecule(source={
                        'source_id': ac,
                        'source_db': 'alphafold',
                        'file_format': 'PDB',
                        'sequence': {
                            'fix_conflicts': True,
                            'fix_engineered_mutations': True,
                            'include_terminal_loops': False
                        }
                    }, molid=0)
        au=m.asymmetric_unit
        self.assertEqual(len(au.atoms),3230)
        self.assertEqual(len(au.residues),397)
        objs=au.objmanager
        topol=objs.get('topol',{})
        seq=objs.get('seq',{})
        os.remove(f'{ac}.pdb')
        os.remove(f'{ac}.json')
        self.assertEqual(len(topol),2)
        self.assertEqual(len(seq),4)
        self.assertTrue('terminals' in seq)
        self.assertEqual(len(seq['terminals']),1)


    def test_molecule_links(self):
        ac='4zmj'
        pipeline = PipelineContext()
        task = FetchTask(specs={'source': 'rcsb', 'sourceID': ac})
        task.provision(packet=dict(pipeline=pipeline))
        self.assertIsInstance(task, FetchTask)
        task.execute()
        m=Molecule(source={
                        'source_id': ac,
                        'source_db': 'rcsb',
                        'file_format': 'PDB',
                        'sequence': {
                            'fix_conflicts': True,
                            'fix_engineered_mutations': True,
                            'include_terminal_loops': False
                        }
                    }, molid=0)
        au=m.asymmetric_unit
        self.assertEqual(len(au.residues),659)
        links=au.objmanager['topol']['links']
        l=links[0]
        self.assertEqual(l.residue1.segtype,'protein')
        self.assertEqual(l.residue1.resname,'ASN')
        self.assertEqual(l.residue1.chainID,'G')
        self.assertEqual(l.residue1.resid,ResID(156))
        self.assertEqual(l.atom1.name,'ND2')
        self.assertEqual(l.atom1.altloc,'')
        self.assertEqual(l.residue2.segtype,'glycan')
        self.assertEqual(l.residue2.resname,'NAG')
        self.assertEqual(l.residue2.chainID,'E')
        self.assertEqual(l.residue2.resid,ResID(615))
        self.assertEqual(l.atom2.name,'C1')
        self.assertEqual(l.atom2.altloc,'')
        self.assertTrue(l.residue2 in l.residue1.down)
        self.assertTrue(l.residue1 in l.residue2.up)

    def test_molecule_bioassemb_4zmj(self):
        ac='4zmj'
        pipeline = PipelineContext()
        task = FetchTask(specs={'source': 'rcsb', 'sourceID': ac})
        task.provision(packet=dict(pipeline=pipeline))
        self.assertIsInstance(task, FetchTask)
        task.execute()
        m=Molecule(source={
                        'source_id': ac,
                        'source_db': 'rcsb',
                        'file_format': 'PDB',
                        'biological_assembly': 1,
                        'sequence': {
                            'fix_conflicts': True,
                            'fix_engineered_mutations': True,
                            'include_terminal_loops': False
                        }
                    }, molid=0)
        m.activate_biological_assembly(1)
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
        ac='4zmj'
        pipeline = PipelineContext()
        task = FetchTask(specs={'source': 'rcsb', 'sourceID': ac})
        task.provision(packet=dict(pipeline=pipeline))
        self.assertIsInstance(task, FetchTask)
        task.execute()
        m=Molecule(source={
                        'source_id': ac,
                        'source_db': 'rcsb',
                        'file_format': 'PDB',
                        'biological_assembly': 1,
                        'sequence': {
                            'fix_conflicts': True,
                            'fix_engineered_mutations': True,
                            'include_terminal_loops': False
                        }
                    }, molid=0)
        au=m.asymmetric_unit
        auao=au.parent_molecule
        self.assertEqual(auao,m)
        for s in au.segments:
            sao=s.parent_molecule
            self.assertEqual(sao,m)
            self.assertEqual(sao.molid,0)

    def test_molecule_adjust_serials(self):
        ac='4zmj'
        pipeline = PipelineContext()
        task = FetchTask(specs={'source': 'rcsb', 'sourceID': ac})
        task.provision(packet=dict(pipeline=pipeline))
        self.assertIsInstance(task, FetchTask)
        task.execute()
        m=Molecule(source={
                        'source_id': ac,
                        'source_db': 'rcsb',
                        'file_format': 'PDB',
                        'biological_assembly': 1,
                        'sequence': {
                            'fix_conflicts': True,
                            'fix_engineered_mutations': True,
                            'include_terminal_loops': False
                        }
                    }, molid=0)
        au=m.asymmetric_unit
        ters=au.objmanager.get('seq',{}).get('terminals',[])
        self.assertEqual(len(ters),2)
        atom_serials=[x.serial for x in au.atoms.data]
        orig_atom_serials=[]
        for a in au.atoms.data:
            if len(a.ORIGINAL_ATTRIBUTES)>0:
                orig_atom_serials.append(a.ORIGINAL_ATTRIBUTES['serial'])
            else:
                orig_atom_serials.append(a.serial)
        self.assertEqual(len(atom_serials),4856)
        self.assertEqual(len(orig_atom_serials),4856)
        self.assertFalse(all([x==y for x,y in zip(atom_serials,orig_atom_serials)]))    
        genuine_atoms=[au.atoms.get(lambda y: y.serial == x.serial) for x in ters]
        self.assertEqual(len(genuine_atoms),2)
        self.assertEqual(atom_serials[-1],4856)

    def test_molecule_cif(self):
        directive=self.get_source_dict('4zmj')
        directive["source"]["biological_assembly"]=1
        mol={}
        au={}
        atoms={}
        segts={}
        for fmt in ['pdb','cif']:
            ac='4zmj'
            pipeline = PipelineContext()
            task = FetchTask(specs={'source': 'rcsb', 'sourceID': ac, 'source_format': fmt})
            task.provision(packet=dict(pipeline=pipeline))
            self.assertIsInstance(task, FetchTask)
            task.execute()
            mol[fmt]=Molecule(source={
                        'source_id': ac,
                        'source_db': 'rcsb',
                        'file_format': fmt,
                        'biological_assembly': 1,
                        'sequence': {
                            'fix_conflicts': True,
                            'fix_engineered_mutations': True,
                            'include_terminal_loops': False
                        }
                    }, molid=0)
            au[fmt]=mol[fmt].asymmetric_unit
            mol[fmt].activate_biological_assembly(directive["source"]["biological_assembly"])
            atoms[fmt]=au[fmt].atoms
            segts[fmt]=set([x.segtype for x in au[fmt].residues])
        self.assertEqual(len(atoms['pdb']),len(atoms['cif']))
        atom_mismatches=[]
        eps=1.e-5
        for pa,ca in zip(atoms['pdb'], atoms['cif']):
            if abs(pa.x-ca.x)>eps or abs(pa.y-ca.y)>eps or abs(pa.z-ca.z)>eps:
                atom_mismatches.append([pa,ca])
        self.assertEqual(len(atom_mismatches),0)
        for pa,ca in zip(atoms['pdb'], atoms['cif']):
            if pa.name!=ca.name:
                atom_mismatches.append([pa,ca])
        self.assertEqual(len(atom_mismatches),0)
        for pa,ca in zip(atoms['pdb'], atoms['cif']):
            if pa.resname!=ca.resname:
                atom_mismatches.append([pa,ca])
        for pa,ca in zip(atoms['pdb'], atoms['cif']):
            if pa.resid!=ResID(ca.auth_seq_id) and pa.resid!=ResID(ca.auth_seq_id,ca.pdbx_pdb_ins_code):
                atom_mismatches.append([pa,ca])
        self.assertEqual(len(atom_mismatches),0)
        for pa,ca in zip(atoms['pdb'], atoms['cif']):
            if pa.chainID!=ca.chainID and pa.chainID!=ca.auth_asym_id:
                atom_mismatches.append([pa,ca])
        # msg=''
        # if len(atom_mismatches)>0:
        #     fm=atom_mismatches[0]
        #     msg=f'{fm[0].chainID}_{fm[0].resname}{fm[0].resseqnum}_{fm[0].name} != [{fm[1].chainID}|{fm[1].auth_asym_id}]_{fm[1].resname}[{fm[1].resseqnum}|{fm[1].auth_seq_id}]_{fm[1].name}'
        # self.assertEqual(len(atom_mismatches),0,msg=msg)


        self.assertEqual(segts['pdb'],segts['cif'])
        res_pdb={}
        res_count_pdb={}
        for st in segts['pdb']:
            res_pdb[st]=[x for x in au['pdb'].residues if x.segtype==st]
            res_count_pdb[st]=len(res_pdb[st])
        self.assertEqual(res_count_pdb['protein'],634)
        self.assertEqual(res_count_pdb['glycan'],25)
        res_count_checks={}
        res_cif={}
        res_count_cif={}
        for st in ['protein', 'glycan']:
            res_cif[st]=[x for x in au['cif'].residues if x.segtype==st]
            res_count_cif[st]=len(res_cif[st])
            res_count_checks[st]={}
            res_count_checks[st]['check']=res_count_cif[st]==res_count_pdb[st]
            res_count_checks[st]['msg']=f'mismatch {st} residues PDB({res_count_pdb[st]}) mmCIF({res_count_cif[st]})'
        check_good=all([res_count_checks[st]['check'] for st in segts['pdb']])
        check_msg='; '.join(res_count_checks[st]['msg'] for st in segts['pdb'] if not res_count_checks[st]['check'])
        if not check_good:
            for pr in res_pdb['protein']:
                logger.debug(f'PDB PROTEIN RESIDUE {pr.chainID}_{pr.resname}{pr.resseqnum}{pr.insertion} resolved {pr.resolved}')
            for cr in res_cif['protein']:
                logger.debug(f'CIF PROTEIN RESIDUE {cr.chainID}_{cr.resname}{cr.resseqnum} [auth {cr.auth_asym_id}_{cr.auth_comp_id}{cr.auth_seq_id}{cr.insertion}] resolved {cr.resolved}')
        self.assertTrue(check_good,msg=check_msg)

    def test_molecule_existing(self):

        dest_psf = Path('existing.psf')
        if dest_psf.exists():
            dest_psf.unlink()
        dest_pdb = Path('existing.pdb')
        if dest_pdb.exists():
            dest_pdb.unlink()
        os.symlink(self.inputs_dir / 'existing.psf', dest_psf)
        os.symlink(self.inputs_dir / 'existing.pdb', dest_pdb)

        m=Molecule(source={
            'prebuilt':{
                'psf':'existing.psf',
                'pdb':'existing.pdb'
            }
            }, molid=0)

        au=m.asymmetric_unit
        ssbonds=au.objmanager.get('topol',{}).get('ssbonds',SSBondList([]))
        self.assertEqual(len(ssbonds),3)
        links=au.objmanager.get('topol',{}).get('links',LinkList([]))
        self.assertEqual(len(au.segments),3)
        self.assertEqual(len(au.segments.filter(lambda x: x.segtype=='protein')),1)
        Path('existing.psf').unlink()
        Path('existing.pdb').unlink()
