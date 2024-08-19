from pestifer.config import ResourceManager, Config
from pestifer.controller import Controller
from pestifer.charmmtop import getResis, makeBondGraph, getMasses, CharmmMasses
from pestifer.tasks import PsfgenTask
import glob
import os
import logging
import gzip
import networkx as nx
import yaml
logger=logging.getLogger(__name__)
from collections import UserList, UserDict

from rdkit.Chem import ForwardSDMolSupplier
from rdkit import Chem
from rdkit.Chem import AllChem

def make_lipid_database(self):
    r=ResourceManager()
    lmsdsdf=os.path.join(r['lmsd'],'structures.sdf.gz')
    self.assertTrue(os.path.exists(lmsdsdf))
    suppl=ForwardSDMolSupplier(gzip.open(lmsdsdf))
    lmid_idx=[]
    for mol in suppl:
        if mol is not None:
            mH=Chem.AddHs(mol)
            lmid_idx.append(dict(lm_id=mol.GetProp('LM_ID'),numatom=mH.GetNumAtoms(),mol=mH,mheavy=mol))
    toppardir=r['charmmff']['toppar']
    M=[]
    for f in glob.glob(os.path.join(toppardir,'top_*.rtf')):
        M.extend(getMasses(f))
    lipidstreamdir=os.path.join(toppardir,'stream/lipid')
    strlipids=glob.glob(os.path.join(lipidstreamdir,'*.str'))
    idx=['list' in x for x in strlipids].index(True)
    strlipids.remove(strlipids[idx])
    # logger.debug(f'strlipids in {", ".join([os.path.basename(x) for x in strlipids])}')
    tops=[os.path.join(toppardir,'top_all36_lipid.rtf')]+glob.glob(os.path.join(lipidstreamdir,'*.str'))
    for f in tops:
        M.extend(getMasses(f))
    M=CharmmMasses(M)
    resis=[]
    for t in tops:
        sublist=getResis(t,M)
        # logger.debug(f'{os.path.basename(t)} {len(sublist)} RESIs')
        resis.extend(sublist)
    logger.debug(f'{len(resis)} RESIs')
    mappings=[]
    for r in resis:
        r.candidates=[]
        for lm in lmid_idx:
            if lm['numatom']==r.num_atoms():
                mol=lm['mheavy']
                if mol.HasProp('FORMULA'):
                    if mol.GetProp('FORMULA')==r.formula():
                        r.candidates.append(lm)
                else:
                    logger.debug(f'no formula in lm mol; pre-matching to {r.resid} by atom count only ({r.numatoms()})')
                    r.candidates.append(lm)
                    # else:
                    #     if r.resid=='PSM': logger.debug(f'{mol.GetProp("FORMULA")} no match for {r.formula()}')
        # if r.resid=='PSM': logger.debug(f'{r.resid} {r.formula()} matches {len(candidates)} lm candidates')#, {r.synonym}, {r.num_atoms()}, {os.path.basename(r.topfile)}'')
        resig=r.to_graph(includeH=False)
        # logger.debug(f'{resig}')
        gmatch=None
        for candidate in r.candidates:
            mol=candidate['mheavy']
            molg,mol_atoms=makeBondGraph(mol)
            # logger.debug(f'  cf {molg}')
            rdict=nx.vf2pp_isomorphism(molg,resig)
            # logger.debug(f'{rdict}')
            if rdict:
                gmatch=candidate
                break
        if gmatch is not None:
            mol=gmatch['mol']
            props=mol.GetPropNames()
            p=['LM_ID','FORMULA','EXACT_MASS','SMILES']
            mp={k:mol.GetProp(k) for k in p if k in props}
            # logger.debug(f'{r.resid}, {r.formula()}, {r.mass():.6f} matches {mp}')
            mappings.append(dict(lmsd_id=mp['LM_ID'],charmm_resid=r.resid,charmm_topfile=r.topfile.split(toppardir)[1][1:],smiles=mp.get('SMILES','')))
        # else:
        #     logger.debug(f'{r.resid} has no lm matches out of {len(candidates)} candidates')
    logger.debug(f'{len(mappings)} successful matches of charmm lipid resids (out of {len(resis)}) to lmsd database')
    with open('lmsd_charmm.yaml','w') as f:
        yaml.dump(mappings,f)

class LMSDDatabase(UserList):
    def __init__(self):
        self.data=[]
        self.resources=ResourceManager()
        yaml_infile=os.path.join(self.resources['lmsd'],'lmsd_charmm.yaml')
        assert os.path.exists(yaml_infile),f'{yaml_infile}: not fould.'
        if yaml_infile:
            with open(yaml_infile,'r') as f:
                self.data=yaml.safe_load(f)
        self.D=UserDict({a['charmm_resid']:a for a in self.data})
        toppardir=self.resources['charmmff']['toppar']
        M=[]
        for f in glob.glob(os.path.join(toppardir,'top_*.rtf')):
            M.extend(getMasses(f))
        lipidstreamdir=os.path.join(toppardir,'stream/lipid')
        strlipids=glob.glob(os.path.join(lipidstreamdir,'*.str'))
        idx=['list' in x for x in strlipids].index(True)
        strlipids.remove(strlipids[idx])
        tops=[os.path.join(toppardir,'top_all36_lipid.rtf')]+glob.glob(os.path.join(lipidstreamdir,'*.str'))
        for f in tops:
            M.extend(getMasses(f))
        self.M=CharmmMasses(M)
        for t in tops:
            sublist=getResis(t,self.M)
            for resi in sublist:
                if resi.resname in self.D:
                    self.D[resi.resname]['charmmtop']=resi
                else:
                    self.D[resi.resname]={}
                    self.D[resi.resname]['charmmtop']=resi
            # logger.debug(f'{os.path.basename(t)} {len(sublist)} RESIs')
        self.charmm_resnames=list(sorted(list(self.D.keys())))

    def get_LMID(self,charmm_resid):
        elem=self.D.get(charmm_resid,{})
        if elem:
            return elem['lmsd_id']
        return ''
    
    def get_smiles(self,charmm_resid):
        elem=self.D.get(charmm_resid,{})
        if elem:
            return elem['smiles']
        return ''
    
    def get_topo(self,charmm_resid):
        elem=self.D.get(charmm_resid,{})
        if elem:
            return elem['charmmtop']
        return None
    
    def yaml_names(self):
        with open('lipid_charmm_names.yaml','w') as f:
            yaml.dump(self.charmm_resnames,f)
    
def make_charmm_pdb(resid,DB):
    smiles=DB.get_smiles(resid)
    topo=DB.get_topo(resid)
    if smiles:
        mol=Chem.MolFromSmiles(smiles)
        if mol is not None and topo is not None:
            molH=Chem.AddHs(mol)
            mol_graph,mol_atoms=makeBondGraph(molH)
            top_graph=topo.to_graph()
            map=nx.vf2pp_isomorphism(mol_graph,top_graph)
            assert not len(map)==0,f'Failure in graph isomorphism'
            params=AllChem.ETKDGv3()
            params.randomSeed=0xf00d
            result=AllChem.EmbedMolecule(molH,params)
            logger.debug(f'embed result {result}')
            if result==-1:
                cids=AllChem.EmbedMultipleConfs(molH,numConfs=30)
                logger.debug(f'#cids {len(cids)}')
                if len(cids)==0:
                    logger.error(f'Could not embed {resid}')
            if map:
                new_names=[map[x] for x in mol_atoms]
                for a,n in zip(molH.GetAtoms(),new_names):
                    resinfo=Chem.rdchem.AtomPDBResidueInfo()
                    resinfo.SetName(n)
                    resinfo.SetResidueName(resid)
                    resinfo.SetResidueNumber(1)
                    resinfo.SetChainId('A')
                    a.SetPDBResidueInfo(resinfo)
                w=Chem.PDBWriter(f'{resid}.pdb')
                w.write(molH,confId=0)
                w.flush()
                w.close()
                badresid=resid[:3]+' '
                logger.debug(f'badresid [{badresid}]')
                with open(f'{resid}.pdb','r') as f:
                    dat=f.read()
                    logger.debug(f'# badresid {dat.count(badresid)}')
                    dat=dat.replace(badresid,resid)
                with open(f'{resid}.pdb','w') as f:
                    f.write(dat)
                return 0
    return -1

def do_psfgen(resid,DB):
    if make_charmm_pdb(resid,DB)==0:
        logger.debug(f'{resid} {DB.D[resid]}')
        orient_par=DB.D[resid].get('align_along_z',[])
        C=Controller(userspecs={'tasks':[
                        {'psfgen':{
                            'source':{
                                'id':resid,
                                'file_format':'PDB',
                                'biological_assembly':0
                                },
                            'mods':{}
                            }},
                        {'md': {'ensemble':'minimize','nsteps':0,'minimize':1000,'dcdfreq':0,'xstfreq':0,'temperature':100}},
                        {'manipulate':{'mods':{'orient':[f'z,{orient_par[0]},{orient_par[1]}']}}},
                        {'md': {'ensemble':'NVT','nsteps':10000,'dcdfreq':100,'xstfreq':100,'temperature':100,'other_parameters':{'cylindricalbc':'on','cylindricalbcaxis':"z",'cylindricalbccenter':"0.0,0.0,0.0",'cylindricalbcl1':100,'cylindricalbcr1':0.1,'cylindricalbck1':0.1}}},
                ]
            }
        )
        C.do_tasks()
    return 0