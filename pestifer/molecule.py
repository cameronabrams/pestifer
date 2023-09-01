"""

.. module:: molecule
   :synopsis: Manages all molecules
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
from pidibble.pdbparse import PDBParser
from .cifutil import CIFload
from .basemod import AncestorAwareMod
from .asymmetricunit import AsymmetricUnit
from .bioassemb import BioAssembList,BioAssemb
from .scriptwriters import Psfgen
from .chainids import ChainIDManager
from .mods import MutationList, apply_psf_info
logger=logging.getLogger(__name__)

class Molecule(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['config','molid','chainIDmanager','source','asymmetric_unit','biological_assemblies','parsed_struct','rcsb_file_format']
    opt_attr=AncestorAwareMod.opt_attr+['active_biological_assembly']
    _molcounter=0
    def __init__(self,source=None,config=None,chainIDmanager=None,excludes={},**options):
        reset=options.get('reset_counter',False)
        if reset:
            Molecule._molcounter=0
        else:
            Molecule._molcounter+=1
        if not config:
            logger.warning('Molecule initialized without config.')
        if not source:
            logger.debug('Molecule initialized without source.')
            p_struct=None
        else:
            if config['rcsb_file_format']=='PDB':
                p_struct=PDBParser(PDBcode=source).parse().parsed
            elif config['rcsb_file_format']=='mmCIF':
                logger.debug(f'CIF source {source}')
                p_struct=CIFload(source)
                logger.debug(f'p_struct type {type(p_struct)}')
        use_psf=options.get('use_psf',None)
        if use_psf:
            apply_psf_info(p_struct,f'{source}.psf')
        # if os.path.exists(f'{source}.psf'):
        #     apply_psf_info(p_struct,f'{source}.psf')
        if chainIDmanager==None:
            logger.debug(f'Molecule instantiating its own ChainIDManager')
            chainIDmanager=ChainIDManager()
        input_dict={
            'config': config,
            'chainIDmanager':chainIDmanager,
            'rcsb_file_format': config['rcsb_file_format'],
            'molid': Molecule._molcounter,
            'source': source,
            'parsed_struct': p_struct,
            'asymmetric_unit': AsymmetricUnit(p_struct,config,chainIDmanager,excludes),
            'biological_assemblies': BioAssembList(p_struct)
        }
        super().__init__(input_dict)
        self.asymmetric_unit.claim_descendants(self,0)
        self.biological_assemblies.claim_descendants(self,0)
    
    def activate_biological_assembly(self,index):
        if index==0 or len(self.biological_assemblies)==0: # we will use the unadulterated A.U. as the B.A.
            self.active_biological_assembly=BioAssemb(self.asymmetric_unit)
            if index!=0:
                logger.warning(f'No biological assemblies specified in input structure.  Using A.U.; ignoring your request of B.A. "{index}"')
        else:
            self.active_biological_assembly=self.biological_assemblies.get(index=index)
            assert self.activate_biological_assembly!=None,f'No biological assembly "{index:d}" found.'
            logger.info(f'Activating biological assembly {self.active_biological_assembly.name} (idx {index})')
        self.active_biological_assembly.activate(self.asymmetric_unit,self.chainIDmanager)
        return self
    
    def write_TcL(self,W:Psfgen,user_mods):
        au=self.asymmetric_unit
        ba=self.active_biological_assembly
        config=self.config
        allmods=user_mods.copy()
        if not 'Mutations' in allmods:
            allmods['Mutations']=MutationList([])
        if not 'Conflicts' in allmods:
            allmods['Conflicts']=MutationList([])
        for k in ba.transforms[0].chainIDmap.keys():
            if au.Mutations:
                allmods['Mutations'].extend(au.Mutations.filter(chainID=k))
            if au.Conflicts:
                allmods['Conflicts'].extend(au.Conflicts.filter(chainID=k))
        user_mutations=user_mods.get('Mutations',MutationList([]))
        au.SSBonds.prune_mutations(user_mutations)
        au.Links.prune_mutations(user_mutations,au.Segments)
        if config['Fix_engineered_mutations']:
            au.SSBonds.prune_mutations(au.Mutations)
            au.Links.prune_mutations(au.Mutations,au.Segments)
            # au.Links.prune_mutations(au.Conflicts) # problematic?
        if config['Fix_conflicts']:
            au.SSBonds.prune_mutations(au.Conflicts)
            au.Links.prune_mutations(au.Conflicts,au.Segments)
            # au.Links.prune_mutations(au.Conflicts) # problematic?
        for transform in ba.transforms:
            W.banner(f'TRANSFORM {transform.index} BEGINS')
            W.banner('The following mappings of A.U. asym ids is used:')
            for k,v in transform.chainIDmap.items():
                W.comment(f'{k}: {v}')
            W.banner('SEGMENTS FOLLOW')
            au.Segments.write_TcL(W,transform,allmods)
            W.banner('DISU PATCHES FOLLOW')
            au.SSBonds.write_TcL(W,transform,allmods)
            W.banner('LINK PATCHES FOLLOW')
            au.Links.write_TcL(W,transform,allmods)
            W.banner(f'TRANSFORM {transform.index} ENDS')
    
    def get_chainmaps(self):
        ba=self.active_biological_assembly
        maps={}
        for transform in ba.transforms:
            for oc,mc in transform.chainIDmap.items():
                if not oc in maps:
                    maps[oc]=[]
                maps[oc].append({'transform':transform.index,'mappedchain':mc})
        return maps
    
    def has_loops(self,min_length=4):
        nloops=0
        au=self.asymmetric_unit
        for S in au.Segments:
            # chainID=S.chainID
            if S.segtype=='PROTEIN':
                for b in S.subsegments:
                    if b.state=='MISSING' and b.num_residues()>=min_length:
                        nloops+=1
        return nloops

    def write_loop_lines(self,writer,cycles=100,min_length=4):
        ba=self.active_biological_assembly
        au=self.asymmetric_unit
        for S in au.Segments:
            # chainID=S.chainID
            if S.segtype=='PROTEIN':
                asymm_segname=S.segname
                for b in S.subsegments:
                    if b.state=='MISSING':
                        if b.num_residues()>=min_length:
                            reslist=[f'{r.resseqnum}{r.insertion}' for r in S.residues[b.bounds[0]:b.bounds[1]+1]]
                            tcllist='[list '+' '.join(reslist)+']'
                            for transform in ba.transforms:
                                cm=transform.chainIDmap
                                act_segID=cm.get(asymm_segname,asymm_segname)
                                writer.addline(f'declash_loop $mLL {act_segID} {tcllist} {cycles}')
    
    def write_gaps(self,writer,min_length=4):
        ba=self.active_biological_assembly
        au=self.asymmetric_unit
        writer.addline('# fields: segname loop-begin-res loop-end-res connect-to-res')
        for S in au.Segments:
            # chainID=S.chainID
            if S.segtype=='PROTEIN':
                asymm_segname=S.segname
                for i,b in enumerate(S.subsegments):
                    if b.state=='MISSING':
                        if b.num_residues()>=min_length and i<(len(S.subsegments)-1):
                            reslist=[f'{r.resseqnum}{r.insertion}' for r in S.residues[b.bounds[0]:b.bounds[1]+1]]
                            bpp=S.subsegments[i+1]
                            nreslist=[f'{r.resseqnum}{r.insertion}' for r in S.residues[bpp.bounds[0]:bpp.bounds[1]+1]]
                            assert bpp.state=='RESOLVED'
                            for transform in ba.transforms:
                                cm=transform.chainIDmap
                                act_segID=cm.get(asymm_segname,asymm_segname)
                                writer.addline(f'{act_segID} {reslist[0]} {reslist[-1]} {nreslist[0]}')

    def write_connect_patches(self,writer,min_length=4):
        ba=self.active_biological_assembly
        au=self.asymmetric_unit
        for S in au.Segments:
            # chainID=S.chainID
            if S.segtype=='PROTEIN':
                asymm_segname=S.segname
                for i,b in enumerate(S.subsegments):
                    if b.state=='MISSING':
                        if b.num_residues()>=min_length and i<(len(S.subsegments)-1):
                            llres=S.residues[b.bounds[1]-1]
                            lres=S.residues[b.bounds[1]]
                            nextb=S.subsegments[i+1]
                            rres=S.residues[nextb.bounds[0]]
                            rrres=S.residues[nextb.bounds[0]+1]
                            for transform in ba.transforms:
                                cm=transform.chainIDmap
                                act_segID=cm.get(asymm_segname,asymm_segname)
                                writer.addline(f'patch HEAL {act_segID}:{llres.resseqnum}{llres.insertion} {act_segID}:{lres.resseqnum}{lres.insertion} {act_segID}:{rres.resseqnum}{rres.insertion} {act_segID}:{rrres.resseqnum}{rrres.insertion}')
