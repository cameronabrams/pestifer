"""

.. module:: molecule
   :synopsis: Manages all molecules
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
from pidibble.pdbparse import PDBParser
from .basemod import AncestorAwareMod
from .asymmetricunit import AsymmetricUnit
from .bioassemb import BioAssembList,BioAssemb
from .stringthings import ByteCollector
from .mods import MutationList, apply_psf_info
logger=logging.getLogger(__name__)

class Molecule(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['config','molid','source','asymmetric_unit','biological_assemblies','parsed_struct']
    opt_attr=AncestorAwareMod.opt_attr+['active_biological_assembly']
    _molcounter=0
    def __init__(self,source=None,config=None,**options):
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
            p_struct=PDBParser(PDBcode=source).parse().parsed
        use_psf=options.get('use_psf',None)
        if use_psf:
            apply_psf_info(p_struct,f'{source}.psf')
        # if os.path.exists(f'{source}.psf'):
        #     apply_psf_info(p_struct,f'{source}.psf')
        input_dict={
            'config': config,
            'molid': Molecule._molcounter,
            'source': source,
            'parsed_struct': p_struct,
            'asymmetric_unit': AsymmetricUnit(p_struct,config),
            'biological_assemblies': BioAssembList(p_struct)#,reset=True)
        }
        super().__init__(input_dict)
        self.asymmetric_unit.claim_descendants(self,0)
        self.biological_assemblies.claim_descendants(self,0)
    
    def activate_biological_assembly(self,index,chainIDmanager):
        if index==0 or len(self.biological_assemblies)==0: # we will use the unadulterated A.U. as the B.A.
            self.active_biological_assembly=BioAssemb(self.asymmetric_unit)
            if index!=0:
                logger.warning(f'No biological assemblies specified in input structure.  Using A.U.; ignoring your request of B.A. "{index}"')
        else:
            biological_assembly=[x for x in self.biological_assemblies if x.index==index]
            assert biological_assembly!=[],f'No biological assembly "{index:d}" found.'
            assert len(biological_assembly)==1,f'No unique biological assembly "{index}" found.'
            self.active_biological_assembly=biological_assembly[0]
            logger.info(f'Activating biological assembly {self.active_biological_assembly.name} (idx {index})')
        ba=self.active_biological_assembly
        auChainIDs=self.asymmetric_unit.chainIDs
        ba.biomt[0].chainIDmap={c:c for c in auChainIDs}
        for biomt in ba.biomt[1:]:
            biomt.chainIDmap=chainIDmanager.generate_next_map(auChainIDs)            
        return self
    
    def write_TcL(self,B:ByteCollector,user_mods,file_collector=None):
        au=self.asymmetric_unit
        ba=self.active_biological_assembly
        config=self.config
        allmods=user_mods.copy()
        if not 'Mutations' in allmods:
            allmods['Mutations']=MutationList([])
        if not 'Conflicts' in allmods:
            allmods['Conflicts']=MutationList([])
        for k in ba.biomt[0].chainIDmap.keys():
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
        for biomt in ba.biomt:
            B.banner(f'TRANSFORM {biomt.index} BEGINS')
            B.banner('The following mappings of A.U. chains is used:')
            for k,v in biomt.chainIDmap.items():
                B.comment(f'{k}: {v}')
            B.banner('SEGMENTS FOLLOW')
            au.Segments.write_TcL(B,biomt,allmods,file_collector=file_collector)
            B.banner('DISU PATCHES FOLLOW')
            au.SSBonds.write_TcL(B,biomt,allmods)
            B.banner('LINK PATCHES FOLLOW')
            au.Links.write_TcL(B,biomt,allmods)
            B.banner(f'TRANSFORM {biomt.index} ENDS')

