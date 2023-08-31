import logging
logger=logging.getLogger(__name__)
from .basemod import AncestorAwareMod, AncestorAwareModList
from .mods import MutationList
# from .config import ConfigGetParam
from .residue import Residue,ResidueList
from .util import isidentity, reduce_intlist
from .scriptwriters import Psfgen

class Segment(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['segtype','segname','chainID','residues','subsegments','parent_chain','resname_charmify']
    opt_attr=AncestorAwareMod.opt_attr+['mutations','deletions','grafts','attachments','psfgen_segname','config_params']
    config_attr=['Include_terminal_loops','Fix_engineered_mutations','Fix_conflicts','Sacrificial_residue','PDB_to_CHARMM_Resnames']

    def __init__(self,parm_dict,input_obj,segname):
        assert all([x in parm_dict for x in self.config_attr])
        if type(input_obj)==dict:
            input_dict=input_obj
        elif type(input_obj)==Residue:
            res=input_obj
            apparent_chainID=res.chainID
            apparent_segtype=res.segtype
            myRes=ResidueList([])
            myRes.append(res.clone())
            subsegments=myRes.state_bounds(lambda x: 'RESOLVED' if len(x.atoms)>0 else 'MISSING')
            # olc=ConfigGetParam('Segname_chars').get(apparent_segtype)
            input_dict={
                'segtype': apparent_segtype,
                'segname': segname,
                'chainID': apparent_chainID,
                'residues': myRes,
                'subsegments':subsegments,
                'parent_chain': apparent_chainID,
                'resname_charmify':parm_dict['PDB_to_CHARMM_Resnames']
            }
        elif type(input_obj)==ResidueList:
            Residues=input_obj
            apparent_chainID=Residues[0].chainID
            apparent_segtype=Residues[0].segtype
            myRes=Residues.clone()
            if apparent_segtype=='PROTEIN':
                # a protein segment must have unique residue numbers
                assert myRes.puniq(['resseqnum','insertion'])
                # a protein segment may not have more than one protein chain
                assert all([x.chainID==myRes[0].chainID for x in myRes])
            else:
                myRes.puniquify(fields=['resseqnum','insertion'],make_common=['chainID'])
            myRes.sort()
            subsegments=myRes.state_bounds(lambda x: 'RESOLVED' if len(x.atoms)>0 else 'MISSING')
            # olc=ConfigGetParam('Segname_chars').get(apparent_segtype)
            input_dict={
                'segtype': apparent_segtype,
                'segname': segname,
                'chainID': apparent_chainID,
                'residues': myRes,
                'subsegments':subsegments,
                'parent_chain': apparent_chainID,
                'resname_charmify':parm_dict['PDB_to_CHARMM_Resnames']
            }
        else:
            logger.error(f'Cannot initialize {self.__class__} from object type {type(input_obj)}')
        input_dict['config_params']=parm_dict
        super().__init__(input_dict)

    def has_graft(self,modlist):
        if not modlist:
            return False
        for g in modlist:
            if g.segname==self.psfgen_segname:
                return True
        return False
    
    def __str__(self):
        return f'{self.segname}: type {self.segtype} chain {self.chainID} with {len(self.residues)} residues'

    def write_TcL(self,W:Psfgen,transform,mods):
        if self.segtype=='PROTEIN':
            self.protein_stanza(W,transform,mods)
        elif self.segtype=='GLYCAN':
            self.glycan_stanza(W,transform,mods)
        else:
            self.generic_stanza(W,transform,mods)
    
    def glycan_stanza(self,W,transform,mods):
        if self.has_graft(mods.get('Grafts',[])):
            W.comment('Glycan grafts not yet implemented.')
        else:
            self.generic_stanza(W,transform,mods)
    
    def generic_stanza(self,W,transform,mods):
        assert len(self.subsegments)==1,'No missing atoms allowed in generic stanza'
        parent_molecule=self.ancestor_obj
        chainIDmap=transform.chainIDmap
        seglabel=self.segname
        image_seglabel=chainIDmap.get(seglabel,seglabel)
        not_asymm=image_seglabel!=seglabel
        # if rep_chainID!=asym_chainID:
        #     seglabel=f'{rep_chainID}{seglabel[1:]}' # first byte is chain ID
        #     logger.info(f'relabeling chain {the_chainID} to {rep_chainID} in segment {self.psfgen_segname} -> new segment {seglabel}')
        transform.register_mapping(self.segtype,image_seglabel,seglabel)

        W.banner(f'BEGIN SEGMENT {image_seglabel}')
        W.comment(f'This segment is referenced as chain {seglabel} in the asymmetric unit (possibly after segtype-processing by pestifer)')
        serial_list=self.residues.atom_serials(as_type=int)
        vmd_red_list=reduce_intlist(serial_list)
        pdb=f'GENERIC_{image_seglabel}.pdb'
        selname=image_seglabel
        W.addfile(pdb)
        W.addline(f'set {selname} [atomselect ${parent_molecule.molid_varname} "serial {vmd_red_list}"]')
        W.addline(f'${selname} set segname {image_seglabel}')
        if not_asymm:
            W.backup_selection(selname)
            W.addline(f'${selname} set chain {image_seglabel}')                 
            W.addline(f'${selname} move {transform.write_TcL()}')
        W.addline(f'${selname} writepdb {pdb}')
        W.addline(f'segment {image_seglabel} '+'{')
        W.addline(f'    pdb {pdb}')
        W.addline('}')
        W.addline(f'coordpdb {pdb} {image_seglabel}')
        if not_asymm:
            W.banner(f'Restoring A.U. state for {seglabel}')
            W.restore_selection(selname)
        W.banner(f'END SEGMENT {image_seglabel}')

    def protein_stanza(self,W:Psfgen,transform,mods):
        parent_molecule=self.ancestor_obj
        chainIDmap=transform.chainIDmap
        seglabel=self.segname
        image_seglabel=chainIDmap.get(seglabel,seglabel)
        not_asymm=image_seglabel!=seglabel
        # if rep_chainID!=asym_chainID:
        #     seglabel=f'{rep_chainID}{seglabel[1:]}' # first byte is chain ID
        #     logger.info(f'relabeling chain {the_chainID} to {rep_chainID} in segment {self.psfgen_segname} -> new segment {seglabel}')
        transform.register_mapping(self.segtype,image_seglabel,seglabel)

        sac_rd=self.config_params['Sacrificial_residue']
        sac_r=sac_rd['name']
        sac_n=sac_rd['min_loop_length']
        # intra-segment, no need to use image_seglabel
        Mutations=mods.get('Mutations',MutationList([])).filter(chainID=seglabel)
        Conflicts=mods.get('Conflicts',MutationList([])).filter(chainID=seglabel)

        W.banner(f'BEGIN SEGMENT {image_seglabel}')
        for i,b in enumerate(self.subsegments):
            if b.state=='RESOLVED':
                """ for a resolved subsegment, generate its pdb file """
                b.selname=f'{image_seglabel}{i:02d}'
                run=ResidueList(self.residues[b.bounds[0]:b.bounds[1]+1])
                b.pdb=f'PROTEIN_{image_seglabel}_{run[0].resseqnum}{run[0].insertion}_to_{run[-1].resseqnum}{run[-1].insertion}.pdb'
                W.addfile(b.pdb)
                # TODO: account for any deletions?  Maybe better at the residue stage...
                serial_list=run.atom_serials(as_type=int)
                at=parent_molecule.asymmetric_unit.Atoms.get(serial=serial_list[-1])
                assert at.resseqnum==run[-1].resseqnum
                vmd_red_list=reduce_intlist(serial_list)
                W.addline(f'set {b.selname} [atomselect ${parent_molecule.molid_varname} "serial {vmd_red_list}"]')
                W.addline(f'${b.selname} set segname {image_seglabel}')
                if hasattr(at,'_ORIGINAL_'):
                    if at._ORIGINAL_["serial"]!=at.serial:
                        W.banner(f'Atom with serial {at._ORIGINAL_["serial"]} in PDB needs serial {at.serial} for VMD')
                """ Relabel chain ID and request coordinate transformation """
                if not_asymm:
                    W.backup_selection(b.selname)
                    W.addline(f'${b.selname} set chain {image_seglabel}')                 
                    W.addline(f'${b.selname} move {transform.write_TcL()}')
                W.addline(f'${b.selname} writepdb {b.pdb}')
        W.addline(f'segment {image_seglabel} '+'{')
        # print(len(self.subsegments),type(self.subsegments))
        if not self.config_params['Include_terminal_loops']:
            if self.subsegments[0].state=='MISSING':
                Nterminal_missing_subsegment=self.subsegments.pop(0)
                # logger.info(f'Since terminal loops are not included, ignoring {str(Nterminal_missing_subsegment)}')
            if self.subsegments[-1].state=='MISSING':
                Cterminal_missing_subsegment=self.subsegments.pop(-1)
                # logger.info(f'Since terminal loops are not included, ignoring {str(Cterminal_missing_subsegment)}')
        for b in self.subsegments:
            if b.state=='RESOLVED':
                W.addline(f'    pdb {b.pdb}')
            elif b.state=='MISSING':
                for r in self.residues[b.bounds[0]:b.bounds[1]+1]:
                    rname=self.resname_charmify.get(r.name,r.name)
                    W.addline(f'    residue {r.resseqnum}{r.insertion} {rname} {image_seglabel}')
                if (b.bounds[1]-b.bounds[0])>(sac_n-1):
                    lrr=self.residues[b.bounds[1]]
                    sac_resseqnum=lrr.resseqnum
                    sac_insertion='A' if lrr.insertion in [' ',''] else chr(ord(lrr.insertion)+1)
                    assert sac_insertion<='Z',f'Residue {lrr.resseqnum} of chain {seglabel} already has too many insertion instances (last: {lrr.insertion}) to permit insertion of a sacrificial {sac_r}'
                    b.sacres=Residue({'name':sac_r,'resseqnum':sac_resseqnum,'insertion':sac_insertion,'chainID':seglabel,'segtype':'PROTEIN'})
                    W.addline(f'    residue {sac_resseqnum}{sac_insertion} {sac_r} {image_seglabel}')
        if self.config_params['Fix_engineered_mutations']:
            for m in Mutations:
                if m.source=='SEQADV':
                    if self.config_params['Fix_engineered_mutations']:
                        W.comment('Below reverts an engineered mutation:')
                        W.addline(m.write_TcL())
                else:
                    W.comment('Below is a user-specified mutation:')
                    W.addline(m.write_TcL())
        if self.config_params['Fix_conflicts']:
            for m in Conflicts:
                W.comment('Below reverts a database conflict:')
                W.addline(m.write_TcL())
        W.addline('}')
        W.banner(f'End segment {image_seglabel}')
        W.banner('Coordinate-specification commands')
        for b in self.subsegments:
            if b.state=='RESOLVED':
                W.comment(f'Subsegment {[self.subsegments.index(b)]} is a resolved run')
                W.addline(f'coordpdb {b.pdb} {image_seglabel}')
        for b in self.subsegments:
            if b.state=='MISSING':
                # will issue the atom-reorienting command to join the C-terminus of prior run to N-terminus of this one, which is model-built using guesscoord
                if (self.subsegments.index(b)==0 or self.subsegments.index(b)==(len(self.subsegments)-1)) and not self.config_params['Include_terminal_loops']:
                    # this is a terminal loop and we are not including terminal loops
                    pass
                else:
                    # not terminal OR we ARE including terminal loops
                    if self.subsegments.index(b)>0:
                        # this is either interior OR C-terminal to be included
                        W.comment(f'Subsegment {[self.subsegments.index(b)]}/{len(self.subsegments)} is a missing loop')
                        this_run=ResidueList(self.residues[b.bounds[0]:b.bounds[1]+1])
                        prior_b=self.subsegments[self.subsegments.index(b)-1]
                        W.comment(f'...attached to subsegment {self.subsegments.index(prior_b)}')
                        prior_run=ResidueList(self.residues[prior_b.bounds[0]:prior_b.bounds[1]+1])
                        W.comment(f'Seeding orientation of model-built loop starting at {str(this_run[0])} from {str(prior_run[-1])}')
                        W.addline(f'{this_run.caco_str(prior_run,image_seglabel,parent_molecule.molid_varname)}')
        W.banner('Intra-segmental terminal patches')
        for i,b in enumerate(self.subsegments):
            if b.state=='MISSING' and i>0 and hasattr(b,'sacres'):
                Cterm=self.residues[b.bounds[1]]
                W.addline(f'patch CTER {image_seglabel}:{Cterm.resseqnum}{Cterm.insertion}')
                nextb=self.subsegments[i+1]
                Nterm=self.residues[nextb.bounds[0]]
                patchname='NTER'
                if Nterm.name=='PRO':
                    patchname='PROP'
                elif Nterm.name=='GLY':
                    patchname='GLYP'
                W.addline(f'patch {patchname} {image_seglabel}:{Nterm.resseqnum}{Nterm.insertion}')
                W.addline(f'delatom {image_seglabel} {b.sacres.resseqnum}{b.sacres.insertion}')
        W.banner('Restoring A.U. state for all resolved subsegments')
        for b in self.subsegments:
            if b.state=='RESOLVED':
                if not_asymm:
                    W.restore_selection(b.selname)
        W.banner(f'END SEGMENT {image_seglabel}')

class SegmentList(AncestorAwareModList):

    # counters_by_segtype={}
    # def append(self,item:Segment):
    #     # olc=self.config['Segname_chars'][item.segtype]
    #     # itemkey=f'{item.chainID}{olc}'
    #     if not itemkey in self.counters_by_segtype:
    #         self.counters_by_segtype[itemkey]=0
    #     self.counters_by_segtype[itemkey]+=1
    #     if item.segtype=='PROTEIN':
    #         item.psfgen_segname=itemkey
    #     else:
    #         item.psfgen_segname=f'{itemkey}{self.counters_by_segtype[itemkey]:1d}'
    #     assert item.psfgen_segname not in self.segnames
    #     self.segnames.append(item.psfgen_segname)
    #     item.resname_charmify=self.config['PDB_to_CHARMM_Resnames']
    #     super().append(item)

    def __init__(self,*objs):
        if len(objs)==1 and objs[0]==[]:
            super().__init__([])
        elif len(objs)==3:
            config=objs[0]
            input_obj=objs[1]
            chainIDmanager=objs[2]
            self.config={k:config.get(k,'no-config') for k in Segment.config_attr}
            self.counters_by_segtype={}
            self.segnames=[]
            self.daughters={}
            if type(input_obj)==ResidueList:
                super().__init__([])
                residues=input_obj
                for stype in config['Segtypes']:
                    self.counters_by_segtype[stype]=0
                    res=residues.filter(segtype=stype)
                    for chainID in res.unique_chainIDs():
                        this_chainID=chainID
                        c_res=residues.filter(segtype=stype,chainID=this_chainID)
                        if chainID in self.segnames:
                            logger.debug(f'chain {chainID} has already been used as a segname')
                            this_chainID=chainIDmanager.next_chain()
                            if not chainID in self.daughters:
                                self.daughters[chainID]=[]
                            self.daughters[chainID].append(this_chainID)
                            c_res.set(chainID=this_chainID)
                            for r in c_res:
                                r.atoms.set(chainID=this_chainID)
                        thisSeg=Segment(self.config,c_res,segname=this_chainID)
                        self.append(thisSeg)
                        logger.debug(f'Making segments: stype {stype} chainID {chainID} segname {thisSeg.segname}')
                        self.segnames.append(thisSeg.segname)
                        self.counters_by_segtype[stype]+=1
            else:
                logger.error(f'Cannot initialize {self.__class__} from objects {objs}')

    def collect_working_files(self):
        working_files=[]
        for seg in self:
            for subseg in seg.subsegments:
                if hasattr(subseg,'pdb'):
                    logger.debug(f'wf: {subseg.pdb}')
                    working_files.append(subseg.pdb)
        return working_files
    
    def write_TcL(self,W:Psfgen,transform,mods):
        for seg in self:
            W.comment(f'Writing seg {seg.segname}')
            seg.write_TcL(W,transform,mods)
            W.comment(f'Done with seg {seg.segname}')
    
    def get_segment_of_residue(self,residue):
        for S in self:
            if residue in S.residues:
                return S
        return None
