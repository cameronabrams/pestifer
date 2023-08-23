import logging
logger=logging.getLogger(__name__)
from .basemod import AncestorAwareMod, AncestorAwareModList
from .mods import MutationList
# from .config import ConfigGetParam
from .residue import Residue,ResidueList
import pestifer.sel as sel
from .util import isidentity, reduce_intlist
from .stringthings import ByteCollector

class Segment(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['segtype','segname','chainID','residues','subsegments']
    opt_attr=AncestorAwareMod.opt_attr+['mutations','deletions','grafts','attachments','psfgen_segname','config_params']
    config_attr=['Include_terminal_loops','Fix_engineered_mutations','Fix_conflicts','Sacrificial_residue','Segname_chars','PDB_to_CHARMM_Resnames']

    def __init__(self,parm_dict,input_obj):
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
                'segname': apparent_chainID,
                'chainID': apparent_chainID,
                'residues': myRes,
                'subsegments':subsegments
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
                'segname': apparent_chainID,
                'chainID': apparent_chainID,
                'residues': myRes,
                'subsegments':subsegments
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

    def write_TcL(self,B:ByteCollector,transform,mods,file_collector=None):
        if self.segtype=='PROTEIN':
            self.protein_stanza(B,transform,mods,file_collector=file_collector)
        elif self.segtype=='GLYCAN':
            self.glycan_stanza(B,transform,mods,file_collector=file_collector)
        else:
            self.generic_stanza(B,transform,mods,file_collector=file_collector)
    
    def glycan_stanza(self,B,transform,mods,file_collector=None):
        if self.has_graft(mods.get('Grafts',[])):
            B.comment('Glycan grafts not yet implemented.')
        else:
            self.generic_stanza(B,transform,mods,file_collector=file_collector)
    
    def generic_stanza(self,B,transform,mods,file_collector=None):
        assert len(self.subsegments)==1,'No missing atoms allowed in generic stanza'
        parent_molecule=self.ancestor_obj
        the_chainID=self.residues[0].chainID
        chainIDmap=transform.chainIDmap
        rep_chainID=chainIDmap.get(the_chainID,the_chainID)
        seglabel=self.psfgen_segname
        if rep_chainID!=the_chainID:
            seglabel=f'{rep_chainID}{seglabel[1:]}' # first byte is chain ID
            logger.info(f'relabeling chain {the_chainID} to {rep_chainID} in {self.psfgen_segname} -> {seglabel}')
        transform.register_mapping(self.segtype,the_chainID,seglabel)

        B.banner(f'BEGIN SEGMENT {seglabel}')
        B.comment(f'This segment is referenced as chain {the_chainID} in the input structure')
        serial_list=self.residues.atom_serials(as_type=int)
        vmd_red_list=reduce_intlist(serial_list)
        pdb=f'GENERIC_{seglabel}.pdb'
        selname=seglabel
        file_collector.append(pdb)
        B.addline(f'set {selname} [atomselect ${parent_molecule.molid_varname} "serial {vmd_red_list}"]')
        if rep_chainID!=the_chainID:
            B.addline(sel.backup(f'{selname}'))
            B.addline(f'${selname} set chain {rep_chainID}')                 
            B.addline(f'${selname} move {transform.write_TcL()}')
        B.addline(f'${selname} writepdb {pdb}')
        B.addline(f'segment {seglabel} '+'{')
        B.addline(f'    pdb {pdb}')
        B.addline('}')
        B.addline(f'coordpdb {pdb} {seglabel}')
        if rep_chainID!=the_chainID:
            B.banner(f'Restoring A.U. state for {seglabel}')
            B.addline(sel.restore(f'{selname}'))
        B.banner(f'END SEGMENT {seglabel}')

    def protein_stanza(self,B:ByteCollector,transform,mods,file_collector=None):
        parent_molecule=self.ancestor_obj
        the_chainID=self.residues[0].chainID
        chainIDmap=transform.chainIDmap
        tmat=transform.tmat
        rep_chainID=chainIDmap.get(the_chainID,the_chainID)
        seglabel=self.psfgen_segname # protein, so should just be letter
        if rep_chainID!=the_chainID:
            logger.info(f'relabeling chain {the_chainID} to {rep_chainID} in {seglabel}')
            seglabel=f'{rep_chainID}'# protein segment names are chainIDs
        transform.register_mapping(self.segtype,the_chainID,seglabel)

        sac_rd=self.config_params['Sacrificial_residue']
        sac_r=sac_rd['name']
        sac_n=sac_rd['min_loop_length']
        Mutations=mods.get('Mutations',MutationList([])).filter(chainID=the_chainID)
        Conflicts=mods.get('Conflicts',MutationList([])).filter(chainID=the_chainID)

        B.banner(f'BEGIN SEGMENT {seglabel}')
        for i,b in enumerate(self.subsegments):
            if b.state=='RESOLVED':
                """ for a resolved subsegment, generate its pdb file """
                b.selname=f'{seglabel}{i:02d}'
                run=ResidueList(self.residues[b.bounds[0]:b.bounds[1]+1])
                b.pdb=f'PROTEIN_{rep_chainID}_{run[0].resseqnum}{run[0].insertion}_to_{run[-1].resseqnum}{run[-1].insertion}.pdb'
                file_collector.append(b.pdb)
                # TODO: account for any deletions?  Maybe better at the residue stage...
                serial_list=run.atom_serials(as_type=int)
                at=parent_molecule.asymmetric_unit.Atoms.get(serial=serial_list[-1])
                assert at.resseqnum==run[-1].resseqnum
                vmd_red_list=reduce_intlist(serial_list)
                B.addline(f'set {b.selname} [atomselect ${parent_molecule.molid_varname} "serial {vmd_red_list}"]')
                if hasattr(at,'_ORIGINAL_'):
                    if at._ORIGINAL_["serial"]!=at.serial:
                        B.banner(f'Atom with serial {at._ORIGINAL_["serial"]} in PDB needs serial {at.serial} for VMD')
                """ Relabel chain ID and request coordinate transformation """
                if rep_chainID!=the_chainID:
                    B.addline(sel.backup(f'{b.selname}'))
                    B.addline(f'${b.selname} set chain {rep_chainID}')                 
                    B.addline(f'${b.selname} move {transform.write_TcL()}')
                B.addline(f'${b.selname} writepdb {b.pdb}')
        B.addline(f'segment {seglabel} '+'{')
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
                B.addline(f'    pdb {b.pdb}')
            elif b.state=='MISSING':
                for r in self.residues[b.bounds[0]:b.bounds[1]+1]:
                    rname=self.resname_charmify.get(r.name,r.name)
                    B.addline(f'    residue {r.resseqnum}{r.insertion} {rname} {seglabel}')
                if (b.bounds[1]-b.bounds[0])>(sac_n-1):
                    lrr=self.residues[b.bounds[1]]
                    sac_resseqnum=lrr.resseqnum
                    sac_insertion='A' if lrr.insertion in [' ',''] else chr(ord(lrr.insertion)+1)
                    assert sac_insertion<='Z',f'Residue {lrr.resseqnum} of chain {the_chainID} already has too many insertion instances (last: {lrr.insertion}) to permit insertion of a sacrificial {sac_r}'
                    b.sacres=Residue({'name':sac_r,'resseqnum':sac_resseqnum,'insertion':sac_insertion,'chainID':seglabel,'segtype':'PROTEIN'})
                    B.addline(f'    residue {sac_resseqnum}{sac_insertion} {sac_r} {seglabel}')
        if self.config_params['Fix_engineered_mutations']:
            for m in Mutations:
                if m.source=='SEQADV':
                    if self.config_params['Fix_engineered_mutations']:
                        B.comment('Below reverts an engineered mutation:')
                        B.addline(m.write_TcL())
                else:
                    B.comment('Below is a user-specified mutation:')
                    B.addline(m.write_TcL())
        if self.config_params['Fix_conflicts']:
            for m in Conflicts:
                B.comment('Below reverts a database conflict:')
                B.addline(m.write_TcL())
        B.addline('}')
        B.banner('Coordinate-specification commands')
        for b in self.subsegments:
            if b.state=='RESOLVED':
                B.comment(f'Subsegment {[self.subsegments.index(b)]} is a resolved run')
                B.addline(f'coordpdb {b.pdb} {seglabel}')
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
                        B.comment(f'Subsegment {[self.subsegments.index(b)]}/{len(self.subsegments)} is a missing loop')
                        this_run=ResidueList(self.residues[b.bounds[0]:b.bounds[1]+1])
                        prior_b=self.subsegments[self.subsegments.index(b)-1]
                        B.comment(f'...attached to subsegment {self.subsegments.index(prior_b)}')
                        prior_run=ResidueList(self.residues[prior_b.bounds[0]:prior_b.bounds[1]+1])
                        B.comment(f'Seeding orientation of model-built loop starting at {str(this_run[0])} from {str(prior_run[-1])}')
                        B.addline(f'{this_run.caco_str(prior_run,seglabel,parent_molecule.molid_varname)})')
        B.banner('Intra-segmental terminal patches')
        for i,b in enumerate(self.subsegments):
            if b.state=='MISSING' and i>0 and hasattr(b,'sacres'):
                Cterm=self.residues[b.bounds[1]]
                B.addline(f'patch CTER {seglabel}:{Cterm.resseqnum}{Cterm.insertion}')
                nextb=self.subsegments[i+1]
                Nterm=self.residues[nextb.bounds[0]]
                patchname='NTER'
                if Nterm.name=='PRO':
                    patchname='PROP'
                elif Nterm.name=='GLY':
                    patchname='GLYP'
                B.addline(f'patch {patchname} {seglabel}:{Nterm.resseqnum}{Nterm.insertion}')
                B.addline(f'delatom {seglabel} {b.sacres.resseqnum}{b.sacres.insertion}')
        B.banner('Restoring A.U. state for all resolved subsegments')
        for b in self.subsegments:
            if b.state=='RESOLVED':
                if not isidentity(tmat) or rep_chainID!=the_chainID:
                    B.addline(sel.restore(f'{b.selname}'))
        B.banner(f'END SEGMENT {seglabel}')

class SegmentList(AncestorAwareModList):

    # counters_by_segtype={}
    def append(self,item:Segment):
        olc=self.config['Segname_chars'][item.segtype]
        itemkey=f'{item.chainID}{olc}'
        if not itemkey in self.counters_by_segtype:
            self.counters_by_segtype[itemkey]=1
        self.counters_by_segtype[itemkey]+=1
        if item.segtype=='PROTEIN':
            item.psfgen_segname=itemkey
        else:
            item.psfgen_segname=f'{itemkey}{self.counters_by_segtype[itemkey]:02d}'
        item.resname_charmify=self.config['PDB_to_CHARMM_Resnames']
        super().append(item)

    def __init__(self,config,input_obj):
        self.config={k:config.get(k,'no-config') for k in Segment.config_attr}
        self.counters_by_segtype={}
        if type(input_obj)==ResidueList:
            super().__init__([])
            residues=input_obj
            for stype in ['PROTEIN','WATER','ION','GLYCAN','LIGAND','OTHER']:
                res=residues.filter(segtype=stype)
                for chainID in res.unique_chainIDs():
                    c_res=residues.get(segtype=stype,chainID=chainID)
                    thisSeg=Segment(self.config,c_res)
                    self.append(thisSeg)            
        else:
            logger.error(f'Cannot initialize {self.__class__} from object type {type(input_obj)}')

    def collect_working_files(self):
        working_files=[]
        for seg in self:
            for subseg in seg.subsegments:
                if hasattr(subseg,'pdb'):
                    logger.debug(f'wf: {subseg.pdb}')
                    working_files.append(subseg.pdb)
        return working_files
    
    def write_TcL(self,B:ByteCollector,transform,mods,file_collector=None):
        for seg in self:
            B.comment(f'Writing seg {seg.segname}')
            seg.write_TcL(B,transform,mods,file_collector=file_collector)
            B.comment(f'Done with seg {seg.segname}')
    
    def get_segment_of_residue(self,residue):
        for S in self:
            if residue in S.residues:
                return S
        return None
