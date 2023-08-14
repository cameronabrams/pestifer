import logging
logger=logging.getLogger(__name__)
from .basemod import AncestorAwareMod, AncestorAwareModList
from .mods import MutationList
from .config import ConfigGetParam
from .residue import Residue,ResidueList
import pestifer.sel as sel
from .util import isidentity, reduce_intlist
from .stringthings import ByteCollector

class Segment(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['segtype','segname','chainID','residues','subsegments']
    opt_attr=AncestorAwareMod.opt_attr+['mutations','deletions','grafts','attachments','psfgen_segname','config_params']

    def __init__(self,input_obj):
        parm_dict={}
        for p in ['Include_terminal_loops','Fix_engineered_mutations','Fix_conflicts','Sacrificial_residue']:
            parm_dict[p]=ConfigGetParam(p)
        if type(input_obj)==dict:
            input_dict=input_obj
        elif type(input_obj)==Residue:
            res=input_obj
            apparent_chainID=res.chainID
            apparent_segtype=res.segtype
            myRes=ResidueList([])
            myRes.append(res.clone())
            subsegments=myRes.state_bounds(lambda x: 'RESOLVED' if len(x.atoms)>0 else 'MISSING')
            olc=ConfigGetParam('Segname_chars').get(apparent_segtype)
            input_dict={
                'segtype': apparent_segtype,
                'segname': f'{apparent_chainID}{olc}',
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
            olc=ConfigGetParam('Segname_chars').get(apparent_segtype)
            input_dict={
                'segtype': apparent_segtype,
                'segname': f'{apparent_chainID}{olc}',
                'chainID': apparent_chainID,
                'residues': myRes,
                'subsegments':subsegments
            }
        else:
            logger.error(f'Cannot initialize {self.__class__} from object type {type(input_obj)}')
        input_dict['config_params']=parm_dict
        super().__init__(input_dict)

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
        B.banner(f'Not generating segment {self.segname} ({self.segtype})')
        B.banner('GLYCANS ARE NOT YET IMPLEMENTED')
    
    def generic_stanza(self,B,transform,mods,file_collector=None):
        B.banner(f'Not generating generic segment {self.segname} ({self.segtype})')

    def protein_stanza(self,B:ByteCollector,transform,mods,file_collector=None):
        parent_molecule=self.ancestor_obj
        the_chainID=self.residues[0].chainID
        chainIDmap=transform.chainIDmap
        tmat=transform.tmat
        rep_chainID=chainIDmap.get(the_chainID,the_chainID)
        seglabel=self.psfgen_segname # protein, so should just be letter
        if rep_chainID!=the_chainID:
            logger.info(f'relabeling chain {the_chainID} to {rep_chainID} in {seglabel}')
            seglabel=f'{rep_chainID}' # first byte is chain ID
        transform.segnamemap[self.psfgen_segname]=seglabel
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
                logger.info(f'Since terminal loops are not included, ignoring {str(Nterminal_missing_subsegment)}')
            if self.subsegments[-1].state=='MISSING':
                Cterminal_missing_subsegment=self.subsegments.pop(-1)
                logger.info(f'Since terminal loops are not included, ignoring {str(Cterminal_missing_subsegment)}')
        for b in self.subsegments:
            if b.state=='RESOLVED':
                B.addline(f'    pdb {b.pdb}')
            elif b.state=='MISSING':
                for r in self.residues[b.bounds[0]:b.bounds[1]+1]:
                    B.addline(f'    residue {r.resseqnum}{r.insertion} {r.resname_charmify()} {seglabel}')
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
                    B.addline(m.write_TL())
        if self.config_params['Fix_conflicts']:
            for m in Conflicts:
                B.comment('Below reverts a database conflict:')
                B.addline(m.write_TcL())
        B.addline('}')
        B.banner('Coordinate-specification commands')
        for b in self.subsegments:
            if b.state=='RESOLVED':
                B.addline(f'coordpdb {b.pdb} {seglabel}')
        for b in self.subsegments:
            if b.state=='MISSING':
                # will issue the atom-reorienting command to join the C-terminus of prior run to N-terminus of this one, which is model-built using guesscoord
                if (self.subsegments.index(b)==0 or self.subsegments.index(b)==(len(self.subsegments)-1)) and not ConfigGetParam('Include_terminal_loops'):
                    # this is a terminal loop and we are not including terminal loops
                    pass
                else:
                    # not terminal OR we ARE including terminal loops
                    if self.subsegments.index(b)>0:
                        # this is either interior OR C-terminal to be included
                        this_run=ResidueList(self.residues[b.bounds[0]:b.bounds[1]+1])
                        prior_b=self.subsegments[self.subsegments.index(b)-1]
                        prior_run=ResidueList(self.residues[prior_b.bounds[0]:prior_b.bounds[1]+1])
                        B.banner(f'Seeding orientation of model-built loop starting at {str(this_run[0])}')
                        B.addline(f'{this_run.caco_str(prior_run,seglabel)}')
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
    counters_by_segtype={}
    def append(self,item:Segment):
        olc=ConfigGetParam('Segname_chars').get(item.segtype)
        itemkey=f'{item.chainID}{olc}'
        if not itemkey in self.counters_by_segtype:
            self.counters_by_segtype[itemkey]=0
        self.counters_by_segtype[itemkey]+=1
        if item.segtype=='PROTEIN':
            item.psfgen_segname=item.segname
        else:
            item.psfgen_segname=f'{itemkey}{self.counters_by_segtype[itemkey]:03d}'
        super().append(item)

    def __init__(self,input_obj):
        if type(input_obj)==ResidueList:
            super().__init__([])
            residues=input_obj
            tmp_counters={}
            for stype in ['PROTEIN','WATER','ION','GLYCAN','LIGAND','OTHER']:
                res=residues.get(segtype=stype)
                for chainID in res.unique_chainIDs():
                    c_res=residues.get(segtype=stype,chainID=chainID)
                    thisSeg=Segment(c_res)
                    olc=ConfigGetParam('Segname_chars').get(thisSeg.segtype)
                    if thisSeg.segtype=='PROTEIN':
                        thisSeg.psfgen_segname=thisSeg.segname
                    else:
                        thisSeg.psfgen_segname=f'{itemkey}{tmp_counters[itemkey]:03d}'
                    itemkey=f'{thisSeg.chainID}{olc}'
                    if not itemkey in tmp_counters:
                        tmp_counters[itemkey]=0
                    tmp_counters[itemkey]+=1
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
    
    # def __init__(self,r,parent_chain=None,subcounter=''):
    #     """Initializes a segment instance by passing in first residue of segment"""
    #     self.parent_chain=parent_chain
    #     self.segtype=_seg_typedict_byresname_[r.name]
    #     if self.segtype=='PROTEIN':
    #         self.segname=r.chainID
    #     else:
    #         self.segname=r.chainID+_segname_second_character_[_seg_typedict_byresname_[r.name]]+subcounter
    #     self.residues=[r]
    #     self.mutations=[]
    #     self.deletions=[]
    #     self.graft=''
    #     self.rootres=''
    #     self.attach=''
    #     self.pdbfiles=[]
    #     self.Runs=[]
    #     if self.segtype=='GLYCAN':
    #         if len(r.up)==0:
    #             print('ERROR: {}:{}{} has no uplink in PDB file.  You may need to add one to a user-link input.'.format(self.segname,r.name,r.resseqnum))
    #             exit()
    #         self.rootres=r.up[0]
    #     r.segname=self.segname
    #     for a in r.atoms:
    #         a.segname=self.segname

class OldSegment:
    def get_chainID(self):
        return self.parent_chain.source_chainID
    def get_molecule(self):
#        print(self.parent_chain.parent_molecule)
        return self.parent_chain.parent_molecule
    def __str__(self):
        retbase='({}){}[{}] {} - {}'.format(self.parent_chain.chainID,self.segname,self.segtype,self.residues[0].ri(),self.residues[-1].ri())
        if self.rootres!='':
            r=self.rootres
            retbase='('+str(r)+')->'+retbase
        if self.graft!='':
            retbase+='::'+self.graft.source_pdb+':'+str(self.graft.source_segment)
        return retbase
    def add_residue(self,r):
        self.residues.append(r)
        r.segname=self.segname
        for a in r.atoms:
            a.segname=self.segname
    def apply_graft(self,g):
        self.graft=g

    def doRuns(self,replica_chainID):
        self.Runs=[]
        Loops=[]
        for r in self.residues:
            if len(self.Runs)==0:
                self.Runs.append(Run(r,replica_chainID))
                ri=0
                self.Runs[ri].term='N'
                if len(r.atoms)==0:
                    self.Runs[ri].typ=='LOOP'
                else:
                    self.Runs[ri].typ=='FRAGMENT'
            else:
                if self.Runs[ri].typ=='FRAGMENT':
                    if len(r.atoms)>0:
                        self.Runs[ri].add_residue(r)
                    else:
                        self.Runs.append(Run(r,replica_chainID))
                        self.Runs[ri].next=self.Runs[ri+1]
                        self.Runs[ri+1].previous=self.Runs[ri]
                        ri+=1
                elif self.Runs[ri].typ=='LOOP':
                    if len(r.atoms)>0:
                        Loops.append(self.Runs[ri])
                        self.Runs.append(Run(r,replica_chainID))
                        self.Runs[ri].next=self.Runs[ri+1]
                        self.Runs[ri+1].previous=self.Runs[ri]
                        ri+=1
                    else:
                        self.Runs[ri].add_residue(r)
        self.Runs[ri].term='C'
        return Loops

    def psfgen_stanza(self,includeTerminalLoops=False,tmat=None):
        stanza=[]
        my_chainID=self.get_chainID()
        rep_chainID=tmat.get_replica_chainID(my_chainID)
        rep_segname=self.segname
        if tmat==None:
            print('ERROR: psfgen_stanza needs a tmat!')
            exit()
        if self.segtype=='PROTEIN':
            if self.graft!='' or self.attach!='':
                print('ERROR: Protein grafts and attachment segments are not yet implemented')
                return 'ERROR',[]
#            Loops=self.subsegments(my_chainID,rep_chainID)
            Loops=self.doRuns(rep_chainID)
            ''' PART 1:  Process selections and write PDB files for fragments '''
            for ss in self.Runs:
                if ss.typ=='FRAGMENT':
                    r0=ss.residues[0]
                    r1=ss.residues[-1]
                    delstr=''
                    for d in self.deletions:
                        if r0.ri() <= d.resseqnumi <= r1.ri():
                            delstr+=' and not resid {}'.format(d.resseqnumi)
                    stanza.append('set mysel [atomselect ${} "protein and chain {} and resid {} to {} {}"]\n'.format(self.get_molecule().molid_varname,
                                self.parent_chain.source_chainID,r0.ri(),r1.ri(),delstr))

                    if not isidentity(tmat):
                         stanza.extend(sel.backup('mysel'))
                         stanza.append('$mysel move {}'.format(tmat.write_TcL()))
                    ''' tmat transformation '''
                    stanza.append('$mysel writepdb {}'.format(ss.pdb_str()))
                    if not isidentity(tmat):
                         stanza.extend(sel.restore('mysel'))
                    self.pdbfiles.append(ss.pdb_str())
            ''' PART 2:  Build segment stanza '''
            stanza.append('segment {} {{\n'.format(rep_segname))
            for i,ss in enumerate(self.Runs):
                if ss.typ=='FRAGMENT':
                    stanza.append('   pdb {}\n'.format(ss.pdb_str()))
                elif ss.typ=='LOOP':
                    if (i==0 or i==(len(self.Runs)-1)) and not includeTerminalLoops:
                        ''' shunt this if this is a terminal loop and includeTerminalLoops is False '''
                        pass
                    else:
                        ''' this is either NOT a terminal loop, or if it is, includeTerminalLoops is True '''
                        for rr in ss.residues:
                            take_it=True
                            for d in self.deletions:
                                if d.resseqnumi == rr.ri():
                                    take_it=False
                                    break
                            if take_it:
                                nm=ResnameCharmify(rr.name)
                                stanza.append('   residue {}{} {} {}\n'.format(rr.resseqnum,rr.insertion,nm,tmat.get_replica_chainID(rr.chainID)))
                        ss.sacrins='0'
                        if len(ss.residues)>3:
                            ''' modeled-in loops longer than three residues need extra processing to relax long bonds. 
                                this starts with inserting a sacrificial glycine so that a pair of CTER/NTER patches
                                can be used to slice the bonding continuity of the segment '''
                            rr=ss.residues[-1]
                            ss.sacrins='A' if rr.insertion == '' or rr.insertion == ' ' else chr(ord(rr.insertion)+1)
                            stanza.append('   residue {}{} {} {}'.format(rr.resseqnum,ss.sacrins,'GLY',tmat.get_replica_chainID(rr.chainID)))
            ''' PART 2.1:  Include mutations '''
            #print('### {} mutations'.format(len(self.mutations)))
            for m in self.mutations:
                stanza.append(m.psfgen_segment_str())
            stanza.append('}') # end of segment

            ''' PART 3:  Issue coordinate-setting commands '''
            ''' coordpdb calls '''
            for ss in self.Runs:
                if ss.typ=='FRAGMENT':
                    stanza.append('coordpdb {} {}'.format(ss.pdb_str(),rep_segname))
            ''' caco calls '''
            for i in range(0,len(self.Runs)):
                ss=self.Runs[i]
                if ss.typ=='LOOP':
                    if (i==0 or i==(len(self.Runs)-1)) and not includeTerminalLoops:
                        pass
                    else:
                        stanza.append(ss.caco_str())
            ''' slice calls '''
            for i in range(0,len(self.Runs)):
                ss=self.Runs[i]
                if ss.typ=='LOOP':
                    if (i==0 or i==(len(self.Runs)-1)) and not includeTerminalLoops:
                        pass
                    else:
                        if (ss.sacrins!='0' and i>0 and i<(len(self.Runs)-1)):
                            nter='NTER'
                            nextres=ss.next.residues[0]
                            rname=nextres.name
                            if rname=='PRO':
                                nter='PROP'
                            if rname=='GLY':
                                nter='GLYP'
                            stanza.append('## Slicing segment at sacrificial glycine to make CTER-NTER')
                            stanza.append('patch CTER {}:{}{}'.format(rep_segname,ss.residues[-1].resseqnum,ss.residues[-1].insertion))
                            stanza.append('## making residue {}{} an nter with patch {}\n'.format(rname,nextres.ri(),nter))
                            stanza.append('patch {} {}:{}\n'.format(nter,rep_segname,nextres.ri()))
                            stanza.append('delatom {} {}{}\n'.format(rep_segname,ss.residues[-1].resseqnum,ss.sacrins))
            return stanza,Loops
        elif self.segtype=='GLYCAN':
            # TODO: stanzstr->stanza[]
            stanza=[]
            if self.graft!='':
                ''' this is a segment with an associated graft '''
                g=self.graft
                g.ingraft_segname=self.segname # graft inherits segname
                g.ingraft_chainID=self.parent_chain.chainID
                stanza.extend(g.transform(self.parent_chain.parent_molecule))
                self.pdbfiles.append(g.transformed_pdb)
                ''' g.transformed_pdb is now the file with inputs! '''
                stanza.append(r'segment {} {{'.format(rep_segname))
                stanza.append('     pdb {}'.format(g.transformed_pdb))
                stanza.append('}')
                stanza.append('coordpdb {} {}'.format(g.transformed_pdb,rep_segname))
            elif self.attach!='':
                ''' this is a segment from another molecule to be attached '''
                logger.info('Attachments are not yet implemented.')
                pass
            else:
                ''' this is an existing glycan segment that will not be overwritten with a graft '''
                f=Fragment(my_chainID,rep_chainID,self.residues[0].resseqnum,self.residues[0].insertion,self.residues[-1].resseqnum,self.residues[-1].insertion)
                pdb=f.pdb_str()
                self.pdbfiles.append(pdb)
                #print('#### segment {}'.format(self.segname))
                stanza.append('set mysel [atomselect ${} "chain {} and resid {} to {}"]'.format(self.get_molecule().molid_varname,self.parent_chain.source_chainID,self.residues[0].resseqnum,self.residues[-1].resseqnum))
                if not tmat.isidentity():
                     stanza.extend(sel.backup('mysel'))
                     stanza.append('$mysel move {}'.format(tmat.write_TcL()))
                stanza.extend(sel.charmm_namify('mysel'))
                stanza.append('$mysel writepdb {}'.format(pdb))
                if not tmat.isidentity():
                     stanza.append(sel.restore('mysel'))
                stanza.append('segment {} {{'.format(rep_segname))
                stanza.append('    pdb {}'.format(pdb))
                stanza.append('}')
                stanza.append('coordpdb {} {}'.format(pdb,rep_segname))
            return stanza,[]
#        elif self.segtype == 'LIGAND':
#            # working here
#            pass
        elif self.segtype in ['WATER','ION','LIGAND','OTHER']:
            f=Fragment(my_chainID,tmat.get_replica_chainID(my_chainID),self.residues[0].resseqnum,self.residues[0].insertion,self.residues[-1].resseqnum,self.residues[-1].insertion)
            pdb=f.pdb_str()
            self.pdbfiles.append(pdb)
            stanzastr+='set mysel [atomselect ${} "chain {} and resid {}{} to {}{}"]\n'.format(self.get_molecule().molid_varname,my_chainID,self.residues[0].resseqnum,self.residues[0].insertion,self.residues[-1].resseqnum,self.residues[-1].insertion)
            stanzastr+=sel.charmm_namify('mysel')
            if not tmat.isidentity():
                 stanzastr+=sel.backup('mysel')
                 stanzastr+='$mysel move {}\n'.format(tmat.write_TcL())
            stanzastr+='$mysel writepdb {}\n'.format(pdb)
            if not tmat.isidentity():
                 stanzastr+=sel.restore('mysel')
            stanzastr+='segment {} {{\n'.format(rep_segname)
            stanzastr+='    pdb {}\n'.format(pdb)
            stanzastr+='}\n'
            stanzastr+='coordpdb {} {}\n'.format(pdb,rep_segname)
            return stanzastr,[] 
        else:
            print('ERROR: Unrecognized segment type {}'.format(self.segtype))
            return 'ERROR',[]

