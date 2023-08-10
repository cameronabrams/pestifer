import logging
logger=logging.getLogger(__name__)
from .basemod import AncestorAwareMod, AncestorAwareModList
from .config import ConfigGetParam
from .residue import ResidueList
import pestifer.sel as sel
from .util import isidentity, reduce_intlist


class Segment(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['segtype','segname','chainID','residues','subsegments']
    opt_attr=AncestorAwareMod.opt_attr+['mutations','deletions','grafts','attachments','psfgen_segname']

    @classmethod
    def from_residue_list(cls,Residues):
        if len(Residues)==0:
            return None
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
        inst=cls(input_dict)
        return inst
    
    def __str__(self):
        return f'{self.segname}: type {self.segtype} chain {self.chainID} with {len(self.residues)} residues'

    def psfgen_stanza(self,tmat,replica_chains={}):
        if self.segtype=='PROTEIN':
            return self.protein_stanza(tmat,replica_chains)
        elif self.segtype=='GLYCAN':
            return self.glycan_stanza(tmat,replica_chains)
        else:
            return self.generic_stanza(tmat,replica_chains)
        
    def protein_stanza(self,tmat,replica_chains={}):
        parent_molecule=self.ancestor_obj
        the_chainID=replica_chains.get(self.residues[0].chainID,self.residues[0].chainID)
        sac_rd=ConfigGetParam('Sacrificial_residue',default={})
        include_terminal_loops=ConfigGetParam('Include_terminal_loops',default=False)
        fix_engineered_muations=ConfigGetParam('Fix_engineered_mutations',default=False)
        fix_conflicts=ConfigGetParam('Fix_conflicts',False)
        sac_r=sac_rd['name']
        sac_n=sac_rd['min_loop_length']
        stanza=f'### BEGIN SEGMENT {self.segname} ###\n'
        for i,b in enumerate(self.subsegments):
            if b.state=='RESOLVED':
                selname=f'{self.psfgen_segname}'
                run=ResidueList(self.residues[b.bounds[0]:b.bounds[1]+1])
                b.pdb=f'PROTEIN_{the_chainID}_{run[0].resseqnum}{run[0].insertion}_to_{run[-1].resseqnum}{run[-1].insertion}.pdb'
                # TODO: account for any deletions
                serial_list=run.atom_serials(as_type=int)
                vmd_red_list=reduce_intlist(serial_list)
                stanza+=f'set {selname} [atomselect {parent_molecule.molid} "serial {vmd_red_list}"]\n'
                # stanza+=f'set {selname} [atomselect {parent_molecule.molid} "serial '+' '.join(serial_list)+'"]\n'
                if replica_chains:
                    chainIDs=[the_chainID for x in serial_list]
                    stanza+=f'${selname} set chain [list {" ".join(chainIDs)}]\n'                        
                if not isidentity(tmat):
                    stanza+=sel.backup(f'{selname}')
                    stanza+=f'${selname} move {tmat.OneLiner()}\n'
                    stanza+=sel.restore(f'{selname}')
                stanza+=f'${selname} writepdb {b.pdb}\n'
        stanza+=f'segment {self.segname} '+'{\n'
        print(len(self.subsegments),type(self.subsegments))
        if not include_terminal_loops:
            if self.subsegments[0].state=='MISSING':
                Nterminal_missing_subsegment=self.subsegments.pop(0)
                logger.info(f'Since terminal loops are not included, ignoring {str(Nterminal_missing_subsegment)}')
            if self.subsegments[-1].state=='MISSING':
                Cterminal_missing_subsegment=self.subsegments.pop(-1)
                logger.info(f'Since terminal loops are not included, ignoring {str(Cterminal_missing_subsegment)}')
        for b in self.subsegments:
            if b.state=='RESOLVED':
                stanza+=f'    pdb {b.pdb}\n'
            elif b.state=='MISSING':
                for r in self.residues[b.bounds[0]:b.bounds[1]+1]:
                    stanza+=f'    residue {r.resseqnum}{r.insertion} {r.resname_charmify()} {the_chainID}\n'
                if (b.bounds[1]-b.bounds[0])>(sac_n-1):
                    lrr=self.residues[b.bounds[1]]
                    sac_resseqnum=lrr.resseqnum
                    sac_insertion='A' if lrr.insertion in [' ',''] else chr(ord(lrr.insertion)+1)
                    assert sac_insertion<='Z',f'Residue {lrr.resseqnum} of chain {the_chainID} already has too many insertion instances (last: {lrr.insertion}) to permit insertion of a sacrificial {sac_r}'
                    stanza+=f'    residue {sac_resseqnum}{sac_insertion} {sac_r} {the_chainID}\n'
        # TODO: mutations
        stanza+='}\n'
        for b in self.subsegments:
            if b.state=='RESOLVED':
                stanza+=f'coordpdb {b.pdb} {self.segname}\n'
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
                        stanza+=f'{this_run.caco_str(prior_run)}\n'
        stanza+=f'### END SEGMENT {self.psfgen_segname} ###\n'
        return stanza

    # def psfgen_atomselect(self):
    #     serial_list=self.residues.atom_serials()
    #     restr=f'set {self.segname}_sel [atomselect {parent_molecule.molid} serial '+' '.join(serial_list)+' ]'
    #     return restr

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

    @classmethod
    def from_residues(cls,residues:ResidueList):
        tmp_counters={}
        Slist=[]
        protein_res=residues.get(segtype='PROTEIN')
        protein_chains=[]
        for r in protein_res:
            if not r.chainID in protein_chains:
                protein_chains.append(r.chainID)
        
        for pc in protein_chains:
            pc_res=residues.get(segtype='PROTEIN',chainID=pc)
            thisSeg=Segment.from_residue_list(pc_res)
            olc=ConfigGetParam('Segname_chars').get(thisSeg.segtype)
            itemkey=f'{thisSeg.chainID}{olc}'
            if not itemkey in tmp_counters:
                tmp_counters[itemkey]=0
            tmp_counters[itemkey]+=1
            if thisSeg.segtype=='PROTEIN':
                thisSeg.psfgen_segname=thisSeg.segname
            else:
                thisSeg.psfgen_segname=f'{itemkey}{tmp_counters[itemkey]:03d}'
            Slist.append(thisSeg)

        # TODO: Water, ligand, glycan, etc

        return cls(Slist)

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
                         stanza.append('$mysel move {}'.format(tmat.OneLiner()))
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
                     stanza.append('$mysel move {}'.format(tmat.OneLiner()))
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
                 stanzastr+='$mysel move {}\n'.format(tmat.OneLiner())
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

