import logging
logger=logging.getLogger(__name__)
from .basemod import CloneableMod, CloneableModList
from .config import ConfigGetParam
#class SubsegmentBounds:
#    def __init__(self,l=-1,r=-1,typ='NONE',d=''):
#        self.l=l
#        self.r=r
#        self.typ=typ
#        self.d=d

# class Run:
#     def __init__(self,r,replica_chainID,previous=None,next=None):
#         self.chainID=r.chainID
#         self.replica_chainID=replica_chainID
#         self.segname=r.segname
#         self.residues=[r]
#         self.term='None' # 'N' if r[0] is an N-terminus, 'C' if r[-1] is a C-terminus
#         self.previous=previous
#         self.next=next
#         self.typ='FRAGMENT' if len(r.atoms)>0 else 'LOOP'
#     def add_residue(self,r):
#         if r.chainID==self.chainID and r.segname==self.segname:
#             self.residues.append(r)
#         else:
#             print(f'Error: cannot add residue with chainID {r.chainID} and segname {r.segname} to Run initialized with {self.chainID} and {self.segname}')
#     def pdb_str(self):
#         r0=self.residues[0]
#         r1=self.residues[-1]
#         ins0='' if r0.insertion ==' ' else r0.insertion
#         ins1='' if r1.insertion ==' ' else r1.insertion
#         return '{}_{}{}_to_{}{}.pdb'.format(self.replica_chainID,r0.resseqnum,ins0,r1.resseqnum,ins1)
#     def __str__(self):
#         r0=self.residues[0]
#         r1=self.residues[-1]
#         ins0='' if r0.insertion ==' ' else r0.insertion
#         ins1='' if r1.insertion ==' ' else r1.insertion
#         return '{} {} {}{} to {}{}'.format(self.typ,self.replica_chainID,r0.resseqnum,ins0,r1.resseqnum,ins1)
#     def caco_str(self):
#         return 'coord {} {}{} N [cacoIn_nOut {}{} {} 0]'.format(self.replica_chainID,self.residues[0].resseqnum,self.residues[0].insertion,
#         self.previous.residues[-1].resseqnum,self.previous.residues[-1].insertion,self.replica_chainID)
#     def heal_str(self):
#         rll=self.residues[-2]
#         rl=self.residues[-1]
#         rr=self.next.residues[0]
#         rrr=self.next.residues[1]
#         return 'patch HEAL {c}:{ll} {c}:{l} {c}:{r} {c}:{rr}\n'.format(c=self.replica_chainID,
#                             ll=rll.ri(),l=rl.ri(),r=rr.ri(),rr=rrr.ri())
#     def input_str(self):
#         rl=self.residues[-1]
#         rr=self.next.residues[0]
#         return '{} {} {}\n'.format(self.replica_chainID,rl.ri(),rr.ri())

# class Fragment:
#     ''' a set of contiguous residues with no gaps that can be loaded into a psfgen
#         segment stanza by invoking a pdb file 
#     '''
#     def __init__(self,source_chainID,replica_chainID,resseqnum1,insertion1,resseqnum2,insertion2):
#         self.source_chainID=source_chainID
#         self.replica_chainID=replica_chainID
#         self.resseqnum1=resseqnum1
#         self.resseqnum2=resseqnum2
#         self.insertion1=insertion1
#         self.insertion2=insertion2
#     def pdb_str(self):
#         ins1='' if self.insertion1==' ' else self.insertion1
#         ins2='' if self.insertion2==' ' else self.insertion2
#         return '{}_{}{}_to_{}{}.pdb'.format(self.replica_chainID,self.resseqnum1,ins1,self.resseqnum2,ins2)
#     def __str__(self):
#         ins1='' if self.insertion1==' ' else self.insertion1
#         ins2='' if self.insertion2==' ' else self.insertion2
#         return 'FRAGMENT: {} {}{} to {}{}'.format(self.replica_chainID,self.resseqnum1,ins1,self.resseqnum2,ins2)

#class Loop:
#    ''' a set of contiguous residues from REMARK 465 pdb entries; i.e., they
#        are missing from the base pdb file but present in the construct.  They
#        must be included in a psfgen segment stanza via the 'residue' invocation.
#    '''
#    def __init__(self,source_chainID,replica_chainID,resseqnum0,insertion0,r):
#        self.source_chainID=source_chainID
#        self.replica_chainID=replica_chainID
#        self.resseqnum0=resseqnum0
#        self.insertion0=insertion0
#        self.residues=[r]
#        self.term='UNSET'
#    def add_residue(self,r):
#        self.residues.append(r)
#    def __str__(self):
#        return 'LOOP: {} ({}{})-[{}{} to {}{}]'.format(self.replica_chainID,self.resseqnum0,self.insertion0,self.residues[0].resseqnum,self.residues[0].insertion,self.residues[-1].resseqnum,self.residues[-1].insertion)
#    def caco_str(self):
#        return 'coord {} {} N [cacoIn_nOut {} {} 0]\n'.format(self.replica_chainID,self.residues[0].resseqnum,self.resseqnum0,self.replica_chainID)

class Segment(CloneableMod):
    req_attr=['segtype','segname','residues','parent_molecule','subsegments']
    opt_attr=['mutations','deletions','grafts','attachments']

    @classmethod
    def from_residue_list(cls,Residues,parent_molecule):
        if len(Residues)==0:
            return None
        apparent_chainID=Residues[0].chainID
        apparent_segtype=Residues[0].segtype
        myRes=Residues.clone()
        if apparent_segtype=='PROTEIN':
            assert myRes.puniq(['resseqnum','insertion'])
        else:
            myRes.uniquify(['resseqnum','insertion'])
        myRes.sort()
        subsegments=myRes.binary_bounds(lambda x: 'RESOLVED' if len(x.atoms)>0 else 'MISSING')
        olc=ConfigGetParam('Segname_chars').get(apparent_segtype)
        input_dict={
            'segtype': apparent_segtype,
            'segname': f'{apparent_chainID}{olc}',
            'residues': Residues,
            'parent_molecule': parent_molecule,
            'subsegments':subsegments
        }
        inst=cls(input_dict)
        return inst
    
    def psfgen_atomselect(self):
        serial_list=self.residues.atom_serials()
        restr=f'set {self.segname}_sel [atomselect {self.parent_molecule.molid} serial '+' '.join(serial_list)+' ]'
        return restr

class SegmentList(CloneableModList):
    pass
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

                    if not tmat.isidentity():
                         stanza.extend(sel.backup('mysel'))
                         stanza.append('$mysel move {}'.format(tmat.OneLiner()))
                    ''' tmat transformation '''
                    stanza.append('$mysel writepdb {}'.format(ss.pdb_str()))
                    if not tmat.isidentity():
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

