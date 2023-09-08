# Author: Cameron F. Abrams <cfa22@drexel.edu
""" Modifications to sequence, topology, and coordinates

    This module defines several types of modifications (or "mods")
    and their respective list collections.  

    Mods are classified into four types:
    1. 'seqmod's are modifications to the amino-acid sequence of a 
    protein.  These mods are typically applied while pestifer
    is building its internal representation of the molecule.
    2. 'topomod's are modifications to the bonding topology, 
    including disulfides, covalent bonds linking sugars
    or other non-protein residues, etc.
    3. 'coormod's are mods to the coordinates of a molecule,
    such as rotation of specific dihedral angles
    4. 'generic' -- anything not in the other three categories.

    seqmods
    -------
    * Ter -- This is a representation of a "TER" record from an
    old-style PDB file. These are tracked in order to modify
    atom serial numbers, since a TER acquires a unique atom
    serial number even though it is not an atom.
    * Seqadv -- A reprentation of any single-point discrepancies
    between the sequence in the actual structure and a database
    sequence. Engineered mutations, conflicts, etc.
    * Mutation -- A point mutation of a particular residue 
    in a particular chain to a new residue name
    * Deletion -- deletion of one or more contiguous residues
    * Substitution -- substitution of one or more contiguous
    residues with a different block of residues of possibly
    different length

    topomods
    --------
    * SSBond -- a disulfide
    * Link -- a non-disulfide covalent bond, typically involving 
    non-protein residues

    coormods
    --------
    * Crot -- any of a variety of rotations around certain dihedrals

"""
import logging
logger=logging.getLogger(__name__)
from pidibble.pdbrecord import PDBRecord
from .cifutil import CIFdict
from .basemod import AncestorAwareMod,AncestorAwareModList
from .stringthings import split_ri
from .scriptwriters import Psfgen
from .config import res_123
from functools import singledispatchmethod

ModTypes=['seqmod','topomod','coormod','generic']

class Ter(AncestorAwareMod):
    """A class for handing TER records in PDB files
    
    Why is this necessary?  TER records are used to indicate
    'breaks' in the list of ATOM records in a PDB file to signify
    physical discontinuity of the amino acid chain(s).  The problem
    we have to solve is that TER records get a unique atom serial
    number, even though they are not atoms. So in a list of atoms
    the serial numbers can be discontinuous.  This is not a problem
    per se, *but* when VMD reads in a PDB file it ignores both
    TER records *and* the explicit atom serial numbers.  Pidibble
    by default does not, so to make the atom records VMD-compliant,
    we have to adjust atom serial numbers, and to do that we
    have to know where the TER records are.
    
    Attributes
    ----------
    req_attr: list
        * serial: (int) atom serial number
        * resname: (str) residue name
        * chainID: (str) chain identifier
        * resseqnum: (int) residue number
        * insertion: (chr) residue insertion code
    
    yaml_header: (str)
        label used for yaml format input/output
    
    modtype: (str)
        type of mod from among ModTypes

    """
    req_attr=['serial','resname','chainID','resseqnum','insertion']
    yaml_header='terminals'
    modtype='seqmod'
    PDB_keyword='TER'

    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)
    @__init__.register(PDBRecord)
    def _from_pdbrecord(self,pdbrecord):
        input_dict={
            'serial':pdbrecord.serial,
            'resname':pdbrecord.residue.resName,
            'chainID':pdbrecord.residue.chainID,
            'resseqnum':pdbrecord.residue.seqNum,
            'insertion':pdbrecord.residue.iCode
        }
        super().__init__(input_dict)
    
class TerList(AncestorAwareModList):
    pass

class Seqadv(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['idCode','resname','chainID','resseqnum','insertion','typekey']
    opt_attr=AncestorAwareMod.opt_attr+['database','dbAccession','dbRes','dbSeq','pdbx_pdb_strand_id','pdbx_auth_seq_num','pdbx_ordinal','residue']
    attr_choices=AncestorAwareMod.attr_choices.copy()
    attr_choices.update({'typekey':['conflict','cloning','expression','typekey','engineered','variant','insertion','deletion','microheterogeneity','chromophore','user','_other_']})
    yaml_header='seqadvs'
    modtype='seqmod'
    PDB_keyword='SEQADV'
    mmCIF_name='struct_ref_seq_dif'
    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)
    @__init__.register(PDBRecord)
    def _from_pdbrecord(self,pdbrecord):
        input_dict={
            'idCode':pdbrecord.idCode,
            'resname':pdbrecord.residue.resName,
            'chainID':pdbrecord.residue.chainID,
            'resseqnum':pdbrecord.residue.seqNum,
            'insertion':pdbrecord.residue.iCode,
            'database':pdbrecord.database,
            'dbAccession':pdbrecord.dbAccession,
            'dbRes':pdbrecord.dbRes,
            'dbSeq':pdbrecord.dbSeq,
            'typekey':self.seqadv_details_keyword(pdbrecord.conflict),
            'residue':None
        }
        super().__init__(input_dict)
    @__init__.register(CIFdict)
    def _from_cifdict(self,cd):
        input_dict={
            'idCode':cd['pdbx_pdb_id_code'],
            'resname':cd['mon_id'],
            'chainID':'UNSET', # author
            'resseqnum':int(cd['seq_num']), # mmcif
            'insertion':cd['pdbx_pdb_ins_code'], # author
            'database':cd['pdbx_seq_db_name'],
            'dbAccession':cd['pdbx_seq_db_accession_code'],
            'dbRes':cd['db_mon_id'],
            'dbSeq':cd['pdbx_seq_db_seq_num'], # author
            'typekey':self.seqadv_details_keyword(cd['details']),
            'pdbx_auth_seq_num':int(cd['pdbx_auth_seq_num']), # author
            'pdbx_ordinal':cd['pdbx_ordinal'],
            'pdbx_pdb_strand_id':cd['pdbx_pdb_strand_id'],
            'residue':None
        }
        super().__init__(input_dict)

    def seqadv_details_keyword(self,text):
        text=text.lower()
        keyword_list=self.__class__.attr_choices['typekey']
        for keyword in keyword_list:
            if keyword in text:
                return keyword
        return '_other_'

    def pdb_line(self):
        return f'SEQADV {self.idCode:3s} {self.resname:>3s} {self.chainID:1s} {self.resseqnum:>4d}{self.insertion:1s} {self.database:>4s} {self.dbAccession:9s} {self.dbRes:3s} {self.dbSeq:>5d} {self.typekey:21s}          '
    
    def assign_residue(self,Residues):
        # logger.debug(f'Searching {len(Residues)} Residues for auth_chain {self.pdbx_pdb_strand_id} auth_seq {self.pdbx_auth_seq_num}')
        assert self.residue==None
        if hasattr(self,'pdbx_pdb_strand_id') and hasattr(self,'pdbx_auth_seq_num'):
            self.assign_obj_to_attr('residue',Residues,
                            auth_asym_id='pdbx_pdb_strand_id',
                            auth_seq_id='pdbx_auth_seq_num',
                            insertion='insertion')
            if self.residue!=None:
                self.chainID=self.residue.chainID
        else:
            self.assign_obj_to_attr('residue',Residues,chainID='chainID',resseqnum='resseqnum',insertion='insertion')

    def update_from_residue(self):
        assert self.residue!=None
        self.chainID=self.residue.chainID
        # else:
        #     logger.debug(f'...seqadv {self.typekey} auth {self.pdbx_pdb_strand_id}:{self.pdbx_auth_seq_num} cannot be resolved from current set of residues')
        # we'll assume that if this residue is not found, then this seqadv is never used anyway

class SeqadvList(AncestorAwareModList):
    def assign_residues(self,Residues):
        delete_us=[]
        for s in self:
            s.assign_residue(Residues)
        delete_us=self.__class__([s for s in self if s.residue==None])
        for s in delete_us:
            self.remove(s)
        return delete_us
    
    def update_from_residues(self):
        for s in self:
            s.update_from_residue()

class Mutation(AncestorAwareMod):
    """A class for handling single-residue mutations
    
    Single-residue mutations are handled in psfgen by including a "mutate" directive within
    a protein segment.  The mutate directive needs only the resid and 3-letter desired
    residue name at that resid position.  Mutate directives come after all pdb and residue
    directives in segment, typically.

    Mutations can be inferred directly from coordinate-file metadata or supplied by the
    user.
    
    """
    req_attr=AncestorAwareMod.req_attr+['chainID','origresname','resseqnum','insertion','newresname','typekey']
    opt_attr=AncestorAwareMod.opt_attr+['pdbx_auth_seq_num']
    yaml_header='mutations'
    modtype='seqmod'

    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)
        if not hasattr(self,'insertion'):
            self.insertion=''
        if not hasattr(self,'typekey'):
            self.typekey='unknown'

    @__init__.register(Seqadv)
    def _from_seqadv(self,sq):
        input_dict={
            'chainID':sq.chainID,
            'origresname':sq.resname,
            'newresname':sq.dbRes,
            'resseqnum':sq.resseqnum,
            'insertion':sq.insertion,
            'typekey':sq.typekey
        }
        if hasattr(sq,'pdbx_auth_seq_num'):
            input_dict['pdbx_auth_seq_num']=sq.__dict__['pdbx_auth_seq_num']
        super().__init__(input_dict)

    @__init__.register(str)
    def _from_shortcode(self,shortcode):
        ### shortcode format: c:nnn,rrr,mmm or c:A###B
        ### c: chainID
        ### nnn: old resname, 1-byte rescode allowed
        ### rrr: resseqnum
        ### mmm: new resname, 1-byte rescode allowed
        s1=shortcode.split(':')
        chainID=s1[0]
        s2=s1[1].split(',')
        codetype='3' if len(s2)==3 else '1'
        if codetype=='3':
            ri=s2[1]
            r,i=split_ri(ri)
            orn=s2[0]
            nrn=s2[2]
        else:
            s2=s2[0]
            orn=res_123[s2[0]]
            nrn=res_123[s2[-1]]
            ri=s2[1:-1]
            r,i=split_ri(ri)
        input_dict={
            'chainID':chainID, # assume author
            'origresname':orn,
            'resseqnum':int(r), # assume author!
            'insertion':i,
            'newresname':nrn,
            'typekey':'user' # shortcodes are how users indicate mutations
        }
        logger.debug(f'user mutation {input_dict}')
        super().__init__(input_dict)

    def __str__(self):
        return f'{self.chainID}:{self.origresname}{self.resseqnum}{self.insertion}{self.newresname}'

    def write_TcL(self):
        if hasattr(self,'pdbx_auth_seq_num'): # mmCIF!
            return f'    mutate {self.resseqnum} {self.newresname}'
        else:
            resseqnumi=f'{self.resseqnum}{self.insertion}'
            return f'    mutate {resseqnumi} {self.newresname}'

class MutationList(AncestorAwareModList):
    pass

class Substitution(AncestorAwareMod):
    """A class for handling substitutions 
    
    A substitution is a user-specified modification in which a sequence of one or
    more contiguous residues are replaced by one or more residues.

    There are four typical cases for substitutions, provided they are not very long:
    1. The substitution is fully contained within a resolved run of residues;
    2. The substitution is fully contained within an unresolved (missing) run of residues;
    3. The substitution starts in a resolved run and ends in an unresolved run;
    4. The substitution starts in an unresolves run and ends in a resolved run.

    Very long substitutions (those replacing a lot of residues) could conceivably 
    span a residue range that contains multiple distinct resolved and unresolved runs,
    but this will be assumed to be very rare.

    We also denote three types of substitution:
    1. The difference between the number of residues in the range where the subsitution
    is made and the length of the substitution is *positive*, meaning residues in the
    original sequence will be deleted;
    2. The difference between the number of residues in the range where the subsitution
    is made and the length of the substitution is *negative*, meaning residues will be
    *added* to the original sequence;
    3. There is no difference between the number of residues in the range where the 
    subsitution is made and the length of the substitution.

    The easiest place to implement a substitution is during the residue/segment
    build (I think).

    """
    req_attr=AncestorAwareMod.req_attr+['chainID','resseqnum1','insertion1','resseqnum2','insertion2','subseq']
    yaml_header='substitutions'
    modtype='seqmod'

    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)

    @__init__.register(str)
    def _from_shortcode(self,shortcode):
        # shortcode format: C:nnn-ccc,abcdef
        # C -- chainID
        # nnn -- N-terminal resid of sequence to be replaced
        # ccc -- C-terminal resid of sequence to be replaced
        # abcdef -- the 1-byte rescode sequence to be substituted
        p1=shortcode.split(':')
        chainID=p1[0]
        p2=p1[1].split(',')
        seqrange=p2[0]
        subseq=p2[1]
        seq=seqrange.split('-')
        r1,i1=split_ri(seq[0])
        r2,i2=split_ri(seq[1])
        input_dict={
            'chainID':chainID,
            'resseqnum1':int(r1),
            'insertion1':i1,
            'resseqnum2':int(r2),
            'insertion2':i2,
            'subseq':subseq
        }
        super().__init__(input_dict)

class SubstitutionList(AncestorAwareModList):
    pass

class Deletion(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['chainID','resseqnum1','insertion1','resseqnum2','insertion2']
    opt_attr=AncestorAwareMod.opt_attr+['model']
    yaml_header='deletions'
    modtype='seqmod'

    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)

    @__init__.register(str)
    def _from_shortcode(self,shortcode):
        # shortcode format: C:nnn-ccc
        # C -- chainID
        # nnn -- N-terminal resid of sequence to be deleted
        # ccc -- C-terminal resid of sequence to be deleted
        p1=shortcode.split(':')
        chainID=p1[0]
        seq=p1[1].split('-')
        r1,i1=split_ri(seq[0])
        r2,i2=split_ri(seq[1])          
        input_dict={
            'chainID':chainID,
            'resseqnum1':int(r1),
            'insertion1':i1,
            'resseqnum2':int(r2),
            'insertion2':i2
        }
        super().__init__(input_dict)

class DeletionList(AncestorAwareModList):
    pass

class Crot(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['angle','degrees']
    opt_attr=AncestorAwareMod.opt_attr+['chainID','resseqnum1','resseqnum2','segname','atom1','atom2','segname1','segname2','segnamei','resseqnumi','atomi','segnamejk','resseqnumj','atomj','resseqnumk','atomk']
    attr_choices=AncestorAwareMod.attr_choices.copy()
    attr_choices.update({'angle':['PHI','PSI','OMEGA','CHI1','CHI2','GLYCAN','LINK','ANGLEIJK']})
    opt_attr_deps=AncestorAwareMod.opt_attr_deps.copy()
    opt_attr_deps.update({
        'PHI':['chainID','resseqnum1','resseqnum2'],
        'PSI':['chainID','resseqnum1','resseqnum2'],
        'OMEGA':['chainID','resseqnum1','resseqnum2'],
        'CHI1':['chainID','resseqnum1'],
        'CHI2':['chainID','resseqnum1'],
        'GLYCAN':['segname','resseqnum1','atom1','resseqnum2','atom2'],
        'LINK':['segname1','segname2','resseqnum1','atom1','resseqnum2','atom2'],
        'ANGLEIJK':['segnamei','resseqnumi','atomi','segnamejk','resseqnumj','atomj','resseqnumk','atomk']
        })
    yaml_header='crotations'
    modtype='coormod'

    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)
    
    @__init__.register(str)
    def _from_shortcode(self,shortcode):
        dat=shortcode.split(',')
        input_dict={}
        input_dict['angle']=dat[0].upper()
        if input_dict['angle']=='PHI' or input_dict['angle']=='PSI' or input_dict['angle']=='OMEGA':
            # this is a backbone torsion, so we need both an owner
            # residue and a residue indicating the end of the 
            # set of residues that will be reoriented by the
            # rotation
            input_dict['chainID']=dat[1]
            input_dict['resseqnum1']=int(dat[2])
            input_dict['resseqnum2']=int(dat[3])
            input_dict['degrees']=float(dat[4])
        elif input_dict['angle']=='CHI1' or input_dict['angle']=='CHI2':
            input_dict['chainID']=dat[1]
            input_dict['resseqnum1']=int(dat[2])
            input_dict['resseqnum2']=-1
            input_dict['degrees']=float(dat[3])
        elif input_dict['angle']=='GLYCAN':
            input_dict['segname']=dat[1]
            input_dict['resseqnum1']=int(dat[2])
            input_dict['atom1']=dat[3]
            input_dict['resseqnum2']=int(dat[4])
            input_dict['atom2']=dat[5]
            input_dict['degrees']=float(dat[6])
        elif input_dict['angle']=='LINK':
            input_dict['segname1']=dat[1]
            input_dict['resseqnum1']=int(dat[2])
            input_dict['atom1']=dat[3]
            input_dict['segname2']=dat[4]
            input_dict['resseqnum2']=int(dat[5])
            input_dict['atom2']=dat[6]
            input_dict['degrees']=float(dat[7])
        elif input_dict['angle']=='ANGLEIJK':
            input_dict['segnamei']=dat[1]
            input_dict['resseqnumi']=int(dat[2])
            input_dict['atomi']=dat[3]
            input_dict['segnamejk']=dat[4]
            input_dict['resseqnumj']=int(dat[5])
            input_dict['atomj']=dat[6]
            input_dict['resseqnumk']=int(dat[7])
            input_dict['atomk']=dat[8]
            input_dict['degrees']=float(dat[9])
        super().__init__(input_dict)
    
    def to_shortcode(self):
        if 'shortcode' in self.__dict__:
            return
        ret=[f'{self.angle}']
        if self.angle=='PHI' or self.angle=='PSI' or self.angle=='OMEGA':
            # this is a backbone torsion, so we need both an owner
            # residue and a residue indicating the end of the 
            # set of residues that will be reoriented by the
            # rotation
            ret.append(f'{self.chainID}')
            ret.append(f'{self.resseqnum1}')
            ret.append(f'{self.resseqnum2}')
            ret.append(f'{self.degrees:.4f}')
        elif self.angle=='CHI1' or self.angle=='CHI2':
            ret.append(f'{self.chainID}')
            ret.append(f'{self.resseqnum1}')
            ret.append(f'-1')
            ret.append(f'{self.degrees:.4f}')
        elif self.angle=='GLYCAN':
            ret.append(f'{self.segname}')
            ret.append(f'{self.resseqnum1}')
            ret.append(f'{self.atom1}')
            ret.append(f'{self.resseqnum2}')
            ret.append(f'{self.atom2}')
            ret.append(f'{self.degrees:.4f}')
        elif self.angle=='LINK':
            ret.append(f'{self.segname1}')
            ret.append(f'{self.resseqnum1}')
            ret.append(f'{self.atom1}')
            ret.append(f'{self.segname2}')
            ret.append(f'{self.resseqnum2}')
            ret.append(f'{self.atom2}')
            ret.append(f'{self.degrees:.4f}')
        elif self.angle=='ANGLEIJK':
            ret.append(f'{self.segnamei}')
            ret.append(f'{self.resseqnumi}')
            ret.append(f'{self.atomi}')
            ret.append(f'{self.segnamejk}')
            ret.append(f'{self.resseqnumj}')
            ret.append(f'{self.atomj}')
            ret.append(f'{self.resseqnumk}')
            ret.append(f'{self.atomk}')
            ret.append(f'{self.degrees:.4f}')
        self.shortcode=','.join(ret)
    
    def __str__(self):
        self.to_shortcode()
        return self.shortcode
        
    def psfgen_lines(self,**kwargs):
        molid=kwargs.get('molid','top')
        endIsCterm=kwargs.get('endIsCterm',False)
        retlines=[]
        if self.angle in ['PHI','PSI','OMEGA']:  # this is a backbone bond
            retlines.append('set r1 [[atomselect {} "chain {} and resid {} and name CA"] get residue]'.format(molid,self.chainID,self.resseqnum1))
            retlines.append('set r2 [[atomselect {} "chain {} and resid {} and name CA"] get residue]'.format(molid,self.chainID,self.resseqnum2))
            if endIsCterm:
                retlines.append('Crot_{}_toCterm $r1 $r2 {} {} {}'.format(self.angle.lower(),self.chainID,molid,self.degrees))
            else:
                retlines.append('Crot_{} $r1 $r2 {} {} {}'.format(self.angle.lower(),self.chainID,molid,self.degrees))
        elif self.angle in ['CHI1','CHI2']:  # this is a side-chain bond
            retlines.append('set r1 [[atomselect {} "chain {} and resid {} and name CA"] get residue]'.format(molid,self.chainID,self.resseqnum1))
            retlines.append('SCrot_{} $r1 {} {} {}'.format(self.angle.lower(),self.chainID,molid,self.degrees))
        elif self.angle=='GLYCAN':  # intra-glycan rotation
            retlines.append('set sel [atomselect {} "segname {}"]'.format(molid,self.segname))
            retlines.append('set i [[atomselect {} "segname {} and resid {} and name {}"] get index]'.format(molid,self.segname,self.resseqnum1,self.atom1))
            retlines.append('set j [[atomselect {} "segname {} and resid {} and name {}"] get index]'.format(molid,self.segname,self.resseqnum2,self.atom2))
            retlines.append('genbondrot {} $sel $i $j {}'.format(molid,self.degrees))
        elif self.angle=='LINK': # ASN-GLYcan rotation
            retlines.append('set sel [atomselect {} "segname {} {}"]'.format(molid,self.segname1,self.segname2))
            retlines.append('set i [[atomselect {} "segname {} and resid {} and name {}"] get index]'.format(molid,self.segname1,self.resseqnum1,self.atom1))
            retlines.append('set j [[atomselect {} "segname {} and resid {} and name {}"] get index]'.format(molid,self.segname2,self.resseqnum2,self.atom2))
            retlines.append('genbondrot {} $sel $i $j {}'.format(molid,self.degrees))
        elif self.angle=='ANGLEIJK':
            retlines.extend([
                'set rotsel [atomselect {} "segname {}"]'.format(molid,self.segnamejk),
                'set ri [lindex [[atomselect {} "segname {} and resid {} and name {}"] get {{x y z}}] 0]'.format(molid,self.segnamei,self.resseqnumi,self.atomi),
                'set rj [lindex [[atomselect {} "segname {} and resid {} and name {}"] get {{x y z}}] 0]'.format(molid,self.segnamejk,self.resseqnumj,self.atomj),
                'set rk [lindex [[atomselect {} "segname {} and resid {} and name {}"] get {{x y z}}] 0]'.format(molid,self.segnamejk,self.resseqnumk,self.atomk),
                'set rij [vecsub $ri $rj]',
                'set rjk [vecsub $rj $rk]',
                'set cijk [veccross $rij $rjk]',
                '$rotsel move [trans center $rj origin $rj axis $cijk {} degrees]'.format(self.degrees)
            ])
        return retlines

class CrotList(AncestorAwareModList):
    pass

class SSBond(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['chainID1','resseqnum1','insertion1','chainID2','resseqnum2','insertion2']
    opt_attr=AncestorAwareMod.opt_attr+['serial_number','residue1','residue2','resname1','resname2','sym1','sym2','length','ptnr1_auth_asym_id','ptnr2_auth_asym_id','ptnr1_auth_seq_id','ptnr2_auth_seq_id']
    yaml_header='ssbonds'
    modtype='topomod'
    PDB_keyword='SSBOND'
    mmCIF_name='struct_conn'

    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)

    @__init__.register(PDBRecord)
    def _from_pdbrecord(self,pdbrecord):
        input_dict={
            'serial_number':pdbrecord.serNum,
            'resname1':pdbrecord.residue1.resName,
            'chainID1':pdbrecord.residue1.chainID,
            'resseqnum1':pdbrecord.residue1.seqNum,
            'insertion1':pdbrecord.residue1.iCode,
            'resname2':pdbrecord.residue2.resName,
            'chainID2':pdbrecord.residue2.chainID,
            'resseqnum2':pdbrecord.residue2.seqNum,
            'insertion2':pdbrecord.residue2.iCode,
            'sym1':pdbrecord.sym1,
            'sym2':pdbrecord.sym2,
            'length':pdbrecord.length,
            'residue1':None,
            'residue2':None
        }
        super().__init__(input_dict)
    
    @__init__.register(CIFdict)
    def _from_cifdict(self,cd):
        input_dict={
            'serial_number':int(cd['id'].strip('disulf')),
            'resname1':'CYS',
            'resname2':'CYS',
            'chainID1':cd['ptnr1_label_asym_id'],
            'chainID2':cd['ptnr2_label_asym_id'],
            'resseqnum1':int(cd['ptnr1_label_seq_id']),
            'resseqnum2':int(cd['ptnr2_label_seq_id']),
            'insertion1':cd['pdbx_ptnr1_pdb_ins_code'],
            'insertion2':cd['pdbx_ptnr2_pdb_ins_code'],
            'sym1':cd['ptnr1_symmetry'],
            'sym2':cd['ptnr2_symmetry'],
            'length':float(cd['pdbx_dist_value']),
            'residue1':None,
            'residue2':None,
            'ptnr1_auth_asym_id':cd['ptnr1_auth_asym_id'],
            'ptnr2_auth_asym_id':cd['ptnr2_auth_asym_id'],
            'ptnr1_auth_seq_id':cd['ptnr1_auth_seq_id'],
            'ptnr2_auth_seq_id':cd['ptnr2_auth_seq_id']
        }
        super().__init__(input_dict)
    
    @__init__.register(str)
    def _from_shortcode(self,shortcode):
        # shortcode format: C_RRR-D_SSS
        # C, D chainIDs
        # RRR, SSS resids
        s1=shortcode.split('-')
        r1=s1[0].split('_')
        r2=s1[1].split('_')
        input_dict={
            'serial_number':0,
            'resname1':'CYS',
            'resname2':'CYS',
            'chainID1':r1[0],
            'chainID2':r2[0],
            'resseqnum1':int(r1[1]),
            'resseqnum2':int(r2[1]),
            'insertion1':'',
            'insertion2':'',
            'length':0.0,
            'sym1':'',
            'sym2':''
        }
        super().__init__(input_dict)

    def __str__(self):
        return f'{self.chainID1}_{self.resseqnum1}{self.insertion1}-{self.chainID2}_{self.resseqnum2}{self.insertion2}'

    def pdb_line(self):
        pdbline='SSBOND'+\
                '{:4d}'.format(self.serial_number)+\
                '{:>4s}'.format(self.resname1)+\
                '{:>2s}'.format(self.chainID1)+\
                '{:5d}'.format(self.resseqnum1)+\
                '{:1s}'.format(self.insertion1)+\
                '{:>6s}'.format(self.resname2)+\
                '{:>2s}'.format(self.chainID2)+\
                '{:5d}'.format(self.resseqnum2)+\
                '{:1s}'.format(self.insertion2)+\
                ' '*23+\
                '{:>6s}'.format(self.sym1)+\
                '{:>7s}'.format(self.sym2)+\
                '{:6.2f}'.format(self.length)
        return pdbline
    
    def write_TcL(self,W:Psfgen,transform):
        chainIDmap=transform.chainIDmap 
        # ok since these are only going to reference protein segments; protein segment names are the chain IDs
        c1=chainIDmap.get(self.chainID1,self.chainID1)
        c2=chainIDmap.get(self.chainID2,self.chainID2)
        r1=self.resseqnum1
        r2=self.resseqnum2
        W.addline(f'patch DISU {c1}:{r1} {c2}:{r2}')

class SSBondList(AncestorAwareModList):
    def assign_residues(self,Residues):
        ignored_by_ptnr1=self.assign_objs_to_attr('residue1',Residues,resseqnum='resseqnum1',chainID='chainID1',insertion='insertion1')
        ignored_by_ptnr2=self.assign_objs_to_attr('residue2',Residues,resseqnum='resseqnum2',chainID='chainID1',insertion='insertion2')
        return self.__class__(ignored_by_ptnr1+ignored_by_ptnr2)
    def write_TcL(self,W:Psfgen,transform):
        for s in self:
            # if mods.ssbondsdelete.is_deleted(s,transform):
            #     W.comment(f'Deleted ssbond: {str(s)}')
            #     continue
            s.write_TcL(W,transform)
        # user may have added some ssbonds
        # for s in mods.ssbonds:
        #     W.comment(f'Below is a user-added ssbond:')
        #     s.write_TcL(W,transform)
    def prune_mutations(self,Mutations):
        pruned=self.__class__([])
        for m in Mutations:
            left=self.get(chainID1=m.chainID,resseqnum1=m.resseqnum,insertion1=m.insertion)
            if left:
                pruned.append(self.remove(left))
            right=self.get(chainID2=m.chainID,resseqnum2=m.resseqnum,insertion2=m.insertion)
            if right:
                pruned.append(self.remove(right))
        return pruned
    
class SSBondDelete(SSBond):
    yaml_header='ssbondsdelete'
    modtype='topomod'

class SSBondDeleteList(SSBondList):
    def is_deleted(self,a_SSBond):
        if self.get(
            chainID1=a_SSBond.chainID1,
            chainID2=a_SSBond.chainID2,
            resseqnum1=a_SSBond.resseqnum1,
            resseqnum2=a_SSBond.resseqnum2):
            return True
        if self.get(
            chainID2=a_SSBond.chainID1,
            chainID1=a_SSBond.chainID2,
            resseqnum2=a_SSBond.resseqnum1,
            resseqnum1=a_SSBond.resseqnum2):
            return True
        return False

class Link(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['name1','chainID1','resseqnum1','insertion1','name2','chainID2','resseqnum2','insertion2']
    opt_attr=AncestorAwareMod.opt_attr+['altloc1','altloc2','resname1','resname2','sym1','sym2','link_distance','segname1','segname2','residue1','residue2','atom1','atom2','empty','segtype1','segtype2','ptnr1_auth_asym_id','ptnr2_auth_asym_id','ptnr1_auth_seq_id','ptnr2_auth_seq_id','ptnr1_auth_comp_id','ptnr2_auth_comp_id']    
    yaml_header='links'
    PDB_keyword='LINK'
    modtype='topomod'

    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)
    
    @__init__.register(PDBRecord)
    def _from_pdbrecord(self,pdbrecord):
        input_dict={
            'name1': pdbrecord.name1,
            'resname1':pdbrecord.residue1.resName,
            'chainID1':pdbrecord.residue1.chainID,
            'resseqnum1':pdbrecord.residue1.seqNum,
            'insertion1':pdbrecord.residue1.iCode,
            'name2': pdbrecord.name2,
            'resname2':pdbrecord.residue2.resName,
            'chainID2':pdbrecord.residue2.chainID,
            'resseqnum2':pdbrecord.residue2.seqNum,
            'insertion2':pdbrecord.residue2.iCode,
            'altloc2':pdbrecord.altLoc2,
            'altloc1':pdbrecord.altLoc1,
            'sym1':pdbrecord.sym1,
            'sym2':pdbrecord.sym2,
            'link_distance':pdbrecord.length,
            'segname1':pdbrecord.residue1.chainID,
            'segname2':pdbrecord.residue2.chainID,
            }
        input_dict.update({
                'residue1':None,
                'residue2':None,
                'atom1':None,
                'atom2':None,
                'empty':False,
                'segtype1':'UNSET',
                'segtype2':'UNSET',
            })
        super().__init__(input_dict)
    
    @__init__.register(CIFdict)
    def _from_cifdict(self,cd):
        p1_seq_id=cd['ptnr1_label_seq_id']
        if p1_seq_id=='':
            p1_seq_id=cd['ptnr1_auth_seq_id']
        p2_seq_id=cd['ptnr2_label_seq_id']
        if p2_seq_id=='':
            p2_seq_id=cd['ptnr2_auth_seq_id']
        input_dict={}
        input_dict['name1']=cd['ptnr1_label_atom_id']
        input_dict['altloc1']=cd['pdbx_ptnr1_label_alt_id']
        input_dict['resname1']=cd['ptnr1_label_comp_id']
        input_dict['chainID1']=cd['ptnr1_label_asym_id']
        input_dict['resseqnum1']=int(p1_seq_id)
        input_dict['insertion1']=cd['pdbx_ptnr1_pdb_ins_code']
        input_dict['name2']=cd['ptnr2_label_atom_id']
        input_dict['altloc2']=cd['pdbx_ptnr2_label_alt_id']
        input_dict['resname2']=cd['ptnr2_label_comp_id']
        input_dict['chainID2']=cd['ptnr2_label_asym_id']
        input_dict['resseqnum2']=int(p2_seq_id)
        input_dict['insertion2']=cd['pdbx_ptnr2_pdb_ins_code']
        input_dict['sym1']=cd['ptnr1_symmetry']
        input_dict['sym2']=cd['ptnr2_symmetry']
        input_dict['link_distance']=float(cd['pdbx_dist_value'])
        input_dict['segname1']=input_dict['chainID1']
        input_dict['segname2']=input_dict['chainID2']
        input_dict.update({'ptnr1_auth_asym_id':cd['ptnr1_auth_asym_id'],
            'ptnr2_auth_asym_id':cd['ptnr2_auth_asym_id'],
            'ptnr1_auth_seq_id':cd['ptnr1_auth_seq_id'],
            'ptnr2_auth_seq_id':cd['ptnr2_auth_seq_id']})
        input_dict.update({
                'residue1':None,
                'residue2':None,
                'atom1':None,
                'atom2':None,
                'empty':False,
                'segtype1':'UNSET',
                'segtype2':'UNSET',
            })
        super().__init__(input_dict)

    def write_TcL(self,W:Psfgen,transform):
        chainIDmap=transform.chainIDmap
        seg1=self.residue1.chainID
        seg1=chainIDmap.get(seg1,seg1)
        seg2=self.residue2.chainID
        seg2=chainIDmap.get(seg2,seg2)
        if self.resname1=='ASN' and self.segtype2=='glycan':
            W.addline(f'patch NGLB {seg1}:{self.resseqnum1}{self.insertion1} {seg2}:{self.resseqnum2}{self.insertion2}')
        else:
            # this is likely an intra-glycan linkage
            if self.name2=='C1' and self.segtype1=='glycan':
                W.addline(f'set cn {self.name1[1]}')
                W.addline(f'set abi [axeq {self.resseqnum2} 0 {seg2} {self.name2} {self.resseqnum1}]')
                W.addline(f'set abj [axeq {self.resseqnum1} 0 {seg1} {self.name1} -1]')
                if self.name1=='O6':
                    W.addline('if { $abi == "a" } { set abi A }')
                    W.addline('if { $abi == "b" } { set abi B; set abj T }')
                    W.addline('if { $abj == "b" } { set abj T }')
                W.addline('set pres "1$cn$abi$abj"')
                W.addline(f'patch $pres {seg1}:{self.resseqnum1}{self.insertion1} {seg2}:{self.resseqnum2}{self.insertion2}')
            elif self.name1=='C1' and self.segtype2=='glycan':
                cmdj=f'[axeq {self.resseqnum2} 0 {seg2} {self.name2} {self.resseqnum1}]'
                cmdi=f'[axeq {self.resseqnum1} 0 {seg1} {self.name1} -1]'
                W.addline(f'patch 1{self.name2[1]:1s}{cmdi}{cmdj} {seg2}:{self.resseqnum2}{self.insertion2} {seg1}:{self.resseqnum1}{self.insertion1}')
            elif self.name1=='O6' and self.name2=='C2':
                W.addline(f'patch SA26AT {seg1}:{self.resseqnum1}{self.insertion1} {seg2}:{self.resseqnum2}{self.insertion2}')
            else:
                logger.warning(f'Could not identify patch for link: {str(self)}')
                W.comment(f'No patch found for {str(self)}')
        
    def __str__(self):
        return f'{self.chainID1}{self.resname1}{self.resseqnum1}{self.insertion1}-{self.chainID2}{self.resname2}{self.resseqnum2}{self.insertion2}'

class LinkList(AncestorAwareModList):
    def assign_residues(self,Residues):
        ignored_by_ptnr1=self.assign_objs_to_attr('residue1',Residues,resseqnum='resseqnum1',chainID='chainID1',insertion='insertion1')
        ignored_by_ptnr2=self.assign_objs_to_attr('residue2',Residues,resseqnum='resseqnum2',chainID='chainID2',insertion='insertion2')
        for link in self:
            link.residue1.linkTo(link.residue2,link)
            link.atom1=link.residue1.atoms.get(name=link.name1,altloc=link.altloc1)
            link.atom2=link.residue2.atoms.get(name=link.name2,altloc=link.altloc2)
            link.segtype1=link.residue1.segtype
            link.segtype2=link.residue2.segtype
        # do cross-assignment to find true orphan links and dangling links
        orphan_1=ignored_by_ptnr1.assign_objs_to_attr('residue2',Residues,resseqnum='resseqnum2',chainID='chainID2',insertion='insertion2')
        orphan_2=ignored_by_ptnr2.assign_objs_to_attr('residue1',Residues,resseqnum='resseqnum1',chainID='chainID1',insertion='insertion1')
        orphans=orphan_1+orphan_2
        rlist=[]
        for link in ignored_by_ptnr1:
            rlist,list=link.residue2.get_down_group()
            rlist.insert(0,link.residue2)
            for r in rlist:
                Residues.remove(r)
        return Residues.__class__(rlist),self.__class__(ignored_by_ptnr1+ignored_by_ptnr2)
    
    def write_TcL(self,W:Psfgen,transform):
        for l in self:
            l.write_TcL(W,transform)

    def remove_orphan_residues(self,Links,Residues):
        for dl in self:
            reslist,lnlist=dl.residue2.get_down_group()
            reslist.insert(0,dl.residue2)
            for r in reslist:
                Residues.remove(r)
            for l in lnlist:
                Links.remove(l)
        return reslist,lnlist

    def prune_mutations(self,Mutations,Segments):
        pruned={'residues':[],'links':self.__class__([]),'segments':[]}
        for m in Mutations:
            rlist,llist=[],[]
            left=self.get(chainID1=m.chainID,resseqnum1=m.resseqnum,insertion1=m.insertion)
            right=self.get(chainID2=m.chainID,resseqnum2=m.resseqnum,insertion2=m.insertion)
            if left:  # this is a link in which this mutation is the left member
                self.remove(left) # get rid of this link
                # we need to remove residue2 and everything downstream
                # remove downstream residues!
                rlist,llist=left.residue2.get_down_group()
                rlist.insert(0,left.residue2)
            elif right: # this is a link in which this mutation is the right member (should be very rare)
                self.remove(right)
                rlist,llist=right.residue2.get_down_group()
                rlist.insert(0,right.residue2)
            if rlist and llist:
                # logger.debug(f'Deleting residues down from and including {str(rlist[0])} due to a mutation')
                S=Segments.get_segment_of_residue(rlist[0])
                for r in rlist:
                    # logger.debug(f'...{str(r)}')
                    pruned['residues'].append(S.residues.remove(r))
                if len(S.residues)==0:
                    # logger.debug(f'All residues of {S.psfgen_segname} are deleted; {S.psfgen_segname} is deleted')
                    pruned['segments'].append(Segments.remove(S))
                for l in llist:
                    pruned['links'].append(self.remove(l))
        return pruned

    def apply_segtypes(self,map):
        self.map_attr('segtype1','resname1',map)
        self.map_attr('segtype2','resname2',map)

class Graft(AncestorAwareMod):
    yaml_header='grafts'
    modtype='topomod'
    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)

class GraftList(AncestorAwareModList):
    pass

class Insertion(AncestorAwareMod):
    yaml_header='insertions'
    modtype='seqmod'
    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)
    pass

class InsertionList(AncestorAwareModList):
    pass

class Cleavage(AncestorAwareMod):
    yaml_header='cleavages'
    modtype='seqmod'
    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)
    pass

class CleavageList(AncestorAwareModList):
    pass


# def apply_psf_info(p_struct,psf):
#     if os.path.exists(psf):
#         with open(psf,'r') as f:
#             psf_lines=f.read().split('\n')
#     metadata_lines=[x for x in psf_lines if x.startswith('REMARKS')]
#     for ml in metadata_lines:
#         words=ml.split()
#         if ' '.join(words[1:3])=='patch DISU':
#             c1,r1i=words[3].split(':')
#             c2,r2i=words[4].split(':')
#             r1,i1=split_ri(r1i)
#             r2,i2=split_ri(r2i)
#             p_struct['SSBONDS'].append(SSBond({'chainID1':c1,'chainID2':c2,'resseqnum1':r1,'resseqnum2':r2,'insertion1':i1,'insertion2':i2}))
        # elif ' '.join(words[1:3])=='patch '
            # add a new SSBOND record to p_struct
        # elif ...