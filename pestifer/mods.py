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
    * Cfusion -- fuse residues from named residue range of named 
    protein chain of named input coordinate file to C-terminus of 
    named segment of base molecule

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
from .scriptwriters import Psfgen,VMD
from .config import res_123
from functools import singledispatchmethod
from .coord import measure_dihedral, ic_reference_closest

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
    """A class for handling SEQADV/seq_dif records in input structure files 
    
    SEQADV/seq_dif records declare differences between the actual sequence of the molecule and the sequence
    of its database reference.  This is how authors report things like accidental mutations (conflicts),
    engineered mutations, and other differences.  These records are residue-specific.

    Attributes
    ----------
    req_attr: list
        * idCode: identification code of the seqadv record
        * resname: residue name
        * chainID: chain identifier
        * resseqnum: sequence number of residue in chain
        * insertion: insertion code
        * typekey: reports the type of the seqadv/seq_dif based on the description provided by the author
    opt_attr: list
        * database: name of database
        * dbAccession: accession code of this molecule's sequence in that database
        * dbRes: name of residue in the database
        * dbSeq: sequence position of residue in the database
        * pdbx_pdb_strand_id: (mmCIF format) chain identifier according to author
        * pdbx_auth_seq_num: (mmCIF format) residue sequence position according to author
        * pdbx_ordinal: some integer, i don't know
        * residue: holder for actual residue this seqadv refers to

    Methods
    -------
    seqadv_details_keyword(text)
        returns the typekey for the given text
    
    pdb_line()
        write the seqadv as it would appear in a PDB file
    
    assign_residue(Residues)
        given the list of residues 'Residues', finds
        the residue that matches the specs of the caller
        assigns it to the 'residue ' element
    
    update_residue()
        assigns the chainID attribute of the caller to
        the value of the chainID of the caller's 
        residue attribute
    
    """
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
    """A class for handling lists of Seqadvs
    
    Methods
    -------
    assign_residues(Residues):
        calls assign_residue for each member element
        and returns a list of residues that were
        not assigned (out of the list Residues)
    
    update_from_residues()
        calls update_from_residue for each member element
    
    """
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
    
    Attributes
    ----------
    req_attr: list
        * chainID: chain identifier
        * origresname: original (wild-type) residue name
        * resseqnum: residue sequence position
        * insertion: residue insertion code
        * newresname: mutated residue name
        * typekey: indicates what type of mutation this is
    opt_attr: list
        * pdb_auth_seq_num: author's sequence number of the residue (value used in PDB file)

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
        """Returns the string to be written in a psfgen input file within a segment """
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
    more contiguous residues are replaced by a new sequence of one or more residues.

    Attributes
    ----------
    req_attr: list
        * chainID: chain identifier
        * resseqnum1: N-terminal resid of sequence to be replaced
        * insertion1: insertion code of N-terminal residue
        * resseqnum2: C-terminal resid of sequence to be replaced
        * insertion2: insertion code of the C-terminal residue
        * subseq: 1-byte rescode sequence to be substituted in

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
    """A class for handling deletions 
    
    A deletion is a user-specified modification in which a sequence of one or
    more contiguous residues are deleted.

    Attributes
    ----------
    req_attr: list
        * chainID: chain identifier
        * resseqnum1: N-terminal resid of sequence to be deleted
        * insertion1: insertion code of N-terminal residue
        * resseqnum2: C-terminal resid of sequence to be deleted
        * insertion2: insertion code of the C-terminal residue

    """
    req_attr=AncestorAwareMod.req_attr+['chainID','resseqnum1','insertion1','resseqnum2','insertion2']
    # opt_attr=AncestorAwareMod.opt_attr+['model']
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

class Cfusion(AncestorAwareMod):
    """A class for handling fusions of residues represented by an existing 
    coordinate file to the C-termini of base-molecule segments
    
    A Cfusion is a user-specified modification which fuses residues from a
    named residue range of a named protein chain of a named input coordinate 
    file to the C-terminus of a named segment of base molecule.

    Attributes
    ----------
    req_attr: list
        * sourcefile: name of source coordinate PDB file of fusion
        * chainID: chainID of the fusion sequence in sourcefile
        * resseqnum1: N-terminal resid of fusion sequence
        * insertion1: insertion code of N-terminal residue
        * resseqnum2: C-terminal resid of fusion sequence
        * insertion2: insertion code of the C-terminal residue
        * tosegment: name of segment to which fusion is made

    """
    req_attr=AncestorAwareMod.req_attr+['sourcefile','chainID','resseqnum1','insertion1','resseqnum2','insertion2','tosegment','id']
    yaml_header='Cfusions'
    modtype='seqmod'
    _Cfusion_counter=0
    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)

    @__init__.register(str)
    def _from_shortcode(self,shortcode):
        # shortcode format: filename:C:nnn-ccc,S
        # filename -- fusion coordinate filename
        # C -- chainID
        # nnn -- N-terminal resid of fusion sequence
        # ccc -- C-terminal resid of fusion sequence
        # S -- tosegment, segment in base-molecule the fusion is fused to
        p1=shortcode.split(':')
        sourcefile=p1[0]
        chainID=p1[1]
        p2=p1[2].split(',')
        seqrange=p2[0]
        tosegment=p2[1]
        seq=seqrange.split('-')
        r1,i1=split_ri(seq[0])
        r2,i2=split_ri(seq[1])
        input_dict={
            'sourcefile':sourcefile,
            'chainID':chainID,
            'resseqnum1':int(r1),
            'insertion1':i1,
            'resseqnum2':int(r2),
            'insertion2':i2,
            'tosegment':tosegment,
            'id':Cfusion._Cfusion_counter
        }
        Cfusion._Cfusion_counter+=1
        super().__init__(input_dict)    

    def write_pre_segment(self,W:Psfgen):
        W.addline(f'set topid [molinfo top get id]')
        W.addline(f'mol new {self.sourcefile}')
        W.addline(f'set cfusid [molinfo top get id]')
        W.addline(f'mol top $topid')
        W.addline(f'set fusres [atomselect $cfusid "protein and chain {self.chainID} and resid {self.resseqnum1}{self.insertion1} to {self.resseqnum2}{self.insertion2}"]')
        self.segfile=f'Cfusion{self.id}_{self.chainID}_{self.resseqnum1}{self.insertion1}_to_{self.resseqnum2}{self.insertion2}.pdb'
        W.addline(f'$fusres writepdb {self.segfile}')
        W.addline(f'delete $cfusid')
    def write_in_segment(self,W:Psfgen):
        W.addline (f'    pdb {self.segfile}')
    def write_post_segment(self,W:Psfgen):
        W.addline(f'coordpdb {self.segfile} {self.tosegment}')

class CfusionList(AncestorAwareModList):
    pass

class Crot(AncestorAwareMod):
    """A class for managing so-called 'C-rotations'
    
    A C-rotation is a transformation in which atoms are rotated around a given bond by a given amount.  The "C" 
    designation means that only the "downstream" atoms of the bond are moved; upstream atoms, along with the atoms
    of the bond itself, are not moved.  The set of upstream atoms and the set of downstream atoms must naturally
    have no topological connection *other* than the bond itself.  Typically, this can be used to execute rotation
    of a backbone look in a C-terminal loop, a side-chain angle, or a glycan angle, usually in the service
    of reducing steric clashses.  The primary job of this class is to translate the C-rotation shortcodes
    specified by the user into TcL commands to be incorporated in a psfgen script.

    NOTE: This is currently implemented in the cfapdbparse (2020) format, and has not been thoroughly tested.

    """
    req_attr=AncestorAwareMod.req_attr+['angle']
    opt_attr=AncestorAwareMod.opt_attr+['chainID','resseqnum1','resseqnum2','segname','atom1','atom2','segname1','segname2','segnamei','resseqnumi','atomi','segnamejk','resseqnumj','atomj','resseqnumk','atomk','degrees']
    attr_choices=AncestorAwareMod.attr_choices.copy()
    attr_choices.update({'angle':['PHI','PSI','OMEGA','CHI1','CHI2','GLYCAN','LINK','ANGLEIJK','ALPHA']})
    opt_attr_deps=AncestorAwareMod.opt_attr_deps.copy()
    opt_attr_deps.update({
        'PHI':['chainID','resseqnum1','resseqnum2'],
        'PSI':['chainID','resseqnum1','resseqnum2'],
        'OMEGA':['chainID','resseqnum1','resseqnum2'],
        'CHI1':['chainID','resseqnum1'],
        'CHI2':['chainID','resseqnum1'],
        'GLYCAN':['segname','resseqnum1','atom1','resseqnum2','atom2'],
        'LINK':['segname1','segname2','resseqnum1','atom1','resseqnum2','atom2'],
        'ANGLEIJK':['segnamei','resseqnumi','atomi','segnamejk','resseqnumj','atomj','resseqnumk','atomk'],
        'ALPHA':['chainID','resseqnum1','resseqnum2']
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
        elif input_dict['angle']=='ALPHA':
            input_dict['chainID']=dat[1]
            input_dict['resseqnum1']=int(dat[2])
            input_dict['resseqnum2']=int(dat[3])

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
        elif self.angle=='ALPHA':
            ret.append(f'{self.chainID}')
            ret.append(f'{self.resseqnum1}')
            ret.append(f'{self.resseqnum2}')
        self.shortcode=','.join(ret)
    
    def __str__(self):
        self.to_shortcode()
        return self.shortcode
    
    def write_TcL(self,W:VMD,transform,**kwargs):
        chainIDmap=transform.chainIDmap
        the_chainID=chainIDmap.get(self.chainID,self.chainID)
        molid_varname=W.molid_varname
        molid=f'${molid_varname}'
        endIsCterm=kwargs.get('endIsCterm',True)
        if self.angle in ['PHI','PSI','OMEGA']:
            W.addline('set r1 [[atomselect {} "chain {} and resid {} and name CA"] get residue]'.format(molid,the_chainID,self.resseqnum1))
            W.addline('set r2 [[atomselect {} "chain {} and resid {} and name CA"] get residue]'.format(molid,the_chainID,self.resseqnum2))
            if endIsCterm:
                W.addline('Crot_{}_toCterm $r1 $r2 {} {} {}'.format(self.angle.lower(),the_chainID,molid,self.degrees))
            else:
                W.addline('Crot_{} $r1 $r2 {} {} {}'.format(self.angle.lower(),the_chainID,molid,self.degrees))
        elif self.angle in ['CHI1','CHI2']:  # this is a side-chain bond
            W.addline('set r1 [[atomselect {} "chain {} and resid {} and name CA"] get residue]'.format(molid,the_chainID,self.resseqnum1))
            W.addline('SCrot_{} $r1 {} {} {}'.format(self.angle.lower(),the_chainID,molid,self.degrees))
        elif self.angle=='GLYCAN':  # intra-glycan rotation
            W.addline('set sel [atomselect {} "segname {}"]'.format(molid,self.segname))
            W.addline('set i [[atomselect {} "segname {} and resid {} and name {}"] get index]'.format(molid,self.segname,self.resseqnum1,self.atom1))
            W.addline('set j [[atomselect {} "segname {} and resid {} and name {}"] get index]'.format(molid,self.segname,self.resseqnum2,self.atom2))
            W.addline('genbondrot {} $sel $i $j {}'.format(molid,self.degrees))
        elif self.angle=='LINK': # ASN-GLYcan rotation
            W.addline('set sel [atomselect {} "segname {} {}"]'.format(molid,self.segname1,self.segname2))
            W.addline('set i [[atomselect {} "segname {} and resid {} and name {}"] get index]'.format(molid,self.segname1,self.resseqnum1,self.atom1))
            W.addline('set j [[atomselect {} "segname {} and resid {} and name {}"] get index]'.format(molid,self.segname2,self.resseqnum2,self.atom2))
            W.addline('genbondrot {} $sel $i $j {}'.format(molid,self.degrees))
        elif self.angle=='ANGLEIJK':
            W.addline('set rotsel [atomselect {} "segname {}"]'.format(molid,self.segnamejk))
            W.addline('set ri [lindex [[atomselect {} "segname {} and resid {} and name {}"] get {{x y z}}] 0]'.format(molid,self.segnamei,self.resseqnumi,self.atomi))
            W.addline('set rj [lindex [[atomselect {} "segname {} and resid {} and name {}"] get {{x y z}}] 0]'.format(molid,self.segnamejk,self.resseqnumj,self.atomj))
            W.addline('set rk [lindex [[atomselect {} "segname {} and resid {} and name {}"] get {{x y z}}] 0]'.format(molid,self.segnamejk,self.resseqnumk,self.atomk))
            W.addline('set rij [vecsub $ri $rj]')
            W.addline('set rjk [vecsub $rj $rk]')
            W.addline('set cijk [veccross $rij $rjk]')
            W.addline('$rotsel move [trans center $rj origin $rj axis $cijk {} degrees]'.format(self.degrees))
        elif self.angle=='ALPHA':
            W.addline('set r1 [[atomselect {} "chain {} and resid {} and name CA"] get residue]'.format(molid,the_chainID,self.resseqnum1))
            W.addline('set r2 [[atomselect {} "chain {} and resid {} and name CA"] get residue]'.format(molid,the_chainID,self.resseqnum2))
            if endIsCterm:
                W.addline('fold_alpha $r1 $r2 {} {}'.format(the_chainID,molid))

class CrotList(AncestorAwareModList):
    def write_TcL(self,W:Psfgen,transform,**kwargs):
        for c in self:
            c.write_TcL(W,transform,**kwargs)

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
        logger.debug(f'parsed {self}')
    
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
            'sym2':'',
            'residue1':None,
            'residue2':None
        }
        super().__init__(input_dict)

    def __str__(self):
        c1=self.chainID1
        r1=self.resseqnum1
        i1=self.insertion1
        if self.residue1:
            c1=self.residue1.chainID
            r1=self.residue1.resseqnum
            i1=self.residue1.insertion
        c2=self.chainID2
        r2=self.resseqnum2
        i2=self.insertion2
        if self.residue2:
            c2=self.residue2.chainID
            r2=self.residue2.resseqnum
            i2=self.residue2.insertion
        # return f'{self.chainID1}_{self.resseqnum1}{self.insertion1}-{self.chainID2}_{self.resseqnum2}{self.insertion2}'
        return f'{c1}_{r1}{i1}-{c2}_{r2}{i2}'

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
        logger.debug(f'writing patch for {self}')
        c1=chainIDmap.get(self.chainID1,self.chainID1)
        c2=chainIDmap.get(self.chainID2,self.chainID2)
        r1=self.resseqnum1
        r2=self.resseqnum2
        W.addline(f'patch DISU {c1}:{r1} {c2}:{r2}')

class SSBondList(AncestorAwareModList):
    def assign_residues(self,Residues):
        ignored_by_ptnr1=self.assign_objs_to_attr('residue1',Residues,resseqnum='resseqnum1',chainID='chainID1',insertion='insertion1')
        ignored_by_ptnr2=self.assign_objs_to_attr('residue2',Residues,resseqnum='resseqnum2',chainID='chainID2',insertion='insertion2')
        return self.__class__(ignored_by_ptnr1+ignored_by_ptnr2)
    def write_TcL(self,W:Psfgen,transform):
        for s in self:
            s.write_TcL(W,transform)
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
    """A class for handling covalent bonds between residues where at least one residue is non-protein
    
    Attributes
    ----------

    req_attr: list
    * name1: name of partner-1 residue
    * chainID1: chain ID of partner-1 residue
    * resseqnum1: residue sequence number of partner-1 residue
    * insertion1: insertion code of partner-1 residue
    * name2: name of partner-2 residue
    * chainID2: chain ID of partner-2 residue
    * resseqnum2: residue sequence number of partner-2 residue
    * insertion2: insertion code of partner-2 residue
    
    """
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
        logger.debug(f'parsed {self}')
    
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

    def set_patchname(self):
        """Determine the charmff patch residue name for this link based on
        residue names, atom names, and when necessary, 3D geometry of a 
        particular dihedral angle around the link's bond
        
        New Attributes
        --------------
        patchname: str
           charmff pres name
        
        patchorder: list
            [1,2] if residue 1 appears first in the pres listing
            [2,1] if residue 2 appears first in the pres listing
        """
        self.patchname=''
        self.patchorder=[1,2]
        logger.debug(f'patch assignment for link {str(self)}')
        if not self.residue1 and not self.residue2:
            logger.debug(f'missing residue')
            logger.debug(f'1 {self.residue1}')
            logger.debug(f'2 {self.residue2}')
            return
        my_res12=[self.residue1,self.residue2]
        if self.resname1=='ASN' and self.segtype2=='glycan':
            # N-linked glycosylation site (toppar_all36_carb_glycopeptide)
            ICmap=[
                {'ICatomnames':['1CG','1ND2','2C1','2O5'],
                 'mapping':{'NGLA':168.99,'NGLB':-70.91}}
            ]
            self.patchname=ic_reference_closest(my_res12,ICmap)
        elif self.resname1=='SER' and self.segtype2=='glycan':
            # O-linked to serine (toppar_all36_carb_glycopeptide)
            ICmap=[
                {'ICatomnames':['1CB','1OG','2C1','2O5'],
                 'mapping':{'SGPA':45.37,'SGPB':19.87}}
            ]
            self.patchname=ic_reference_closest(my_res12,ICmap)
        elif self.resname1=='THR' and self.segtype2=='glycan':
            # O-linked to serine (toppar_all36_carb_glycopeptide)
            ICmap=[
                {'ICatomnames':['1CB','1OG1','2C1','2O5'],
                 'mapping':{'SGPA':69.9,'SGPB':33.16}}
            ]
            self.patchname=ic_reference_closest(my_res12,ICmap)
        elif self.name2=='C1' and self.segtype2=='glycan' and self.segtype1=='glycan':
            # all taken from 1xyz pres in top_all36_carb.rtf
            # including PHI angles of ICs with atoms in both residues
            if self.name1=='O1': # 1->1 link
                ICmap=[
                    {'ICatomnames':'1O5  1C1  1O1  2C1'.split(),
                     'mapping':{'11aa':103.46,'11ab':121.75,'11bb':-56.58}
                    },
                    {'ICatomnames':'1C1  1O1  2C1  2O5'.split(),
                     'mapping':{'11aa':103.54,'11ab':51.80,'11bb':-79.64}
                    },
                    {'atomnames':'1O1  2C1  2O5  2C5'.split(),
                     'mapping':{'11aa':64.56,'11ab':167.51,'11bb':172.18}}
                ]
                self.patchname=ic_reference_closest(my_res12,ICmap)
            elif self.name1=='O2': # 1->2 link
                ICmap=[
                    {'ICatomnames':'1C1  1C2  1O2  2C1'.split(),
                     'mapping':{'12aa':-132.81,'12ab':115.32,'12ba':-133.78,'12bb':117.14}
                    },
                    {'ICatomnames':'1C2  1O2  2C1  2O5'.split(),
                     'mapping':{'12aa':47.16,'12ab':86.93,'12ba':168.07,'12bb':-168.07}
                    }
                ]
                self.patchname=ic_reference_closest(my_res12,ICmap)
            elif self.name1=='O3': # 1->3 link
                ICmap=[
                    {'ICatomnames':'1C2  1C3  1O3  2C1'.split(),
                     'mapping':{'13aa':113.19,'13ab':-141.32,'13ba':-131.68,'13bb':-141.32}
                    },
                    {'ICatomnames':'1C3  1O3  2C1  2O5'.split(),
                     'mapping':{'13aa':65.46,'13ab':65.46,'13ba':-100.16,'13bb':-130.16}
                    }
                ]
                self.patchname=ic_reference_closest(my_res12,ICmap)
            elif self.name1=='O4': # 1->4 link
                ICmap=[
                    {'ICatomnames':'1C3  1C4  1O4  2C1'.split(),
                     'mapping':{'14aa':-86.29,'14ab':72.71,'14ba':-86.3,'14bb':81.86}},
                    {'ICatomnames':'1C4  1O4  2C1  2O5'.split(),
                     'mapping':{'14aa':133.57,'14ab':48.64,'14ba':-130.97,'14bb':-130.97}}
                ]
                self.patchname=ic_reference_closest(my_res12,ICmap)
            elif self.name1=='O6': # 1->6 link
                ICmap=[
                    {'ICatomnames':'1C6  1O6  2C1  2O5'.split(),
                     'mapping':{'16AT':71.24,'16BT':-63.49}}
                ]
                self.patchname=ic_reference_closest(my_res12,ICmap)
        elif self.name2=='C2' and self.segtype2=='glycan' and self.segtype1=='glycan':
            if self.name1=='O6':
                self.patchname='SA26AT'
            elif self.name1=='O8':
                self.patchname='SA28AA'
            elif self.name1=='O9':
                self.patchname='SA29AT'

        elif self.name1=='O6' and self.name2=='C2':
                self.patchname='SA26AT' 
        elif 'ZN' in self.resname1 and 'HIS' in self.resname2:
                self.patchname='ZNHD'
                self.patchorder=[2,1]
        else:
            logger.warning(f'Could not identify patch for link: {str(self)}')
            self.patchname='UNFOUND'

    def write_TcL(self,W:Psfgen,transform):
        """Insert the appropriate TcL commands to add this link in a psfgen script
        
        Assumes that if one of the two residues is an asparagine and the other is a 
        from a glycan, this requires an CHARMM NGLB patch

        A 2->6 intraglycan linkage requies the SA26T patch

        Others' patches are named 1xij where x is the carbon number of the upstream monomer,
        and i and j indicate whether the bond is axial or equatorial.  This is currently
        determined using the custom 'axeq' TcL procedure.

        Parameters
        ----------
        W: Psfgen
            the psfgen scriptwriter object
        transform: BiomT
            the designated transform under which this link is operational; used for its chainIDmap
        """
        chainIDmap=transform.chainIDmap
        seg1=self.residue1.chainID
        seg1=chainIDmap.get(seg1,seg1)
        seg2=self.residue2.chainID
        seg2=chainIDmap.get(seg2,seg2)
        if not self.patchname=='UNFOUND':
            if self.patchorder==[1,2]:
                W.addline(f'patch {self.patchname} {seg1}:{self.resseqnum1}{self.insertion1} {seg2}:{self.resseqnum2}{self.insertion2}')
            elif self.patchorder==[2,1]:
                W.addline(f'patch {self.patchname} {seg2}:{self.resseqnum2}{self.insertion2} {seg1}:{self.resseqnum1}{self.insertion1}')
        else:
            logger.warning(f'Could not identify patch for link: {str(self)}')
            W.comment(f'No patch found for {str(self)}')
    
    def __str__(self):
        return f'{self.chainID1}_{self.resname1}{self.resseqnum1}{self.insertion1}-{self.chainID2}_{self.resname2}{self.resseqnum2}{self.insertion2}'

class LinkList(AncestorAwareModList):
    """A class for handling lists of Links

    Methods
    -------
    assign_resiudes(Residues)
        scans the provided list of residues to assign actual residues to the 'residue1' and 'residue2'
        attributes of each link; sets the 'atom1' and 'atom2' attributes to point to the actual
        atoms; sets up the 'up' and 'down' pointers for every link; sets the segtype1 and segtype2
        attributes of every link; returns list of residues not assigned to links and links to which
        no residues were assigned.

    write_TcL(W,transform)
        calls write_TcL for each link

    prune_mutations(Mutations,Segments)
        Given a list of mutations, removes any links that were declared such that any mutation
        would make the link chemically impossible.  E.g., mutating an ASN of an N-linked
        glycosylation site results in loss of the ASN-BGLNA link.
    
    apply_segtypes(map)
        using the resname-to-segtype map provided, set the 'segtype1' and 'seqtype2' attributes
        of every element.
                
    """
    def assign_residues(self,Residues):
        """Assigns residue and atom pointers to each link; sets up the up and down links of both
        residues so that linked residue objects can reference one another; flags residues from
        list of residues passed in that are not assigned to any links

        Arguments
        ---------
        Residues: ResidueList
            List of residues to search in order to make residue assignments
        
        Returns
        -------
        ResidueList: list of residues from Residues that are not used for any assignments

        """
        ignored_by_ptnr1=self.assign_objs_to_attr('residue1',Residues,resseqnum='resseqnum1',chainID='chainID1',insertion='insertion1')
        ignored_by_ptnr2=self.assign_objs_to_attr('residue2',Residues,resseqnum='resseqnum2',chainID='chainID2',insertion='insertion2')
        for link in self:
            link.residue1.linkTo(link.residue2,link)
            link.atom1=link.residue1.atoms.get(name=link.name1,altloc=link.altloc1)
            link.atom2=link.residue2.atoms.get(name=link.name2,altloc=link.altloc2)
            link.segtype1=link.residue1.segtype
            link.segtype2=link.residue2.segtype
            # TODO: determine patch residue
            link.set_patchname()
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

    # def remove_orphan_residues(self,Links,Residues):
    #     for dl in self:
    #         reslist,lnlist=dl.residue2.get_down_group()
    #         reslist.insert(0,dl.residue2)
    #         for r in reslist:
    #             Residues.remove(r)
    #         for l in lnlist:
    #             Links.remove(l)
    #     return reslist,lnlist

    def prune_mutations(self,Mutations,Segments):
        """Prune off any links and associated objects as a result of mutations
        
        Arguments
        ---------
        Mutations: MutationList
            list of mutations
            
        Segments: SegmentList
            list of assembled segments; might need to be modified if pruning gets rid
            of a whole segment's worth of residues
        
        """
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
                    Segments.remove(S)
                    pruned['segments'].append(S)
                for l in llist:
                    pruned['links'].append(self.remove(l))
        return pruned

    def apply_segtypes(self,map):
        """Apply segtype values to each of the two residues using the map
        
        Parameters
        ----------
        map: dict
            map of segtypes for given resnames
        """
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