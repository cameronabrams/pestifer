# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import logging
logger=logging.getLogger(__name__)
from functools import singledispatchmethod

from .seqadv import Seqadv

from ..baseobj import AncestorAwareObj, AncestorAwareObjList
from ..config import res_123
from ..stringthings import split_ri

class Mutation(AncestorAwareObj):
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
    req_attr=AncestorAwareObj.req_attr+['chainID','origresname','resseqnum','insertion','newresname','typekey']
    opt_attr=AncestorAwareObj.opt_attr+['pdbx_auth_seq_num']
    yaml_header='mutations'
    objcat='seq'

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
        if hasattr(sq,'label_asym_id'):
            input_dict['chainID']=sq.__dict__['label_aym_id']
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

class MutationList(AncestorAwareObjList):
    pass
