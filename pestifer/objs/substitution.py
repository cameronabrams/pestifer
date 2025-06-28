# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import logging
logger=logging.getLogger(__name__)
from functools import singledispatchmethod

from ..core.baseobj import AncestorAwareObj, AncestorAwareObjList
from ..core.stringthings import split_ri

class Substitution(AncestorAwareObj):
    """A class for handling substitutions 
    
    A substitution is a user-specified modification in which a sequence of one or
    more contiguous residues are replaced by a new sequence of one or more residues.

    Attributes
    ----------
    req_attr : list
        * chainID : str
            chain ID of the segment to which substitution is made
        * resseqnum1 : int
            N-terminal resid of sequence to be replaced
        * insertion1 : str
            insertion code of N-terminal residue
        * resseqnum2 : int
            C-terminal resid of sequence to be replaced
        * insertion2 : str
            insertion code of the C-terminal residue
        * subseq : str
            amino acid sequence to substitute expressed using one-letter codes
    This class represents a substitution modification in a protein structure, allowing for the
    specification of a range of residues to be replaced and the new sequence to be inserted.
    The substitution is defined by the chain ID, the N-terminal and C-terminal residue numbers
    (including insertion codes), and the sequence of residues to be substituted in.
    The class provides methods to initialize from a shortcode, which is a compact representation
    of the substitution, and to handle the attributes required for the substitution.
    The shortcode format is `C:nnn-ccc,abcdef`, where `C` is the chain ID, `nnn` is the N-terminal residue number,
    `ccc` is the C-terminal residue number, and `abcdef` is the sequence of residues to be substituted.
    The `subseq` attribute is a string of one-letter amino acid codes representing the new sequence to be inserted.
    The class inherits from `AncestorAwareObj`, which provides basic functionality for handling attributes
    and initialization from various input types.
    The `yaml_header` and `objcat` attributes are used for YAML serialization and categorization of the object.
    The `req_attr` attribute lists the required attributes for the substitution object.
    The class is designed to be flexible and can be initialized from a string shortcode or a dictionary
    representation of the substitution.

    """
    req_attr=AncestorAwareObj.req_attr+['chainID','resseqnum1','insertion1','resseqnum2','insertion2','subseq']
    yaml_header='substitutions'
    objcat='seq'

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

class SubstitutionList(AncestorAwareObjList):
    pass