# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import logging
logger=logging.getLogger(__name__)
from functools import singledispatchmethod
from pidibble.pdbrecord import PDBRecord

from ..baseobj import AncestorAwareObj, AncestorAwareObjList

class Ter(AncestorAwareObj):
    """A class for handing TER records in PDB files
    
    Why is this necessary?  TER records are used to indicate
    'breaks' in the list of ATOM records in a PDB file to signify
    physical discontinuity of the amino acid chain(s).  The problem
    we have to solve is that sometimes a TER record gets a unique 
    atom serial number, even though they are not atoms. So in a list 
    of atoms the sequence of serial numbers can be discontinuous.
    This is not a problem per se, *but* when VMD reads in a PDB file 
    it ignores both TER records *and* the explicit atom serial numbers.
    Pidibble by default does not, so to make the atom records 
    VMD-compliant, we have to adjust atom serial numbers, and to do 
    that we have to know where any non-empty TER records are.
    
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
    
    objcat: (str)
        type of mod from among objcats

    """
    req_attr=['serial','resname','chainID','resseqnum','insertion']
    yaml_header='terminals'
    objcat='seq'
    PDB_keyword='TER'

    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)
    @__init__.register(PDBRecord)
    def _from_pdbrecord(self,pdbrecord):
        if not pdbrecord.serial:
            # this TER is empty
            input_dict={
                'serial':None,
                'resname':None,
                'chainID':None,
                'resseqnum':None,
                'insertion':None
            }
        else:
            input_dict={
                'serial':pdbrecord.serial,
                'resname':pdbrecord.residue.resName,
                'chainID':pdbrecord.residue.chainID,
                'resseqnum':pdbrecord.residue.seqNum,
                'insertion':pdbrecord.residue.iCode
            }
        super().__init__(input_dict)
    
class TerList(AncestorAwareObjList):
    pass
