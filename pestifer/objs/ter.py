# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A class for handling TER records in PDB files.     
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
"""

import logging
logger=logging.getLogger(__name__)
from functools import singledispatchmethod
from pidibble.pdbrecord import PDBRecord

from ..core.baseobj import AncestorAwareObj, AncestorAwareObjList

class Ter(AncestorAwareObj):
    """
    A class for handing TER records in PDB files
    """

    req_attr=['serial','resname','chainID','resseqnum','insertion']
    """
    Required attributes for a Ter object.
    These attributes must be provided when creating a Ter object.
    
    - ``serial``: The serial number of the TER record.
    - ``resname``: The residue name of the TER record.
    - ``chainID``: The chain ID of the TER record.
    - ``resseqnum``: The residue sequence number of the TER record.
    - ``insertion``: The insertion code for the TER record.
    """

    yaml_header='terminals'
    """ 
    YAML header for Ter objects.
    This header is used to identify Ter objects in YAML files.
    """
    
    objcat='seq'
    """
    Category of the Ter object.
    This categorization is used to group Ter objects in the object manager.
    """
    
    PDB_keyword='TER'
    """
    The PDB keyword for TER records.
    This keyword is used to identify TER records in PDB files.
    """
    
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
