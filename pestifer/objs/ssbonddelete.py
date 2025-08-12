# Author: Cameron F. Abrams

"""
One modification a user might want to make is reducing a disulfide bond
to free cysteine residues.  This is handled by the SSBondDelete class.
"""

from .ssbond import SSBond, SSBondList
from typing import ClassVar

class SSBondDelete(SSBond):
    """
    A class for handling deletion of SSBonds.
    """

    _yaml_header: ClassVar[str] = 'ssbondsdelete'
    """
    YAML header for SSBondDelete objects.
    This header is used to identify SSBondDelete objects in YAML files.
    """

class SSBondDeleteList(SSBondList):
    """
    A class for handling a list of deleted SSBonds.
    """
    def is_deleted(self, a_SSBond: SSBond) -> bool:
        """
        Check if a given SSBond is deleted in this list.
        
        Parameters
        ----------
        a_SSBond : SSBond
            The SSBond to check for deletion.
        
        Returns
        -------
        bool
            True if the SSBond is deleted, False otherwise.
        """
        if self.get(lambda x: 
            x.chainID1 == a_SSBond.chainID1 and 
            x.chainID2 == a_SSBond.chainID2 and
            x.resid1 == a_SSBond.resid1 and
            x.resid2 == a_SSBond.resid2):
            return True
        if self.get(lambda x: 
            x.chainID2 == a_SSBond.chainID1 and
            x.chainID1 == a_SSBond.chainID2 and
            x.resid2 == a_SSBond.resid1 and
            x.resid1 == a_SSBond.resid2):
            return True
        return False