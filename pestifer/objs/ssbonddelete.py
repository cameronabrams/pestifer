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

    def describe(self):
        return f"SSBondDelete(chainID1={self.chainID1}, resseqnum1={self.resseqnum1}, insertion1={self.insertion1}, chainID2={self.chainID2}, resseqnum2={self.resseqnum2}, insertion2={self.insertion2})"

    _yaml_header: ClassVar[str] = 'ssbondsdelete'
    """
    YAML header for SSBondDelete objects.
    This header is used to identify SSBondDelete objects in YAML files.
    """

    _objcat: ClassVar[str] = 'topol'
    """
    Category of the SSBondDelete object.
    This categorization is used to group SSBondDelete objects in the object manager.
    """

class SSBondDeleteList(SSBondList):
    """
    A class for handling a list of deleted SSBonds.
    """
    def is_deleted(self,a_SSBond):
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