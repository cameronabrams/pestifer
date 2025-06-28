# Author: Cameron F. Abrams
from .ssbond import SSBond, SSBondList

class SSBondDelete(SSBond):
    """A class for handling deletion of SSBonds."""
    yaml_header='ssbondsdelete'
    objcat='topol'

class SSBondDeleteList(SSBondList):
    """A class for handling a list of deleted SSBonds."""
    def is_deleted(self,a_SSBond):
        """Check if a given SSBond is deleted in this list.
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