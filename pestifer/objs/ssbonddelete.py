# Author: Cameron F. Abrams
from .ssbond import SSBond, SSBondList

class SSBondDelete(SSBond):
    yaml_header='ssbondsdelete'
    objcat='topol'

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