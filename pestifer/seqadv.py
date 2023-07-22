"""

.. module:: seqdv
   :synopsis: Manages SEQADV records in PDB files
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
class Seqadv:
    def __init__(self,pdbrecord):
        if len(pdbrecord)>0:
            self.pdbrecord=pdbrecord
            self.record_name=pdbrecord[0:6]
            self.idCode=pdbrecord[7:11]
            self.rawresName=pdbrecord[12:15]
            self.resName=self.rawresName.strip()
            self.chainID=pdbrecord[16:17]
            self.seqNum=int(pdbrecord[18:22])
            self.iCode=pdbrecord[22:23]
            self.database=pdbrecord[24:28]
            self.dbAccession=pdbrecord[29:38]
            self.dbRes=pdbrecord[39:42].strip()
            self.dbSeqraw=pdbrecord[43:48]
            if self.dbSeqraw.strip().isdigit()==True:
                self.dbSeq=int(self.dbSeqraw)
            else:
                self.dbSeq=''
            self.conflict=pdbrecord[49:70].strip()
    def pdb_line(self):
        return f'{self.record_name:6s} {self.idCode:3s} {self.rawresName:>3s} {self.chainID:1s} {self.seqNum:>4d}{self.iCode:1s} {self.database:>4s} {self.dbAccession:9s} {self.dbRes:3s} {self.dbSeqraw:5s} {self.conflict:21s}          '
    def clone(self,chain):
        newSeqadv=Seqadv(pdbrecord=self.pdb_line())
        newSeqadv.chainID=chain
        newSeqadv.pdbrecord=newSeqadv.pdb_line()
        return newSeqadv
    def __eq__(self,other):
        result=True
        result&=self.pdbrecord==other.pdbrecord
        result&=self.record_name==other.record_name
        result&=self.idCode==other.idCode
        result&=self.rawresName==other.rawresName
        result&=self.resName==other.resName
        result&=self.chainID==other.chainID
        result&=self.seqNum==other.seqNum
        result&=self.iCode==other.iCode
        result&=self.database==other.database
        result&=self.dbAccession==other.dbAccession
        result&=self.dbRes==other.dbRes
        result&=self.dbSeqraw==other.dbSeqraw
        result&=self.dbSeq==other.dbSeq
        result&=self.conflict==other.conflict
        return result

    def printshort(self):
        return '{}-{}{}{}'.format(self.chainID,self.resName,self.seqNum,self.dbRes)
    def __str__(self):
        retstr='{}\n'+\
               ' idCode      {}\n'+\
               ' resName     {}\n'+\
               ' chainID     {}\n'+\
               ' seqNum      {}\n'+\
               ' iCode       {}\n'+\
               ' database    {}\n'+\
               ' dbAccession {}\n'+\
               ' dbRes       {}\n'+\
               ' dbSeq       {}\n'+\
               ' conflict    {}\n'
        return retstr.format(self.record_name,
                             self.idCode,
                             self.resName,
                             self.chainID,
                             self.seqNum,
                             self.iCode,
                             self.database,
                             self.dbAccession,
                             self.dbRes,
                             self.dbSeq,
                             self.conflict)

