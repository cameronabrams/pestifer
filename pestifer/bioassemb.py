

import numpy as np
from pidibble.pdbparse import get_symm_ops
import logging
logger=logging.getLogger(__name__)

class BioAssemb:
    ''' Container for handling info for "REMARK 350 BIOMOLECULE: #" stanzas in RCSB PDB files
        or _pdbx_struct blocks in mmCIF files '''
    def __init__(self):
        self.chains_depot=[]
        self.biomt=[]
        self.pdbx_struct={}
    
    @classmethod
    def from_asymmetricUnit(cls,parent_molecule=''):
        inst=cls()
        inst.name='Asymmetric Unit'
        inst.parent_molecule=parent_molecule
        inst.biomt=[BiomT()]
        return inst
    
    @classmethod
    def from_PDBRecord(cls,pdbrecord,parent_molecule=''):
        inst=cls()
        inst.parent_molecule=parent_molecule
        inst.parsePDBRecord(pdbrecord)
        return inst
            
    def parsePDBRecord(self,pdbrecord):
        rs=pdbrecord.key.split('.')
        self.name='.'.join(rs[1:3])
        self.chainIDs=pdbrecord.header
        M,T=get_symm_ops(pdbrecord)
        self.biomt=[]
        for i,(m,t) in enumerate(zip(M,T)):
            b=BiomT().parse_rot_trans(m,t)
            # b.apply_to_chainIDs=self.apply_to_chainIDs
            b.index=i
            self.biomt.append(b)

        self.pdbx_struct['author_biol_unit']=pdbrecord.get_token('AUTHOR DETERMINED BIOLOGICAL UNIT')
        self.pdbx_struct['software_biol_unit']=pdbrecord.get_token('SOFTWARE DETERMINED QUATERNARY STRUCTURE')
        self.pdbx_struct['software_used']=pdbrecord.get_token('SOFTWARE USED')
        tok=pdbrecord.get_token('TOTAL BURIED SURFACE AREA')
        if tok:
            self.pdbx_struct['total_buried_surface_area'],self.pdbx_struct['surface_area_units']=tok.split()
        self.pdbx_struct['surface_area_complex']=pdbrecord.get_token('SURFACE AREA OF THE COMPLEX')
        tok=pdbrecord.get_token('CHANGE IN SOLVENT FREE ENERGY:')
        if tok:
            self.pdbx_struct['change_solv_fe'],self.pdbx_struct['fe_units']=tok.split()

    # def show(self,isActive=False,indent='    '):
    #     activeLabel='**ACTIVE**' if isActive else ''
    #     print('{}{} {} {:d} {}'.format(indent,activeLabel,self.name,self.index,activeLabel))
    #     if len(self.pdbx_struct)>0:
    #         print(indent*2,self.pdbx_struct)
    #     for b in self.biomt:
    #         b.show(indent*2)
    
    # def CIFBiomT(self,cifdict):
    #     self.biomt.append(BiomT(index=len(self.biomt)))
    #     self.biomt[-1].CIFBiomT(cifdict)

#     def inheritConstructs(self,au,cid_depot):
#         for b in self.biomt:
#             aumd=au.biomt[0].md
#             ''' clone '''
# #            for bcid in self.apply_to_chainIDs:
# #               if bcid not in au.apply_to_chainIDs:
# #                    print(f'PDB indicates Biomolecule {b.index} applies its transformations to chain {bcid} which is not in the AU')
#             self.apply_to_chainIDs_as_read=self.apply_to_chainIDs[:]
#             self.apply_to_chainIDs=au.apply_to_chainIDs[:]
#             b.apply_to_chainIDs=self.apply_to_chainIDs
#             b.apply_to_chainIDs_as_read=self.apply_to_chainIDs_as_read
#             if b.isidentity():
#                 ''' this biomolecular assembly has one biomt that is the 
#                     asymmetric unit '''
#                 b.md=au.biomt[0].md
#             else:
#                 for aucid in au.apply_to_chainIDs:
#                     try:
#                         localcid=cid_depot.pop(0)
#                     except:
#                         print('Error: not enough alphabetical chain IDs!')
#                         exit(-1)
#                     b.mapChainIDs(aucid,localcid)
                
#                 b.md=aumd.Clone(chainmap=b.replicachainID_from_sourcechainID,invchainmap=b.sourcechainID_from_replicachainID)
#                 #print(f'inheritConstructs: {len(b.md.Mutations)} mutations.')
#         return cid_depot                    

class BiomT:
    def __init__(self,index=0):
        # self.md=None
        # self.apply_to_chainIDs=[]
        # self.apply_to_chainIDs_as_read=[]
#        self.ownChainIDs=[]
        self.index=index
        self.tmat=np.array([[1, 0, 0, 0],[0, 1, 0, 0],[0, 0, 1, 0]],dtype=float)
        # self.replicachainID_from_sourcechainID={}
        # self.sourcechainID_from_replicachainID={}
    def parse_rot_trans(self,M,T):
        for i in range(3):
            for j in range(3):
                self.tmat[i][j]=M[i][j]
            self.tmat[i][3]=T[i]
        return self
    # def parseBIOMT(self,ax,words):
    #     if self.index==-1:
    #         self.index=int(words[3])
    #     vals=[]
    #     for w in words[4:]:
    #        vals.append(float(w))
    #     self.tmat[ax-1]=vals
    # def CIFBiomT(self,cifdict):
    #     self.index=int(cifdict['id'])
    #     for i in range(3):
    #         for j in range(3):
    #             self.tmat[i][j]=float(cifdict['matrix[{}][{}]'.format(i+1,j+1)])
    #         self.tmat[i][3]=float(cifdict['vector[{}]'.format(i+1)])

    # def show(self,indent='    '):
    #     print('{}BIOMT {:d} operates on chains {:s}'.format(indent,self.index,', '.join(self.apply_to_chainIDs)))
    #     if len(self.apply_to_chainIDs_as_read)>len(self.apply_to_chainIDs):
    #         print('{}{}(Note: This may not reflect the chains listed for this biomolecule in the PDB.'.format(indent,indent))
    #         print('{}{}{}Original chains are {:s}'.format(indent,indent,indent,', '.join(self.apply_to_chainIDs_as_read)))
    #         print('{}{}One or more of these chains may have been deleted after processing mutations/deletions.)'.format(indent,indent))
    #     if not self.isidentity():
    #         print('{}    TMAT'.format(indent),self.tmat)
    #         if len(self.replicachainID_from_sourcechainID)>0:
    #             print('{}    REPC'.format(indent),self.replicachainID_from_sourcechainID)
    #     else:
    #         print('{}    IDENTITY'.format(indent))
    # def isidentity(self):
    #     t=self.tmat
    #     if t[0][0]==1.0 and t[1][1]==1.0 and t[2][2]==1.0:
    #         return True
    #     else:
    #         return False 
    # def mapChainIDs(self,c,newc):
    #     self.replicachainID_from_sourcechainID[c]=newc
    #     self.sourcechainID_from_replicachainID[newc]=c
    # def get_replica_chainID(self,c):
    #     if c in self.replicachainID_from_sourcechainID:
    #        return self.replicachainID_from_sourcechainID[c]
    #     else:
    #         return c
    # def report_chain_replicas(self):
    #     for k,v in self.replicachainID_from_sourcechainID.items():
    #         print('#### {} -> {}'.format(k,v))
    # def get_base_chainID(self,newc):
    #     if newc in self.sourcechainID_from_replicachainID:
    #         return self.sourcechainID_from_replicachainID[newc]
    #     else:
    #         return newc
    def OneLiner(self):
        retstr=r'{ '
        for i in range(3):
            retstr+=r'{ '
            for j in range(4):
               retstr+='{} '.format(self.tmat[i][j])
            retstr+=r' } '
        retstr+='{ 0 0 0 1 } }'
        return retstr

