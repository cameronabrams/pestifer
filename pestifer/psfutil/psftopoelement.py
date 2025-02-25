# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" Defines the PSFTopoElement class and many descendants
    Defines the PSFContent class for parsing PSF files
"""
import numpy as np
import logging

from argparse import Namespace
from collections import UserList
from copy import deepcopy
from itertools import product

from ..util.coord import mic_shift
from ..util.util import countTime

logger=logging.getLogger(__name__)

class LineList(UserList):
    def __init__(self,data):
        super().__init__(data)

class PSFTopoElement(Namespace):
    def __init__(self,idx_list):
        a_dict={}
        for i,idx in enumerate(idx_list):
            a_dict.update({f'serial{i+1}':idx})
        super().__init__(**a_dict)
        self.idx_list=idx_list.copy()

    def __len__(self):
        return len(self.idx_list)
    
    def __str__(self):
        return '-'.join([str(x) for x in self.idx_list])

    def calculate_stuff(self):
        self.COM=np.mean(self.P,axis=0)
        self.Radii=self.P-self.COM
        self.radii=np.array([n/d for n,d in zip(self.Radii,np.linalg.norm(self.Radii,axis=1))])
        self.B=np.array([p[1]-p[0] for p in product(self.P,repeat=2)]).reshape((len(self.P),len(self.P),3))
        self.b=np.linalg.norm(self.B.reshape((len(self.P)**2,3)),axis=1).reshape(len(self.P),len(self.P))

    def injest_coordinates(self,A,pos_key=['posX','posY','posZ'],meta_key=[]):
        self.P=np.array(A.loc[self.idx_list][pos_key].values)
        if np.any(np.isnan(self.P)):
            logger.debug(f'nan detected in {self.P}')
        if len(self.P)>1:
            self.calculate_stuff()
        for k in meta_key:
            self.__dict__[k]=A.loc[self.idx_list[0],k]

    def assign_cell_index(self,LC):
        self.linkcell_idx=LC.ldx_of_cellndx(LC.cellndx_of_point(self.COM))

    def validate_image(self,box):
        boxlen=np.array([box[i][i] for i in [0,1,2]])
        in_same_img=np.linalg.norm(self.P-self.COM)<(boxlen/2)
        self.same_image=np.all(in_same_img)
        return self.same_image
    
    def mic_shift(self,reference_point,box):
        holder=self.copy()
        for i in range(len(self.P)):
            holder.P[i],bl=mic_shift(self.P[i],reference_point,box)
        if np.any(bl): holder.calculate_stuff()
        return holder
    
    def copy(self):
        return deepcopy(self)

    def match(self,**crits):
        return all([self.__dict__[k]==v for k,v in crits.items()])

class PSFTopoElementList(UserList):
    def __init__(self,data):
        super().__init__(data)

    @countTime
    def injest_coordinates(self,A,pos_key=['posX','posY','posZ'],meta_key=[]):
        for b in self:
            b.injest_coordinates(A,pos_key=pos_key,meta_key=meta_key)

    @countTime
    def assign_cell_indexes(self,LC):
        for b in self:
            b.assign_cell_index(LC)

    def validate_images(self,box):
        return all([b.validate_image(box) for b in self])

    def get(self,**crits):
        result=[]
        for e in self:
            if e.match(**crits):
                result.append(e)
        return type(self)(result)
    



    
