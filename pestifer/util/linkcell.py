# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" Custom link-cell algorithm for pair-wise searches in 3D space

The Linkcell object is used by the RingCheck algorithm to detect pierced rings

"""
import numpy as np
import logging
from itertools import product
logger=logging.getLogger(__name__)

class Linkcell:
    """ Handles the link-cell algorithm """
    def __init__(self,corners=np.array([[0,0,0],[100,100,100]],dtype=float),cutoff=10.0,atidxlabel='serial',atposlabel=['x','y','z']):
        """__init__ Constructor for an empty Linkcell object

        :param corners: lower-left and upper-right corner coordinates for the box
        :type corners: np.array
        :param cutoff: cutoff distance, defaults to None
        :type cutoff: float, optional
        :param wrap_point: function used on any 3-D point to wrap it into the central image
        :type cutoff: bool
        """

        self.cutoff=cutoff
        self.lower_left_corner,self.upper_right_corner=corners
        self.sidelengths=self.upper_right_corner-self.lower_left_corner
        self.atidxlabel=atidxlabel
        self.atposlabel=atposlabel
        # number of cells along x, y, and z directions
        self.cells_per_dim=np.floor(self.sidelengths/self.cutoff).astype(int)
        self.ncells=np.prod(self.cells_per_dim)
        # dimensions of one cell
        self.celldim=self.sidelengths/self.cells_per_dim
        logger.debug(f'Linkcell structure: {self.ncells} cells ({self.cells_per_dim}) size {self.celldim} A^3')

    def wrap_point(self,ri):
        """wrap_point wraps point ri into the base periodic image

        :param ri: a point
        :type ri: np.ndarray(3,float)
        :return: a tuple containing (1) the wrapped point and (2) number of box lengths required to wrap this point, per dimension
        :rtype: tuple
        """
        R=ri.copy()
        box_lengths=np.array([0,0,0],dtype=int)
        for i in range(3):
            while R[i]<self.lower_left_corner[i]:
                R[i]+=self.sidelengths[i]
                box_lengths[i]+=1
            while R[i]>=self.upper_right_corner[i]:
                R[i]-=self.sidelengths[i]
                box_lengths[i]-=1
        if not np.all(R==R):
            logger.debug(f'{ri} {self.sidelengths} -> {R}')
            exit(-1)
        return R,box_lengths

    def cellndx_of_point(self,R):
        """cellndx_of_point returns the (i,j,k) cell index of point R
        """
        wrapR,bl=self.wrap_point(R)
        # logger.debug(f'{R} {wrapR}')
        wrapR-=self.lower_left_corner
        C=np.floor(wrapR*np.reciprocal(self.celldim)).astype(int)
        lowdim=(C<np.zeros(3).astype(int)).astype(int) # will never happen if R is wrapped
        hidim=(C>=self.cells_per_dim).astype(int) # could happen if exactly there
        if (any(lowdim) or any(hidim)):
            logger.warning(f'Warning: point {R} maps to out-of-bounds-cell {C} ({self.cells_per_dim})')
            logger.warning(f'lower-left: {self.lower_left_corner}')
            logger.warning(f'upper-right: {self.upper_right_corner}')
        C+=lowdim
        C-=hidim
        return C

    def ldx_of_cellndx(self,C):
        """ldx_of_cellndx returns scalar index of cell with (i,j,k)-index C
        """
        nc=self.cells_per_dim
        # logger.debug(f'{C} {nc}')
        xc=C[2]*nc[0]*nc[1]+C[1]*nc[1]+C[0]
        return xc

    def cellndx_of_ldx(self,i):
        """cellndx_of_ldx returns (i,j,k)-index of cell with scalar index i
        """
        nc=self.cells_per_dim
        t1=i//(nc[0]*nc[1])
        r1=i-t1*(nc[0]*nc[1])
        t2=r1//nc[1]
        t3=r1-t2*nc[1]
        return np.array([t3,t2,t1])

    def ldx_searchlist_of_ldx(self,i):
        """searchlist_of_ldx returns the list of cells to search around and including central cell i """
        # assert i!=-1
        C=self.cellndx_of_ldx(i)
        nc=self.cells_per_dim
        d=np.array([-1,0,1])
        searchlist=[self.ldx_of_cellndx(np.mod(C+p,nc)) for p in product(d,d,d)]
        assert i in searchlist
        return searchlist
