# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import numpy as np
import pandas as pd
import logging
from itertools import product
from multiprocessing import Pool
from functools import partial
from .stringthings import my_logger
from .util import countTime
import cProfile, pstats, io
from pstats import SortKey
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
        """

        self.cutoff=cutoff
        # self.box=box # box cell dimensions
        # self.origin=origin # by NAMD convention, origin is the CENTER of the box
        self.lower_left_corner,self.upper_right_corner=corners
        self.sidelengths=self.upper_right_corner-self.lower_left_corner
        self.atidxlabel=atidxlabel
        self.atposlabel=atposlabel
        # number of cells along x, y, and z directions
        self.cells_per_dim=np.floor(self.sidelengths/self.cutoff).astype(int)
        self.ncells=np.prod(self.cells_per_dim)
        # logger.debug(f'{self.ncells}')
        # dimensions of one cell
        self.celldim=self.sidelengths/self.cells_per_dim
        # 3-d array of lower left corner as a 3-space point, indexed by i,j,k
        # initialized to all zeros, calculated below
        # self.cells=np.zeros((*self.ncells,3))
        # # 1-d array of (i,j,k) indices indexed by linear cell index (0...ncells-1)
        # self.cellndx=np.array(list(product(*[np.arange(x) for x in self.ncells])))
        # # logger.debug(f'{self.cellndx[0]} to {self.cellndx[-1]}')
        # # logger.debug(f'{self.ldx_of_cellndx(self.cellndx[0])} to {self.ldx_of_cellndx(self.cellndx[-1])}')
        # # logger.debug(f'{self.celldim} {len(self.cells)} {len(self.cellndx)}')
        # # 3-d array of lower left corner as a 3-space point, indexed by i,j,k
        # for t in self.cellndx:
        #     i,j,k=t
        #     self.cells[i,j,k]=self.celldim*np.array([i,j,k])+self.lower_left_corner
        # # logger.debug(f'{self.cellndx[-1]}')
        # # set up neighbor lists using linear indices
        # self.make_neighborlists()
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

    # @countTime
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

    # def point_in_cellndx(self,R,C):
    #     """point_in_cellndx returns True if point R is located in cell with (i,j,k) index C

    #     :param R: 3-space point
    #     :type R: numpy.ndarray
    #     :param C: (i,j,k) cell index
    #     :type C: (int,int,int)
    #     :return: 
    #         - True: R in C
    #         - False: R not in C
    #     :rtype: boolean
    #     """
    #     LL,UU=self.corners_of_cellndx(C)
    #     wrapR,bl=self.wrap_point(R)
    #     return all(wrapR<UU) and all(wrapR>=LL)

    # def corners_of_cellndx(self,C):
    #     """corners_of_cellndx returns the lower-left and upper-right corners of cell with (i,j,k) index C, as an array of 3-space points

    #     :param C: (i,j,k) index
    #     :type C: (int,int,int)
    #     :return: 2x3 array of lower-left and upper-right corner coordinates
    #     :rtype: numpy.ndarray
    #     """
    #     LL=self.cells[C[0],C[1],C[2]]
    #     UU=LL+self.celldim
    #     return np.array([LL,UU])

    # def cellndx_in_structure(self,C):
    #     """cellndx_in_structure test to see if (i,j,k) index given is within the established linkcell structure

    #     :param C: (i,j,k) index
    #     :type C: (int,int,int)
    #     :return:  - True: C is in the linkcell structure
    #               - False: C is not in the linkcell structure
    #     :rtype: boolean
    #     """
    #     return all(np.zeros(3)<=C) and all(C<self.ncells)

    def ldx_of_cellndx(self,C):
        """ldx_of_cellndx returns scalar index of cell with (i,j,k)-index C
        """
        nc=self.cells_per_dim
        # logger.debug(f'{C} {nc}')
        xc=C[2]*nc[0]*nc[1]+C[1]*nc[1]+C[0]
        # xc=C[0]*nc[1]*nc[2]+C[1]*nc[1]+C[2]
        return xc

    def cellndx_of_ldx(self,i):
        """cellndx_of_ldx returns (i,j,k)-index of cell with scalar index i
        """
        nc=self.cells_per_dim
        t1=i//(nc[0]*nc[1])
        r1=i-t1*(nc[0]*nc[1])
        t2=r1//nc[1]
        t3=r1-t2*nc[1]
        # r=t1+t2+t3
        # assert(r==self.cellndx[i])
        return np.array([t3,t2,t1]) # self.cellndx[i]

    @countTime
    def populate_par(self,adf):
        """populate_par populate the linkcell structure by setting the "linkcell_idx" attribute of each atom in the coordinates dataframe adf

        :param adf: gromacs-format coordinate dataframe
        :type adf: pandas.DataFrame
        :return: modified dataframe
        :rtype: pandas.DataFrame
        """
        sp=adf[self.atposlabel]
        # logger.debug(f'\n{sp.head()}')
        for i,srow in sp.iterrows():
            ri=srow.values
            if not np.any(ri==ri):
                logger.debug(f'{i} {ri}')
            C=self.cellndx_of_point(ri)
            idx=self.ldx_of_cellndx(C)
            adf.loc[i,'linkcell_idx']=idx
        return adf

    def _return_list_lens(self,idx_list,mlists):
        """_return_list_lens returns a list of lengths of lists in list mlists

        :param idx_list: indicies of mlists that will be tallied
        :type idx_list: list of ints
        :param mlists: list of lists, each a list of integer atom global indices
        :type mlists: list of list of ints
        :return: list of lengths of lists
        :rtype: list of ints
        """
        return [len(mlists[i]) for i in idx_list]

    @countTime
    def populate(self,adf,ncpu=1,profile=False,report_info=True):
        """populate Populates linkcell structure
        """
        if profile:
            pr=cProfile.Profile()
            pr.enable()
        N=adf.shape[0]
        logger.debug(f'Assigning {N} atoms to {len(self.cellndx)} cells')
        adf['linkcell_idx']=-1*np.ones(N).astype(int)
        ess='s' if ncpu>1 else ''
        logger.debug(f'Cell assignment will use {ncpu} processor{ess}')
        if ncpu>1:
            p=Pool(processes=ncpu)
            adf_split=np.array_split(adf,ncpu)
            result=p.map(partial(self.populate_par),adf_split)
            p.close()
            p.join()
            rdf=pd.DataFrame()
            for a in result:
                rdf=pd.concat((rdf,a))
        else:
            rdf=self.populate_par(adf)
        my_logger(rdf,logger.debug,dfoutmode='info',just='<',fill='')
        adf['linkcell_idx']=rdf['linkcell_idx'].copy()

    # @countTime
    # def make_memberlists(self,adf,ncpu=1,report_info=True):
        # self.memberlists=[[] for _ in range(self.cellndx.shape[0])]

        # if self.atidxlabel==None:
        #     for i,row in adf.iterrows():
        #         lc_idx=row['linkcell_idx']
        #         try:
        #             self.memberlists[lc_idx].append(i)
        #         except:
        #             logger.debug(f'Linear linkcell index {lc_idx} of atom {i} is out of range.\ncellndx.shape[0] is {self.cellndx.shape[0]}\nThis is a bug.')
        #             raise Exception
        # else:
        #     idx_list=adf[self.atidxlabel].to_list()
        #     for i in idx_list:
        #         lc_idx=adf[adf[self.atidxlabel]==i]['linkcell_idx'].values[0]
        #         try:
        #             self.memberlists[lc_idx].append(i)
        #         except:
        #             logger.debug(f'Linear linkcell index {lc_idx} of atom {i} is out of range.\ncellndx.shape[0] is {self.cellndx.shape[0]}\nThis is a bug.')
        #             raise Exception

        # if report_info:
        counts=adf.groupby('linkcell_idx').count().iloc[:,1]

        # idx_list=list(range(len(self.memberlists)))
        # p=Pool(processes=ncpu)
        # idx_list_split=np.array_split(idx_list,ncpu)
        # result=p.map(partial(self._return_list_lens,mlists=self.memberlists),idx_list_split)
        # p.close()
        # p.join()
        # result=np.array([item for sublist in result for item in sublist])
        self.avg_cell_pop=counts.iloc[:].mean()
        self.min_cell_pop=counts.iloc[:].min()
        self.max_cell_pop=counts.iloc[:].max()
        logger.debug(f'Avg/min/max cell pop: {self.avg_cell_pop:>8.3f}/{self.min_cell_pop:>8d}/{self.max_cell_pop:>8d}')
        # if profile:
        #     pr.disable()
        #     s = io.StringIO()
        #     sortby = SortKey.CUMULATIVE
        #     ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        #     ps.print_stats()
        #     msg='\n'.join(s.getvalue().split('\n')[:20])
        #     my_logger(msg,logger.debug,just='<',fill='')

        # logger.debug(f'Linkcell.populate() ends.')

    # def make_neighborlists(self):
    #     """make_neighborlists populates the neighborlist member, one element per cell; each element is the list of neighbors of that cell
    #     """
    #     self.neighborlists=[[] for _ in range(self.cellndx.shape[0])]
    #     # logger.debug(f'{self.cellndx[-1]}')
    #     for C in self.cellndx:
    #         idx=self.ldx_of_cellndx(C)
    #         # if any(C>self.cellndx[-1]):
    #         # logger.debug(f'{C} {self.cellndx[-1]} {C>self.cellndx[-1]} {idx}')
    #         for D in self.neighbors_of_cellndx(C):
    #             if self.ldx_of_cellndx(D)!=idx:
    #                 self.neighborlists[idx].append(self.ldx_of_cellndx(D))
    # def make_memberlists(self,cdf):
    #     """make_memberlists populates the memberlists member, one element per cell; each element is the list of atom indices in that cell 

    #     :param cdf: coordinates data frame
    #     :type cdf: pd.DataFrame
    #     """
    #     self.memberlists=[[] for _ in range(self.cellndx.shape[0])]
    #     rdf=cdf[cdf['linkcell_idx']!=-1]
    #     # logger.debug(f'Generated {len(self.memberlists)} empty memberlists.')
    #     for i,r in rdf.iterrows():
    #         cidx=r['linkcell_idx']
    #         idx=r['globalIdx']
    #         self.memberlists[cidx].append(idx)
    #     rl=np.array([len(self.memberlists[i]) for i in range(self.cellndx.shape[0])])
    #     assert int(rl.sum())==rdf.shape[0] # check to make sure all atoms are counted
    #     avg_cell_pop=rl.mean()
    #     min_cell_pop=int(rl.min())
    #     max_cell_pop=int(rl.max())
    #     logger.debug(f'Avg/min/max cell pop: {avg_cell_pop:>8.3f}/{min_cell_pop:>8d}/{max_cell_pop:>8d}')

    # def neighbors_of_cellndx(self,Ci):
    #     """neighbors_of_cellndx returns the list of neighbors of cell Ci by their (i,j,k) indices

    #     :param Ci: _description_
    #     :type Ci: _type_
    #     :return: _description_
    #     :rtype: _type_
    #     """
    #     assert self.cellndx_in_structure(Ci),f'Error: cell {Ci} outside of cell structure {self.ncells}'
    #     retlist=[]
    #     dd=np.array([-1,0,1])
    #     S=list(product(dd,dd,dd))
    #     S.remove((0,0,0))
    #     for s in S:
    #         nCi=Ci+np.array(s)
    #         for d in range(3):
    #             if nCi[d]==self.ncells[d]:
    #                 nCi[d]=0
    #             elif nCi[d]==-1:
    #                 nCi[d]=self.ncells[d]-1
    #         retlist.append(nCi)
    #     assert len(retlist)==(len(dd)**3)-1
    #     return retlist

    def ldx_searchlist_of_ldx(self,i):
        """searchlist_of_ldx returns the list of cells to search around and including central cell i """
        # assert i!=-1
        C=self.cellndx_of_ldx(i)
        nc=self.cells_per_dim
        d=np.array([-1,0,1])
        searchlist=[self.ldx_of_cellndx(np.mod(C+p,nc)) for p in product(d,d,d)]
        assert i in searchlist
        return searchlist
        # :
        #     nc_ndx=
        # retlist=[]
        # C=self.cellndx[i]
        # for c in self.neighbors_of_cellndx(C):
        #     retlist.append(self.ldx_of_cellndx(c))
        # assert len(retlist)==26,f'Error: not counting enough neighbor cells; found {len(retlist)}'
        # return retlist

    # def are_cellndx_neighbors(self,Ci,Cj):
    #     """are_cellndx_neighbors returns True if cells with (i,j,k) indices Ci and Cj are neighbors

    #     :param Ci: cell index
    #     :type Ci: np.ndarray(3,int)
    #     :param Cj: cell index
    #     :type Cj: np.ndarray(3,int)
    #     :return: True if cells are neighbors, False otherwise
    #     :rtype: bool
    #     """
    #     assert self.cellndx_in_structure(Ci),f'Error: cell {Ci} outside of cell structure {self.ncells}'
    #     assert self.cellndx_in_structure(Cj),f'Error: cell {Cj} outside of cell structure {self.ncells}'
    #     dij=Ci-Cj
    #     low=(dij<-self.ncells/2).astype(int)
    #     hi=(dij>self.ncells/2).astype(int)
    #     dij+=(low-hi)*self.ncells
    #     return all([x in [-1,0,1] for x in dij])

    # def are_ldx_neighbors(self,ildx,jldx):
    #     """are_ldx_neighbors returns True if cells with scalar indices ildx and jldx are neighbors

    #     :param ildx: scalar cell index
    #     :type ildx: int
    #     :param jldx: scalar cell index
    #     :type jldx: int
    #     :return: True if cells are neighbors, False otherwise
    #     :rtype: bool
    #     """
    #     # should never call this for atoms with unset lc indices
    #     assert ildx!=-1
    #     assert jldx!=-1
    #     return jldx in self.neighborlists[ildx]
