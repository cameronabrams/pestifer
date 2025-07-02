# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
Defines the :class:`PSFTopoElement` and :class:`PSFTopoElementList` classes for handling PSF topology elements.

These classes provide a framework for representing and manipulating the various elements found in PSF files, including atoms, bonds, angles, and dihedrals.
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
    """
    A class for handling a list of lines, typically used for reading and writing PSF topology elements.
    This class inherits from :class:`collections.UserList` and provides additional functionality for managing a list of lines
    that represent PSF topology elements."""
    def __init__(self,data):
        super().__init__(data)

class PSFTopoElement(Namespace):
    """
    A class representing a single PSF topology element, such as an atom, bond, angle, or dihedral.
    This class inherits from :class:`argparse.Namespace`, allowing it to store attributes as if they were command-line arguments.
    It provides methods for ingesting coordinates, calculating properties, and validating images.
    The class also includes methods for mic-shifting coordinates, copying the element, and matching against specific criteria.
    
    Parameters
    ----------
    idx_list : list
        A list of indices representing the serial numbers of the atoms involved in this topology element.
        The indices are used to uniquely identify the atoms in the PSF file.
    """
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
        """
        Calculate various properties of the topology element based on its coordinates.
        This method computes the center of mass (COM), radii, normalized radii, and the distance matrix between points.
        It also calculates the vector differences between all pairs of points and their corresponding distances.
        
        Attributes
        ----------
        COM : numpy.ndarray
            The center of mass of the points in the topology element.
        Radii : numpy.ndarray
            The vector differences of the points from the center of mass.
        radii : numpy.ndarray
            The normalized radii of the points, calculated as the vector differences divided by their magnitudes.
        B : numpy.ndarray
            The vector differences between all pairs of points in the topology element.
        b : numpy.ndarray
            The distances between all pairs of points, calculated as the magnitudes of the vector differences.
        """
        self.COM=np.mean(self.P,axis=0)
        self.Radii=self.P-self.COM
        self.radii=np.array([n/d for n,d in zip(self.Radii,np.linalg.norm(self.Radii,axis=1))])
        self.B=np.array([p[1]-p[0] for p in product(self.P,repeat=2)]).reshape((len(self.P),len(self.P),3))
        self.b=np.linalg.norm(self.B.reshape((len(self.P)**2,3)),axis=1).reshape(len(self.P),len(self.P))

    def ingest_coordinates(self,A,pos_key=['posX','posY','posZ'],meta_key=[]):
        """
        Ingest coordinates from a DataFrame and store them in the topology element.
        This method extracts the coordinates of the points associated with the indices in `self.idx_list`
        and stores them in the `self.P` attribute. It also calculates properties if there are multiple points.
        
        Parameters
        ----------
        A : pandas.DataFrame
            A DataFrame containing the coordinates and metadata for the topology element.
        pos_key : list, optional
            A list of column names in the DataFrame that contain the position coordinates (default is ['posX', 'posY', 'posZ']).
        meta_key : list, optional
            A list of column names in the DataFrame that contain metadata for the topology element (default is an empty list).
            
        Attributes
        ----------
        P : numpy.ndarray
            An array of shape (n, 3) containing the coordinates of the points in the topology element.
        same_image : bool
            A boolean indicating whether all points are in the same image based on their coordinates.
        """
        self.P=np.array(A.loc[self.idx_list][pos_key].values)
        if np.any(np.isnan(self.P)):
            logger.debug(f'nan detected in {self.P}')
        if len(self.P)>1:
            self.calculate_stuff()
        for k in meta_key:
            self.__dict__[k]=A.loc[self.idx_list[0],k]

    def assign_cell_index(self,LC):
        """
        Assign the link cell index to the topology element based on its center of mass.
        This method uses the `linkcell` module to determine the link cell index of the center of mass (COM) of the topology element.
        
        Parameters  
        ----------
        LC : LinkCell
            An instance of the :class:`LinkCell <..util.linkcell.LinkCell>` class that provides methods for determining link cell indices.
        """
        self.linkcell_idx=LC.ldx_of_cellndx(LC.cellndx_of_point(self.COM))

    def validate_image(self,box):
        """
        Validate whether all points in the topology element are within the same image based on the provided box dimensions.
        This method checks if the distance from the center of mass (COM) to each point is less than half the box length in each dimension.
        
        Parameters
        ----------
        box : numpy.ndarray
            A 2D array representing the box dimensions, where each row corresponds to a dimension (x, y, z) and contains the lower and upper bounds of that dimension.

        Returns
        -------
        bool
            A boolean indicating whether all points are in the same image. If all points are within half the box length in each dimension, it returns `True`; otherwise, it returns `False`.
        """
        boxlen=np.array([box[i][i] for i in [0,1,2]])
        in_same_img=np.linalg.norm(self.P-self.COM)<(boxlen/2)
        self.same_image=np.all(in_same_img)
        return self.same_image
    
    def mic_shift(self,reference_point,box):
        """
        Perform a minimum image convention (MIC) shift on the coordinates of the topology element.
        This method shifts the coordinates of the points in the topology element so that they are within the minimum image convention,
        relative to a reference point, and within the specified box dimensions.
        
        Parameters
        ----------
        reference_point : numpy.ndarray
            A 1D array representing the reference point for the MIC shift.
        box : numpy.ndarray
            A 2D array representing the box dimensions, where each row corresponds to a dimension (x, y, z) and contains the lower and upper bounds of that dimension.
        """
        holder=self.copy()
        for i in range(len(self.P)):
            holder.P[i],bl=mic_shift(self.P[i],reference_point,box)
        if np.any(bl): holder.calculate_stuff()
        return holder
    
    def copy(self):
        """
        Create a deep copy of the PSFTopoElement instance.
        This method returns a new instance of the PSFTopoElement class with the same attributes as the original instance.
        """
        return deepcopy(self)

    def match(self,**crits):
        """
        Check if the topology element matches the specified criteria.
        This method compares the attributes of the topology element against the provided keyword arguments (criteria).
        
        Parameters
        ----------
        **crits : dict
            A dictionary of criteria where the keys are attribute names and the values are the expected values for those attributes.

        Returns
        -------
        bool
            A boolean indicating whether the topology element matches the specified criteria. If all criteria are met, it returns `True`; otherwise, it returns `False`.
        """
        return all([self.__dict__[k]==v for k,v in crits.items()])

class PSFTopoElementList(UserList):
    """
    A class representing a list of PSF topology elements, such as atoms, bonds, angles, or dihedrals.
    This class inherits from :class:`collections.UserList` and provides methods for managing a list of PSF topology elements.
    It can be initialized from a list of lines or from another :class:`PSFTopoElementList <.psftopoelement.PSFTopoElementList>` object.

    Parameters
    ----------
    data : list
        A list of PSF topology elements or lines representing the topology elements.
    """
    def __init__(self,data):
        super().__init__(data)

    @countTime
    def ingest_coordinates(self,A,pos_key=['posX','posY','posZ'],meta_key=[]):
        """
        Ingest coordinates from a DataFrame and store them in each topology element in the list.
        This method iterates over each PSF topology element in the list and calls the ``ingest_coordinates`` method
        of each element to extract and store the coordinates from the DataFrame.

        Parameters
        ----------
        A : pandas.DataFrame
            The DataFrame containing the coordinates to be ingested.
        pos_key : list of str
            The list of column names in the DataFrame corresponding to the position coordinates.
        meta_key : list of str
            The list of column names in the DataFrame corresponding to the metadata.
        """
        for b in self:
            b.ingest_coordinates(A,pos_key=pos_key,meta_key=meta_key)

    @countTime
    def assign_cell_indexes(self,LC):
        """
        Assign the link cell index to each topology element in the list based on their center of mass.
        This method iterates over each PSF topology element in the list and calls the ``assign_cell_index`` method
        of each element to determine and assign the link cell index using the provided link cell object.

        Parameters
        ----------
        LC : LinkCell
            The link cell object used to assign the link cell index to each topology element.
        """
        for b in self:
            b.assign_cell_index(LC)

    def validate_images(self,box):
        """
        Validate whether all topology elements in the list are within the same image based on the provided box
        dimensions. 
        This method checks if each topology element's coordinates are within half the box length in each dimension.

        Parameters
        ----------
        box : Box
            The box object representing the simulation box dimensions.
        
        Returns
        -------
        bool
            A boolean indicating whether all topology elements are in the same image. If all elements are within
            half the box length in each dimension, it returns ``True``; otherwise, it returns ``False``.
        """
        return all([b.validate_image(box) for b in self])

    def get(self,**crits):
        result=[]
        for e in self:
            if e.match(**crits):
                result.append(e)
        return type(self)(result)
    



    
