"""

.. module:: test_bioassemb
   :synopsis: tests pestifer.bioassemb
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""

import unittest
from pestifer.bioassemb import BioAssemb, Transform
from pestifer.asymmetricunit import AsymmetricUnit
import numpy as np

class TestBiomT(unittest.TestCase):
    def test_identity(self):
        bi=Transform('identity')
        I=np.array([[1, 0, 0, 0],[0, 1, 0, 0],[0, 0, 1, 0],[0, 0, 0, 1]],dtype=float)
        self.assertTrue(np.array_equal(I,bi.tmat))
    def test_actual(self):
        R=np.array([[1, 1, 1],[0, 1, 0],[0, 0, 1]],dtype=float)
        RT=np.array([[1, 1, 1, 5],[0, 1, 0, 4],[0, 0, 1, 3],[0, 0, 0, 1]],dtype=float)
        T=np.array([5, 4, 3],dtype=float)
        bi=Transform(R,T,[],0)
        self.assertTrue(np.array_equal(RT,bi.tmat))
    def test_au(self):
        au=AsymmetricUnit() # empty
        ba=BioAssemb(au)
        bi=ba.transforms[0]
        I=np.array([[1, 0, 0, 0],[0, 1, 0, 0],[0, 0, 1, 0],[0, 0, 0, 1]],dtype=float)
        self.assertTrue(np.array_equal(I,bi.tmat))

