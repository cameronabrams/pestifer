import unittest
from pestifer.linkcell import Linkcell
from pestifer.ring import RingList
from pestifer.coord import coorddf_from_pdb
from pestifer.util import cell_from_xsc
from pestifer.psf import PSFContents,PSFBondList
import os
import logging
logger=logging.getLogger(__name__)
import numpy as np

class TestRingCheck(unittest.TestCase):
    def test_ring_check_coords(self):
        self.box,self.orig=cell_from_xsc('test.xsc')
        self.coorddf=coorddf_from_pdb('test.pdb')
        self.LC=Linkcell(self.box,10.0,origin=self.orig)
        self.LC.populate(self.coorddf,os.cpu_count())
        self.topol=PSFContents('test.psf',parse_topology=True)
        self.topol.nonsolvent_atoms.injest_coordinates(self.coorddf,idx_key='serial',pos_key=['x','y','z'],meta_key=['linkcell_idx','residue'],box=self.box)
        self.assertTrue(self.topol.nonsolvent_atoms[0].meta['residue'].resName=='LEU')
        self.assertTrue(self.topol.bonds.validate_images(self.box))
        for bond in self.topol.bonds:
            logger.debug(f'bond com {bond.COM}')
            bond.oc_ldx=self.LC.ldx_of_cellndx(self.LC.cellndx_of_point(bond.COM))
        self.Rings=RingList(self.topol.G)
        self.Rings.injest_coordinates(self.coorddf,idx_key='serial',pos_key=['x','y','z'],meta_key=['linkcell_idx','residue'],box=self.box)
        for ring in self.Rings:
            oc=self.LC.ldx_of_cellndx(self.LC.cellndx_of_point(ring.O))
            logger.debug(f'ring center in cell {oc} -- members in {ring.M["linkcell_idx"]}')
            nc=self.LC.searchlist_of_ldx(oc)
            logger.debug(f'...search for bonds in cells {nc} and {oc}')
            search_bonds=PSFBondList([])
            for ldx in nc+[oc]:
                this_sb=[b for b in self.topol.bonds if b.oc_ldx==ldx]
                search_bonds.extend(this_sb)
            logger.debug(f'...search will include {len(search_bonds)} bonds')
            piercings=[]
            for bond in search_bonds:
                p,r=ring.pierced_by(np.array([bond.atom1.r,bond.atom2.r]))
                piercings.append(p)
            if any(piercings):
                logger.debug(f'bond {bond.atom1.serial}-{bond.atom2.serial} pierces this ring')