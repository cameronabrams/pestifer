import logging
from pestifer.ring import RingList
import unittest
from pestifer.psfutil.psfcontents import PSFContents
from pestifer.util.coord import coorddf_from_pdb
from pestifer.util.util import cell_from_xsc

logger=logging.getLogger(__name__)

class TestRings(unittest.TestCase):
    def test_ringlist_count(self):
        psf=PSFContents('test.psf',topology_segtypes=['lipid'],parse_topology=['bonds'])
        Rings=RingList(psf.G)
        logger.debug(f'{len(Rings)} rings detected.')
        nrings=len(Rings)
        self.assertEqual(nrings,560)
        nrings_by_size={}
        for r in Rings:
            s=len(r)
            if not s in nrings_by_size:
                nrings_by_size[s]=0
            nrings_by_size[s]+=1
        self.assertTrue(6 in nrings_by_size)
        self.assertTrue(5 in nrings_by_size)
        self.assertFalse(4 in nrings_by_size)
        self.assertFalse(7 in nrings_by_size)

    def test_ringlist_coords(self):
        box,orig=cell_from_xsc('test.xsc')
        coorddf=coorddf_from_pdb('test.pdb')
        topol=PSFContents('test.psf',parse_topology=['bonds'])
        self.assertEqual(coorddf.shape[0],len(topol.atoms))
        Rings=RingList(topol.G)
        Rings.injest_coordinates(coorddf,pos_key=['x','y','z'])
        Rings.validate_images(box)
        self.assertTrue(all([x.same_image for x in Rings]))