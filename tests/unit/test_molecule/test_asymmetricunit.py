import unittest
from pestifer.core.objmanager import ObjManager
from pestifer.molecule.segment import Segment, SegmentList
from pestifer.molecule.chainidmanager import ChainIDManager
from pestifer.molecule.asymmetricunit import AsymmetricUnit
from pidibble.pdbparse import PDBParser
from pestifer.util.cifutil import CIFload
from mmcif.api.PdbxContainers import DataContainer
from pathlib import Path
import logging
logger = logging.getLogger(__name__)
class TestAsymmetricUnit(unittest.TestCase):
    
    def test_asymmetric_unit_initialization_empty(self):
        AU = AsymmetricUnit()
        self.assertIsInstance(AU, AsymmetricUnit)
        self.assertFalse(hasattr(AU, 'segments'))

    def test_asymmetric_unit_initialization_from_pdb_downloads(self):
        downloads = {
            'PDB': 'fixtures/data/6pti.pdb', 
            'mmCIF': 'fixtures/data/6pti.cif'
        }
        for fmt, filepath in downloads.items():
            logger.debug(f'Testing asymmetric unit initialization from {fmt} file: {filepath}')
            if fmt == 'PDB':
                pstruct = PDBParser(filepath=filepath, input_format=fmt).parse().parsed
                self.assertIsInstance(pstruct, dict)
            elif fmt == 'mmCIF':
                pstruct = CIFload(Path(filepath))
                self.assertIsInstance(pstruct, DataContainer)
            logger.debug(f'Parsed structure for {fmt} format: {type(pstruct)}')
            AU = AsymmetricUnit(parsed=pstruct,chainIDmanager=ChainIDManager(format=fmt),objmanager=ObjManager())
            self.assertIsInstance(AU, AsymmetricUnit)
            self.assertTrue(hasattr(AU, 'segments'))
            self.assertIsInstance(AU.segments, SegmentList)
            # expected number of segments in 6pti is 3 (protein, water, ions?)
            self.assertEqual(len(AU.segments), 3)
            prot_segment = AU.segments[0]
            self.assertIsInstance(prot_segment, Segment)
            self.assertEqual(prot_segment.chainID, 'A')
            self.assertEqual(prot_segment.segtype, 'protein')
            ion_segment = AU.segments[1]
            self.assertIsInstance(ion_segment, Segment)
            self.assertEqual(ion_segment.chainID, 'B')
            self.assertEqual(ion_segment.segtype, 'ion')
            ion_residues = ion_segment.residues
            self.assertEqual(len(ion_residues), 1)
            ion = ion_residues[0]
            self.assertEqual(ion.resname, 'PO4')
            ion_atoms = ion.atoms
            self.assertEqual(len(ion_atoms), 5)
            water_segment = AU.segments[2]
            self.assertIsInstance(water_segment, Segment)
            self.assertEqual(water_segment.chainID, 'C')
            self.assertEqual(water_segment.segtype, 'water')

            obj_manager = AU.objmanager
            self.assertIsInstance(obj_manager, ObjManager)
            if fmt == 'PDB':
                self.assertIn('seq', obj_manager) 
                # from TER records; not present in mmCIF files
            
            self.assertTrue('topol' in obj_manager)
            self.assertTrue('ssbonds' in obj_manager['topol'])
            for ssbond in obj_manager['topol']['ssbonds']:
                logger.debug(f"SSBond: {ssbond.describe()}")
            self.assertEqual(len(obj_manager['topol']['ssbonds']), 3)
            ssbonds = obj_manager['topol']['ssbonds']
            self.assertEqual(ssbonds[0].resseqnum1, 5)
            self.assertEqual(ssbonds[0].resseqnum2, 55)