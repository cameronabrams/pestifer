import unittest
from pathlib import Path

from pidibble.pdbparse import PDBParser

from pestifer.molecule.asymmetricunit import AsymmetricUnit
from pestifer.molecule.chainidmanager import ChainIDManager
from pestifer.molecule.molecule import Molecule
from pestifer.core.objmanager import ObjManager
from pestifer.objs.graft import Graft, GraftList
from pestifer.objs.resid import ResID

INPUTS = Path(__file__).parents[2] / 'inputs'


class TestGraftIntegration(unittest.TestCase):
    """
    Integration test for the graft mechanism.

    Source: 4b7i chain C (8-residue glycan: NAG1, NAG2, BMA3, MAN4-7, FUC8)
    Base:   1gc1 chains G (protein) + A (glycan NAG1→FUC2, linked to ASN G234)

    After grafting NAG C1 from 4b7i onto NAG A1 in 1gc1:
      - chain A glycan (NAG1, FUC2) → renamed to segname 'GG01' (parent chain G)
      - graft donor residues (NAG2, BMA3, MAN4-7, FUC8) → new segname 'GG02'
      - external links (index→donor): segname1='GG01', segname2='GG02'
      - internal links (donor→donor): segname1='GG02', segname2='GG02'
    """

    @classmethod
    def setUpClass(cls):
        # Build source AsymmetricUnit from 4b7i (chain C glycan only,
        # skip_glycan_renaming keeps chain C as segname 'C').
        parsed_4b7i = PDBParser(
            filepath=INPUTS / '4b7i.pdb', input_format='PDB'
        ).parse().parsed
        om_source = ObjManager()
        source_au = AsymmetricUnit(
            parsed=parsed_4b7i,
            sourcespecs={'sequence': {'skip_glycan_renaming': True}},
            objmanager=om_source,
            chainIDmanager=ChainIDManager(format='PDB'),
        )
        # activate() needs mol.asymmetric_unit.segments and mol.objmanager
        cls.source_mol = Molecule()
        cls.source_mol.asymmetric_unit = source_au
        cls.source_mol.objmanager = om_source

        # Parse 1gc1 for use in test methods (absolute path avoids conftest dir change)
        cls.parsed_1gc1 = PDBParser(
            filepath=INPUTS / '1gc1.pdb', input_format='PDB'
        ).parse().parsed

    def test_graft_activation_counts(self):
        g = Graft(
            chainID='A',
            target_root=ResID(1),
            source_pdbid='4b7i',
            source_chainID='C',
            source_root=ResID(1),
        )
        g.activate(self.source_mol)

        # NAG C1 is the sole index residue
        self.assertEqual(len(g.index_residues), 1)
        self.assertEqual(g.index_residues[0].resname, 'NAG')
        self.assertEqual(g.index_residues[0].resid, ResID(1))

        # Remaining 7 chain-C glycan residues are donors
        self.assertEqual(len(g.donor_residues), 7)

        # NAG C1→NAG C2 and NAG C1→FUC C8 are external links
        self.assertEqual(len(g.donor_external_links), 2)

        # NAG C2→BMA C3, BMA C3→MAN C4, BMA C3→MAN C7,
        # MAN C4→MAN C5, MAN C4→MAN C6 are internal links
        self.assertEqual(len(g.donor_internal_links), 5)

    def test_graft_segname_and_link_segnames(self):
        g = Graft(
            chainID='A',
            target_root=ResID(1),
            source_pdbid='4b7i',
            source_chainID='C',
            source_root=ResID(1),
        )
        g.activate(self.source_mol)

        graft_list = GraftList([g])
        om_base = ObjManager()
        om_base.ingest(graft_list)

        # Build the base AU from 1gc1, protein chain G + glycan chain A only.
        # Chain G has several embedded glycans plus chain A (NAG1→FUC2).
        base_au = AsymmetricUnit(
            parsed=self.parsed_1gc1,
            sourcespecs={"include": ["chainID in ('G', 'A')"]},
            objmanager=om_base,
            chainIDmanager=ChainIDManager(format='PDB'),
        )

        # Receiver residue was found (NAG A1) and assigned
        self.assertIsNotNone(g.residues)
        self.assertEqual(len(g.residues), 1)

        # After segment generation, g.chainID was set to the receiver's psfgen segname
        self.assertEqual(g.chainID, g.graft_segname)

        # graft_segname equals the receiver's existing glycan segment
        self.assertIsNotNone(g.graft_segname)
        self.assertIn(g.graft_segname, base_au.segments.segnames)
        self.assertTrue(g.graft_segname.startswith('GG'))

        # g.chainID was set to the glycan segname (for inherit_objs routing)
        self.assertEqual(g.chainID, g.graft_segname)

        # Donor counts unchanged
        self.assertEqual(len(g.donor_residues), 7)
        self.assertEqual(len(g.donor_external_links), 2)
        self.assertEqual(len(g.donor_internal_links), 5)

        # All graft links (external and internal) reference the same glycan segment
        # because donors join the receiver's existing segment.
        for lnk in g.donor_external_links.data:
            self.assertEqual(lnk.segname1, g.graft_segname,
                             f'External link {lnk}: expected segname1={g.graft_segname}')
            self.assertEqual(lnk.segname2, g.graft_segname,
                             f'External link {lnk}: expected segname2={g.graft_segname}')

        for lnk in g.donor_internal_links.data:
            self.assertEqual(lnk.segname1, g.graft_segname,
                             f'Internal link {lnk}: expected segname1={g.graft_segname}')
            self.assertEqual(lnk.segname2, g.graft_segname,
                             f'Internal link {lnk}: expected segname2={g.graft_segname}')
