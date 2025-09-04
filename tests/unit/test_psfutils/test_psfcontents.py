# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import logging
import unittest

from pidibble.pdbparse import PDBParser

from pestifer.core.command import Command
from pestifer.molecule.atom import AtomList
from pestifer.molecule.chainidmanager import ChainIDManager
from pestifer.molecule.molecule import Molecule
from pestifer.molecule.residue import ResidueList
from pestifer.molecule.segment import Segment, SegmentList
from pestifer.objs.link import Link
from pestifer.psfutil.psfcontents import PSFContents, get_toppar_from_psf

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

class TestPSFContents(unittest.TestCase):
    
    def test_psfcontents_initialize(self):
        source = 'test.psf'
        c = Command(r'grep "\!NATOM" '+f'{source}')
        c.run()
        expected_natom = int(c.stdout.strip().split()[0])
        psfcontents = PSFContents(source)
        self.assertIsInstance(psfcontents, PSFContents)
        self.assertEqual(len(psfcontents.atoms), expected_natom)
        self.assertEqual(len(psfcontents.segments), 45)
        c = Command(r'grep -c "DISU" ' + source)
        c.run()
        expected_ndisu = int(c.stdout.strip().split()[0])
        self.assertEqual(len(psfcontents.ssbonds), expected_ndisu)
        self.assertEqual(psfcontents.segments[0].segname, 'A')
        self.assertEqual(psfcontents.segments[0].segtype, 'protein')
        self.assertEqual(psfcontents.segments[2].segname, 'C')
        self.assertEqual(psfcontents.segments[2].segtype, 'glycan')

    def test_psfcontents_atom_include_logic(self):
        source = 'test.psf'
        psfcontents = PSFContents(source)
        atom_include_logic = ['chainID == "A"', 'chainID == "B"']
        ignored_atom_count = psfcontents.apply_atom_logics(atom_include_logic, [])
        self.assertEqual(ignored_atom_count, 25540)
        self.assertEqual(len(psfcontents.segments), 2)
        self.assertEqual(psfcontents.segments[0].segname, 'A')
        self.assertEqual(psfcontents.segments[1].segname, 'B')

    def test_psf_atomcount(self):
        source='test.psf'
        c=Command(r'grep "\!NATOM" '+f'{source}')
        c.run()
        expected_natom=int(c.stdout.strip().split()[0])
        psf=PSFContents(source)
        self.assertEqual(len(psf.atoms),expected_natom)
    
    def test_psf_ssbondcount(self):
        source='test.psf'
        c=Command(f'grep REMARKS {source} | grep -c DISU')
        c.run()
        expected_nssbonds=int(c.stdout)
        psf=PSFContents(source)
        self.assertEqual(len(psf.ssbonds),expected_nssbonds)
        if expected_nssbonds>0:
            c=Command(f'grep REMARKS {source} | grep DISU | head -1')
            c.run()
            p=c.stdout.split()
            self.assertEqual(p[0],'REMARKS')
            self.assertEqual(p[1],'patch')
            self.assertEqual(p[2],'DISU')
            fss=psf.ssbonds[0]
            x=fss.chainID1
            y=fss.resid1
            self.assertEqual(p[3],f'{x}:{y}')
            x=fss.chainID2
            y=fss.resid2
            self.assertEqual(p[4],f'{x}:{y}')

    def test_psf_linkcount(self):
        source='test.psf'
        psf=PSFContents(source)
        for patchtype in Link._patch_atomnames.keys():
            c=Command(f'grep REMARKS {source} | grep patch| grep -c {patchtype}')
            c.run(ignore_codes=[1])
            expected=int(c.stdout)
            if not expected:
                self.assertTrue(patchtype not in psf.patches)
            else:
                patches = list(filter(lambda x: x.patchname == patchtype, psf.links))
                self.assertEqual(len(patches),expected)
                if expected>0:
                    t_links=psf.links.filter(lambda x: x.patchname == patchtype)
                    c=Command(f'grep REMARKS {source} | grep {patchtype} | head -1')
                    c.run()
                    p=c.stdout.split()
                    self.assertEqual(p[0],'REMARKS')
                    self.assertEqual(p[1],'patch')
                    self.assertEqual(p[2],patchtype)
                    fss=t_links[0]
                    x=fss.chainID1
                    y=fss.resid1
                    self.assertEqual(p[3],f'{x}:{y}')
                    x=fss.chainID2
                    y=fss.resid2
                    self.assertEqual(p[4],f'{x}:{y}')

    def test_psf_bondcount(self):
        source='test.psf'
        psf=PSFContents(source,parse_topology=['bonds'])
        c=Command(r'grep "\!NBOND" '+f'{source}')
        c.run()
        expected_nbond=int(c.stdout.strip().split()[0])
        self.assertTrue(expected_nbond>0)
        self.assertEqual(len(psf.bonds),expected_nbond)
        c=Command(r'grep -A1 "\!NBOND" '+f'{source} | tail -1')
        c.run()
        firstbondline=[int(_) for _ in c.stdout.strip().split()]
        self.assertEqual(psf.bonds[0].serial1,firstbondline[0])
        self.assertEqual(psf.bonds[0].serial2,firstbondline[1])
        self.assertEqual(psf.bonds[1].serial1,firstbondline[2])
        self.assertEqual(psf.bonds[1].serial2,firstbondline[3])
        self.assertEqual(psf.bonds[2].serial1,firstbondline[4])
        self.assertEqual(psf.bonds[2].serial2,firstbondline[5])
        self.assertEqual(psf.bonds[3].serial1,firstbondline[6])
        self.assertEqual(psf.bonds[3].serial2,firstbondline[7])

    def test_psf_anglecount(self):
        source='test.psf'
        psf=PSFContents(source,parse_topology=['angles'])
        c=Command(r'grep "\!NTHETA" '+f'{source}')
        c.run()
        expected_nangle=int(c.stdout.strip().split()[0])
        self.assertTrue(expected_nangle>0)
        self.assertEqual(len(psf.angles),expected_nangle)
        c=Command(r'grep -A1 "\!NTHETA" '+f'{source} | tail -1')
        c.run()
        firstangleline=[int(_) for _ in c.stdout.strip().split()]
        self.assertEqual(psf.angles[0].serial1,firstangleline[0])
        self.assertEqual(psf.angles[0].serial2,firstangleline[1])
        self.assertEqual(psf.angles[0].serial3,firstangleline[2])
        self.assertEqual(psf.angles[1].serial1,firstangleline[3])
        self.assertEqual(psf.angles[1].serial2,firstangleline[4])
        self.assertEqual(psf.angles[1].serial3,firstangleline[5])
        self.assertEqual(psf.angles[2].serial1,firstangleline[6])
        self.assertEqual(psf.angles[2].serial2,firstangleline[7])
        self.assertEqual(psf.angles[2].serial3,firstangleline[8])

    def test_psf_dihedralcount(self):
        source='test.psf'
        psf=PSFContents(source,parse_topology=['dihedrals'])
        c=Command(r'grep "\!NPHI" '+f'{source}')
        c.run()
        expected_ndihedral=int(c.stdout.strip().split()[0])
        self.assertTrue(expected_ndihedral>0)
        self.assertEqual(len(psf.dihedrals),expected_ndihedral)
        c=Command(r'grep -A1 "\!NPHI" '+f'{source} | tail -1')
        c.run()
        firstdihedralline=[int(_) for _ in c.stdout.strip().split()]
        self.assertEqual(psf.dihedrals[0].serial1,firstdihedralline[0])
        self.assertEqual(psf.dihedrals[0].serial2,firstdihedralline[1])
        self.assertEqual(psf.dihedrals[0].serial3,firstdihedralline[2])
        self.assertEqual(psf.dihedrals[0].serial4,firstdihedralline[3])
        self.assertEqual(psf.dihedrals[1].serial1,firstdihedralline[4])
        self.assertEqual(psf.dihedrals[1].serial2,firstdihedralline[5])
        self.assertEqual(psf.dihedrals[1].serial3,firstdihedralline[6])
        self.assertEqual(psf.dihedrals[1].serial4,firstdihedralline[7])

    def test_psf_get_toppar_from_psf(self):
        source='test.psf'
        toppars=get_toppar_from_psf(source)
        self.assertEqual(len(toppars),4)
        self.assertTrue('toppar_all36_carb_glycopeptide.str' in toppars)
        self.assertTrue('toppar_all36_prot_modify_res.str' in toppars)
        self.assertTrue('toppar_water_ions.str' in toppars)
        self.assertTrue('toppar_all36_moreions.str' in toppars)

        source='equilibrate.psf'
        toppars=get_toppar_from_psf(source)
        self.assertEqual(len(toppars),8)
        self.assertTrue('toppar_water_ions.str' in toppars)
        self.assertTrue('toppar_all36_carb_glycopeptide.str' in toppars)
        self.assertTrue('toppar_all36_prot_modify_res.str' in toppars)
        self.assertTrue('toppar_all36_moreions.str' in toppars)
        self.assertTrue('toppar_all36_lipid_lps.str' in toppars)
        self.assertTrue('toppar_all36_carb_imlab.str' in toppars)
        self.assertTrue('toppar_all36_lipid_cholesterol.str' in toppars)
        self.assertTrue('toppar_all36_lipid_sphingo.str' in toppars)

    
    def test_psf_pdb_congruency(self):
        logger.debug('Parsing companion.pdb..')
        pdb = PDBParser(filepath='companion.pdb').parse().parsed
        logger.debug(f'Building PDB AtomList...')
        pdbatoms = AtomList.from_pdb(pdb).reserialize()
        logger.debug(f'Building PDB ResidueList...')
        pdbresidues = ResidueList.from_residuegrouped_atomlist(pdbatoms)
        logger.debug(f'Parsing companion.psf...')
        psf = PSFContents('companion.psf')
        logger.debug(f'Checking congruency of {len(pdbatoms)} PDB atoms with {len(psf.atoms)} PSF atoms...')
        self.assertEqual(len(pdbatoms), len(psf.atoms))
        self.assertEqual(len(pdbresidues), len(psf.residues))
        self.assertTrue(all([a1.serial == a2.serial for a1, a2 in zip(pdbatoms, psf.atoms)]))
        self.assertTrue(all([a1.name == a2.atomname for a1, a2 in zip(pdbatoms, psf.atoms)]))
        pdbatoms.apply_psf_attributes(psf.atoms)
        for a1, a2 in zip(pdbatoms, psf.atoms):
            self.assertEqual(a1.resid, a2.resid)
            self.assertEqual(a1.segname, a2.segname)
        pdbresidues = ResidueList.from_residuegrouped_atomlist(pdbatoms)
        # self.assertTrue(all([a1.resid == a2.resid for a1, a2 in zip(pdbresidues, psf.residues)]))
        for r1, r2 in zip(pdbresidues, psf.residues):
            self.assertEqual(r1.resid, r2.resid)
            self.assertEqual(r1.resname, r2.resname)
            self.assertEqual(r1.segname, r2.segname)

        pdbsegments = SegmentList([])
        pdbsegments.chainIDmanager = ChainIDManager()
        pdbsegments.residues = pdbresidues
        pdbsegments.seq_spec = {}
        pdbsegments.psfcompanion = psf.segments
        pdbsegments.build_from_psf_and_pdb_data()

        self.assertEqual(len(pdbsegments), len(psf.segments))
        for seg1, seg2 in zip(pdbsegments, psf.segments):
            self.assertEqual(seg1.segname, seg2.segname)
            self.assertEqual(len(seg1.residues), len(seg2.residues))