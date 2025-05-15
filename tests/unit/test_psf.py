# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import unittest

from pestifer.command import Command
from pestifer.objs.link import Link
from pestifer.psfutil.psfbond import PSFBond
from pestifer.psfutil.psfangle import PSFAngle
from pestifer.psfutil.psfdihedral import PSFDihedral
from pestifer.psfutil.psfcontents import PSFContents, get_toppar_from_psf


class TestPSF(unittest.TestCase):

    def test_psf_bond_eq(self):
        b1=PSFBond([1,2])
        b2=PSFBond([2,1])
        self.assertEqual(b1,b2)

    def test_psf_angle_eq(self):
        a1=PSFAngle([1,2,3])
        a2=PSFAngle([3,2,1])
        self.assertEqual(a1,a2)

    def test_psf_dihedral_eq(self):
        d1=PSFDihedral([1,2,3,4])
        d2=PSFDihedral([4,3,2,1])
        self.assertEqual(d1,d2)

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
            y=f'{fss.resseqnum1}{fss.insertion1}'
            self.assertEqual(p[3],f'{x}:{y}')
            x=fss.chainID2
            y=f'{fss.resseqnum2}{fss.insertion2}'
            self.assertEqual(p[4],f'{x}:{y}')

    def test_psf_linkcount(self):
        source='test.psf'
        psf=PSFContents(source)
        for patchtype in Link.patch_atomnames.keys():
            c=Command(f'grep REMARKS {source} | grep patch| grep -c {patchtype}')
            c.run(ignore_codes=[1])
            expected=int(c.stdout)
            if not expected:
                self.assertTrue(patchtype not in psf.patches)
            else:
                self.assertEqual(len(psf.patches[patchtype]),expected)
                if expected>0:
                    t_links=psf.links.filter(patchname=patchtype)
                    c=Command(f'grep REMARKS {source} | grep {patchtype} | head -1')
                    c.run()
                    p=c.stdout.split()
                    self.assertEqual(p[0],'REMARKS')
                    self.assertEqual(p[1],'patch')
                    self.assertEqual(p[2],patchtype)
                    fss=t_links[0]
                    x=fss.chainID1
                    y=f'{fss.resseqnum1}{fss.insertion1}'
                    self.assertEqual(p[3],f'{x}:{y}')
                    x=fss.chainID2
                    y=f'{fss.resseqnum2}{fss.insertion2}'
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