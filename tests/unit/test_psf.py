# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import unittest
from pestifer.psf import *
from pestifer.command import Command
from pestifer.mods import *

class TestPSF(unittest.TestCase):
    def test_read_psf(self):
        source='test.psf'
        c=Command(f'grep "\!NATOM" {source}')
        c.run()
        expected_natom=int(c.stdout.strip().split()[0])
        psf=PSFContents(source)
        self.assertEqual(len(psf.atoms),expected_natom)
    
        c=Command(f'grep REMARKS {source} | grep -c DISU')
        c.run()
        expected_nssbonds=int(c.stdout)
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

        for patchtype in Link.allowed_patchnames:
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
