import unittest
from pestifer.objs.crot import Crot, CrotList

class TestCrot(unittest.TestCase):
    def test_crot_creation_phi(self):
        crot = Crot(
            angle='PHI',
            chainID="A",
            resseqnum1=5,
            insertion1="",
            resseqnum2=15,
            insertion2="",
            degrees=180.0
        )
        self.assertIsInstance(crot, Crot)
        self.assertEqual(crot.chainID, "A")
        self.assertEqual(crot.resseqnum1, 5)
        self.assertEqual(crot.insertion1, "")
        self.assertEqual(crot.resseqnum2, 15)
        self.assertEqual(crot.insertion2, "")
        self.assertEqual(crot._yaml_header, 'crotations')
        self.assertEqual(crot._objcat, 'coord')
        self.assertEqual(repr(crot), "Crot(angle='PHI', chainID='A', resseqnum1=5, insertion1='', resseqnum2=15, insertion2='', degrees=180.0)")

    def test_crot_creation_psi(self):
        crot = Crot(
            angle='PSI',
            chainID="B",
            resseqnum1=10,
            insertion1="",
            resseqnum2=20,
            insertion2="",
            degrees=120.0
        )
        self.assertIsInstance(crot, Crot)
        self.assertEqual(crot.chainID, "B")
        self.assertEqual(crot.resseqnum1, 10)
        self.assertEqual(crot.insertion1, "")
        self.assertEqual(crot.resseqnum2, 20)
        self.assertEqual(crot.insertion2, "")
        self.assertEqual(crot._yaml_header, 'crotations')
        self.assertEqual(crot._objcat, 'coord')
        self.assertEqual(repr(crot), "Crot(angle='PSI', chainID='B', resseqnum1=10, insertion1='', resseqnum2=20, insertion2='', degrees=120.0)")

    def test_crot_creation_omega(self):
        crot = Crot(
            angle='OMEGA',
            chainID="C",
            resseqnum1=15,
            insertion1="",
            resseqnum2=25,
            insertion2="",
            degrees=240.0
        )
        self.assertIsInstance(crot, Crot)
        self.assertEqual(crot.chainID, "C")
        self.assertEqual(crot.resseqnum1, 15)
        self.assertEqual(crot.insertion1, "")
        self.assertEqual(crot.resseqnum2, 25)
        self.assertEqual(crot.insertion2, "")
        self.assertEqual(crot._yaml_header, 'crotations')
        self.assertEqual(crot._objcat, 'coord')
        self.assertEqual(repr(crot), "Crot(angle='OMEGA', chainID='C', resseqnum1=15, insertion1='', resseqnum2=25, insertion2='', degrees=240.0)")

    def test_crot_creation_chi1(self):
        crot = Crot(
            angle='CHI1',
            chainID="D",
            resseqnum1=1,
            insertion1="X",
            degrees=60.0
        )
        self.assertIsInstance(crot, Crot)
        self.assertEqual(crot.chainID, "D")
        self.assertEqual(crot.resseqnum1, 1)
        self.assertEqual(crot.insertion1, "X")
        self.assertEqual(crot._yaml_header, 'crotations')
        self.assertEqual(crot._objcat, 'coord')
        self.assertEqual(repr(crot), "Crot(angle='CHI1', chainID='D', resseqnum1=1, insertion1='X', degrees=60.0)")

    def test_crot_creation_angleijk(self):
        crot = Crot(
            angle='ANGLEIJK',
            segnamei="E",
            resseqnumi=2,
            insertioni="Y",
            atomi="NC",
            segnamejk="F",
            resseqnumj=3,
            atomj="CN",
            insertionj="Z",
            resseqnumk=4,
            insertionk="W",
            atomk="CA",
            degrees=90.0
        )
        self.assertIsInstance(crot, Crot)
        self.assertEqual(crot.segnamei, "E")
        self.assertEqual(crot.resseqnumi, 2)
        self.assertEqual(crot.insertioni, "Y")
        self.assertEqual(crot.segnamejk, "F")
        self.assertEqual(crot.resseqnumj, 3)
        self.assertEqual(crot.insertionj, "Z")
        self.assertEqual(crot.resseqnumk, 4)
        self.assertEqual(crot.insertionk, "W")
        self.assertEqual(crot.atomk, "CA")
        self.assertEqual(crot._yaml_header, 'crotations')
        self.assertEqual(crot._objcat, 'coord')
        self.maxDiff = None
        self.assertEqual(repr(crot), "Crot(angle='ANGLEIJK', segnamei='E', resseqnumi=2, insertioni='Y', atomi='NC', segnamejk='F', resseqnumj=3, insertionj='Z', atomj='CN', resseqnumk=4, insertionk='W', atomk='CA', degrees=90.0)")

    def test_crot_creation_alpha(self):
        crot = Crot(
            angle='ALPHA',
            chainID="G",
            resseqnum1=8,
            insertion1="",
            resseqnum2=12,
            insertion2="",
            resseqnum3= 16,
            insertion3=""
        )
        self.assertIsInstance(crot, Crot)
        self.assertEqual(crot.chainID, "G")
        self.assertEqual(crot.resseqnum1, 8)
        self.assertEqual(crot.insertion1, "")
        self.assertEqual(crot.resseqnum2, 12)
        self.assertEqual(crot.insertion2, "")
        self.assertEqual(crot._yaml_header, 'crotations')
        self.assertEqual(crot._objcat, 'coord')
        self.assertEqual(repr(crot), "Crot(angle='ALPHA', chainID='G', resseqnum1=8, insertion1='', resseqnum2=12, insertion2='', resseqnum3=16, insertion3='')")

    def test_crot_creation_invalid_angle(self):
        with self.assertRaises(ValueError):
            Crot(
                angle='invalid',
                chainID="D",
                resseqnum1=1,
                insertion1="",
                resseqnum2=2,
                insertion2="",
                degrees=90.0
            )

    def test_crot_from_shortcode_phi(self):
        crot = Crot("PHI,A,5,15,180.0")
        self.assertIsInstance(crot, Crot)
        self.assertEqual(crot.angle, "PHI")
        self.assertEqual(crot.chainID, "A")
        self.assertEqual(crot.resseqnum1, 5)
        self.assertEqual(crot.insertion1, "")
        self.assertEqual(crot.resseqnum2, 15)
        self.assertEqual(crot.insertion2, "")

    def test_crot_from_shortcode_psi(self):
        crot = Crot("PSI,B,10,20,120.0")
        self.assertIsInstance(crot, Crot)
        self.assertEqual(crot.angle, "PSI")
        self.assertEqual(crot.chainID, "B")
        self.assertEqual(crot.resseqnum1, 10)
        self.assertEqual(crot.insertion1, "")
        self.assertEqual(crot.resseqnum2, 20)
        self.assertEqual(crot.insertion2, "")
    
    def test_crot_from_shortcode_omega(self):
        crot = Crot("OMEGA,C,15,25,240.0")
        self.assertIsInstance(crot, Crot)
        self.assertEqual(crot.angle, "OMEGA")
        self.assertEqual(crot.chainID, "C")
        self.assertEqual(crot.resseqnum1, 15)
        self.assertEqual(crot.insertion1, "")
        self.assertEqual(crot.resseqnum2, 25)
        self.assertEqual(crot.insertion2, "")
    
    def test_crot_from_shortcode_chi1(self):
        crot = Crot("CHI1,D,1X,60.0")
        self.assertIsInstance(crot, Crot)
        self.assertEqual(crot.angle, "CHI1")
        self.assertEqual(crot.chainID, "D")
        self.assertEqual(crot.resseqnum1, 1)
        self.assertEqual(crot.insertion1, "X")
        self.assertEqual(crot.degrees, 60.0)

    def test_crot_from_shortcode_angleijk(self):
        crot = Crot("ANGLEIJK,E,2Y,NC,F,3Z,CN,4W,CA,90.0")
        self.assertIsInstance(crot, Crot)
        self.assertEqual(crot.angle, "ANGLEIJK")
        self.assertEqual(crot.segnamei, "E")
        self.assertEqual(crot.resseqnumi, 2)
        self.assertEqual(crot.insertioni, "Y")
        self.assertEqual(crot.segnamejk, "F")
        self.assertEqual(crot.resseqnumj, 3)
        self.assertEqual(crot.insertionj, "Z")
        self.assertEqual(crot.resseqnumk, 4)
        self.assertEqual(crot.insertionk, "W")
        self.assertEqual(crot.atomk, "CA")
        self.assertEqual(crot.degrees, 90.0)
    
    def test_crot_from_shortcode_alpha_full(self):
        crot = Crot("ALPHA,G,8,12,16")
        self.assertIsInstance(crot, Crot)
        self.assertEqual(crot.angle, "ALPHA")
        self.assertEqual(crot.chainID, "G")
        self.assertEqual(crot.resseqnum1, 8)
        self.assertEqual(crot.insertion1, "")
        self.assertEqual(crot.resseqnum2, 12)
        self.assertEqual(crot.insertion2, "")
        self.assertEqual(crot.resseqnum3, 16)
        self.assertEqual(crot.insertion3, "")

    def test_crot_from_shortcode_alpha_short(self):
        crot = Crot("ALPHA,G,8,12")
        self.assertIsInstance(crot, Crot)
        self.assertEqual(crot.angle, "ALPHA")
        self.assertEqual(crot.chainID, "G")
        self.assertEqual(crot.resseqnum1, 8)
        self.assertEqual(crot.insertion1, "")
        self.assertEqual(crot.resseqnum2, 12)
        self.assertEqual(crot.insertion2, "")
        self.assertEqual(crot.resseqnum3, 12)
        self.assertEqual(crot.insertion3, "")

    def test_crot_list_creation(self):
        crot_list = CrotList()
        self.assertIsInstance(crot_list, CrotList)
        self.assertEqual(len(crot_list), 0)