import unittest
from pestifer.objs.crot import Crot, CrotList
from pestifer.objs.resid import ResID

class TestCrot(unittest.TestCase):
    def test_crot_creation_phi(self):
        crot = Crot(
            angle = 'PHI',
            chainID = "A",
            resid1 = ResID(5),
            resid2 = ResID(15),
            degrees = 180.0
        )
        self.assertIsInstance(crot, Crot)
        self.assertEqual(crot.chainID, "A")
        self.assertEqual(crot.resid1.resid, 5)
        self.assertEqual(crot.resid2.resid, 15)
        self.assertEqual(crot._yaml_header, 'crotations')
        self.assertEqual(crot._objcat, 'coord')
        self.assertEqual(repr(crot), "Crot(angle='PHI', chainID='A', resid1=ResID(resseqnum=5), resid2=ResID(resseqnum=15), degrees=180.0)")

    def test_crot_creation_psi(self):
        crot = Crot(
            angle = 'PSI',
            chainID = "B",
            resid1 = ResID(10),
            resid2 = ResID(20),
            degrees =  120.0
        )
        self.assertIsInstance(crot, Crot)
        self.assertEqual(crot.chainID, "B")
        self.assertEqual(crot.resid1.resid, 10)
        self.assertEqual(crot.resid2.resid, 20)
        self.assertEqual(crot._yaml_header, 'crotations')
        self.assertEqual(crot._objcat, 'coord')
        self.assertEqual(repr(crot), "Crot(angle='PSI', chainID='B', resid1=ResID(resseqnum=10), resid2=ResID(resseqnum=20), degrees=120.0)")

    def test_crot_creation_omega(self):
        crot = Crot(
            angle = 'OMEGA',
            chainID = "C",
            resid1 = ResID(15),
            resid2 = ResID(25),
            degrees = 240.0
        )
        self.assertIsInstance(crot, Crot)
        self.assertEqual(crot.chainID, "C")
        self.assertEqual(crot.resid1.resid, 15)
        self.assertEqual(crot.resid2.resid, 25)
        self.assertEqual(crot._yaml_header, 'crotations')
        self.assertEqual(crot._objcat, 'coord')
        self.assertEqual(repr(crot), "Crot(angle='OMEGA', chainID='C', resid1=ResID(resseqnum=15), resid2=ResID(resseqnum=25), degrees=240.0)")

    def test_crot_creation_chi1(self):
        crot = Crot(
            angle = 'CHI1',
            chainID = "D",
            resid1 = ResID('1X'),
            degrees = 60.0
        )
        self.assertIsInstance(crot, Crot)
        self.assertEqual(crot.chainID, "D")
        self.assertEqual(crot.resid1.resid, '1X')
        self.assertEqual(crot._yaml_header, 'crotations')
        self.assertEqual(crot._objcat, 'coord')
        self.assertEqual(repr(crot), "Crot(angle='CHI1', chainID='D', resid1=ResID(resseqnum=1, insertion='X'), degrees=60.0)")

    def test_crot_creation_angleijk(self):
        crot = Crot(
            angle='ANGLEIJK',
            segnamei="E",
            residi = ResID('2Y'),
            atomi="NC",
            segnamejk="F",
            residj = ResID(3),
            atomj="CN",
            residk = ResID(4, 'W'),
            atomk="CA",
            degrees=90.0
        )
        self.assertIsInstance(crot, Crot)
        self.assertEqual(crot.segnamei, "E")
        self.assertEqual(crot.residi.resid, '2Y')
        self.assertEqual(crot.segnamejk, "F")
        self.assertEqual(crot.residj.resid, 3)
        self.assertEqual(crot.residk.resid, '4W')
        self.assertEqual(crot.atomk, "CA")
        self.assertEqual(crot._yaml_header, 'crotations')
        self.assertEqual(crot._objcat, 'coord')
        self.maxDiff = None
        self.assertEqual(repr(crot), "Crot(angle='ANGLEIJK', segnamei='E', residi=ResID(resseqnum=2, insertion='Y'), atomi='NC', segnamejk='F', residj=ResID(resseqnum=3), atomj='CN', residk=ResID(resseqnum=4, insertion='W'), atomk='CA', degrees=90.0)")

    def test_crot_creation_alpha(self):
        crot = Crot(
            angle='ALPHA',
            chainID="G",
            resid1=ResID(8),
            resid2=ResID(12),
            resid3=ResID(16),
        )
        self.assertIsInstance(crot, Crot)
        self.assertEqual(crot.chainID, "G")
        self.assertEqual(crot.resid1.resid, 8)
        self.assertEqual(crot.resid2.resid, 12)
        self.assertEqual(crot._yaml_header, 'crotations')
        self.assertEqual(crot._objcat, 'coord')
        self.assertEqual(repr(crot), "Crot(angle='ALPHA', chainID='G', resid1=ResID(resseqnum=8), resid2=ResID(resseqnum=12), resid3=ResID(resseqnum=16))")

    def test_crot_creation_invalid_angle(self):
        with self.assertRaises(ValueError):
            Crot(
                angle='invalid',
                chainID="D",
                resid1=ResID(1),
                resid2=ResID(2),
                degrees=90.0
            )

    def test_crot_from_shortcode_phi(self):
        crot = Crot("PHI,A,5,15,180.0")
        self.assertIsInstance(crot, Crot)
        self.assertEqual(crot.angle, "PHI")
        self.assertEqual(crot.chainID, "A")
        self.assertEqual(crot.resid1.resid, 5)
        self.assertEqual(crot.resid2.resid, 15)

    def test_crot_from_shortcode_psi(self):
        crot = Crot("PSI,B,10,20,120.0")
        self.assertIsInstance(crot, Crot)
        self.assertEqual(crot.angle, "PSI")
        self.assertEqual(crot.chainID, "B")
        self.assertEqual(crot.resid1.resid, 10)
        self.assertEqual(crot.resid2.resid, 20)

    def test_crot_from_shortcode_omega(self):
        crot = Crot("OMEGA,C,15,25,240.0")
        self.assertIsInstance(crot, Crot)
        self.assertEqual(crot.angle, "OMEGA")
        self.assertEqual(crot.chainID, "C")
        self.assertEqual(crot.resid1.resid, 15)
        self.assertEqual(crot.resid2.resid, 25)

    def test_crot_from_shortcode_chi1(self):
        crot = Crot("CHI1,D,1X,60.0")
        self.assertIsInstance(crot, Crot)
        self.assertEqual(crot.angle, "CHI1")
        self.assertEqual(crot.chainID, "D")
        self.assertEqual(crot.resid1.resid, '1X')
        self.assertEqual(crot.degrees, 60.0)

    def test_crot_from_shortcode_angleijk(self):
        crot = Crot("ANGLEIJK,E,2Y,NC,F,3Z,CN,4W,CA,90.0")
        self.assertIsInstance(crot, Crot)
        self.assertEqual(crot.angle, "ANGLEIJK")
        self.assertEqual(crot.segnamei, "E")
        self.assertEqual(crot.residi.resid, '2Y')
        self.assertEqual(crot.segnamejk, "F")
        self.assertEqual(crot.residj.resid, '3Z')
        self.assertEqual(crot.residk.resid, '4W')
        self.assertEqual(crot.atomk, "CA")
        self.assertEqual(crot.degrees, 90.0)
    
    def test_crot_from_shortcode_alpha_full(self):
        crot = Crot("ALPHA,G,8,12,16")
        self.assertIsInstance(crot, Crot)
        self.assertEqual(crot.angle, "ALPHA")
        self.assertEqual(crot.chainID, "G")
        self.assertEqual(crot.resid1.resid, 8)
        self.assertEqual(crot.resid2.resid, 12)
        self.assertEqual(crot.resid3.resid, 16)

    def test_crot_from_shortcode_alpha_short(self):
        crot = Crot("ALPHA,G,8,12")
        self.assertIsInstance(crot, Crot)
        self.assertEqual(crot.angle, "ALPHA")
        self.assertEqual(crot.chainID, "G")
        self.assertEqual(crot.resid1.resid, 8)
        self.assertEqual(crot.resid2.resid, 12)
        self.assertEqual(crot.resid3.resid, 12)

    def test_crot_list_creation(self):
        crot_list = CrotList()
        self.assertIsInstance(crot_list, CrotList)
        self.assertEqual(len(crot_list), 0)