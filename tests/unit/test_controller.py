import unittest
from pestifer.controller import Controller
from pestifer.mods import *

class TestController(unittest.TestCase):
    def test_controller(self):
        c=Controller('user_config.yaml')
        self.assertEqual(c.config.defs['Title'],'Test User Config')
        self.assertEqual(c.config.defs['SourcePDB'],'4zmj')
        self.assertEqual(len(c.build_steps),1)
        s=c.build_steps[0]
        self.assertEqual(s.solvate,False)
        self.assertTrue('Mutations' in s.mods)
        self.assertEqual(len(s.mods['Mutations']),4)
        self.assertEqual(s.smdclose,True)
        self.assertEqual(s.relax_steps,[1000])
        m=s.mods['Mutations'][0]
        self.assertEqual(type(m),Mutation)
        self.assertEqual(m.chainID,'A')
        m=s.mods['Mutations'][1]
        self.assertEqual(m.origresname,'MET')
        self.assertTrue('Crotations' in s.mods)
        self.assertEqual(len(s.mods['Crotations']),5)
        c=s.mods['Crotations'][0]
        self.assertEqual(c.angle,'PHI')
        self.assertEqual(c.degrees,60)
        c=s.mods['Crotations'][-1]
        self.assertEqual(c.angle,'CHI2')
        self.assertEqual(c.degrees,75)
        self.assertTrue('SSBondsDelete' in s.mods)
        sd=s.mods['SSBondsDelete'][0]
        self.assertEqual(sd.chainID1,'A')
        self.assertEqual(sd.chainID2,'C')
        self.assertEqual(sd.resseqnum1,123)
        self.assertEqual(sd.resseqnum2,321)
        sd=s.mods['SSBondsDelete'][-1]
        self.assertEqual(sd.chainID1,'B')
        self.assertEqual(sd.chainID2,'D')
        self.assertEqual(sd.resseqnum1,345)
        self.assertEqual(sd.resseqnum2,543)

        # self.assertEqual(m['resseqnum'],151)
        # self.assertEqual(m['origresname'],'ALA')
        # self.assertEqual(m['newresname'],'CYS')
        # m=s.mods['Mutations'][1]
        # self.assertTrue('shortcode' in m)

