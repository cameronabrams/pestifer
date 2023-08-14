import unittest
from pestifer.controller import Controller
from pestifer.mods import *

class TestController(unittest.TestCase):
    def test_controller_base(self):
        C=Controller('user_config.yaml').build_molecules()
        self.assertEqual(C.config.defs['Title'],'Test User Config')
        self.assertEqual(C.config.defs['Source']['pdb'],'4zmj')
        self.assertEqual(len(C.build_steps),1)
        s=C.build_steps[0]
        self.assertEqual(s.solvate,False)
        self.assertTrue('Mutations' in s.mods)
        self.assertEqual(len(s.mods['Mutations']),4)
        self.assertTrue(type(s.mods['Mutations']),MutationList)
        sublist=s.mods.get('Mutations',MutationList([]))
        for m in sublist:
            print(str(m))
        res=sublist.filter(chainID='A')
        self.assertEqual(len(res),2)
        self.assertTrue(hasattr(res,'__len__'))
        self.assertEqual(type(sublist),MutationList)
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

        baseMol=C.molecules['4zmj']
        self.assertEqual(baseMol.molid,1)
        self.assertEqual(baseMol.active_biological_assembly.index,1)
        b=baseMol.active_biological_assembly
        self.assertEqual(b.index,1)
        self.assertEqual(len(b.biomt),3)
        t1=b.biomt[0]
        self.assertEqual(t1.chainIDmap,{'A': 'A', 'B': 'B', 'C': 'C', 'D': 'D', 'G': 'G'})
        t2=b.biomt[1]
        self.assertEqual(t2.chainIDmap,{'A': 'E', 'B': 'F', 'C': 'H', 'D': 'I', 'G': 'J'})
        t3=b.biomt[2]
        self.assertEqual(t3.chainIDmap,{'A': 'K', 'B': 'L', 'C': 'M', 'D': 'N', 'G': 'O'})
        # with open('molout.tcl','w') as f:
        #     f.write(baseMol.write_TcL())
        # self.assertEqual(m['resseqnum'],151)
        # self.assertEqual(m['origresname'],'ALA')
        # self.assertEqual(m['newresname'],'CYS')
        # m=s.mods['Mutations'][1]
        # self.assertTrue('shortcode' in m)
    def test_controller_do(self):
        C=Controller('user_config.yaml')
        C.do(clean_up=True)
