"""

.. module:: test_objlist
   :synopsis: tests of the objlist
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
from pestifer.baseobj import ObjList, BaseObj, CloneableObj, CloneableObjList
from pestifer.objs.mutation import MutationList,Mutation
import unittest
import random

class aObj(BaseObj):
    req_attr=['a','b']
    opt_attr=['c']
class bObj(BaseObj):
    req_attr=['x','y','A']

class aList(ObjList):
    pass

class bList(ObjList):
    pass

class cObj(CloneableObj):
    req_attr=['c','d']

class cList(CloneableObjList):
    pass

class TestObjList(unittest.TestCase):
    def test_init(self):
        L=aList([])
        self.assertEqual(len(L),0)
    def test_init_item(self):
        a=aObj({'a':1,'b':2})
        L=aList([a])
        self.assertEqual(len(L),1)
        for item in L:
            self.assertEqual(item.a,1)
    def test_extend(self):
        t=[]
        for i in range(20):
            t.append(aObj({'a':i,'b':25-i}))
        L=aList(t)
        u=[]
        for i in range(10):
            u.append(aObj({'a':i,'b':i**2}))
        L.extend(u)
        self.assertEqual(len(L),30)
    def test_getitem(self):
        L=aList([])
        for i in range(20):
            L.append(aObj({'a':i,'b':i**2}))
        item=L[5]
        self.assertEqual(item.a,5)
    def test_index(self):
        L=aList([])
        for i in range(20):
            L.append(aObj({'a':i,'b':i**2}))
        item=aObj({'a':5,'b':25})
        i=L.index(item)
        self.assertEqual(i,5)
    def test_get(self):
        L=aList([])
        for i in range(20):
            L.append(aObj({'a':i,'b':i**2}))
        item=L.get(a=5,b=25)
        self.assertEqual(L.index(item),5)    
    def test_getitem(self):
        L=aList([])
        for i in range(20):
            L.append(aObj({'a':i,'b':i**2}))
        self.assertEqual(L[5].a,5)
        self.assertEqual(L[5].b,25)
    def test_add(self):
        L=aList([])
        for i in range(20):
            L.append(aObj({'a':i,'b':i**2}))
        M=aList([])
        for i in range(20):
            M.append(aObj({'a':i,'b':i**2}))
        D=L+M
        self.assertEqual(len(D),40)
    def test_del(self):
        L=aList([])
        for i in range(20):
            L.append(aObj({'a':i,'b':i**2}))
        self.assertEqual(L[5].a,5)
        self.assertEqual(len(L),20)
        del L[5]
        self.assertEqual(L[5].a,6)
        self.assertEqual(len(L),19)
    def test_remove(self):
        L=aList([])
        for i in range(20):
            L.append(aObj({'a':i,'b':i**2}))
        self.assertEqual(L[5].a,5)
        self.assertEqual(len(L),20)
        L.remove(aObj({'a':5,'b':25}))
        self.assertEqual(L[5].a,6)
        self.assertEqual(len(L),19)
    def test_filter(self):
        L=aList([])
        for i in range(20):
            L.append(aObj({'a':i%4,'b':i**2}))
        subL=L.filter(a=0)
        self.assertEqual(len(subL),5)
    def test_nested(self):
        L=bList([])
        for i in range(20):
            A=aList([])
            for j in range(10):
                A.append(aObj({'a':j%3,'b':3*j,'c':'foundme'}))
            L.append(bObj({'x':i%5,'y':5*i,'A':A}))
        item=L.get(x=1,y=5)
        self.assertEqual(len(item.A),10)
        a_item=L.get_attr(('A',{'a':1,'b':3}),x=1,y=5)
        self.assertEqual(a_item.c,'foundme')

    def test_clone(self):
        L=cList([])
        for i in range(20):
            L.append(cObj({'c':i//2,'d':i**3}))
        M=L.clone()
        self.assertEqual(len(M),len(L))

    def test_clone_woptions(self):
        L=cList([])
        for i in range(20):
            L.append(cObj({'c':i//2,'d':i**3}))
        item=L[0]
        self.assertEqual(item.d,0)
        M=L.clone()#d='replace')
        self.assertEqual(len(M),len(L))
        # print(M[0].__dict__)
        item=M[0]
        self.assertEqual(item.d,0)

    def test_sort_strict1(self):
        L=cList([])
        I=list(range(20))
        random.shuffle(I)
        for i in I:
            L.append(cObj({'c':i,'d':(i)**3}))
        L.sort()
        self.assertTrue(all([x.c<=y.c for x,y in zip(L[:-1],L[1:])]))

    def test_sort_strict2(self):
        L=cList([])
        I=list(range(20,0,-1))
        # make an unsortable list because < always gives False
        for i in I:
            L.append(cObj({'c':i,'d':(19-i)**3}))
        L0_presort=L[0]
        L.sort() # will NOT sort
        L0=L[0]
        self.assertTrue(all([x.c==y for x,y in zip(L,I)]))

    def test_sort_nonstrict(self):
        L=cList([])
        I=list(range(20,0,-1))
        # make an strictly unsortable list because < always gives False
        for i in I:
            L.append(cObj({'c':i,'d':(19-i)**3}))
        L.sort(by=['c']) # will sort by 'c' attr only
        self.assertTrue(all([x.c<=y.c for x,y in zip(L[:-1],L[1:])]))
        # list will NOT be sorted by ['d']
        self.assertFalse(all([x.d<=y.d for x,y in zip(L[:-1],L[1:])]))
        L.sort(by=['d'])
        self.assertTrue(all([x.d<=y.d for x,y in zip(L[:-1],L[1:])]))
        self.assertFalse(all([x.c<=y.c for x,y in zip(L[:-1],L[1:])]))

    def test_shuffle(self):
        L=cList([])
        I=list(range(20))
        for i in I:
            L.append(cObj({'c':i,'d':(19-i)**3}))
        self.assertTrue(all([x.c<=y.c for x,y in zip(L[:-1],L[1:])]))
        random.shuffle(L)
        # can we assume the list must change order?
        self.assertFalse(all([x.c<=y.c for x,y in zip(L[:-1],L[1:])]))

class TestMutationList(unittest.TestCase):
    def test_mutation_list(self):
        L=MutationList([])
        c=L.filter(chainID='A')
        self.assertEqual(len(c),0)
        for i in range(10):
            L.append(Mutation(f'{chr(ord("A")+i)}:PHE,{123+i},TYR'))
        self.assertEqual(len(L),10)
        self.assertTrue(type(L),MutationList)
        m=L[0]
        self.assertEqual(m.chainID,'A')
        m=L[-1]
        self.assertEqual(m.chainID,'J')
        subl=L.filter(chainID='J')
        self.assertEqual(hasattr(subl,'__len__'),True)
        self.assertEqual(type(subl),MutationList)

