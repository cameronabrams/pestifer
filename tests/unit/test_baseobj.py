import unittest
import pytest
from pestifer.baseobj import BaseObj, CloneableObj, AncestorAwareObj, ObjList, StateInterval, StateIntervalList
from argparse import Namespace

class TestBaseObj(unittest.TestCase):
    def test_baseobj_empty(self):
        a=BaseObj({})
        self.assertEqual(type(a),BaseObj)
    def test_baseobj_undifferentiated_fromdict(self):
        b=BaseObj({'a':1,'b':2})
        self.assertEqual(type(b),BaseObj)
        self.assertEqual(b.a,1)
        self.assertEqual(b.b,2)
    def test_baseobj_fromnamespace(self):
        z=Namespace(a=1,b=2)
        b=BaseObj(z)
        self.assertEqual(type(b),BaseObj)
        self.assertEqual(b.a,1)
        self.assertEqual(b.b,2)

    def test_baseobj_catch_mutuals_notrequired(self):
        # required attributes are not in mutual exclusives
        class tbm(BaseObj):
            req_attr=['R1','R2','R3']
            opt_attr=['O1','O2']
            alt_attr=[('R1','O2')] # wrong!
        with pytest.raises(AssertionError) as E:
            a=tbm({'R1':1,'R2':2,'R3':3})
            print(a.__dict__)
        self.assertEqual(E.type,AssertionError)

    def test_baseobj_catch_required(self):
        # all required attributes have values
        class tbm(BaseObj):
            req_attr=['R1','R2','R3']
        with pytest.raises(AssertionError) as E:
            a=tbm({'R1':1,'R2':2})
        self.assertEqual(E.type,AssertionError)

    def test_baseobj_mutuals_satisfied(self):
        # all mutual exclusivity requirements are met
        class tbm(BaseObj):
            req_attr=['R1','R2','R3']
            opt_attr=['O1','O2']
            alt_attr=[('O1','O2')]
        with pytest.raises(AssertionError) as E:
            a=tbm({'R1':1,'R2':2,'R3':3,'O1':4,'O2':5})
        self.assertEqual(E.type,AssertionError)

    def test_baseobj_choices_met(self):
        # all attributes with limited set of possible values have valid values
        class tbm(BaseObj):
            req_attr=['R1','R2','R3']
            opt_attr=['O1','O2']
            alt_attr=[('O1','O2')]
            attr_choices={'R1':[7,8,9]}
        with pytest.raises(AssertionError) as E:
            a=tbm({'R1':1,'R2':2,'R3':3})
        self.assertEqual(E.type,AssertionError)

    def test_baseobj_optional_dependents_met(self):
        # for each optional attribute present, all required dependent attributes are also present
        class tbm(BaseObj):
            req_attr=['R1','R2','R3']
            opt_attr=['O1','O2']
            opt_attr_deps={'O1':['O2']}
        with pytest.raises(AssertionError) as E:
            a=tbm({'R1':1,'R2':2,'R3':3,'O1':4})
        self.assertEqual(E.type,AssertionError)

    def test_baseobj_equality(self):
        class tbm(BaseObj):
            req_attr=['a','b']
        id1={'a':1,'b':2}
        bm1=tbm(id1)
        bm2=BaseObj({})
        self.assertFalse(bm1==bm2)
        bm2=tbm(id1)
        self.assertTrue(bm1==bm2)

    def test_baseobj_relational(self):
        class tbm(BaseObj):
            req_attr=['a','b']
        bm1=tbm({'a':1,'b':2})
        bm2=tbm({'a':0,'b':2})
        self.assertTrue(bm2<bm1)
        bm1=tbm({'a':1,'b':2})
        bm2=tbm({'a':0,'b':3})
        self.assertFalse(bm2<bm1)
        self.assertTrue(bm2.weak_lt(bm1,['a']))
        self.assertFalse(bm2.weak_lt(bm1,['b']))
        self.assertFalse(bm2.weak_lt(bm1,[]))
        bm2=tbm({'a':1,'b':2})
        self.assertTrue(bm1==bm2)
        self.assertFalse(bm2<bm1)
        self.assertFalse(bm1<bm2)
        class xbm(BaseObj):
            req_attr=['a','b']
        bm3=xbm({'a':1,'b':2})
        self.assertFalse(bm3==bm1)
        self.assertFalse(bm3<bm1)
        self.assertFalse(bm1<bm3)

    def test_baseobj_set(self):
        class tbm(BaseObj):
            req_attr=['a','x']
        class sbm(BaseObj):
            req_attr=['b','x']
        L=ObjList([sbm({'b':x,'x':x+1}) for x in range(1,10)])
        A=tbm({'a':1,'x':2})
        A.sub=sbm({'b':'sub','x':99})
        A.sublist=L
        A.set(x=77)
        self.assertEqual(A.x,77)
        self.assertEqual(A.sub.x,77)
        self.assertTrue(all([item.x==77 for item in A.sublist]))
        L=ObjList([sbm({'b':x,'x':x+1}) for x in range(1,10)])
        A=tbm({'a':1,'x':2})
        A.sub=sbm({'b':'sub','x':99})
        A.sublist=L
        A.set(shallow=True,x=88)
        self.assertEqual(A.x,88)
        self.assertEqual(A.sub.x,99) # not set
        self.assertEqual([item.x for item in A.sublist],list(range(2,11)))
    
    def test_allneg(self):
        class tbm(BaseObj):
            req_attr=['a','x']
        i=tbm({'a':1,'x':2})
        # j=tbm({'a':2,'x':1})
        self.assertTrue(i.allneg(a=3,x=1))
        self.assertFalse(i.allneg(a=3,x=2))
        L=ObjList([tbm({'a':'H','x':x+1}) for x in range(1,10)])+ObjList([tbm({'a':'C','x':x+1}) for x in range(1,10)])
        M=L.negfilter(a='H')
        self.assertEqual(len(L),18)
        self.assertEqual(len(M),9)
        self.assertEqual(M[0].a,'C')
        self.assertTrue(all([x.a=='C' for x in M]))
        
    def test_swap_attr(self):
        class tbm(BaseObj):
            req_attr=['a','x']
        i=tbm({'a':1,'x':2})
        i.swap_attr('a','x')
        self.assertEqual(i.a,2)
        self.assertEqual(i.x,1)
        self.assertEqual(i.swapped['a'],'x')

    def test_copy_attr(self):
        class tbm(BaseObj):
            req_attr=['a','x']
        i=tbm({'a':1,'x':2})
        i.copy_attr('a','x')
        self.assertEqual(i.a,2)
        self.assertEqual(i.x,2)
        self.assertEqual(i.copied_over['a'],1)

class TestStateInterval(unittest.TestCase):
    def test_stateintervals(self):
        s=StateInterval({'state':'no state','bounds':[0,0]})
        self.assertEqual(s.state,'no state')
        self.assertEqual(s.bounds,[0,0])
        t=StateInterval({'state':'sleeping','bounds':[0,10]})
        self.assertEqual(t.state,'sleeping')
        self.assertEqual(t.bounds,[0,10])
        self.assertEqual(len(StateInterval.req_attr),2)

class TestCloneableObj(unittest.TestCase):
    def test_cloneablemod(self):
        class MyCloneable(CloneableObj):
            req_attr=CloneableObj.req_attr.copy()
            req_attr.extend(['a','b'])
        c1=MyCloneable({'a':1,'b':2})
        self.assertEqual(len(MyCloneable.req_attr),2)
        self.assertEqual(len(MyCloneable.opt_attr),1)
        self.assertTrue('clone_of' in MyCloneable.opt_attr)
        self.assertFalse('clone_of' in c1.__dict__)
        self.assertEqual(c1.a,1)
        self.assertEqual(c1.b,2)
        c2=c1.clone()
        self.assertEqual(c2.a,1)
        self.assertEqual(c2.b,2)
        self.assertTrue(c2.is_clone())
        self.assertEqual(c2.clone_of,c1)
        self.assertEqual(c2.get_original(),c1)
        self.assertFalse(c1.is_clone())
        c3=c1.clone(a=10)
        self.assertEqual(c1.a,1)
        self.assertEqual(c3.a,10)
        self.assertEqual(c3.b,2)
        self.assertTrue(c3.is_clone())
        self.assertEqual(c3.clone_of,c1)
        self.assertEqual(c3.get_original(),c1)

class TestAncestorAwareObj(unittest.TestCase):
    def test_ancestorawaremod_init(self):
        self.assertEqual(AncestorAwareObj.req_attr,[])
        self.assertEqual(AncestorAwareObj.opt_attr,['ancestor_obj'])
        class MyAncestorAwareObj(AncestorAwareObj):
            req_attr=AncestorAwareObj.req_attr.copy()
            req_attr.extend(['a','b'])
        self.assertEqual(MyAncestorAwareObj.req_attr,['a','b'])
        self.assertEqual(MyAncestorAwareObj.opt_attr,['ancestor_obj'])
        a=MyAncestorAwareObj({'a':1,'b':2})
        self.assertTrue('a' in a.__dict__)
        self.assertFalse('ancestor_obj' in a.__dict__) # it is optional
        self.assertTrue('ancestor_obj' in a.__class__.opt_attr)
        b=MyAncestorAwareObj({'a':3,'b':4,'ancestor_obj':a})
        self.assertEqual(b.ancestor_obj,a)
        b=MyAncestorAwareObj({'a':7,'b':5})

        class MyNewAncestorAwareObj(AncestorAwareObj):
            req_attr=AncestorAwareObj.req_attr.copy()
            req_attr.extend(['obj1','obj2'])

        self.assertEqual(MyNewAncestorAwareObj.req_attr,['obj1','obj2'])
        self.assertEqual(MyNewAncestorAwareObj.opt_attr,['ancestor_obj'])
        
        c=MyNewAncestorAwareObj({'obj1':a,'obj2':b})
        self.assertEqual(c.obj1,a)
        self.assertEqual(c.obj2,b)
        self.assertFalse('ancestor_obj' in a.__dict__) # no one has claimed him yet!
        self.assertFalse('ancestor_obj' in b.__dict__) # no one has claimed her yet!
        c.claim_descendants(c)
        self.assertEqual(a.ancestor_obj,c)
        self.assertEqual(b.ancestor_obj,c)

class TestBaseObjList(unittest.TestCase):

    def test_uniquify(self):
        L=ObjList([])
        class tbm(BaseObj):
            req_attr=['a','b'] # by default, BaseObj.req_attr is [], so this is ok
        a_mod=tbm({'a':1,'b':1})
        L.append(a_mod) # 0
        L.append(tbm({'a':1,'b':1})) # 1
        L.append(tbm({'a':2,'b':1})) # 2
        L.append(tbm({'a':3,'b':1})) # 3
        L.append(tbm({'a':4,'b':1})) # 4
        L.append(tbm({'a':4,'b':1})) # 5
        L.append(tbm({'a':4,'b':1})) # 6
        L.append(tbm({'a':1,'b':1})) # 7
        L.puniquify(['a'])
        self.assertEqual(L[0].a,1)
        self.assertTrue(not hasattr(L[0],'_ORIGINAL_'))
        self.assertEqual(L[1].a,5)
        self.assertTrue(hasattr(L[1],'_ORIGINAL_'))
        self.assertEqual(L[1]._ORIGINAL_['a'],1)        
        self.assertEqual(L[5].a,6)
        self.assertTrue(hasattr(L[5],'_ORIGINAL_'))
        self.assertEqual(L[5]._ORIGINAL_['a'],4)        
        self.assertEqual(L[6].a,7)
        self.assertTrue(hasattr(L[6],'_ORIGINAL_'))
        self.assertEqual(L[6]._ORIGINAL_['a'],4)        
        self.assertEqual(L[7].a,8)
        self.assertTrue(hasattr(L[7],'_ORIGINAL_'))
        self.assertEqual(L[7]._ORIGINAL_['a'],1)
        self.assertTrue(a_mod in L)

    def test_uniquify_commonize(self):
        L=ObjList([])
        class tbm(BaseObj):
            req_attr=['a','b']# by default, BaseObj.req_attr is [], so this is ok
        L.append(tbm({'a':1,'b':0})) # 0
        L.append(tbm({'a':1,'b':2})) # 1
        L.append(tbm({'a':2,'b':4})) # 2
        L.append(tbm({'a':3,'b':6})) # 3
        L.append(tbm({'a':4,'b':7})) # 4
        L.append(tbm({'a':4,'b':9})) # 5
        L.append(tbm({'a':4,'b':11})) # 6
        L.append(tbm({'a':1,'b':13})) # 7
        L.puniquify(['a'],make_common=['b'])
        self.assertTrue(all([x.b==0 for x in L]))
        self.assertEqual(L[1]._ORIGINAL_['b'],2)
        self.assertEqual(L[2]._ORIGINAL_['b'],4)
        self.assertEqual(L[3]._ORIGINAL_['b'],6)

    def test_map_attr(self):
        map={
            1:'okImappedyou.',
            2:'okImappedyou..',
            3:'okImappedyou...',
            5:'okImappedyou......',
            7:'okImappedyou.......',
            9:'okImappedyou........'
        }
        class tbm(BaseObj):
            req_attr=['a','b','c']
        L=ObjList([])
        L.append(tbm({'a':1,'b':2,'c':'mapme1'}))
        L.append(tbm({'a':2,'b':2,'c':'mapme3'}))
        L.append(tbm({'a':3,'b':4,'c':'mapme7'}))
        L.append(tbm({'a':5,'b':2,'c':'mapme8'}))
        L.append(tbm({'a':7,'b':6,'c':'mapme2'}))
        L.append(tbm({'a':9,'b':2,'c':'mapme1'}))
        L.map_attr('c','a',map)
        self.assertEqual(L[0].c,'okImappedyou.')


    def test_state_interval_computation(self):
        L=ObjList([])
        class tbm(BaseObj):
            req_attr=['a','b'] # by default, BaseObj.req_attr is [], so this is ok
        L.append(tbm({'a':[],'b':1})) # 0
        L.append(tbm({'a':[],'b':1})) # 1
        L.append(tbm({'a':[2,3,4],'b':1})) # 2
        L.append(tbm({'a':[5,6],'b':1})) # 3
        L.append(tbm({'a':[7,8,9,10],'b':1})) # 4
        L.append(tbm({'a':[],'b':1})) # 5
        L.append(tbm({'a':[11,12],'b':1})) # 6
        L.append(tbm({'a':[13],'b':1})) # 7
        b=L.state_bounds(lambda x: 'RESOLVED' if len(x.a)>0 else 'MISSING')
        self.assertEqual(type(b),StateIntervalList)
        self.assertEqual(b[0],StateInterval({'state':'MISSING','bounds':[0,1]}))    
        self.assertEqual(b[1],StateInterval({'state':'RESOLVED','bounds':[2,4]}))   
        self.assertEqual(b[2],StateInterval({'state':'MISSING','bounds':[5,5]}))
        self.assertEqual(b[3],StateInterval({'state':'RESOLVED','bounds':[6,7]}))
        r=[x.pstr() for x in b]
        self.assertEqual(r,['MISSING(2)','RESOLVED(3)','MISSING(1)','RESOLVED(2)'])
        c=b.clone()
        self.assertTrue(c is not b)
        self.assertEqual(type(c),StateIntervalList)
        self.assertEqual(c[0],StateInterval({'state':'MISSING','bounds':[0,1]}))    
        self.assertEqual(c[1],StateInterval({'state':'RESOLVED','bounds':[2,4]}))   
        self.assertEqual(c[2],StateInterval({'state':'MISSING','bounds':[5,5]}))
        self.assertEqual(c[3],StateInterval({'state':'RESOLVED','bounds':[6,7]}))
        r=[x.pstr() for x in b]
        self.assertEqual(r,['MISSING(2)','RESOLVED(3)','MISSING(1)','RESOLVED(2)'])

    # def test_ancestor_aware(self):
    #     class tam(AncestorAwareObj):
    #         def __init__(self,input_dict):
    #             tam.req_attr.append('a','b')
    #             super().__init__(input_dict)
