# Author: Cameorn F. Abrams, <cfa22@drexel.edu>

from __future__ import annotations
import unittest
from pestifer.core.baseobj import BaseObj, BaseObjList
from pestifer.objs.resid import ResID
from argparse import Namespace
from pydantic import Field, ValidationError, ConfigDict
from typing import Optional, List, Dict, Annotated
import logging

logger = logging.getLogger(__name__)
class TestBaseObj(unittest.TestCase):
    def test_baseobj_initialization(self):
        # empty instantiation is ok
        a = BaseObj()
        self.assertIsInstance(a, BaseObj)
        self.assertEqual(a._required_fields, set())
        self.assertEqual(a._optional_fields, set())
        self.assertEqual(a._mutually_exclusive, set())
        self.assertEqual(a._attr_choices, {})
        self.assertEqual(repr(a), "BaseObj()")
        
        # Missing required fields should raise ValidationError
        with self.assertRaises(ValidationError):
            b = BaseObj(name='Test Object')

        with self.assertRaises(ValidationError):
            c = BaseObj(Namespace(name='Namespace Object'))

        with self.assertRaises(ValidationError):
            d = BaseObj(name='Direct Object')

    def test_baseobj_subclassing(self):
        # even though a field is declared, it must be in _required_fields or _optional_fields
        class ConcreteObj(BaseObj):
            name: str = Field(..., description="Name of the object")
        with self.assertRaises(ValidationError):
            obj = ConcreteObj(name='Test Object')
        
        class ConcreteObj(BaseObj):
            _required_fields = {'name'}
            name: str = Field(..., description="Name of the object")

        obj = ConcreteObj(name='Test Object')
        self.assertIsInstance(obj, ConcreteObj)
        self.assertEqual(repr(obj), "ConcreteObj(name='Test Object')")

        with self.assertRaises(ValidationError):
            # Missing required 'name' field
            bad_obj = ConcreteObj(size=10)

        class ConcreteObj(BaseObj):
            _optional_fields = {'description'}
            description: str | None = Field(None, description="Description of the object")

        # empty obj is ok since there are no required fields
        obj = ConcreteObj()
        self.assertIsInstance(obj, ConcreteObj)
        # optional fields with None defaults and having None values are not reported
        self.assertEqual(repr(obj), "ConcreteObj()")
        obj = ConcreteObj(description="This is a test object")

    def test_baseobj_subclass_with_mutually_exclusive_fields(self):
        class ConcreteObj(BaseObj):
            _optional_fields = {'field_a', 'field_b'}
            _mutually_exclusive = {frozenset({'field_a', 'field_b'})}
            field_a: str | None = Field(None, description="Field A")
            field_b: str | None = Field(None, description="Field B")

        obj = ConcreteObj(field_a="A")
        self.assertEqual(repr(obj), "ConcreteObj(field_a='A')")
        obj2 = ConcreteObj(field_b="B")
        self.assertEqual(repr(obj2), "ConcreteObj(field_b='B')")
        with self.assertRaises(ValidationError):
            # Both fields set, should raise error
            bad_obj = ConcreteObj(field_a="A", field_b="B")

    def test_baseobj_subclass_with_attr_choices(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'status'}
            _attr_choices = {'status': {'active', 'inactive'}}
            status: str = Field(..., description="Status of the object")

        obj = ConcreteObj(status="active")
        self.assertEqual(repr(obj), "ConcreteObj(status='active')")
        with self.assertRaises(ValidationError):
            # Invalid status
            bad_obj = ConcreteObj(status="unknown")

    def test_baseobj_subclass_with_attr_dependencies(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'required_field'}
            _optional_fields = {'optional_field'}
            _attr_dependencies = {'required_field': {'Trigger me to require optional field': {'optional_field'}}}
            required_field: str = Field(..., description="Required field")
            optional_field: Optional[str] = Field(None, description="Optional field")

        # with self.assertRaises(ValidationError) as context:
        #     ConcreteObj(required_field="Must be set")
        #     logger.debug(context.exception)
        obj = ConcreteObj(required_field="Must be set")
        self.assertEqual(repr(obj), "ConcreteObj(required_field='Must be set')")
        with self.assertRaises(ValidationError):
            obj2 = ConcreteObj(required_field="Trigger me to require optional field")
        obj2 = ConcreteObj(required_field="Trigger me to require optional field", optional_field="Here I am")
        self.assertEqual(repr(obj2), "ConcreteObj(required_field='Trigger me to require optional field', optional_field='Here I am')")
        with self.assertRaises(ValidationError):
            # Missing required 'required_field'
            bad_obj = ConcreteObj(optional_field="Should not work")

    def test_baseobj_subclass_comparison(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'name', 'number'}
            name: str = Field(..., description="Name of the object")
            number: int = Field(..., description="A number")

        obj1 = ConcreteObj(name="Object 1", number=1)
        obj2 = ConcreteObj(name="Object 2", number=1)
        self.assertTrue(obj1 < obj2)  # Assuming string comparison on name
        self.assertTrue(obj1 <= obj2)
        self.assertFalse(obj2 < obj1)
        obj1 = ConcreteObj(name="Object 1", number=1)
        obj2 = ConcreteObj(name="Object 1", number=1)
        self.assertTrue(obj1 == obj2)  # Assuming string comparison on name
        self.assertTrue(obj1 <= obj2)
        self.assertFalse(obj2 < obj1)
        obj1 = ConcreteObj(name="Object 1", number=2)
        obj2 = ConcreteObj(name="Object 1", number=1)
        self.assertFalse(obj1 < obj2)  # Different numbers
        self.assertFalse(obj1 <= obj2)
        self.assertTrue(obj2 < obj1)

    def test_baseobj_subclass_strhash(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'name', 'number'}
            name: str = Field(..., description="Name of the object")
            number: int = Field(..., description="A number")

        obj = ConcreteObj(name="Object 1", number=1)
        self.assertEqual(obj.strhash(), "name=Object 1-number=1")
        obj2 = ConcreteObj(name="Object 2", number=2)
        self.assertNotEqual(obj.strhash(), obj2.strhash())
        self.assertEqual(obj.strhash(fields=['name']), "name=Object 1")
        self.assertEqual(obj.strhash(fields=['number']), "number=1")
        self.assertEqual(obj.strhash(fields=['name', 'number']), "name=Object 1-number=1")
        self.assertEqual(obj.strhash(digest="md5"), "8def7f1adbc0ea4a563ff5ac2ca02bd8")
        self.assertEqual(obj.strhash(digest="sha256"), "260bc656cb73a249da50424002251df5050ce2769d11771b21ca7a3513112faa")

    def test_baseobj_subclass_matches(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'name', 'number'}
            name: str = Field(..., description="Name of the object")
            number: int = Field(..., description="A number")

        obj1 = ConcreteObj(name="Object 1", number=1)
        # wildmatch tests if key NAMES are substrings of the model attribute names
        self.assertTrue(obj1.wildmatch(nam="Object 1", umber=1))
        obj2 = ConcreteObj(name="Object 2", number=2)
        self.assertTrue(obj1.weak_lt(obj2, attr=['name']))

    def test_baseobj_subclass_nested_wildmatch(self):
        class SubObj(BaseObj):
            _required_fields = {'sub_name'}
            sub_name: str = Field(..., description="Name of the sub-object")

        class RestrictiveObj(SubObj):
            def __eq__(self, other: RestrictiveObj | dict) -> bool:
                return self.sub_name == other.sub_name

        class SuperObj(BaseObj):
            _required_fields = {'super_name', 'sub_object', 'rsub_object'}
            super_name: str = Field(..., description="Name of the super-object")
            sub_object: SubObj = Field(..., description="A sub-object")
            rsub_object: RestrictiveObj = Field(..., description="A restrictive sub-object")

        Sub = SubObj(sub_name="Sub Object 1")
        OtherSub = SubObj(sub_name="Sub Object 1")
        ResSub = RestrictiveObj(sub_name="Sub Object 1")
        self.assertTrue(Sub == OtherSub)
        self.assertTrue(Sub == ResSub)
        Sup = SuperObj(super_name="Super Object 1", sub_object=Sub, rsub_object=ResSub)
        
        # wildmatch means a keyword can be a substring of the model attribute name
        # here we use "sub_obj" instead of "sub_object"
        # using a dict version of a subobject should work thanks pydantic's full serialization of Sup
        self.assertTrue(Sup.wildmatch(sub_obj={"sub_name": "Sub Object 1"}))
        # wildmatch will still work since BaseObj __eq__ allows for other to be a dict serialization or a BaseObj instance
        self.assertTrue(Sup.wildmatch(sub_obj=OtherSub))
        # returns false -- no key matching
        self.assertFalse(Sup.wildmatch(sub_xxx="Sub Object 1"))

        # returns true, since we are using a dict and model_dump() as used
        # in wildmatch() serializes the caller
        self.assertTrue(Sup.wildmatch(rsub_obj={"sub_name": "Sub Object 1"}))

        with self.assertRaises(AttributeError):
            # fails because RestrictiveObj overrides __eq__ and only permits self to be a RestrictiveObj 
            self.assertTrue(Sup.wildmatch(rsub_obj=ResSub))

    def test_baseobj_subclass_yamldump(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'name', 'number'}
            name: str = Field(..., description="Name of the object")
            number: int = Field(..., description="A number")

        obj = ConcreteObj(name="Object 1", number=1)
        self.assertEqual(obj.dump(), "name: Object 1\nnumber: 1\n")
        obj2 = ConcreteObj(name="Object 2", number=2)
        self.assertNotEqual(obj.dump(), obj2.dump())

    def test_baseobj_subclass_inlist(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'name', 'number'}
            name: str = Field(..., description="Name of the object")
            number: int = Field(..., description="A number")

        obj1 = ConcreteObj(name="Object 1", number=1)
        obj2 = ConcreteObj(name="Object 2", number=2)
        obj_list = [obj1, obj2]
        self.assertIn(obj1, obj_list)
        self.assertIn(obj2, obj_list)
        self.assertNotIn(ConcreteObj(name="Object 3", number=3), obj_list)

    def test_baseobj_subclass_map_attr(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'name', 'number'}
            name: str = Field(..., description="Name of the object")
            number: int = Field(..., description="A number")

        obj = ConcreteObj(name="Object 1", number=2)
        obj.map_attr('number','name', {'Object 1':99})
        self.assertEqual(obj.number, 99)

    def test_baseobj_subclass_swap_attributes(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'name', 'number'}
            name: str = Field(..., description="Name of the object")
            number: int = Field(..., description="A number")

        obj = ConcreteObj(name="Object 1", number=2)
        obj.swap_attr('number','name')
        self.assertEqual(obj.name, 2)
        self.assertEqual(obj.number, "Object 1")

    def test_baseobj_subclass_copy_attributes(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'name', 'number'}
            name: str = Field(..., description="Name of the object")
            number: int = Field(..., description="A number")

        obj1 = ConcreteObj(name="Object 1", number=2)
        obj1.copy_attr('name', 'number')
        self.assertEqual(obj1.name, 2)
        self.assertEqual(obj1.number, 2)

    def test_baseobj_subclass_nested(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'name', 'number'}
            _optional_fields = {'spouse', 'children', 'dict_of_pets'}
            name: str = Field(..., description="Name of the object")
            number: int = Field(..., description="A number")
            spouse: Optional["ConcreteObj"] = None
            children: Optional[List["ConcreteObj"]] = Field(default_factory=list, description="Child objects")
            dict_of_pets: Optional[Dict[str, "ConcreteObj"]] = Field(default_factory=dict, description="Dictionary of pet objects")

        grandkids = [ConcreteObj(name="GrandKid 1", number=5), ConcreteObj(name="GrandKid 2", number=6)]
        stepmom = ConcreteObj(name="StepMom", number=1)
        child1 = ConcreteObj(name="Child 1", number=3, children=grandkids)
        child2 = ConcreteObj(name="Child 2", number=4)
        dog = ConcreteObj(name="Dog", number=1)
        cat = ConcreteObj(name="Cat", number=2)
        dad = ConcreteObj(name="Daddy", number=2, spouse=stepmom, children=[child1, child2], dict_of_pets={"older": dog, "younger": cat})

        self.assertEqual(type(dad), ConcreteObj)
        self.assertEqual(type(dad.spouse), ConcreteObj)

        dad.set(name="Daddy Updated",shallow=True)
        self.assertEqual(dad.name, "Daddy Updated")
        self.assertEqual(dad.spouse.name, "StepMom")
        self.assertEqual(dad.children[0].name, "Child 1")
        self.assertEqual(dad.children[1].name, "Child 2")
        self.assertEqual(dad.dict_of_pets["older"].name, "Dog")
        self.assertEqual(dad.dict_of_pets["younger"].name, "Cat")
        self.assertEqual(dad.children[0].children[0].name, "GrandKid 1")
        self.assertEqual(dad.children[0].children[1].name, "GrandKid 2")
        # change everyone's name to Gumby recursively (default behavior of set())
        dad.set(name="Gumby")
        self.assertEqual(dad.name, "Gumby")
        self.assertEqual(dad.spouse.name, "Gumby")
        self.assertEqual(dad.children[0].name, "Gumby")
        self.assertEqual(dad.children[1].name, "Gumby")
        self.assertEqual(dad.dict_of_pets["older"].name, "Gumby")
        self.assertEqual(dad.dict_of_pets["younger"].name, "Gumby")
        self.assertEqual(dad.children[0].children[0].name, "Gumby")
        self.assertEqual(dad.children[0].children[1].name, "Gumby")

    def test_baseobj_subclass_update_attr_from_obj_attr(self):
        class BasicObj(BaseObj):
            _required_fields = {'name'}
            name: str = Field(..., description="Name of the object")
            
        class CompoundObj(BaseObj):
            _required_fields = {'name', 'number', 'subobject'}
            name: str = Field(..., description="Name of the object")
            number: int = Field(..., description="A number")
            subobject: BasicObj = Field(..., description="A subobject")

        obj1 = BasicObj(name="Howdy Doody")
        obj2 = CompoundObj(name="Object 2", number=3, subobject=obj1)
        obj2.update_attr_from_obj_attr('name', 'subobject', 'name')
        self.assertEqual(obj2.name, "Howdy Doody")
        self.assertEqual(obj2.number, 3)

    def test_baseobj_subclass_update_attr_from_objlist_elem_attr(self):
        class BasicObj(BaseObj):
            _required_fields = {'name'}
            name: str = Field(..., description="Name of the object")

        class CompoundObj(BaseObj):
            _required_fields = {'name', 'number', 'subobject_list'}
            name: str = Field(..., description="Name of the object")
            number: int = Field(..., description="A number")
            subobject_list: List[BasicObj] = Field(..., description="A list of subobjects")

        obj1 = BasicObj(name="Howdy Doody")
        obj3 = BasicObj(name="Bumgy")
        obj2 = CompoundObj(name="Object 2", number=3, subobject_list=[obj1,obj3])
        obj2.update_attr_from_objlist_elem_attr('name', 'subobject_list', 0, 'name')
        self.assertEqual(obj2.name, "Howdy Doody")
        obj2.update_attr_from_objlist_elem_attr('name', 'subobject_list', 1, 'name')
        self.assertEqual(obj2.name, "Bumgy")

    def test_baseobj_subclass_new_adapt(self):
        class DataObj:
            def __init__(self, **kwargs):
                self.__dict__.update(kwargs)

        class ConcreteObj(BaseObj):
            _required_fields = {'name', 'number'}
            name: str = Field(..., description="Name of the object")
            number: int = Field(..., description="A number")
            @classmethod
            def _adapt(cls, *args, **kwargs) -> dict:
                if args:
                    if isinstance(args[0], DataObj):
                        return dict(name=args[0].name, number=args[0].number)
                return super()._adapt(*args, **kwargs)

        d = DataObj(name='Data Name', number=777)
        c = ConcreteObj(d)
        self.assertEqual(c.name, 'Data Name')
        self.assertEqual(c.number, 777)

class TestBaseObjList(unittest.TestCase):

    def test_baseobjlist_is_abstract(self):
        with self.assertRaises(TypeError):
            AbstractList = BaseObjList()

    def test_baseobjlist_subclass(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'name'}
            name: str = Field(..., description="Name of the object")

        class ConcreteObjList(BaseObjList[ConcreteObj]):
            def describe(self) -> str:
                return f"Concrete Object List with {len(self)} items."

        obj1 = ConcreteObj(name="Object 1")
        obj2 = ConcreteObj(name="Object 2")
        obj_list = ConcreteObjList([obj1, obj2])
        self.assertEqual(len(obj_list), 2)
        self.assertIn(obj1, obj_list)
        self.assertIn(obj2, obj_list)
        self.assertEqual(obj_list[0].name, "Object 1")
        self.assertEqual(obj_list[1].name, "Object 2")

    def test_baseobj_list_append(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'name'}
            name: str = Field(..., description="Name of the object")

        class ConcreteObjList(BaseObjList[ConcreteObj]):
            def describe(self) -> str:
                return f"Concrete Object List with {len(self)} items."
        
        obj1 = ConcreteObj(name="Object 1")
        obj2 = ConcreteObj(name="Object 2")
        obj_list = ConcreteObjList([obj1])
        self.assertEqual(len(obj_list), 1)
        self.assertIn(obj1, obj_list)
        obj_list.append(obj2)
        self.assertEqual(len(obj_list), 2)
        self.assertIn(obj2, obj_list)

    def test_baseobj_list_extend(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'name'}
            name: str = Field(..., description="Name of the object")

        class ConcreteObjList(BaseObjList[ConcreteObj]):
            def describe(self) -> str:
                return f"Concrete Object List with {len(self)} items."
        
        obj1 = ConcreteObj(name="Object 1")
        obj2 = ConcreteObj(name="Object 2")
        obj3 = ConcreteObj(name="Object 3")
        obj_list = ConcreteObjList([obj1])
        self.assertEqual(len(obj_list), 1)
        self.assertIn(obj1, obj_list)
        obj_list.extend([obj2, obj3])
        self.assertEqual(len(obj_list), 3)
        self.assertIn(obj2, obj_list)
        self.assertIn(obj3, obj_list)

    def test_baseobj_list_insert(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'name'}
            name: str = Field(..., description="Name of the object")

        class ConcreteObjList(BaseObjList[ConcreteObj]):
            def describe(self) -> str:
                return f"Concrete Object List with {len(self)} items."
        
        obj1 = ConcreteObj(name="Object 1")
        obj2 = ConcreteObj(name="Object 2")
        obj_list = ConcreteObjList([obj1])
        self.assertEqual(len(obj_list), 1)
        self.assertIn(obj1, obj_list)
        obj_list.insert(0, obj2)
        self.assertEqual(len(obj_list), 2)
        self.assertIn(obj2, obj_list)
        self.assertEqual(obj_list[0].name, "Object 2")
        self.assertEqual(obj_list[1].name, "Object 1") 

    def test_baseobj_list_remove(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'name'}
            name: str = Field(..., description="Name of the object")
        
        class ConcreteObjList(BaseObjList[ConcreteObj]):
            def describe(self) -> str:
                return f"Concrete Object List with {len(self)} items."
        
        obj1 = ConcreteObj(name="Object 1")
        obj2 = ConcreteObj(name="Object 2")
        obj_list = ConcreteObjList([obj1, obj2])
        self.assertEqual(len(obj_list), 2)
        self.assertIn(obj1, obj_list)
        self.assertIn(obj2, obj_list)
        obj_list.remove(obj1)
        self.assertEqual(len(obj_list), 1)
        self.assertNotIn(obj1, obj_list)
        self.assertIn(obj2, obj_list)

    def test_baseobj_list_add_lists(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'name'}
            name: str = Field(..., description="Name of the object")

        class ConcreteObjList(BaseObjList[ConcreteObj]):
            def describe(self) -> str:
                return f"Concrete Object List with {len(self)} items."
        
        obj1 = ConcreteObj(name="Object 1")
        obj2 = ConcreteObj(name="Object 2")
        obj_list1 = ConcreteObjList([obj1])
        obj_list2 = ConcreteObjList([obj2])
        self.assertEqual(len(obj_list1), 1)
        self.assertIn(obj1, obj_list1)
        self.assertEqual(len(obj_list2), 1)
        self.assertIn(obj2, obj_list2)
        obj_list3 = obj_list1 + obj_list2
        self.assertIsInstance(obj_list3, ConcreteObjList)
        self.assertEqual(len(obj_list3), 2)
        self.assertIn(obj1, obj_list3)
        self.assertIn(obj2, obj_list3)

        class WobblyObj(BaseObj):
            _required_fields = {'name'}
            name: str = Field(..., description="Name of the object")

        class WobblyObjList(BaseObjList[WobblyObj]):
            def describe(self) -> str:
                return f"Wobbly Object List with {len(self)} items."
        
        wob1 = WobblyObj(name="Wobbly 1")
        wob2 = WobblyObj(name="Wobbly 2")
        wob_list1 = WobblyObjList([wob1,wob2])
        with self.assertRaises(TypeError):
            # Trying to add a WobblyObjList to a ConcreteObjList should raise TypeError
            hybrid_list = obj_list1 + wob_list1
        with self.assertRaises(TypeError):
            # Trying to add a ConcreteObjList to a WobblyObjList should raise TypeError
            wob_list1 += obj_list1

    def test_baseobj_list_getitem(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'name'}
            name: str = Field(..., description="Name of the object")

        class ConcreteObjList(BaseObjList[ConcreteObj]):
            def describe(self) -> str:
                return f"Concrete Object List with {len(self)} items."
        
        obj1 = ConcreteObj(name="Object 1")
        obj2 = ConcreteObj(name="Object 2")
        obj_list = ConcreteObjList([obj1, obj2])
        self.assertEqual(obj_list[0].name, "Object 1")
        self.assertEqual(obj_list[1].name, "Object 2")
        with self.assertRaises(IndexError):
            _ = obj_list[2]

    def test_baseobj_list_validation_check(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'name'}
            name: str = Field(..., description="Name of the object")

        class ConcreteObjList(BaseObjList[ConcreteObj]):
            def describe(self) -> str:
                return f"Concrete Object List with {len(self)} items."
        
        class WobblyObj(BaseObj):
            _required_fields = {'name'}
            name: str = Field(..., description="Name of the object")

    
        obj1 = ConcreteObj(name="Object 1")
        obj2 = ConcreteObj(name="Object 2")
        wob1 = WobblyObj(name="Wobbly 1")

        obj_list = ConcreteObjList([obj1, obj2])
        self.assertEqual(len(obj_list), 2)
        self.assertIn(obj1, obj_list)
        self.assertIn(obj2, obj_list)
        with self.assertRaises(TypeError):
            obj_list.append(wob1)

    def test_baseobj_list_iterate(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'name'}
            name: str = Field(..., description="Name of the object")

        class ConcreteObjList(BaseObjList[ConcreteObj]):
            def describe(self) -> str:
                return f"Concrete Object List with {len(self)} items."
        
        obj1 = ConcreteObj(name="Object 1")
        obj2 = ConcreteObj(name="Object 2")
        obj_list = ConcreteObjList([obj1, obj2])
        names = [obj.name for obj in obj_list]
        self.assertEqual(names, ["Object 1", "Object 2"])

    def test_baseobj_list_setitem(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'name'}
            name: str = Field(..., description="Name of the object")

        class ConcreteObjList(BaseObjList[ConcreteObj]):
            def describe(self) -> str:
                return f"Concrete Object List with {len(self)} items."

        obj1 = ConcreteObj(name="Object 1")
        obj2 = ConcreteObj(name="Object 2")
        obj_list = ConcreteObjList([obj1, obj2])
        obj3 = ConcreteObj(name="Object 3")
        obj_list[0] = obj3
        self.assertEqual(obj_list[0].name, "Object 3")

    def test_baseobj_list_filter(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'name', 'color', 'flavor', 'shape'}
            name: str = Field(..., description="Name of the object")
            color: str = Field(..., description="Color of the object")
            flavor: str = Field(..., description="Flavor of the object")
            shape: str = Field(..., description="Shape of the object")

        class ConcreteObjList(BaseObjList[ConcreteObj]):
            def describe(self) -> str:
                return f"Concrete Object List with {len(self)} items."

        obj1 = ConcreteObj(name="Object 1", color="Red", flavor="Sweet", shape="Round")
        obj2 = ConcreteObj(name="Object 2", color="Green", flavor="Sweet", shape="Round")
        obj3 = ConcreteObj(name="Object 3", color="Blue", flavor="Bitter", shape="Triangle")
        obj4 = ConcreteObj(name="Object 4", color="Red", flavor="Bitter", shape="Square")
        obj_list = ConcreteObjList([obj1, obj2, obj3, obj4])
        condition = lambda x: x.color == "Red"
        filtered_list = obj_list.filter(condition)
        self.assertIsInstance(filtered_list, ConcreteObjList)
        self.assertEqual(len(filtered_list), 2)
        self.assertEqual(filtered_list[0].name, "Object 1")
        self.assertEqual(filtered_list[1].name, "Object 4")

        condition = lambda x: x.flavor == "Bitter" and x.color == "Blue"
        filtered_list = obj_list.filter(condition)
        self.assertIsInstance(filtered_list, ConcreteObjList)
        self.assertEqual(len(filtered_list), 1)
        self.assertEqual(filtered_list[0].name, "Object 3")
        self.assertEqual(filtered_list[0].flavor, "Bitter")
        self.assertEqual(filtered_list[0].color, "Blue")

        impossible_condition = lambda x: x.flavor == "Bitter" and x.color == "Green"
        filtered_list = obj_list.filter(impossible_condition)
        self.assertIsInstance(filtered_list, ConcreteObjList)
        self.assertEqual(len(filtered_list), 0)

        # class WobblyObj(BaseObj):
        #     _required_fields = {'name', 'rank', 'serialno'}
        #     name: str = Field(..., description="Name of the object")
        #     rank: int = Field(..., description="Rank of the object")
        #     serialno: int = Field(..., description="Serial number of the object")

        # class WobblyObjList(BaseObjList[WobblyObj]):
        #     def describe(self) -> str:
        #         return f"Wobbly Object List with {len(self)} items."
        
        # wob1 = WobblyObj(name="Wobbly", rank=1, serialno=123)
        # wob2 = WobblyObj(name="Wobbly", rank=2, serialno=456)
        # wob3 = WobblyObj(name="XXXXXX", rank=2, serialno=654)
        # wob4 = WobblyObj(name="Wobbly", rank=1, serialno=789)
        # wob_list = WobblyObjList([wob1, wob2, wob3, wob4])
        # filtered_wob_list = wob_list.filter(name="Wobbly", rank=1)
        # self.assertEqual(len(filtered_wob_list), 2)
        # self.assertEqual(filtered_wob_list[0].name, "Wobbly")
        # self.assertEqual(filtered_wob_list[0].rank, 1)
        # self.assertEqual(filtered_wob_list[0].serialno, 123)
        # self.assertEqual(filtered_wob_list[1].name, "Wobbly")
        # self.assertEqual(filtered_wob_list[1].rank, 1)
        # self.assertEqual(filtered_wob_list[1].serialno, 789)
        # another_filtered_wob_list = wob_list.filter(name="Wobbly", rank=2)
        # self.assertEqual(len(another_filtered_wob_list), 1)
        # self.assertEqual(another_filtered_wob_list[0].name, "Wobbly")
        # self.assertEqual(another_filtered_wob_list[0].rank, 2)
        # self.assertEqual(another_filtered_wob_list[0].serialno, 456)
        # an_empty_filtered_wob_list = wob_list.filter(name="Nonexistent", rank=99)
        # self.assertEqual(len(an_empty_filtered_wob_list), 0)

    # def test_baseobj_list_negfilter(self):
    #     class ConcreteObj(BaseObj):
    #         _required_fields = {'name', 'color', 'flavor', 'shape'}
    #         name: str = Field(..., description="Name of the object")
    #         color: str = Field(..., description="Color of the object")
    #         flavor: str = Field(..., description="Flavor of the object")
    #         shape: str = Field(..., description="Shape of the object")

    #     class ConcreteObjList(BaseObjList[ConcreteObj]):
    #         def describe(self) -> str:
    #             return f"Concrete Object List with {len(self)} items."

    #     obj1 = ConcreteObj(name="Object 1", color="Red", flavor="Sweet", shape="Round")
    #     obj2 = ConcreteObj(name="Object 2", color="Green", flavor="Sweet", shape="Round")
    #     obj3 = ConcreteObj(name="Object 3", color="Blue", flavor="Bitter", shape="Triangle")
    #     obj4 = ConcreteObj(name="Object 4", color="Red", flavor="Bitter", shape="Square")
    #     obj_list = ConcreteObjList([obj1, obj2, obj3, obj4])
    #     condition = lambda x: x.color == "Red"
    #     surviving_list = obj_list.negfilter(condition)
    #     self.assertEqual(len(surviving_list), 2)
    #     self.assertIn(obj2, surviving_list)
    #     self.assertIn(obj3, surviving_list)
    #     self.assertNotIn(obj1, surviving_list)
    #     self.assertNotIn(obj4, surviving_list)

    #     condition = lambda x: x.flavor == "Sweet" and x.shape == "Round"
    #     surviving_list = obj_list.negfilter(condition)
    #     self.assertEqual(len(surviving_list), 2)
    #     self.assertIn(obj3, surviving_list)
    #     self.assertNotIn(obj1, surviving_list)
    #     self.assertNotIn(obj2, surviving_list)
    #     self.assertIn(obj4, surviving_list)

    def test_baseobj_list_get(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'name', 'rank'}
            name: str = Field(..., description="Name of the object")
            rank: int = Field(..., description="Rank of the object")

        class ConcreteObjList(BaseObjList[ConcreteObj]):
            def describe(self) -> str:
                return f"Concrete Object List with {len(self)} items."

        obj1 = ConcreteObj(name="Object 1", rank=1)
        obj2 = ConcreteObj(name="Object 2", rank=2)
        obj_list = ConcreteObjList([obj1, obj2])
        self.assertEqual(obj_list.get(lambda x: x.name == "Object 1").name, "Object 1")
        self.assertEqual(obj_list.get(lambda x: x.name == "Object 2").name, "Object 2")
        self.assertIsNone(obj_list.get(lambda x: x.name == "Nonexistent Object"))
        obj3 = ConcreteObj(name="Object 3", rank=2)
        obj_list.append(obj3)
        self.assertEqual(obj_list.get(lambda x: x.name == "Object 3").name, "Object 3")
        self.assertEqual(obj_list.get(lambda x: x.name == "Object 3").rank, 2)
        rank2 = obj_list.get(lambda x: x.rank == 2)
        self.assertEqual(len(rank2), 2)
        self.assertEqual(rank2[0].name, "Object 2")
        self.assertEqual(rank2[1].name, "Object 3")
    
    def test_baseobj_list_iget(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'name', 'rank'}
            name: str = Field(..., description="Name of the object")
            rank: int = Field(..., description="Rank of the object")

        class ConcreteObjList(BaseObjList[ConcreteObj]):
            def describe(self) -> str:
                return f"Concrete Object List with {len(self)} items."

        obj1 = ConcreteObj(name="Object 1", rank=1)
        obj2 = ConcreteObj(name="Object 2", rank=2)
        obj_list = ConcreteObjList([obj1, obj2])
        self.assertEqual(obj_list.iget(lambda x: x.name == "Object 1"), 0)
        self.assertEqual(obj_list.iget(lambda x: x.name == "Object 2"), 1)
        self.assertIsNone(obj_list.iget(lambda x: x.name == "Nonexistent Object"))
        obj3 = ConcreteObj(name="Object 3", rank=2)
        obj_list.append(obj3)
        self.assertEqual(obj_list.iget(lambda x: x.name == "Object 3"), 2)
        rank2 = obj_list.iget(lambda x: x.rank == 2)
        self.assertEqual(len(rank2), 2)
        self.assertEqual(rank2[0], 1)
        self.assertEqual(rank2[1], 2)

    def test_baseobj_list_set(self):
        class PersonObj(BaseObj):
            _required_fields = {'name', 'number'}
            _optional_fields = {'spouse', 'children', 'dict_of_pets'}
            name: str = Field(..., description="Name of the object")
            number: int = Field(..., description="A number")
            spouse: Optional["PersonObj"] = None
            children: Optional[List["PersonObj"]] = Field(default_factory=list, description="Child objects")
            dict_of_pets: Optional[Dict[str, "PersonObj"]] = Field(default_factory=dict, description="Dictionary of pet objects")

        class FamilyObj(BaseObj):
            _required_fields = {'root', 'surname1', 'address', 'number'}
            _optional_fields = {'surname2', 'gamenight'}
            root: PersonObj = Field(..., description="Root person of the family")
            surname1: str = Field(..., description="Primary surname of the family")
            surname2: Optional[str] = Field(None, description="Secondary surname of the family")
            address: str = Field(..., description="Family address")
            number: int = Field(..., description="Family number")
            gamenight: Optional[bool] = Field(False, description="Is there a game night?")

        class FamilyObjList(BaseObjList[FamilyObj]):
            def describe(self) -> str:
                return f"Concrete Object List with {len(self)} items."
            
        families=FamilyObjList()

        root1 = PersonObj(name="Root 1", number=1)
        family1 = FamilyObj(root=root1, surname1="Smith", address="123 Elm St", number=100)
        families.append(family1)

        spouse2= PersonObj(name="Spouse 2", number=3)
        root2= PersonObj(name="Root 2", number=2, spouse=spouse2)
        family2 = FamilyObj(root=root2, surname1="Johnson", address="456 Oak St", number=101)
        families.append(family2)

        child3 = PersonObj(name="Child 3", number=4)
        spouse3 = PersonObj(name="Spouse 3", number=5)
        root3 = PersonObj(name="Root 3", number=6, spouse=spouse3, children=[child3])
        family3 = FamilyObj(root=root3, surname1="Williams", surname2="Jones", address="789 Pine St", number=102, gamenight=True)
        families.append(family3)

        grandkids = [PersonObj(name="GrandKid 1", number=7), PersonObj(name="GrandKid 2", number=8)]
        children4 = [PersonObj(name="Child 4", number=9, children=grandkids), PersonObj(name="Child 5", number=10)]
        spouse4 = PersonObj(name="StepMom 4", number=11)
        root4 = PersonObj(name="Root 4", number=12, spouse=spouse4, children=children4)
        family4 = FamilyObj(root=root4, surname1="Brown", surname2="Green", address="101 Maple St", gamenight=True, number=103)
        families.append(family4)
        self.assertEqual(len(families), 4)

        families.set(number=666,shallow=True)
        for family in families:
            self.assertEqual(family.number, 666)
            self.assertNotEqual(family.root.number, 666)  # Shallow set does not change nested objects
        families.set(number=777, shallow=False)
        for family in families:
            self.assertEqual(family.number, 777)
            self.assertEqual(family.root.number, 777)
            if family.root.spouse:
                self.assertEqual(family.root.spouse.number, 777)
            for child in family.root.children:
                self.assertEqual(child.number, 777)
                if child.children:
                    for grandchild in child.children:
                        self.assertEqual(grandchild.number, 777)

    def test_baseobj_list_prune_exclusions(self):
        exclusions = lambda x: (x.restype in ['X', 'Y', 'Z'] or x.resname in ['POOPYBUTT', 'FARTYPOO'])
        class Residue(BaseObj):
            _required_fields = {'restype','resname','resid'}
            restype: str = Field(..., description="Type of the residue")
            resname: str = Field(..., description="Name of the residue")
            resid: ResID = Field(..., description="Residue ID")
        
        restype_options=['A','B','C','D','E','F','G','H','I','J','X','Y','Z']
        resname_options=['ALA','GLY','SER','THR','LEU','ILE','VAL','PRO','PHE','TRP','POOPYBUTT','FARTYPOO']

        class ResidueList(BaseObjList[Residue]):
            def describe(self) -> str:
                return f"Residue List with {len(self)} items."
                
        R = ResidueList()
        for _ in range(100):
            rtidx=_%len(restype_options)
            rnidx=_%len(resname_options)
            R.append(Residue(resid=ResID(_), restype=restype_options[rtidx], resname=resname_options[rnidx]))

        self.assertEqual(len(R), 100)
        BadResidues = R.prune_exclusions(exclusions)
        # logger.debug(f"BadResidues: {BadResidues}")
        # logger.debug(f"R.filter(Y): {R.filter(lambda x: x.restype == 'Y')}")
        self.assertEqual(len(BadResidues), 34)
        self.assertEqual(len(R), 66)
        for b in BadResidues:
            self.assertTrue((b.restype in ['X', 'Y', 'Z']) or (b.resname in ['POOPYBUTT', 'FARTYPOO']))
        for g in R:
            self.assertFalse((g.restype in ['X', 'Y', 'Z']) or (g.resname in ['POOPYBUTT', 'FARTYPOO']))

        R = ResidueList()
        for _ in range(100):
            rtidx=_%len(restype_options)
            rnidx=_%len(resname_options)
            R.append(Residue(resid=ResID(_), restype=restype_options[rtidx], resname=resname_options[rnidx]))

        self.assertEqual(len(R), 100)
        exclusions = lambda x: (x.restype in ['X', 'Y', 'Z'] and x.resname in ['POOPYBUTT', 'FARTYPOO'])
        BadResidues = R.prune_exclusions(exclusions)
        # logger.debug(f"BadResidues: {BadResidues}")
        # logger.debug(f"R.filter(Y): {R.filter(lambda x: x.restype == 'Y')}")
        self.assertEqual(len(BadResidues), 3)
        self.assertEqual(len(R), 97)
        for b in BadResidues:
            self.assertTrue((b.restype in ['X', 'Y', 'Z']) and (b.resname in ['POOPYBUTT', 'FARTYPOO']))
        for g in R:
            self.assertFalse((g.restype in ['X', 'Y', 'Z']) and (g.resname in ['POOPYBUTT', 'FARTYPOO']))

    def test_baseobj_list_dict_to_condition(self):
        exclusions = {'restype':['X','Y','Z'],'resname':['POOPYBUTT','FARTYPOO']}
        class Residue(BaseObj):
            _required_fields = {'restype','resname'}
            restype: str = Field(..., description="Type of the residue")
            resname: str = Field(..., description="Name of the residue")
        
        restype_options=['A','B','C','D','E','F','G','H','I','J','X','Y','Z']
        resname_options=['ALA','GLY','SER','THR','LEU','ILE','VAL','PRO','PHE','TRP','POOPYBUTT','FARTYPOO']

        class ResidueList(BaseObjList[Residue]):
            def describe(self) -> str:
                return f"Residue List with {len(self)} items."
                
        R = ResidueList()
        for _ in range(100):
            rtidx=_%len(restype_options)
            rnidx=_%len(resname_options)
            R.append(Residue(restype=restype_options[rtidx], resname=resname_options[rnidx]))

        self.assertEqual(len(R), 100)
        condition = R.dict_to_condition(exclusions, conjunction='or')
        BadResidues = R.prune_exclusions(condition)
        self.assertEqual(len(BadResidues), 34)
        self.assertEqual(len(R), 66)
        for b in BadResidues:
            self.assertTrue((b.restype in exclusions['restype']) or (b.resname in exclusions['resname']))
        for g in R:
            self.assertFalse((g.restype in exclusions['restype']) or (g.resname in exclusions['resname']))

    def test_baseobj_list_uniqattrs(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'name','number'}
            name: str = Field(..., description="Name of the object")
            number: int = Field(..., description="A number")

        class ConcreteObjList(BaseObjList[ConcreteObj]):
            def describe(self) -> str:
                return f"Concrete Object List with {len(self)} items."

        obj1 = ConcreteObj(name="Object 1", number=1)
        obj2 = ConcreteObj(name="Object 2", number=2)
        obj3 = ConcreteObj(name="Object 3", number=1)
        obj_list = ConcreteObjList([obj1, obj2, obj3])
        self.assertEqual(len(obj_list), 3)

        unique_attrs = obj_list.uniqattrs(['name'])
        self.assertEqual(len(unique_attrs), 1)
        unique_names=unique_attrs['name']
        self.assertEqual(len(unique_names), 3)
        self.assertIn("Object 1", unique_names)
        self.assertIn("Object 2", unique_names)
        self.assertIn("Object 3", unique_names)

        unique_numbers = obj_list.uniqattrs(['number'])
        self.assertEqual(len(unique_numbers), 1)
        unique_numbers=unique_numbers['number']
        self.assertIn(1, unique_numbers)
        self.assertIn(2, unique_numbers)

    def test_baseobj_list_binnify(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'name','number'}
            name: str = Field(..., description="Name of the object")
            number: int = Field(..., description="A number")

        class ConcreteObjList(BaseObjList[ConcreteObj]):
            def describe(self) -> str:
                return f"Concrete Object List with {len(self)} items."

        L=ConcreteObjList()
        for _ in range(100):
            L.append(ConcreteObj(name=f"H{_%2}", number=_%3))

        bin_attrs = L.binnify(['number'])
        self.assertIn('number=0', bin_attrs)
        self.assertEqual(len(bin_attrs['number=0']), 34)  # 0 appears 34 times
        self.assertIn('number=1', bin_attrs)
        self.assertEqual(len(bin_attrs['number=1']), 33)  # 1 appears 33 times
        self.assertIn('number=2', bin_attrs)
        self.assertEqual(len(bin_attrs['number=2']), 33)  # 2 appears 33 times
        self.assertNotIn('number=3', bin_attrs)  # 3 is not a valid bin for number % 3

        more_bin_attrs=L.binnify(['name'])
        self.assertIn('name=H0', more_bin_attrs)
        self.assertEqual(len(more_bin_attrs['name=H0']), 50)  # H0 appears 50 times
        self.assertIn('name=H1', more_bin_attrs)
        self.assertEqual(len(more_bin_attrs['name=H1']), 50)  # H1 appears 50 times
        self.assertNotIn('name=H2', more_bin_attrs)

        yet_more_bin_attrs=L.binnify()
        self.assertIn('name=H0-number=0', yet_more_bin_attrs)
        self.assertEqual(len(yet_more_bin_attrs['name=H0-number=0']), 17)
        self.assertIn('name=H0-number=1', yet_more_bin_attrs)
        self.assertEqual(len(yet_more_bin_attrs['name=H0-number=1']), 16)
        self.assertIn('name=H0-number=2', yet_more_bin_attrs)
        self.assertEqual(len(yet_more_bin_attrs['name=H0-number=2']), 17)
        self.assertNotIn('name=H0-number=3', yet_more_bin_attrs)  # 3 is not a valid bin for number % 3
        self.assertIn('name=H1-number=0', yet_more_bin_attrs)
        self.assertEqual(len(yet_more_bin_attrs['name=H1-number=0']), 17)
        self.assertIn('name=H1-number=1', yet_more_bin_attrs)
        self.assertEqual(len(yet_more_bin_attrs['name=H1-number=1']), 17)
        self.assertIn('name=H1-number=2', yet_more_bin_attrs)
        self.assertEqual(len(yet_more_bin_attrs['name=H1-number=2']), 16)
        self.assertNotIn('name=H1-number=3', yet_more_bin_attrs)  # 3 is not a valid bin for number % 3
    
    def test_baseobj_list_puniq(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'name', 'number'}
            name: str = Field(..., description="Name of the object")
            number: int = Field(..., description="A number")

        class ConcreteObjList(BaseObjList[ConcreteObj]):
            def describe(self) -> str:
                return f"Concrete Object List with {len(self)} items."

        L=ConcreteObjList()
        for _ in range(100):
            L.append(ConcreteObj(name=f"H{_}", number=_))

        self.assertTrue(L.puniq(['number']))
        self.assertTrue(L.puniq(['name']))
        self.assertTrue(L.puniq())

        AL=ConcreteObjList()
        for _ in range(100):
            AL.append(ConcreteObj(name=f"H{_%10}", number=_))

        self.assertTrue(AL.puniq(['number']))
        self.assertFalse(AL.puniq(['name']))
        self.assertTrue(AL.puniq())

        BL=ConcreteObjList()
        for _ in range(100):
            BL.append(ConcreteObj(name=f"H{_%5}", number=_%5))

        self.assertFalse(BL.puniq(['number']))
        self.assertFalse(BL.puniq(['name']))
        self.assertFalse(BL.puniq())

    def test_baseobj_list_with_self_fields(self):
        
        class PersonObj(BaseObj):
            _required_fields = {'name'}
            _optional_fields= {'children'}
            model_config = ConfigDict(arbitrary_types_allowed=True)
            name: str = Field(..., description="Name of the object")
            children: Annotated[PersonList, PersonList.validate] = Field(default_factory=lambda: PersonList())

        class PersonList(BaseObjList[PersonObj]):
            def describe(self) -> str:
                return f"Self Referencing Object List with {len(self)} items."
            
        PersonObj.model_rebuild()

        Person1 = PersonObj(name="Alice")
        Person2 = PersonObj(name="Bob")
        Person3 = PersonObj(name="Charlie", children=PersonList([Person1, Person2]))
        self.assertEqual(Person3.name, "Charlie")
        self.assertEqual(len(Person3.children), 2)
        self.assertEqual(Person3.children[0].name, "Alice")
        self.assertEqual(Person3.children[1].name, "Bob")

    def test_baseobj_list_remove_duplicates(self):
        class ConcreteObj(BaseObj):
            _required_fields = {'name', 'number'}
            name: str = Field(..., description="Name of the object")
            number: int = Field(..., description="A number")

        class ConcreteObjList(BaseObjList[ConcreteObj]):
            def describe(self) -> str:
                return f"Concrete Object List with {len(self)} items."

        obj1 = ConcreteObj(name="Object 1", number=1)
        obj2 = ConcreteObj(name="Object 2", number=2)
        obj3 = ConcreteObj(name="Object 1", number=1)
        obj_list = ConcreteObjList([obj1, obj2, obj3])
        self.assertEqual(len(obj_list), 3)
        obj_list.remove_duplicates()
        self.assertEqual(len(obj_list), 2)

    def test_assign_obj_to_attr(self):
        class ConcreteObj(BaseObj):
            model_config = ConfigDict(arbitrary_types_allowed=True)
            _required_fields = {'name', 'number'}
            _optional_fields = {'child'}
            name: str = Field(..., description="Name of the object")
            number: int = Field(..., description="A number")
            child: SecondaryObj | None = Field(None, description="Child object")

        class SecondaryObj(BaseObj):
            _required_fields = {'number_a','number_b','color'}
            number_a: int = Field(..., description="A number")
            number_b: int = Field(..., description="Another number")
            color: str = Field(..., description="Color of the object")
        
        class SecondaryObjList(BaseObjList[SecondaryObj]):
            def describe(self) -> str:
                return f"Secondary Object List with {len(self)} items."

        ConcreteObj.model_rebuild()

        obj1 = ConcreteObj(name="Object 1", number=1)

        secondary_obj1 = SecondaryObj(number_a=1, number_b=1, color="red")
        secondary_obj2 = SecondaryObj(number_a=1, number_b=4, color="blue")

        secondary_obj_list = SecondaryObjList([secondary_obj1, secondary_obj2])

        # assign the first secondary object to the child attribute of both concrete objects
        obj1.assign_obj_to_attr('child', secondary_obj_list, number_b='number')

        self.assertEqual(obj1.child, secondary_obj1)

    def test_assign_objs_to_attr(self):
        class ConcreteObj(BaseObj):
            model_config = ConfigDict(arbitrary_types_allowed=True)
            _required_fields = {'name', 'number'}
            _optional_fields = {'child'}
            name: str = Field(..., description="Name of the object")
            number: int = Field(..., description="A number")
            child: AnotherConObj | None = Field(None, description="Child object")

        class ConcreteObjList(BaseObjList[ConcreteObj]):
            def describe(self) -> str:
                return f"Concrete Object List with {len(self)} items."

        class AnotherConObj(BaseObj):
            _required_fields = {'number_a','number_b','color'}
            number_a: int = Field(..., description="A number")
            number_b: int = Field(..., description="Another number")
            color: str = Field(..., description="Color of the object")

        class SecondaryObjList(BaseObjList[AnotherConObj]):
            def describe(self) -> str:
                return f"Secondary Object List with {len(self)} items."

        ConcreteObj.model_rebuild()

        obj1 = ConcreteObj(name="Object 1", number=1) # this one will be removed
        obj2 = ConcreteObj(name="Object 2", number=2)
        obj3 = ConcreteObj(name="Object 3", number=3)
        obj_list = ConcreteObjList([obj1, obj2, obj3])

        secondary_obj1 = AnotherConObj(number_a=11, number_b=3, color="red")
        secondary_obj2 = AnotherConObj(number_a=21, number_b=2, color="blue")
        secondary_obj3 = AnotherConObj(number_a=31, number_b=4, color="green")

        secondary_obj_list = SecondaryObjList([secondary_obj1, secondary_obj2, secondary_obj3])

        deleted_objs = obj_list.assign_objs_to_attr('child', secondary_obj_list, number_b='number')

        self.assertEqual(obj1.child, None)
        self.assertEqual(obj2.child, secondary_obj2)
        self.assertEqual(obj3.child, secondary_obj1)

        self.assertEqual(len(deleted_objs), 1)
        self.assertIn(obj1, deleted_objs)

    def test_update_attr_from_obj_attr(self):
        class ConcreteObj(BaseObj):
            model_config = ConfigDict(arbitrary_types_allowed=True)
            _required_fields = {'name', 'number'}
            _optional_fields = {'child'}
            name: str = Field(..., description="Name of the object")
            number: int = Field(..., description="A number")
            child: SecondaryObj | None = Field(None, description="Child object")

        class ConcreteObjList(BaseObjList[ConcreteObj]):
            def describe(self) -> str:
                return f"Concrete Object List with {len(self)} items."

        class SecondaryObj(BaseObj):
            _required_fields = {'number_a','number_b','color'}
            number_a: int = Field(..., description="A number")
            number_b: int = Field(..., description="Another number")
            color: str = Field(..., description="Color of the object")

        ConcreteObj.model_rebuild()

        secondary_obj1 = SecondaryObj(number_a=11, number_b=3, color="red")
        secondary_obj2 = SecondaryObj(number_a=21, number_b=2, color="blue")
        secondary_obj3 = SecondaryObj(number_a=31, number_b=4, color="green")
        obj1 = ConcreteObj(name="Object 1", number=1, child=secondary_obj1) # this one will be removed
        obj2 = ConcreteObj(name="Object 2", number=2, child=secondary_obj2)
        obj3 = ConcreteObj(name="Object 3", number=3, child=secondary_obj3)
        obj_list = ConcreteObjList([obj1, obj2, obj3])
        # Update the 'number' attribute of each ConcreteObj based on the 'number_b' attribute of its child SecondaryObj
        obj_list.update_attr_from_obj_attr('number', 'child', 'number_b')
        self.assertEqual(obj1.number, 3)  # obj1's child had number_b=3
        self.assertEqual(obj2.number, 2)  # obj2's child had number_b=2
        self.assertEqual(obj3.number, 4)  # obj3's

    def test_update_attr_from_objlist_elem_attr(self):
        class PersonObj(BaseObj):
            _required_fields = {'name', 'favorite_food', 'most_hated_food'}
            _optional_fields = {'children'}
            model_config = ConfigDict(arbitrary_types_allowed=True)
            name: str = Field(..., description="Name of the object")
            favorite_food: str = Field(..., description="Favorite food of the object")
            most_hated_food: str = Field(..., description="Most hated food of the object")
            children: Annotated[PersonList, PersonList.validate] = Field(default_factory=lambda: PersonList())

        class PersonList(BaseObjList[PersonObj]):
            def describe(self) -> str:
                return f"Self Referencing Object List with {len(self)} items."

        PersonObj.model_rebuild()

        Person1_child1 = PersonObj(name="Bob", favorite_food="Sushi", most_hated_food="Spinach")
        Person1_child2 = PersonObj(name="Charlie", favorite_food="Pasta", most_hated_food="Zucchini")
        Person1 = PersonObj(name="Alice", favorite_food="Pizza", most_hated_food="Broccoli", 
                            children=PersonList([Person1_child1, Person1_child2]))
        
        Person2_child1 = PersonObj(name="Dave", favorite_food="Tacos", most_hated_food="Eggplant")
        Person2_child2 = PersonObj(name="Eve", favorite_food="Burgers", most_hated_food="Liver")
        Person2 = PersonObj(name="Frank", favorite_food="Salad", most_hated_food="Fish", 
                            children=PersonList([Person2_child1, Person2_child2]))
    
        persons = PersonList([Person1, Person2])

        persons.update_attr_from_objlist_elem_attr('favorite_food', 'children', 0, 'most_hated_food')

        self.assertEqual(Person1.favorite_food, "Spinach")
        self.assertEqual(Person2.favorite_food, "Eggplant")

    def test_baseobj_list_remove_instance(self):
        # define a concreteobj class
        class ConcreteObj(BaseObj):
            _required_fields = {'name', 'number'}
            name: str = Field(..., description="Name of the object")
            number: int = Field(..., description="Number of the object")
        # and the corresponding ConcreteObjList class
        class ConcreteObjList(BaseObjList[ConcreteObj]):
            def describe(self) -> str:
                return f"Concrete Object List with {len(self)} items."

        # make a list of concrete objects
        obj1 = ConcreteObj(name="Object 1", number=1)
        obj2 = ConcreteObj(name="Object 2", number=2)
        obj3 = ConcreteObj(name="Object 3", number=3)
        obj_list = ConcreteObjList([obj1, obj2, obj3])

        # Remove obj2 and check the list contents
        obj_list.remove_instance(obj2)
        self.assertNotIn(obj2, obj_list)
        self.assertIn(obj1, obj_list)
        self.assertIn(obj3, obj_list)