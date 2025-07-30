"""
Object subpackage for pestifer.

Each module in this subpackage defines a class that represents a specific type of object used in the Pestifer runtime.

Each class is a subclass of the abstract :class:`pestifer.core.baseobj_new.BaseObj`.  :class:`pestifer.core.baseobj_new.BaseObj` is a subclass of :class:`pydantic.BaseModel` and provides additional functionality for Pestifer objects, such as YAML serialization and deserialization, and a unique object ID.

Every class must implement the `_adapt` static method to become concrete.  The job of _adapt is to allow a subclass to convert a raw input (e.g., a string or dictionary) into a dictionary of keyword arguments that can be passed to the class constructor.  This allows for flexible initialization of objects from various input formats.

"""