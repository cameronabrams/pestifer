"""
Object subpackage for pestifer.

Each module in this subpackage defines a class that represents a specific type of object used in the Pestifer runtime.

Each class is a subclass of the abstract :class:`pestifer.core.baseobj.BaseObj`.  :class:`pestifer.core.baseobj.BaseObj` is a subclass of :class:`pydantic.BaseModel` and provides additional functionality for Pestifer objects, such as YAML serialization and deserialization, and a unique object ID.

Subclasses can define their own required and optional fields, as well as overriding the _adapt class method to handle initialization from a specific input data format, such as a shortcode.

"""