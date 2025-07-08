# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
Defines the :class:`Example` class for managing examples in the documentation.
"""
from collections import UserList
import yaml
import logging
logger = logging.getLogger(__name__)

class Example:
    """
    Represents an example in the documentation.

    Attributes
    ----------
    name : str
        The name of the example.  Used to resolve the associate YAML file in the examples directory (typically `pestifer/resources/examples`). Also used to resolve the entry in the toctree in the examples RST file (typically `docs/source/examples.rst`), and the name of the specific RST file for the example (typically `docs/source/examples/{name}.rst`).
    pdbID : str
        The PDB ID associated with the example.
    description : str
        A description of the example.  Used as the title of the example in its RST file.
    index : int
        The index of the example in the list of examples.
    """

    def __init__(self, name: str, pdbID: str = '', description: str = '', index: int = 0):
        self.name = name
        self.pdbID = pdbID
        self.description = description
        self.index = index

    def to_dict(self) -> dict:
        """
        Convert the example to a dictionary representation.
        
        Returns
        -------
        dict
            A dictionary containing the example's attributes.
        """
        return {
            'name': self.name,
            'pdbID': self.pdbID,
            'description': self.description,
            'index': self.index
        }

    def report_line(self,formatter=r'{:>7s}    {:>4s}  {:<30s}    {}') -> str:
        """
        Generate a formatted line for reporting the example.
        
        Parameters
        ----------
        formatter : str, optional
            A format string for the report. Default is a string that formats the index, name, pdbID, and description of the example.
        
        Returns
        -------
        str
            A formatted string representing the example.
        """
        return formatter.format(str(self.index), self.pdbID, self.name, self.description)

    def to_yaml(self) -> str:
        """
        Convert the example to a YAML string.
        
        Returns
        -------
        str
            A YAML representation of the example.
        """
        return yaml.dump(self.to_dict(), default_flow_style=False)

    def new_rst(self) -> str:
        """
        Generate a new reStructuredText file for the example.
        
        Returns
        -------
        str
            The contents of the new reStructuredText file.
        """
        title= f"Example {self.index}: {self.description}"
        yaml_file_name=self.name if self.name.endswith('.yaml') else f"{self.name}.yaml"
        if self.pdbID.startswith('P') and len(self.pdbID) > 4:
            # this is likely an alphafold ID
            url= f"https://alphafold.com/api/prediction/{self.pdbID}"
        else:
            url= f"https://www.rcsb.org/structure/{self.pdbID}"
        return f".. _example {self.name}:\n\n{title}\n" \
               f"{'-' * len(title)}\n\n" \
               f"`PDB ID {self.pdbID} <{url}>`_ is...\n\n" \
               f"This example demonstrates that ...\n\n" \
               f".. literalinclude:: ../../../pestifer/resources/examples/{yaml_file_name}\n" \
               f"    :language: yaml\n\n"

class ExampleList(UserList):
    """
    Represents a list of examples.

    Inherits from `UserList` to provide list-like behavior for managing examples.
    """

    def __init__(self, examples=None):
        super().__init__(examples or [])
        for i, example in enumerate(self.data):
            if not isinstance(example, Example):
                raise TypeError(f'Expected Example instance, got {type(example)}')
            if not hasattr(example, 'index'):
                example.index = i + 1  # Set index if not already set

    @classmethod
    def from_list_of_dicts(cls, examples_list):
        """
        Create an ExampleList from a list of dictionaries.
        
        Parameters
        ----------
        examples_list : list of dict
            A list of dictionaries, each representing an example.
        
        Returns
        -------
        ExampleList
            An instance of ExampleList containing the examples.
        """
        examples = [Example(**example_data) for example_data in examples_list]
        # set index for each example
        for i, example in enumerate(examples):
            example.index = i + 1
        return cls(examples)

    def append(self, example: Example):
        """
        Append an Example instance to the list.
        
        Parameters
        ----------
        example : Example
            The Example instance to append.
        """
        if not isinstance(example, Example):
            raise TypeError(f'Expected Example instance, got {type(example)}')
        super().append(example)
        example.index = len(self.data)

    def insert(self, index: int, example: Example):
        """
        Insert an Example instance at a specific index.
        
        Parameters
        ----------
        index : int
            The 0-based index at which to insert the example.
        example : Example
            The Example instance to insert.
        """
        if not isinstance(example, Example):
            raise TypeError(f'Expected Example instance, got {type(example)}')
        super().insert(index, example)
        # reindex the examples after insertion
        for i in range(index + 1, len(self.data)):
            self.data[i].index = i + 1
    
    def remove(self, example: Example):
        """
        Remove an Example instance from the list.
        
        Parameters
        ----------
        example : Example
            The Example instance to remove.
        """
        if not isinstance(example, Example):
            raise TypeError(f'Expected Example instance, got {type(example)}')
        super().remove(example)
        # reindex all just to be safe
        for i in range(len(self.data)):
            self.data[i].index = i + 1

    def to_yaml(self) -> str:
        """
        Convert the list of examples to a YAML string.
        
        Returns
        -------
        str
            A YAML representation of the examples.
        """
        return yaml.dump([example.to_dict() for example in self.data], default_flow_style=False)

    def to_list_of_dicts(self) -> list:
        """
        Convert the list of examples to a list of dictionaries.
        
        Returns
        -------
        list
            A list of dictionaries, each representing an example.
        """
        return [example.to_dict() for example in self.data]

    def from_yaml(self, yaml_str: str, overwrite: bool = False):
        """
        Load examples from a YAML string.
        
        Parameters
        ----------
        yaml_str : str
            A YAML string containing the examples.
        overwrite : bool
            Whether to overwrite existing examples.
        """
        examples = yaml.safe_load(yaml_str)
        if not overwrite:
            self.clear()
        for example_data in examples:
            example = Example(**example_data)
            self.append(example)
        for i, example in enumerate(self.data):
            example.index = i + 1

    @classmethod
    def read_yaml(cls, yaml_str: str):
        """
        Load examples from a YAML string.
        
        Parameters
        ----------
        yaml_str : str
            A YAML string containing the examples.
        """
        inst=cls([])
        examples = yaml.safe_load(yaml_str)
        for example_data in examples:
            example = Example(**example_data)
            inst.append(example)
        for i, example in enumerate(inst.data):
            example.index = i + 1
        return inst