# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
Defines the :class:`Example` class for managing examples.
"""
import yaml
import logging
from pathlib import Path
from pydantic import Field
from typing import ClassVar
from .baseobj import BaseObj, BaseObjList

logger = logging.getLogger(__name__)

class Example(BaseObj):
    """
    Represents an example.
    """

    _required_fields = {'example_id', 'title', 'shortname'}
    _optional_fields = {'db_id', 'author_name', 'author_email', 'input', 'auxiliary_inputs', 'outputs'}

    example_id: int = Field(...,title='The unique identifier for the example.')
    title: str = Field(..., title='A single-line descriptive title for the example.')
    shortname: str = Field(..., title='A short name for the example; <shortname>.yaml is assumed to be the name of the main input script.')

    db_id: str | None = Field(None, title='The database ID associated with the example, if applicable.')
    author_name: str | None = Field(None, title='The name of the author of the example.')
    author_email: str | None = Field(None, title='The email of the author of the example.')
    input: str | Path | None = Field(None, title='The input script name for the example, if different from <shortname>.yaml.')
    auxiliary_inputs: list[str | Path] | None = Field(None, title='A list of auxiliary inputs for the example.')
    outputs: list[str | Path] | None = Field(None, title='A list of outputs for the example.')

    inputs_subdir: ClassVar[str] = "inputs"
    outputs_subdir: ClassVar[str] = "outputs"

    folder_name_format: ClassVar[str] = 'ex{example_id:02d}'
    """ format for the name of the root folder of each example """

    def to_dict(self, ignore_none: bool = False) -> dict:
        """
        Convert the example to a dictionary representation.
        
        Returns
        -------
        dict
            A dictionary containing the example's attributes.
        """
        res_dict = {
            'example_id': self.example_id,
            'title': self.title,
            'shortname': self.shortname
            }
        for of in self._optional_fields:
            if not ignore_none:
                res_dict[of] = getattr(self, of)
            elif getattr(self, of) is not None:
                res_dict[of] = getattr(self, of)
        return res_dict

    @property
    def scriptname(self) -> str:
        return f'{self.shortname}.yaml'

    @property
    def rootfolderpath(self) -> Path:
        return Path(self.folder_name_format.format(example_id=self.example_id))
    
    @property
    def inputspath(self) -> Path:
        return Path(f'{self.rootfolderpath}/{self.inputs_subdir}')
    
    @property
    def outputspath(self) -> Path:
        return Path(f'{self.rootfolderpath}/{self.outputs_subdir}')

    @property
    def scriptpath(self) -> Path:
        return Path(f'{self.inputspath}/{self.scriptname}')

    @staticmethod
    def _get_implied_metadata(script: str) -> tuple[str, str | None]:
        """
        Get the implied title and database ID from Pestifer input script.

        Parameters
        ----------
        script : str
            The path to the pestifer script file.

        Returns
        -------
        tuple[str, str | None]
            A tuple containing the implied title and database ID.
        """
        with open(script, 'r') as f:
            main_input: dict = yaml.safe_load(f)
            title = main_input.get('title', 'No title provided')
            tasks = main_input.get('tasks', [])
            if len(tasks) > 0 and 'fetch' in tasks[0]:
                db_id = tasks[0]['fetch']['sourceID']
            else:
                db_id = None
        return title, db_id
    
    def report_line(self, formatter=r'{:>7s}    {:>4s}  {:<30s}    {}') -> str:
        """
        Generate a formatted line for reporting the example.
        
        Parameters
        ----------
        formatter : str, optional
            A format string for the report. Default is a string that formats the index, name, pdbID, and title of the example.
        
        Returns
        -------
        str
            A formatted string representing the example.
        """
        dbidstr = self.db_id if self.db_id else 'N/A'
        return formatter.format(str(self.example_id), dbidstr, self.shortname, self.title)

    def to_yaml(self) -> str:
        """
        Convert the example to a YAML string.
        
        Returns
        -------
        str
            A YAML representation of the example.
        """
        return yaml.dump(self.to_dict(ignore_none=True), default_flow_style=False)

    def update_options(self, author_name: str = None, author_email: str = None, auxiliary_inputs: list[str] = None, outputs: list[str] = None) -> 'Example':
        """
        Update the optional attributes of the example in place.  Return the Example instance itself.
        
        Parameters
        ----------
        author_name : str, optional
            The new name of the author of the example.
        author_email : str, optional
            The new email of the author of the example.
        auxiliary_inputs : list[str], optional
            A new list of auxiliary inputs for the example.
        outputs : list[str], optional
            A new list of outputs for the example.

        Returns
        -------
        Example
            The updated Example instance.
        """
        if author_name is not None:
            self.author_name = author_name
        if author_email is not None:
            self.author_email = author_email
        if auxiliary_inputs is not None:
            self.auxiliary_inputs = auxiliary_inputs
        if outputs is not None:
            self.outputs = outputs
        return self

class ExampleList(BaseObjList[Example]):
    """
    Represents a list of examples.

    Inherits from `UserList` to provide list-like behavior for managing examples.
    """

    def describe(self) -> str:  
        return f'ExampleList with {len(self.data)} examples'

    @classmethod
    def from_list_of_dicts(cls, examples_list) -> 'ExampleList':
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
        return cls(examples)

    def append(self, example: Example):
        """
        Append an Example instance to the list.
        
        Parameters
        ----------
        example : Example
            The Example instance to append.
        """
        assert example.example_id not in [e.example_id for e in self.data], f'Example with ID {example.example_id} already exists'
        super().append(example)
    
    def to_yaml(self) -> str:
        """
        Convert the list of examples to a YAML string.
        
        Returns
        -------
        str
            A YAML representation of the examples.
        """
        return yaml.dump([example.to_dict(ignore_none=True) for example in self.data], default_flow_style=False)

    def to_list_of_dicts(self) -> list:
        """
        Convert the list of examples to a list of dictionaries.
        
        Returns
        -------
        list
            A list of dictionaries, each representing an example.
        """
        return [example.to_dict(ignore_none=True) for example in self.data]

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

    @classmethod
    def read_yaml(cls, yaml_str: str):
        """
        Load examples from a YAML string.
        
        Parameters
        ----------
        yaml_str : str
            A YAML string containing the examples.
        """
        inst = cls([])
        examples = yaml.safe_load(yaml_str)
        for example_data in examples:
            example = Example(**example_data)
            inst.append(example)
        return inst

    def get_example_by_shortname(self, shortname: str) -> Example:
        """
        Get an example by its shortname.

        Parameters
        ----------
        name : str
            The name of the example to retrieve.
        
        Returns
        -------
        Example
            The Example instance with the specified name or None if not found.
        
        """
        query_result = self.get(lambda x: x.shortname == shortname)
        if query_result:
            if isinstance(query_result, list) and len(query_result) > 1:
                logger.warning(f'Multiple examples found with shortname {shortname}. Returning the first one.')
                return query_result[0]
            return query_result
        return None
    
    def get_example_by_example_id(self, example_id: int) -> Example:
        query_result = self.get(lambda x: x.example_id == example_id)
        if query_result:
            if isinstance(query_result, list) and len(query_result) > 1:
                logger.warning(f'Multiple examples found with example_id {example_id}. Returning the first one.')
                return query_result[0]
            return query_result
        return None