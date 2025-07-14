# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
Defines the :class:`Example` class for managing examples in the documentation.
"""
from collections import UserList
import yaml
import logging
import re
import os
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

    def __init__(self, name: str, pdbID: str = '', description: str = '', index: int = 0, author_name: str = '', author_email: str = '', companion_files=[]):
        """
        Initialize an Example instance.

        Parameters
        ----------
        name : str
            The name of the example.
        pdbID : str, optional
            The PDB ID associated with the example.
        description : str, optional
            A description of the example.
        index : int, optional
            The index of the example in the list of examples.
        author_name : str, optional
            The name of the author of the example.
        author_email : str, optional
            The email of the author of the example.
        companion_files : list, optional
            A list of companion files associated with the example.
        """
        self.name = name
        self.pdbID = pdbID
        self.description = description
        self.index = index
        self.author_name = author_name
        self.author_email = author_email
        self.companion_files = companion_files

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
            'index': self.index,
            'author_name': self.author_name,
            'author_email': self.author_email,
            'companion_files': self.companion_files
        }

    @classmethod
    def from_yaml(cls,yaml_file_name:str,description:str='',pdbID:str='',author_name:str='',author_email:str='',companion_files: list = []):
        """
        Create a new Example object based on the provided YAML file and optional parameters.  This provides the convenience of setting the necessary attributes of the Example object based on the contents of the YAML file, while also allowing for customization of the description, pdbID, author name, author email, and companion files.

        Parameters
        ----------
        yaml_file_name : str
            The name of the example file to grab information from. This should be a valid file path.
        description : str, optional
            A description of the example; if not provided, extracts the ``title`` field from the YAML file.
        pdbID : str, optional
            The PDB ID associated with the example; if not provided, extracts the ``id`` field from the ``psfgen`` task of the ``tasks`` list in the YAML file.
        author_name : str, optional
            The name of the author of the example; if not provided, defaults to an empty string.
        author_email : str, optional
            The email of the author of the example; if not provided, defaults to an empty string.
        companion_files : list, optional
            A list of companion files associated with the example; defaults to an empty list.

        Returns
        -------
        Example
            An Example object containing the name, description, pdbID, and companion_files of the example.

        Raises
        ------
        ValueError
            If the name of the YAML file is not provided. 
        FileNotFoundError
            If the example file does not exist at the specified path.
        """
        if not yaml_file_name:
            raise ValueError('Name must be provided')
        with open(yaml_file_name, 'r') as f:
            try:
                new_config=yaml.safe_load(f)
            except yaml.YAMLError as e:
                raise ValueError(f'Invalid YAML file {yaml_file_name}: {e}')
        if not description:
            if 'title' in new_config:
                description = new_config['title']
            else:
                description = yaml_file_name.replace('.yaml', '').replace('_', ' ').title()
        if not pdbID:
            if 'id' in new_config.get('tasks', [{}])[0].get('psfgen', {}).get('source', {}):
                pdbID = new_config['tasks'][0]['psfgen']['source']['id']
            elif 'alphafold' in new_config.get('tasks', [{}])[0].get('psfgen', {}).get('source', {}):
                pdbID = new_config['tasks'][0]['psfgen']['source']['alphafold']
            else:
                raise ValueError(f'PDB ID must be provided or found in the YAML file {yaml_file_name}')
        if not author_name and not author_email:
            with open(yaml_file_name, 'r') as f:
                lines= f.readlines()
            for line in lines:
                """ Recognized format: list of strings that start with `# Author:` and contain the author's name and email.
                The format is flexible, but the following rules apply:
                    - the line must start with `# Author:`
                    - a single string contains the '@' symbol, which is used to identify the email address
                    - the name and email must be separated by a comma or whitespace
                    - the name can contain spaces, internal commas, and periods
                    - a terminal comma in the name can be ignored
                    - any nonalphanumeric characters in email string (e.g., the '<' and '>' in '<email-address>') can be deleted from values
                """
                if line.startswith('# Author:'):
                    author_info = line.replace('# Author:', '').strip()
                    tokens = author_info.split()
                    idx_of_email_tokens = [i for i, token in enumerate(tokens) if '@' in token]
                    if not idx_of_email_tokens:
                        raise ValueError(f'No email found in author line: {line}')
                    if len(idx_of_email_tokens) > 1:
                        raise ValueError(f'Multiple emails found in author line: {line}')
                    email_idx = idx_of_email_tokens[0]
                    raw_email = tokens.pop(email_idx)
                    # remove any non-alphanumeric characters from the beginning and end of the email
                    author_email = re.sub(r'^\W+|\W+$', '', raw_email)
                    author_name = ' '.join(tokens).strip(',;')
                    break
        
        input_dict= {
            # the name attributed is assigned the basename of the YAML file without the extension
            'name': os.path.splitext(os.path.basename(yaml_file_name))[0],
            'description': description,
            'pdbID': pdbID,
            'author_name': author_name,
            'author_email': author_email,
            'companion_files': companion_files or []
        }
        return cls(**input_dict)

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

    
    def update_in_place(self, name: str = None, pdbID: str = None, description: str = None, author_name: str = None, author_email: str = None, companion_files: list = None):
        """
        Update the attributes of the example in place.  Return the Example instance itself.
        
        Parameters
        ----------
        name : str, optional
            The new name of the example.
        pdbID : str, optional
            The new PDB ID associated with the example.
        description : str, optional
            The new description of the example.
        author_name : str, optional
            The new name of the author of the example.
        author_email : str, optional
            The new email of the author of the example.
        companion_files : list, optional
            A new list of companion files associated with the example; note, these are appended if they are unique.

        Returns
        -------
        Example
            The updated Example instance.
        """
        if name is not None:
            self.name = name
        if pdbID is not None:
            self.pdbID = pdbID
        if description is not None:
            self.description = description
        if author_name is not None:
            self.author_name = author_name
        if author_email is not None:
            self.author_email = author_email
        if companion_files is not None:
            for f in companion_files:
                if f not in self.companion_files:
                    self.companion_files.append(f)
        return self

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
    
    def get_example_by_name(self, name: str) -> Example:
        """
        Get an example by its name.
        
        Parameters
        ----------
        name : str
            The name of the example to retrieve.
        
        Returns
        -------
        Example
            The Example instance with the specified name or None if not found.
        
        """
        for example in self.data:
            if example.name == name:
                return example
        return None