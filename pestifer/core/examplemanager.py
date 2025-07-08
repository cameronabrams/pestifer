# Author: Cameron F. Abrams <cfa22@drexel.edu>

""" 
A module for managing example input files in pestifer.
This module provides the :class:`ExampleManager` class, which allows users to check out example YAML files,
report the list of examples, and manage example resources.
"""

import logging
import os
import shutil
import yaml
logger=logging.getLogger(__name__)

from .example import Example, ExampleList
from ..sphinxext.sphinx_examplemanager import SphinxExampleManager

class ExampleManager:
    """
    A class for managing example input files in pestifer.
    This class provides methods to check out example YAML files, report the list of examples, and manage example resources.
    
    Parameters
    ----------
    example_path : str
        The path to the directory containing example input files. This directory should contain an ``info.yaml`` file
        that describes the examples available in that directory.
    docs_path : str
        The path to the directory containing documentation files. This is used to link examples to their documentation
        and is not strictly required for the functionality of the class, but it is useful for generating documentation
        or reports that include examples.
    """

    def __init__(self,example_path,docs_source_path=None):
        if not os.path.isdir(example_path):
            raise FileNotFoundError(f'Directory "{example_path}" containing example YAML files does not exist or is not a directory')
        if docs_source_path is not None and not os.path.isdir(docs_source_path):
            raise FileNotFoundError(f'Docs source path "{docs_source_path}" does not exist or is not a directory')
        self.path=os.path.abspath(example_path)
        self.docs_source_path=os.path.abspath(docs_source_path) if docs_source_path else None
        self._read_info()  # read the info.yaml file to populate the examples list
        if self.docs_source_path:
            self.sphinx_example_manager=SphinxExampleManager(docs_source_path=self.docs_source_path)
        else:
            self.sphinx_example_manager=None

    def _read_info(self):
        """
        Read the info.yaml file and update the examples list.
        
        This method is called internally to refresh the examples list from the info.yaml file.
        """
        info_file=os.path.join(self.path,'info.yaml')
        if not os.path.isfile(info_file):
            raise FileNotFoundError(f'Example path {self.path} does not contain info.yaml')
        with open(info_file,'r') as f:
            self.info=yaml.safe_load(f)
        if 'examples' not in self.info:
            raise KeyError(f'info.yaml in {self.path} does not contain examples key')
        self.examples_list=ExampleList.from_list_of_dicts(self.info['examples'])

    def _write_info(self):
        """
        Write the current examples list to the info.yaml file.  Since self.examples_list is a list of Example objects, it needs to be converted back to a dictionary format.
        
        This method is called internally to save the current state of the examples list to the info.yaml file.
        """
        info_file=os.path.join(self.path,'info.yaml')
        saveme=dict(examples=self.examples_list.to_list_of_dicts())
        with open(info_file,'w') as f:
            yaml.dump(saveme,f,default_flow_style=False)
        logger.debug(f'Wrote info.yaml to {info_file}')

    def checkout_example_yaml(self,index:int):
        """
        Check out an example YAML file by its index and copy it to the current working directory.
        
        Parameters
        ----------
        index : int
            The index of the example to check out. The index is 1-based, meaning the first example has index 1.
        
        Returns
        -------
        str
            The name of the example YAML file that was checked out.

        Raises
        ------
        IndexError
            If the index is out of range for the examples list.
        FileNotFoundError
            If the example YAML file does not exist in the specified path.
        """
        if index < 1 or index > len(self.examples_list):
            raise IndexError(f'Index {index} is out of range for examples list of length {len(self.examples_list)}')
        example_yaml=self.examples_list[index-1].name+'.yaml'
        example_yaml_path=os.path.join(self.path,example_yaml)
        if not os.path.isfile(example_yaml_path):
            raise FileNotFoundError(f'Example YAML file {example_yaml_path} does not exist')
        shutil.copy(example_yaml_path,os.getcwd())
        logger.info(f'Checked out example {index} from {self.path} to current working directory {os.getcwd()}')
        return example_yaml # return the name of the copied file
    
    def new_example_yaml(self,id='ABCD',build_type='minimal'):
        if len(id)==4: # assume a PDB id
            idtype='PDB'
        elif id.startswith('P'): # assume an alphafold id by uniprot id
            idtype='Alphafold'
        else:
            raise ValueError(f'Invalid id {id} for new example YAML; must be a 4-letter PDB ID or an Alphafold/UNIPROT ID starting with "P"')
        example_yaml=self.examples_list[0].name+'.yaml'
        example_yaml_path=os.path.join(self.path,example_yaml)
        with open(example_yaml_path, 'r') as f:
            try:
                example_config=yaml.safe_load(f)
            except yaml.YAMLError as e:
                raise ValueError(f'Invalid YAML file {example_yaml_path}: {e}')
        if build_type=='minimal':
            psfgen_task=example_config['tasks'][0]
            example_config['tasks']=[psfgen_task]  # keep only the psfgen task
        example_config['title']=f'New template pestifer config for id {id} ({idtype})'
        if idtype=='PDB':
            example_config['tasks'][0]['psfgen']['source']['id']=id
        elif idtype=='Alphafold':
            del example_config['tasks'][0]['psfgen']['source']['id']
            example_config['tasks'][0]['psfgen']['source']['alphafold']=id
        if build_type=='full':
            example_config['tasks'][-1]['terminate']['basename']=f'my_{id.lower()}'
            example_config['tasks'][-1]['terminate']['package']['basename']=f'my_{id.lower()}'
        output_yaml=os.path.join(os.getcwd(),f'{id.lower()}.yaml')
        with open(output_yaml, 'w') as f:
            yaml.dump(example_config, f, default_flow_style=False)
        logger.info(f'Generated new example YAML file {output_yaml} for id {id} ({idtype})')

    def report_examples_list(self,header=False,formatter=r'{:>7s}  {:>8s}  {:<30s}  {}'):
        """
        Generate a report of the available examples in the examples list.
        
        Parameters
        ----------
        header : bool, optional
            If True, include a header in the report. Default is False.
        formatter : str, optional
            A format string for the report. Default is a string that formats the index, name, and description of each example. Default is ``r'{:>7s}    {:<30s}    {}'``.
        
        Returns
        -------
        str
            A report of the available examples in the examples list.
        """
        if not self.examples_list:
            return 'No examples available'
        if header:
            report_lines = [formatter.format('Index','ID','Name','Description')+'\n']
        else:
            report_lines = []
        for i, example in enumerate(self.examples_list):
            report_lines.append(example.report_line(formatter=formatter))
        return ''.join(report_lines)


    def new_example(self,yaml_file_name:str,description:str='',pdbID:str=''):
        """
        Grab the information from an example YAML file along with description and return it as an Example object.  This method allows the caller to specify a description and PDB ID, or it will extract these from the YAML file if not provided.  It also checks that the YAML file is valid and contains the necessary fields.

        Parameters
        ----------
        yaml_file_name : str
            The name of the example file to grab information from. This should be a valid file path.
        description : str, optional
            A description of the example; if not provided, extracts the ``title`` field from the YAML file.
        pdbID : str, optional
            The PDB ID associated with the example; if not provided, extracts the ``id`` field from the ``psfgen`` task of the ``tasks`` list in the YAML file.

        Returns
        -------
        Example
            An Example object containing the name, description, and pdbID of the example.

        Raises
        ------
        ValueError
            If the name, description, or pdbID is not provided. 
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
                description = 'No description provided'
        if not pdbID:
            if 'id' in new_config.get('tasks', [{}])[0].get('psfgen', {}).get('source', {}):
                pdbID = new_config['tasks'][0]['psfgen']['source']['id']
            else:
                pdbID = 'No PDB ID provided'
        input_dict= {
            'name': os.path.splitext(os.path.basename(yaml_file_name))[0],
            'description': description,
            'pdbID': pdbID
        }
        return Example(**input_dict)
    
    def add_example(self,name:str,pdbID:str='',description:str=''):
        """
        Add a new example to the examples list.

        Parameters
        ----------
        name : str
            The name of the example to add.  It should be a valid file path to a YAML file, with or without the `.yaml` extension.
        pdbID : str
            The PDB ID associated with the example; if not provided, extracts the ``id`` field from the ``psfgen`` task of the ``tasks`` list in the YAML file.
        description : str
            A description of the example; if not provided, extracts the ``title`` field from the YAML file.

        Returns
        -------
        int
            The index of the newly added example in the examples list (1-based index).

        Raises
        ------
        ValueError
            If the name, description, or pdbID is not provided. 
        FileNotFoundError
            If the example file does not exist at the specified path.
        """
        return self.insert_example(len(self.examples_list)+1,name,description,pdbID)
    
    def delete_example(self,index:int):
        """
        Delete an example from the examples list by its index.
        
        Parameters
        ----------
        index : int
            The index of the example to delete (1-based).

        Returns
        -------
        int
            The index of the deleted example (1-based).

        Raises
        ------
        IndexError
            If the index is out of range for the examples list.
        FileNotFoundError
            If the example file does not exist in the specified path.
        """
        if index < 1 or index > len(self.examples_list):
            raise IndexError(f'Index {index} is out of range for examples list of length {len(self.examples_list)}')
        real_index = index - 1  # convert to zero-based index
        example=self.examples_list[real_index]
        self.examples_list.remove(example)  # remove the example from the list
        example_file=os.path.join(self.path,example.name+'.yaml')
        if os.path.isfile(example_file):
            logger.debug(f'Deleting example file "{example_file}"')
            os.remove(example_file)
        else:
            raise FileNotFoundError(f'Example file "{example_file}" does not exist')
        self._write_info()
        if self.sphinx_example_manager:
            self.sphinx_example_manager.delete_example(example)
        logger.info(f'Deleted example {index}: {example.name}')
        return example
    
    def insert_example(self,index:int,name:str,description:str='',pdbID:str=''):
        """
        Insert a new example into the examples list at a specified index.
        
        Parameters
        ----------
        index : int
            The index at which to insert the new example (1-based).
        name : str
            The name of the example file to insert. This should be a valid file path.
        description : str
            A description of the example.
        pdbID : str
            The PDB ID associated with the example.

        Returns
        -------
        int
            The index of the newly inserted example in the examples list (1-based).
            
        Raises
        ------
        IndexError
            If the index is out of range for the examples list.
        ValueError
            If the name, description, or pdbID is not provided.
        FileNotFoundError
            If the example file does not exist at the specified path.
        """
        yaml_file_path=name if name.endswith('.yaml') else name + '.yaml'
        new_example=self.new_example(yaml_file_path,description,pdbID)
        new_example.index = index  # set the index for the new example
        real_index = index - 1  # convert to zero-based index
        if real_index < 0 or real_index > len(self.examples_list):
            raise IndexError(f'Index {index} is out of range for examples list of length {len(self.examples_list)}')
        shutil.copy(yaml_file_path, self.path)
        self.examples_list.insert(real_index, new_example)
        self._write_info()
        if self.sphinx_example_manager:
            self.sphinx_example_manager.insert_example(index, new_example)
        logger.info(f'Inserted new example at index {index}: {new_example.name}')
        return index

    def update_example(self,index:int,name:str,description:str='',pdbID:str=''):
        """
        Update an existing example in the examples list by its index.  This will edit the entry at position ``index`` in the toctree of "examples.rst" to have name ``name``, delete the existing example YAML file in resources/examples, delete the rst documentation file for the old example, copy the new example YAML file to resources/examples, and generate a new rst documentation file for the updated example.  The index is 1-based, meaning the first example has index 1.  The operations involving files in docs/source/examples are performed by the ``sphinx_example_manager`` if it is initialized.
        
        Parameters
        ----------
        index : int
            The index of the example to update (1-based).
        name : str
            The name of the example file to update.  This can replace the existing example file, or it can be a new file.
            It should be a valid file path to a YAML file, with or without the `.yaml` extension.
        description : str
            A description of the example.
        pdbID : str
            The PDB ID associated with the example.

        Returns
        -------
        Example
            The updated Example object.

        Raises
        ------
        IndexError
            If the index is out of range for the examples list.
        ValueError
            If the name, description, or pdbID is not provided.
        FileNotFoundError
            If the example file does not exist at the specified path.
        """
        if index < 1 or index > len(self.examples_list):
            raise IndexError(f'Index {index} is out of range for examples list of length {len(self.examples_list)}')
        if not name:
            raise ValueError('Name must be provided')
        yaml_file_path=name if name.endswith('.yaml') else name + '.yaml'
        if not os.path.isfile(yaml_file_path):
            raise FileNotFoundError(f'Example file {yaml_file_path} does not exist')
        real_index = index - 1  # convert to zero-based index
        if real_index >= len(self.examples_list):
            raise IndexError(f'Index {index} is out of range for examples list of length {len(self.examples_list)}')
        elif real_index < 0:
            raise IndexError(f'Index {index} is out of range for examples list of length {len(self.examples_list)}')
        # get the existing example to update
        example=self.examples_list[real_index]
        # delete the existing example file if it exists
        example_file=os.path.join(self.path,example.name+'.yaml')
        if os.path.isfile(example_file):
            logger.debug(f'Deleting existing example file "{example_file}"')
            os.remove(example_file)
        new_example=self.new_example(yaml_file_path,description,pdbID)
        new_example.index = index  # set the index for the new example
        self.examples_list[real_index] = new_example  # replace the old example with the new one
        shutil.copy(yaml_file_path, self.path)  # copy the new example
        self._write_info()
        if self.sphinx_example_manager:
            self.sphinx_example_manager.update_example(index, new_example)
        logger.info(f'Updated example at index {index}: {new_example.name}')
        return new_example