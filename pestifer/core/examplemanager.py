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
            raise FileNotFoundError(f'Directory "{example_path}" containing example folders does not exist or is not a directory')
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
        with open(info_file,'w',encoding='utf-8') as f:
            yaml.dump(saveme,f,default_flow_style=False)
        logger.debug(f'Wrote info.yaml to {info_file}')

    def checkout_example(self,index:int):
        """
        Copy example YAML file and associated companion files by example index to the current working directory.

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
        example= self.examples_list[index-1]  # convert to zero-based index
        example_folder=example.name
        example_folder_path=os.path.join(self.path,example_folder)
        if not os.path.isdir(example_folder_path):
            raise FileNotFoundError(f'Example folder {example_folder} does not exist in {self.path}')
        example_yaml=example.name+'.yaml'
        example_yaml_path=os.path.join(example_folder_path,example_yaml)
        if not os.path.isfile(example_yaml_path):
            raise FileNotFoundError(f'Example YAML file {example_yaml_path} does not exist')
        # copy all files in the example_yaml_path directory to the current working directory
        shutil.copy(example_yaml_path, os.getcwd())
        # copy all companion files to the current working directory
        companion_files=self.examples_list[index-1].companion_files
        if companion_files:
            for companion_file in companion_files:
                companion_file_path=os.path.join(self.path,example_folder,companion_file)
                if os.path.isfile(companion_file_path):
                    shutil.copy(companion_file_path, os.getcwd())
                else:
                    logger.warning(f'Declared companion file "{companion_file}" does not exist in {os.path.join(self.path,example_folder)}')
        logger.info(f'Checked out example {index} from {self.path} to current working directory {os.getcwd()}')
        return example_yaml # return the name of the copied file
    
    def new_example_yaml(self,id='ABCD',build_type='minimal'):
        """
        Generate a new example YAML file based on an existing example template.  The id can be a 4-letter PDB ID or an Alphafold/UNIPROT ID starting with "P".  The build_type can be either 'minimal' or 'full', which determines whether the generated YAML file contains only the psfgen task or all tasks including termination.

        Parameters
        ----------
        id : str, optional
            The ID for the new example YAML file. It can be a 4-letter PDB ID or an Alphafold/UNIPROT ID starting with "P". Default is 'ABCD'.
        build_type : str, optional
            The type of build for the new example YAML file. It can be either 'minimal' or 'full'. Default is 'minimal'.

        """
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
            report_lines.append(example.report_line(formatter=formatter)+'\n')
        return ''.join(report_lines)


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
        example_folder_path=os.path.join(self.path,example.name)
        if os.path.isdir(example_folder_path):
            logger.debug(f'Deleting example folder "{example_folder_path}"')
            shutil.rmtree(example_folder_path)
        else:
            raise FileNotFoundError(f'Example folder "{example_folder_path}" does not exist')
        self._write_info()
        if self.sphinx_example_manager:
            self.sphinx_example_manager.delete_example(example)
        logger.info(f'Deleted example {index}: {example.name}')
        return example
    


    def checkin_example(self,example: Example):
        """
        Check in an Example instance by copying its YAML file and companion files to the appropriate example folder, which is created if it doesn't already exist.  This will overwrite the existing files.  If any file referenced by the example does not exist in the current working directory, a warning is logged but no action taken.

        Parameters
        ----------
        example : Example
            The Example instance to check in.
        """
        example_folder=os.path.join(self.path,example.name)
        user_yaml_file_path=example.name+'.yaml' # correct; the name attribute should never have an extension
        if not os.path.isdir(example_folder):
            logger.debug(f'Creating example folder "{example_folder}"')
            if not os.path.isfile(user_yaml_file_path):
                raise FileNotFoundError(f'Example YAML file {user_yaml_file_path} does not exist in your current working directory {os.getcwd()}')
            os.makedirs(example_folder)
            shutil.copy(user_yaml_file_path, example_folder)
        else:
            if not os.path.isfile(user_yaml_file_path):
                logger.debug(f'Example YAML file {user_yaml_file_path} does not exist in your current working directory {os.getcwd()}')
                logger.debug(f'Pestifer will now check that this YAML file is already in the example folder {example_folder}')
                existing_yaml_file_path = os.path.join(example_folder, user_yaml_file_path)
                if not os.path.isfile(existing_yaml_file_path):
                    raise FileNotFoundError(f'Example YAML file {existing_yaml_file_path} does not exist in the example folder {example_folder}')
            else:
                logger.debug(f'Copying example YAML file "{user_yaml_file_path}" to example folder "{example_folder}"')
                shutil.copy(user_yaml_file_path, example_folder)
        for f in example.companion_files:
            if os.path.isfile(f):
                shutil.copy(f, example_folder)
            else:
                logger.warning(f'Companion file {f} does not exist in {os.getcwd()}')

    def insert_example(self,index:int,name:str,description:str='',pdbID:str='',author_name:str='',author_email:str='',companion_files: list = []):
        """
        Insert a new example into the examples list at a specified index.
        
        Parameters
        ----------
        index : int
            The index at which to insert the new example (1-based).
        name : str
            The name of the example file to insert. This should be a valid file path since a new Example is created.
        description : str
            A description of the example.
        pdbID : str
            The PDB ID associated with the example.
        author_name : str
            The name of the author of the example.
        author_email : str
            The email of the author of the example.
        companion_files : list, optional
            A list of companion files associated with the example; defaults to an empty list.

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
        new_example=Example.from_yaml(yaml_file_path,description=description,pdbID=pdbID,author_name=author_name,author_email=author_email,companion_files=companion_files) 
        new_example.index = index  # set the index for the new example
        real_index = index - 1  # convert to zero-based index
        if real_index < 0 or real_index > len(self.examples_list):
            raise IndexError(f'Index {index} is out of range for examples list of length {len(self.examples_list)}')
        self.checkin_example(new_example)  # check in the new example by copying its YAML file and companion files to the appropriate example folder
        self.examples_list.insert(real_index, new_example)
        self._write_info()
        if self.sphinx_example_manager:
            self.sphinx_example_manager.insert_example(index, new_example)
        logger.info(f'Inserted new example at index {index}: {new_example.name}')
        return index

    def update_example(self,index:int,name:str='',description:str='',pdbID:str='',author_name:str='',author_email:str='',companion_files: list = []):
        """
        Update an existing example in the examples list by its unique index.
        
        * if name is not given, then
           * overwrite attributes ``description``, ``pdbID``, ``author_name``, and ``author_email`` of the existing example at index with the 
             provided values, if they are not empty.  The name of the example is unchanged.
        * otherwise, if name is given, then
           * if name matches name of example at index
              * overwrite attributes ``description``, ``pdbID``, ``author_name``, and ``author_email`` of the existing example at index with the 
                provided values, if they are not empty.
              * if there is a file <name>.yaml in the user's cwd, then
                * check this yaml file in, overwriting the existing YAML file (DANGER!)
              * otherwise if there is NOT a file <name>.yaml in the user's cwd, then
                * congratulate the user on their diligence in providing both the unique index and unique name of the existing example they want to 
                  modify, and proceed with the update.
           * otherwise, if name does not match name of example at index
              * if there is already an example with this name, then raise an error
              * otherwise, if there is not already an example with this name, then
                * if there is a file <name>.yaml in the user's cwd, then
                  * assume the user is trying to create a new example by completely overwriting the existing example at index; i.e.,
                    * rename the existing example at index to <name> and overwrite all attributes ``name``, ``description``, ``pdbID``, ``author_name``, and ``author_email`` of the existing example at index with the provided values, if they are not empty.
                * otherwise, if there is NOT a file <name>.yaml in the user's cwd, then
                  * assume the user is simply trying to rename the existing example at index
                  * rename the example at index; i.e., overwrite all attributes ``name``, ``description``, ``pdbID``, ``author_name``, and ``author_email`` of the existing example at index with the provided values, if they are not empty (``name`` is necessarily not empty in this branch).

        Parameters
        ----------
        index : int
            The index of the example to update (1-based).
        name : str
            The name of the example to update.  It should not have an extension, so if it does, it is stripped off.
        description : str
            A description of the example.
        pdbID : str
            The PDB ID (or Alphafold ID) associated with the example.
        author_name : str
            The name of the author of the example.
        author_email : str
            The email of the author of the example.
        companion_files : list, optional
            A list of companion files associated with the example; defaults to an empty list.
            
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
        # get the name of the example with this index by directly querying the examples list
        current_example=self.examples_list[index-1]
        current_example_name=current_example.name
        current_example_folder=current_example_name
        name=name.replace('.yaml','')  # strip .yaml from the name if it is there
        if name == '' or name == current_example_name:
            # if name is '', then the user just wants to update the existing example using the other parameters; there will be no change to the name of the example folder or the YAML file.
            new_example=current_example.update_in_place(description=description,pdbID=pdbID,author_name=author_name,author_email=author_email,companion_files=companion_files)
            self.checkin_example(new_example)
        else:
            # first, check if there is already an example with this name; if so, raise an error
            for i, ex in enumerate(self.examples_list):
                if ex.name == name:
                    raise ValueError(f'You have named an existing example whose index does not match the index you provided: {i+1} ({ex.name})')
            if os.path.isfile(os.path.join(os.getcwd(),name+'.yaml')):
                # if there is a file <name>.yaml in the user's cwd, then assume the user is trying to create a new example; bail and instruct user to use --example-action add instead.
                logger.error(f'The name you provided matches a YAML file in your CWD, but does not match the name of the example at index {index} ({current_example_name})')
                logger.error(f'If you are trying to create a new example, please use --example-action add instead of --example-action update')
                logger.error(f'If you are trying to rename the existing example at index {index}, please be sure to run pestifer from a directory that does not contain a file named {name}.yaml')
                raise FileExistsError(f'Improper use of --example-action update')
            current_example_folder_path=os.path.join(self.path,current_example_folder)
            current_example.update_in_place(description=description,pdbID=pdbID,author_name=author_name,author_email=author_email,companion_files=companion_files)
            current_example.name=name
            renamed_example_folder_path=os.path.join(self.path,name)
            # rename the folder
            os.rename(current_example_folder_path, renamed_example_folder_path)
            # perform a checkin to overwrite any files if they are in the user's cwd
            self.checkin_example(current_example)
        self._write_info()
        if self.sphinx_example_manager:
            self.sphinx_example_manager.update_example(index, new_example)
        logger.info(f'Updated example at index {index}: {new_example.name}')
        return new_example
    
    def add_example(self,name:str,pdbID:str='',description:str='',author_name:str='',author_email:str='',companion_files: list = []):
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
        author_name : str
            The name of the author of the example; if not provided, defaults to an empty string
        author_email : str
            The email of the author of the example; if not provided, defaults to an empty string
        companion_files : list, optional
            A list of companion files associated with the example; defaults to an empty list.

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
        return self.insert_example(len(self.examples_list)+1,name,description=description,pdbID=pdbID,author_name=author_name,author_email=author_email,companion_files=companion_files)

    def rename_example(self,index:int,new_name:str):
        """
        Rename an example in the examples list by its index.
        
        Parameters
        ----------
        index : int
            The index of the example to rename (1-based).
        new_name : str
            The new name for the example. This should be a valid file path.

        Returns
        -------
        Example
            The renamed Example object.

        Raises
        ------
        IndexError
            If the index is out of range for the examples list.
        ValueError
            If the new name is not provided.
        FileNotFoundError
            If the example file does not exist at the specified path.
        """
        self.update_example(index, name=new_name)  # use update_example to handle renaming logic
        logger.info(f'Renamed example {index} from {self.examples_list[index-1].name} to {new_name}')

    def set_example_author(self,index:int,author_name:str,author_email:str):
        """
        Set the author information for an example in the examples list by its index.
        
        Parameters
        ----------
        index : int
            The index of the example to set the author for (1-based).
        author_name : str
            The name of the author.
        author_email : str
            The email address of the author.

        Returns
        -------
        Example
            The updated Example object with the new author information.

        Raises
        ------
        IndexError
            If the index is out of range for the examples list.
        ValueError
            If the author name or email is not provided.
        """
        self.update_example(index,author_name=author_name, author_email=author_email)
        logger.info(f'Set author for example at index {index}: {self.examples_list[index-1].name} by {author_name} <{author_email}>')