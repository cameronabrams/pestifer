# Author: Cameron F. Abrams <cfa22@drexel.edu>

""" 
A module for managing example input files in pestifer.
This module provides the :class:`ExampleManager` class, which allows users to check out example YAML files,
report the list of examples, and manage example resources.

ExampleManager does not directly manage the example documentation, but it relies on the :class:`SphinxExampleManager` 
class to handle the documentation side of things.
"""

import logging
import os
import shutil
import yaml

from pathlib import Path

from .artifacts import *
from .example import Example, ExampleList
from ..sphinxext.sphinx_examplemanager import SphinxExampleManager
from ..util.formatvalidator import FormatValidator

logger = logging.getLogger(__name__)

class ExampleManager:
    """
    A class for managing example YAML input files in pestifer.
    This class provides methods to check out example YAML files, report the list of examples, and manage example resources.
    Each example has a root folder that lives in `examples_path`.
    
    Parameters
    ----------
    example_path : str
        The path to the directory containing example root folders.
    docs_path : str
        The path to the directory containing documentation files. This is passed directly to the :class:`SphinxExampleManager` class.
        If not provided, the documentation management is not enabled.
    """

    def __init__(self, examples_path: str | Path | None = None, 
                       docs_source_path: str | Path | None = None):
        if not examples_path:
            raise ValueError('You must provide a path to the directory containing example input files')
        if not os.path.isdir(examples_path):
            logger.debug(f'Directory "{examples_path}" does not exist; creating it')
            os.makedirs(examples_path)
        self.root_folder_format_validator = FormatValidator(Example.folder_name_format, greedy=True, flexible_ws=True)
        self.path = Path(examples_path).absolute()
        self.examples: ExampleList = ExampleList([])
        for efolder in self.path.iterdir():
            if efolder.is_dir() and self.root_folder_format_validator.fullmatch(efolder.name):
                example_id = self.root_folder_format_validator.extract(efolder.name).get("example_id")
                inputsfolder = efolder / Example.inputs_subdir
                yaml_inputs = list(inputsfolder.glob('*.yaml')) if inputsfolder.is_dir() else []
                title = 'No title provided'
                db_id = None
                shortname = None
                if len(yaml_inputs) == 0:
                    logger.warning(f'No YAML files found in {inputsfolder}')
                elif len(yaml_inputs) > 1:
                    logger.warning(f'Multiple YAML files found in {inputsfolder}: {yaml_inputs}')
                else:
                    main_input_script = yaml_inputs[0]
                    shortname = main_input_script.name.replace('.yaml', '')
                    title, db_id = Example._get_implied_metadata(main_input_script)
                if shortname:
                    self.examples.append(Example(
                        example_id=example_id,
                        title=title,
                        shortname=shortname,
                        db_id=db_id
                    ))
        if docs_source_path:
            # create the SphinxExampleManager instance if docs_source_path is provided
            self.sphinx_example_manager = SphinxExampleManager(docs_source_path=Path(docs_source_path).absolute())
        else:
            self.sphinx_example_manager = None

    def examplefolderpath(self, example: Example) -> Path:
        return self.path / example.rootfolderpath

    def scriptpath(self, example: Example) -> Path:
        return self.path / example.scriptpath
    
    def inputspath(self, example: Example) -> Path:
        return self.path / example.inputspath
    
    def outputspath(self, example: Example) -> Path:
        return self.path / example.outputspath

    def checkout_example(self, example_id: int) -> Example:
        """
        Copy example YAML file and associated companion files by example ID to the current working directory.

        Parameters
        ----------
        example_id : int
            The ID of the example to check out.

        Returns
        -------
        Example
            The example that was checked out.

        Raises
        ------
        IndexError
            If the index is out of range for the examples list.
        FileNotFoundError
            If the example YAML file does not exist in the specified path.
        """
        example = self.examples.get_example_by_example_id(example_id)
        if not example:
            raise IndexError(f'Example with ID {example_id} not found')
        with open(example.scriptname, 'w') as f:
            f.write(self.scriptpath(example).read_text())
        for aux_path in self.inputspath(example).glob('*'):
            if aux_path.name != example.scriptname and aux_path.is_file():
                shutil.copy(aux_path, os.getcwd())
        logger.info(f'Checked out example {example_id} from {self.path.name} to current working directory {os.getcwd()}')
        return example

    def new_example_yaml(self, db_id: str = 'ABCD', build_type: str = 'minimal', outputfilename: str = None, title: str = ''):
        """
        Generate a new example YAML file based on an existing example template.  The id can be a 4-letter PDB ID or an Alphafold/UNIPROT ID starting with "P".  The build_type can be either 'minimal' or 'full', which determines whether the generated YAML file contains only the psfgen task or all tasks including termination.

        Parameters
        ----------
        db_id : str, optional
            The ID for the new example YAML file. It can be a 4-letter PDB ID or an Alphafold/UNIPROT ID starting with "P". Default is 'ABCD'.
        build_type : str, optional
            The type of build for the new example YAML file. It can be either 'minimal' or 'full'. Default is 'minimal'.

        """
        if len(db_id) == 4:  # assume a PDB id
            idtype = 'PDB'
        elif db_id.startswith('P') or db_id.startswith('O'):  # assume an alphafold id by uniprot id
            idtype = 'Alphafold'
        else:
            raise ValueError(f'Invalid id {db_id} for new example YAML; must be a 4-letter PDB ID or an Alphafold/UNIPROT ID starting with "P"')
        example_yaml_path = self.scriptpath(self.examples[0])
        with open(example_yaml_path, 'r') as f:
            try:
                example_config = yaml.safe_load(f)
            except yaml.YAMLError as e:
                raise ValueError(f'Invalid YAML file {example_yaml_path}: {e}')
        example_config['title'] = title if title else f'New template pestifer config for id {db_id} ({idtype})'
        if build_type == 'minimal':
            fetch_task = example_config['tasks'][0]
            example_config['tasks'] = [fetch_task, example_config['tasks'][1]]  # keep only the fetch task
        if idtype == 'PDB' or idtype == 'Alphafold':
            example_config['tasks'][0]['fetch']['sourceID'] = db_id
        if build_type == 'full':
            example_config['tasks'][-1]['terminate']['basename'] = f'my_{db_id.lower()}'
            example_config['tasks'][-1]['terminate']['package']['basename'] = f'my_{db_id.lower()}'
        if outputfilename:
            output_yaml = outputfilename
        else:
            output_yaml = f'{db_id.lower()}.yaml'
        with open(output_yaml, 'w') as f:
            yaml.dump(example_config, f, default_flow_style=False, sort_keys=False)
        logger.info(f'Generated new example YAML file {output_yaml} for id {db_id} ({idtype})')

    def report_examples(self, header=False, formatter=r'{:>7s}  {:>8s}  {:<30s}  {}'):
        """
        Generate a report of the available examples in the examples list.
        
        Parameters
        ----------
        header : bool, optional
            If True, include a header in the report. Default is False.
        formatter : str, optional
            A format string for the report. Default is a string that formats the index, name, and title of each example. Default is ``r'{:>7s}    {:<30s}    {}'``.
        
        Returns
        -------
        str
            A report of the available examples in the examples list.
        """
        if not self.examples:
            return 'No examples available'
        if header:
            report_lines = [formatter.format('ID','DBID','Name','Title')+'\n']
        else:
            report_lines = []
        for i, example in enumerate(self.examples):
            report_lines.append(example.report_line(formatter=formatter)+'\n')
        return ''.join(report_lines)

    def delete_example(self, example_id: int):
        """
        Delete an example from the examples list by its ID and returns the instance.

        Parameters
        ----------
        example_id : int
            The example ID of the example to delete.

        Returns
        -------
        Example
            The deleted example instance.

        Raises
        ------
        IndexError
            If the index is out of range for the examples list.
        FileNotFoundError
            If the example file does not exist in the specified path.
        """

        example = self.examples.get_example_by_example_id(example_id)

        if not example:
            raise IndexError(f'Example ID {example_id} is out of range for examples list of length {len(self.examples)}')

        self.examples.remove(example)  # remove the example from the list
        example_folder_path = os.path.join(self.path, example.rootfolderpath)
        if os.path.isdir(example_folder_path):
            logger.debug(f'Deleting example folder "{example_folder_path}"')
            shutil.rmtree(example_folder_path)
        else:
            raise FileNotFoundError(f'Example folder "{example_folder_path}" does not exist')
        if self.sphinx_example_manager:
            self.sphinx_example_manager.delete_example(example)
        logger.info(f'Deleted example {example_id}: {example.shortname}')
        return example

    def checkin_example(self, example: Example, overwrite: bool = False):
        """
        Check in an Example instance by copying its YAML file and companion files to the appropriate example folder, which is created if it doesn't already exist.  This will overwrite the existing files.  If any file referenced by the example does not exist in the current working directory, a warning is logged but no action taken.

        Parameters
        ----------
        example : Example
            The Example instance defining the destination folder and files to copy.
        overwrite : bool, optional
            If True, overwrite existing files in the example folder. Default is False.
        """
        example_folder = self.examplefolderpath(example)
        example_inputs_subfolder = self.inputspath(example)
        example_outputs_subfolder = self.outputspath(example)
        user_yaml_file_path = Path(example.shortname + '.yaml')  # correct; the name attribute should never have an extension
        if not example_folder.is_dir(): # we are checking in a new example
            logger.debug(f'Creating example folder "{example_folder}"')
            example_folder.mkdir()
            example_inputs_subfolder.mkdir()
            example_outputs_subfolder.mkdir()
            if not user_yaml_file_path.is_file():
                raise FileNotFoundError(f'Example YAML file {user_yaml_file_path.name} does not exist in your current working directory {os.getcwd()}')
            shutil.copy(user_yaml_file_path, example_inputs_subfolder)
            for f in example.auxiliary_inputs:
                if os.path.isfile(f):
                    shutil.copy(f, example_inputs_subfolder)
                else:
                    logger.warning(f'Declared auxiliary input file {f} does not exist in {os.getcwd()}')
            for f in example.outputs:
                if os.path.isfile(f):
                    shutil.copy(f, example_outputs_subfolder)
                else:
                    logger.warning(f'Declared output file {f} does not exist in {os.getcwd()}')
        else:
            existing_yaml_file_path = self.scriptpath(example)
            if overwrite and user_yaml_file_path.is_file():
                shutil.copy(user_yaml_file_path, example_inputs_subfolder)
            for f in example.auxiliary_inputs:
                existing_aux_path = example_inputs_subfolder / f
                if overwrite and existing_aux_path.is_file():
                    os.remove(existing_aux_path)
                    shutil.copy(f, example_inputs_subfolder)
                elif not existing_yaml_file_path.exists():
                    shutil.copy(f, example_inputs_subfolder)
            for f in example.outputs:
                existing_output_path = example_outputs_subfolder / f
                if overwrite and existing_output_path.is_file():
                    os.remove(existing_output_path)
                    shutil.copy(f, example_outputs_subfolder)
                elif not existing_output_path.exists():
                    shutil.copy(f, example_outputs_subfolder)

    def append_example(self, example_id: int, scriptname: str, title: str = '', db_id: str = '', author_name: str = '', author_email: str = '', auxiliary_inputs: list = [], outputs: list = []):
        """
        Create a new Example from keyword arguments and add it to the examples list.
        
        Parameters
        ----------
        example_id : int
            The unique integer ID of the example.
        scriptname : str
            The name of the example script file to insert. This should be a valid file path since a new Example is created.
        title : str
            A title of the example; if not provided, it is extracted from the script.
        db_id : str
            The structure/sequence database ID associated with the example.
        author_name : str
            The name of the author of the example.
        author_email : str
            The email of the author of the example.
        auxiliary_inputs: list[str]
            A list of auxiliary input files associated with the example.
        outputs: list[str]
            A list of output files associated with the example.

        Returns
        -------
        Example:
            The newly created Example object.

        Raises
        ------
        FileNotFoundError
            If the example scriptfile does not exist in the CWD.
        """

        if not os.path.isfile(scriptname):
            raise FileNotFoundError(f'Example scriptfile {scriptname} does not exist in the CWD.')

        new_example = Example(
            example_id=example_id,
            shortname=os.path.splitext(scriptname)[0],
            title=title,
            db_id=db_id,
            author_name=author_name,
            author_email=author_email,
            auxiliary_inputs=auxiliary_inputs,
            outputs=outputs
        )
        self.checkin_example(new_example)  # check in the new example by copying its YAML file and companion files to the appropriate example folder
        self.examples.append(new_example)
        if self.sphinx_example_manager:
            self.sphinx_example_manager.append_example(new_example)
        logger.info(f'Added new example {new_example.example_id}: {new_example.shortname}')
        return new_example

    def update_example(self, example_id: int, shortname: str = '', title: str = '', db_id: str = '', author_name: str = '', author_email: str = '', auxiliary_inputs: list = [], outputs: list = [], skip_sphinx: bool = False):
        """
        Update an existing example in the examples list by its unique index.

        Parameters
        ----------
        example_id : int
            The unique ID of the example to update.
        shortname : str
            The short name of the example to update (without extension).
        title : str
            A title of the example; overrides value of 'title' in the <name>.yaml if it exists.
        pdbID : str
            The PDB ID (or Alphafold ID) associated with the example; overrides value of 'id' or 'alphafold' in the <name>.yaml if it exists.
        author_name : str
            The name of the author of the example; overrides the name in the # Author line of <name>.yaml if it exists.
        author_email : str
            The email of the author of the example; overrides the email in the # Author line of <name>.yaml if it exists.
        companion_files : list, optional
            A list of companion files associated with the example; defaults to an empty list.
        skip_sphinx : bool, optional
            If True, skip updating the Sphinx documentation for this example.

        Returns
        -------
        Example or None
            The updated Example object if the update was successful, or None if the example was not found.

        Raises
        ------

        FileNotFoundError
            If the example file does not exist at the specified path.
        """
        current_example: Example = self.examples.get_example_by_example_id(example_id)
        if not current_example:
            logger.warning(f'No example found with ID {example_id}')
            return None
        if shortname:
            user_script = Path(f'{shortname}.yaml')  # in cwd
            if user_script.is_file():
                # user is providing a new script for this example
                self.scriptpath(current_example).unlink()  # remove the old script
                shutil.copy(user_script, self.inputspath(current_example))
                user_implied_title, user_implied_db_id = Example._get_implied_metadata(user_script)
                current_example.title = user_implied_title
                current_example.db_id = user_implied_db_id
            else:
                # user is not providing a new script, they just want to change the shortname
                self.scriptpath(current_example).rename(self.inputspath(current_example) / f'{shortname}.yaml')
            current_example.shortname = shortname

        if title:
            current_example.title = title
        if db_id:
            current_example.db_id = db_id
        if author_name:
            current_example.author_name = author_name
        if author_email:
            current_example.author_email = author_email
        if auxiliary_inputs:
            if current_example.auxiliary_inputs is not None:
                for f in current_example.auxiliary_inputs:
                    Path(f).unlink(missing_ok=True)  # remove old auxiliary inputs
            current_example.auxiliary_inputs = auxiliary_inputs
            for f in auxiliary_inputs:
                if os.path.isfile(f):
                    shutil.copy(f, self.inputspath(current_example))
        if outputs:
            if current_example.outputs is not None:
                for f in current_example.outputs:
                    Path(f).unlink(missing_ok=True)
            current_example.outputs = outputs
            for f in outputs:
                if os.path.isfile(f):
                    shutil.copy(f, self.outputspath(current_example))
        if self.sphinx_example_manager and not skip_sphinx:
            self.sphinx_example_manager.update_example(example_id, current_example)
        logger.info(f'Updated example {example_id}: {current_example.shortname}')
        return current_example

    def rename_example(self, example_id: int, new_shortname: str):
        """
        Rename an example in the examples list by its ID.

        Parameters
        ----------
        index : int
            The ID of the example to rename.
        new_shortname : str
            The new short name for the example.

        Returns
        -------
        Example
            The renamed Example object.

        """
        return self.update_example(example_id, shortname=new_shortname)  # use update_example to handle renaming logic

    def set_example_author(self, example_id: int, author_name: str, author_email: str):
        """
        Set the author information for an example in the examples list by its ID.

        Parameters
        ----------
        example_id : int
            The ID of the example to set the author for.
        author_name : str
            The name of the author.
        author_email : str
            The email address of the author.

        Returns
        -------
        Example
            The updated Example object with the new author information.
        """
        return self.update_example(example_id, author_name=author_name, author_email=author_email)

    def clone(self, newpath: str | Path) -> 'ExampleManager':
        """
        Spawn a new ExampleManager instance rooted at a directory different from the current ExampleManager's path (newpath).  The new ExampleManager instance will contain copies of all examples in the current ExampleManager's examples list.

        Parameters
        ----------
        subpath : str or Path
            The subdirectory path relative to the current ExampleManager's path.

        Returns
        -------
        ExampleManager
            A new ExampleManager instance rooted at the specified subdirectory.
        """
        new_manager = ExampleManager(examples_path=newpath, docs_source_path=None)
        for example in self.examples:
            new_manager.copy_from(self, example)
        return new_manager

    def copy_from(self, other: 'ExampleManager', example: Example):
        """
        Copy an example from another ExampleManager instance to this instance.

        Parameters
        ----------
        other : ExampleManager
            The other ExampleManager instance to copy the example from.
        example : Example
            The Example instance to copy.
        overwrite : bool, optional
            If True, overwrite existing files in the example folder. Default is False.
        """
        current_ids = [ex.example_id for ex in self.examples.data]
        if example.example_id in current_ids:
            raise ValueError(f'Example ID {example.example_id} already exists in the destination ExampleManager')
        source_example_folder = other.examplefolderpath(example)
        dest_example_folder = self.examplefolderpath(example)
        shutil.copytree(source_example_folder, dest_example_folder, dirs_exist_ok=True)
        # create the outputs subfolder regardless of whether there are any outputs to copy
        dest_outputs_path = self.outputspath(example)
        dest_outputs_path.mkdir(parents=True, exist_ok=True)
        if example.outputs:
            source_outputs_path = other.outputspath(example)
            for output_file in example.outputs:
                source_file_path = source_outputs_path / output_file
                if source_file_path.is_file():
                    shutil.copy(source_file_path, dest_outputs_path)
                else:
                    logger.warning(f'Output file {output_file} does not exist in source ExampleManager; skipping copy')
        self.examples.append(example)