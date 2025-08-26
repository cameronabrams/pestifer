# Author: Cameron F. Abrams <cfa22@drexel.edu>

""" 
A module for managing example input files in pestifer.
This module provides the :class:`ExampleManager` class, which allows users to check out example YAML files,
report the list of examples, and manage example resources.

ExampleManager does not directly manage the example documentation, but it relies on the :class:`SphinxExampleManager` 
class to handle the documentation side of things.
"""

from glob import glob
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
        self.path = os.path.abspath(examples_path)
        self._read_dirtree()
        if docs_source_path:
            # create the SphinxExampleManager instance if docs_source_path is provided
            self.sphinx_example_manager = SphinxExampleManager(docs_source_path=os.path.abspath(docs_source_path))
        else:
            self.sphinx_example_manager = None

    def _read_dirtree(self):
        """
        Read the directory tree to populate the examples list.
        """
        savedir = os.getcwd()
        os.chdir(self.path)
        # do stuff
        self.examples = ExampleList([])
        exdirs = [d for d in os.listdir('.') if os.path.isdir(d) if self.root_folder_format_validator.fullmatch(d)]
        exdirs.sort()
        for exdir in exdirs:
            example_id = self.root_folder_format_validator.extract(exdir).get("example_id")
            inputdir = os.path.join(exdir, 'inputs')
            yaml_inputs = glob(os.path.join(inputdir, '*.yaml'))
            title = 'No title provided'
            db_id = None
            shortname = None
            if len(yaml_inputs) == 0:
                logger.warning(f'No YAML files found in {inputdir}')
            elif len(yaml_inputs) > 1:
                logger.warning(f'Multiple YAML files found in {inputdir}: {yaml_inputs}')
            else:
                main_input_script = yaml_inputs[0]
                shortname = os.path.basename(main_input_script).replace('.yaml', '')
                title, db_id = Example._get_implied_metadata(main_input_script)
            if shortname:
                self.examples.append(Example(
                    example_id=example_id,
                    title=title,
                    shortname=shortname,
                    db_id=db_id
                ))

        os.chdir(savedir)

    def scriptpath(self, example: Example) -> Path:
        return self.path / example.scriptpath

    def checkout_example(self, example_id: int) -> Example:
        """
        Copy example YAML file and associated companion files by example ID to the current working directory.

        Parameters
        ----------
        example_id : int
            The ID of the example to check out.

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
        example = self.examples.get_example_by_example_id(example_id)
        if not example:
            raise IndexError(f'Example with ID {example_id} not found')
        example_folder = f'ex{example_id:02d}'
        example_folder_path = os.path.join(self.path, example_folder)
        if not os.path.isdir(example_folder_path):
            raise FileNotFoundError(f'Example folder {example_folder} does not exist in {self.path}')
        example_inputs_path = os.path.join(example_folder_path, example.inputs_subdir)
        example_yaml = example.shortname + '.yaml'
        example_yaml_path = os.path.join(example_inputs_path, example_yaml)
        all_example_inputs = os.listdir(example_inputs_path)
        auxiliary_example_inputs = [f for f in all_example_inputs if f != example_yaml]
        if not os.path.isfile(example_yaml_path):
            raise FileNotFoundError(f'Example YAML file {example_yaml_path} does not exist')
        # copy all files in the example_yaml_path directory to the current working directory
        shutil.copy(example_yaml_path, os.getcwd())
        # copy all companion files to the current working directory
        for companion_file in auxiliary_example_inputs:
            companion_file_path = os.path.join(self.path, example_folder, companion_file)
            if os.path.isfile(companion_file_path):
                shutil.copy(companion_file_path, os.getcwd())
            else:
                logger.warning(f'Declared companion file "{companion_file}" does not exist in {os.path.join(self.path,example_folder)}')
        logger.info(f'Checked out example {example_id} from {self.path} to current working directory {os.getcwd()}')
        return example

    def new_example_yaml(self, db_id='ABCD', build_type='minimal'):
        """
        Generate a new example YAML file based on an existing example template.  The id can be a 4-letter PDB ID or an Alphafold/UNIPROT ID starting with "P".  The build_type can be either 'minimal' or 'full', which determines whether the generated YAML file contains only the psfgen task or all tasks including termination.

        Parameters
        ----------
        id : str, optional
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
        example_yaml = self.examples[0].shortname + '.yaml'
        example_yaml_path = os.path.join(self.path, self.examples[0].shortname, example_yaml)
        with open(example_yaml_path, 'r') as f:
            try:
                example_config = yaml.safe_load(f)
            except yaml.YAMLError as e:
                raise ValueError(f'Invalid YAML file {example_yaml_path}: {e}')
        if build_type == 'minimal':
            fetch_task = example_config['tasks'][0]
            example_config['tasks'] = [fetch_task, example_config['tasks'][1]]  # keep only the fetch task
        example_config['title'] = f'New template pestifer config for id {id} ({idtype})'
        if idtype == 'PDB' or idtype == 'Alphafold':
            example_config['tasks'][0]['fetch']['sourceID'] = db_id
        if build_type == 'full':
            example_config['tasks'][-1]['terminate']['basename'] = f'my_{db_id.lower()}'
            example_config['tasks'][-1]['terminate']['package']['basename'] = f'my_{db_id.lower()}'
        output_yaml = os.path.join(os.getcwd(), f'{db_id.lower()}.yaml')
        with open(output_yaml, 'w') as f:
            yaml.dump(example_config, f, default_flow_style=False)
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

    def checkin_example(self, example: Example):
        """
        Check in an Example instance by copying its YAML file and companion files to the appropriate example folder, which is created if it doesn't already exist.  This will overwrite the existing files.  If any file referenced by the example does not exist in the current working directory, a warning is logged but no action taken.

        Parameters
        ----------
        example : Example
            The Example instance to check in.
        """
        example_folder = os.path.join(self.path, example.rootfolderpath)
        example_inputs_subfolder = os.path.join(example_folder, example.inputs_subdir)
        example_outputs_subfolder = os.path.join(example_folder, example.outputs_subdir) 
        user_yaml_file_path = example.shortname + '.yaml'  # correct; the name attribute should never have an extension
        if not os.path.isdir(example_folder): # we are checking in a new example
            logger.debug(f'Creating example folder "{example_folder}"')
            if not os.path.isfile(user_yaml_file_path):
                raise FileNotFoundError(f'Example YAML file {user_yaml_file_path} does not exist in your current working directory {os.getcwd()}')
            os.makedirs(example_inputs_subfolder)
            os.makedirs(example_outputs_subfolder)
            shutil.copy(user_yaml_file_path, example_inputs_subfolder)
        else:
            if not os.path.isfile(user_yaml_file_path):
                logger.debug(f'Example YAML file {user_yaml_file_path} does not exist in your current working directory {os.getcwd()}')
                logger.debug(f'Pestifer will now check that this YAML file is already in the example folder {example_folder}')
                existing_yaml_file_path = os.path.join(example_inputs_subfolder, user_yaml_file_path)
                if not os.path.isfile(existing_yaml_file_path):
                    raise FileNotFoundError(f'Example YAML file {existing_yaml_file_path} does not exist in the example folder {example_folder}')
            else:
                logger.debug(f'Copying example YAML file "{user_yaml_file_path}" to example inputs subfolder "{example_inputs_subfolder}"')
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

        new_example = Example(example_id=example_id,shortname=os.path.splitext(scriptname)[0],title=title,db_id=db_id,author_name=author_name,author_email=author_email,auxiliary_inputs=auxiliary_inputs,outputs=outputs)
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
        current_example = self.examples.get_example_by_example_id(example_id)
        if not current_example:
            logger.warning(f'No example found with ID {example_id}')
            return None
        current_example_folder = os.path.join(self.path, current_example.folder_name_format.format(example_id=current_example.example_id))
        current_example_inputs_subfolder = os.path.join(current_example_folder, current_example.inputs_subdir)
        current_example_outputs_subfolder = os.path.join(current_example_folder, current_example.outputs_subdir)
        current_example_scriptpath = os.path.join(current_example_inputs_subfolder, f'{current_example.shortname}.yaml')
        assert os.path.isfile(current_example_scriptpath), f'Example script file {current_example_scriptpath} does not exist in the current example folder {current_example_folder}'

        if shortname:
            user_script = f'{shortname}.yaml' # in cwd
            if os.path.isfile(user_script):
                # user is providing a new script for this example
                shutil.copy(user_script, current_example_scriptpath)
                user_implied_title, user_implied_db_id = Example._get_implied_metadata(user_script)
                current_example.title = user_implied_title
                current_example.db_id = user_implied_db_id
            else:
                # user is not providing a new script, they just want to change the shortname
                new_example_scriptpath = os.path.join(current_example_inputs_subfolder, f'{shortname}.yaml')
                os.rename(current_example_scriptpath, new_example_scriptpath)
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
            for f in current_example.auxiliary_inputs:
                f_path = os.path.join(current_example_inputs_subfolder, f)
                if os.path.isfile(f_path):
                    os.remove(f_path)  # remove old auxiliary inputs
            current_example.auxiliary_inputs = auxiliary_inputs
            for f in auxiliary_inputs:
                if os.path.isfile(f):
                    shutil.copy(f, current_example_inputs_subfolder)
        if outputs:
            current_outputs = os.listdir(current_example_outputs_subfolder)
            for co in current_outputs:
                os.remove(os.path.join(current_example_outputs_subfolder, co))
            current_example.outputs = outputs
            for f in outputs:
                if os.path.isfile(f):
                    shutil.copy(f, current_example_outputs_subfolder)
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
