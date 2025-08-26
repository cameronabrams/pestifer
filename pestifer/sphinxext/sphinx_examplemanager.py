# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
This module provides the `SphinxExampleManager` class, which allows for the management of examples in the documentation.
It includes methods to add, update, and delete examples.  It also provides functionality to modify the TOC tree for the main 
examples RST file.  This is typically the file `docs/source/examples.rst` in the package.

Files pertaining to documentation of specific examples are stored in docs/source/examples/ex##/ where ## is the two-digit 
unique ID for the example.  The RST file in each folder should have the name <shortname>.rst, where <shortname> is the
name of the example (the 'shortname' attribute of the associated :class:~pestifer.core.example.Example instance).

"""
import logging
import os
import shutil

from pidibble.pdbparse import PDBParser

from .toctree_util import modify_toctree, get_name_from_toctree
from ..core.example import Example
from ..util.formatvalidator import FormatValidator
from ..util.stringthings import example_footer

logger = logging.getLogger(__name__)

class SphinxExampleManager:
    """
    Sphinx extension to manage examples in the documentation.

    This extension allows for the inclusion of example input files in the documentation.
    It provides a way to link examples to their documentation and manage example files.

    Attributes
    ----------
    docs_source_path : str
        The path to the source directory (typically PACKAGE/docs/source) containing documentation files parseable by Sphinx.
    """

    def __init__(self, docs_source_path=None, examples_folder_name='examples', examples_rst_name='examples.rst'):
        if docs_source_path is not None and not os.path.isdir(docs_source_path):
            raise FileNotFoundError(f'Docs source path {docs_source_path} does not exist or is not a directory')
        if docs_source_path is None:
            return None
        self.docs_source_path = docs_source_path
        self.example_folder_name_format_validator = FormatValidator(Example.folder_name_format, greedy=True, flexible_ws=True)
        self.examples_basename = examples_rst_name.replace('.rst', '')
        self.examples_folder_name = examples_folder_name
        # determine if this is an empty docs source path and if so, build it
        self.examples_folder_path = os.path.join(docs_source_path, examples_folder_name)
        if not os.path.isdir(self.examples_folder_path):
            logger.debug(f'Examples folder {self.examples_folder_path} under {docs_source_path} does not exist or is not a directory')
            os.makedirs(self.examples_folder_path, exist_ok=True)
        self.examples_rst = os.path.join(docs_source_path, examples_rst_name)
        if not os.path.isfile(self.examples_rst):
            logger.debug(f'Main examples RST file {self.examples_rst} under {docs_source_path} does not exist. Creating a new one.')
            with open(self.examples_rst, 'w') as f:
                f.write(f'.. _{self.examples_basename}:\n\n')
                f.write(f'{self.examples_basename.title()}\n')
                f.write(f'{"=" * len(self.examples_basename)}\n\n')
                f.write('.. toctree::\n   :maxdepth: 1\n\n')
        p = PDBParser()
        self.rcsb_url = p.pdb_format_dict['BASE_URL']
        self.alphafold_url = p.pdb_format_dict['ALPHAFOLD_API_URL']
        self.alphafold_api_key = p.pdb_format_dict['ALPHAFOLD_API_KEY']

    def structure_url(self, example: Example) -> str:
        if example.db_id and example.db_id[0] in 'OP' and len(example.db_id) > 4:
            # this is likely an alphafold ID
            url= f"{self.alphafold_url}/{example.db_id}&key={self.alphafold_api_key}"
        else:
            url= f"{self.rcsb_url}/{example.db_id}"
        return url

    def new_rst(self, example:Example):
        """
        Create a new RST file for the example.

        Parameters
        ----------
        example : Example
            The example to create an RST file for.

        Returns
        -------
        str
            The contents of then new RST file
        """
        title = f"Example {example.example_id}: {example.title}"
        example_folder_name = f'{example.example_id:02d}'
        yaml_script_name = example.shortname if example.shortname.endswith('.yaml') else f"{example.shortname}.yaml"
        url = self.structure_url(example)
        dbname = 'PDB' if 'rcsb' in url else 'Alphafold'
        basestring = \
               f".. _example {example.shortname}:\n\n{title}\n" \
               f"{'-' * len(title)}\n\n" \
               f"`{dbname} ID {example.db_id} <{url}>`_ is...\n\n" \
               f"This example demonstrates that ...\n\n" \
               f".. literalinclude:: ../../../../pestifer/resources/examples/{example_folder_name}/{example.inputs_subdir}/{yaml_script_name}\n" \
               f"    :language: yaml\n\n"
        if example.author_email and example.author_name:
            basestring += example_footer(example.author_name, example.author_email)
        return basestring
    
    def append_example(self, example: Example):
        """
        Add an example to the documentation
        
        Parameters
        ----------

        example : Example
            The example to add.
        """
        example_folder_name = Example.folder_name_format.format(example_id=example.example_id)
        assert not os.path.exists(os.path.join(self.examples_folder_path, example_folder_name))
        os.makedirs(os.path.join(self.examples_folder_path, example_folder_name), exist_ok=True)
        destination_file = os.path.join(self.examples_folder_path, example_folder_name, f'{example.shortname}.rst')
        with open(destination_file, 'w') as f:
            f.write(self.new_rst(example))
        if self.examples_rst:
            modify_toctree(self.examples_rst, action='append', new_entry=f'{example_folder_name}/{example.shortname}', common_prefix=self.examples_folder_name)

    def update_example(self, example_id: int, example: Example):
        """
        Update an existing example's documentation with content in example
        
        Parameters
        ----------
        example_id : int
            The unique ID of the example to update.
        example : Example
            The example whose content will be used to update the existing example
        """
        current_example_folder_name = Example.folder_name_format.format(example_id=example_id)
        current_example_root = os.path.join(self.examples_folder_path, current_example_folder_name)
        if not os.path.isdir(current_example_root):
            raise FileNotFoundError(f'Current example folder {current_example_folder_name} does not exist; cannot update.')
        current_name = get_name_from_toctree(self.examples_rst, current_example_folder_name)
        current_rst_file = os.path.join(current_example_root, f'{current_name}.rst')
        if not os.path.exists(current_rst_file):
            raise FileNotFoundError(f'Current example file {current_rst_file} does not exist. Cannot update non-existing example files.')
        updated_name = example.shortname
        if current_name != updated_name:
            # rename the current example rst file to include the index
            updated_rst_file = os.path.join(current_example_root, f'{updated_name}.rst')
            os.rename(current_rst_file, updated_rst_file)
            logger.debug(f'Renamed {current_rst_file} to {updated_rst_file}')
            if self.examples_rst:
                modify_toctree(self.examples_rst, action='update', target=f'{current_example_folder_name}/{current_name}', new_entry=f'{current_example_folder_name}/{example.shortname}')
            current_rst_file = updated_rst_file
        with open(current_rst_file, 'r') as f:
            rst_lines=f.readlines()
        with open(current_rst_file, 'w') as f:
            for line in rst_lines:
                if updated_name != current_name and current_name in line:
                    # replace the current example name with the new example name
                    line = line.replace(current_name, updated_name)
                elif '<div' in line or '</div' in line:
                    continue
                elif 'Author' in line:
                    # replace the author information with the new example's author information
                    if example.author_name and example.author_email:
                        line = example_footer(example.author_name, example.author_email)
                    else:
                        continue
                elif 'PDB ID' in line or 'Alphafold ID' in line:
                    # replace the PDB ID with the new example's PDB ID
                    if example.db_id:
                        url = self.structure_url(example)
                        if 'rcsb' in url:
                            line = f"`PDB ID {example.db_id} <{url}>`_\n"
                        elif 'alphafold' in url:
                            line = f"`Alphafold ID {example.db_id} <{url}>`_\n"
                f.write(line)

    def delete_example(self, example: Example):
        """
        Delete an example from the examples RST file.  This method first searches the TOC tree for the example's name, then deletes the corresponding file in the examples source path (typically docs/source/examples/name.rst). Since the examples are ordered by their position in the TOC tree, this method will also shift the indices of subsequent examples down by one, by editing the titles of their respective RST files.
        
        Parameters
        ----------
        example : Example
            The example to delete.
        """
        current_example_folder_name = Example.folder_name_format.format(example_id=example.example_id)
        current_example_root = os.path.join(self.examples_folder_path, current_example_folder_name)
        if not os.path.isdir(current_example_root):
            raise FileNotFoundError(f'Current example folder {current_example_folder_name} does not exist; cannot delete.')
        shutil.rmtree(current_example_root)
        if self.examples_rst:
            modify_toctree(self.examples_rst, action='delete', target=f'{current_example_folder_name}/{example.shortname}')

    # def rename_example(self, index:int, new_name:str):
    #     """
    #     Rename an example in the examples RST file.
        
    #     Parameters
    #     ----------
    #     index : int
    #         The 1-based index of the example to rename.
    #     new_name : str
    #         The new name for the example.
    #     """
    #     if self.examples_folder_path:
    #         current_name = get_name_from_toctree(self.examples_rst, index)
    #         current_rst_file = os.path.join(self.examples_folder_path, f'{current_name}.rst')
    #         if not os.path.exists(current_rst_file):
    #             raise FileNotFoundError(f'Current example file {current_rst_file} does not exist. Cannot rename non-existing example files.')
    #         new_rst_file = os.path.join(self.examples_folder_path, f'{new_name}.rst')
    #         os.rename(current_rst_file, new_rst_file)
    #         logger.debug(f'Renamed {current_rst_file} to {new_rst_file}')
    #         # update the explicit yaml include in the example RST file
    #         with open(new_rst_file, 'r') as f:
    #             lines = f.readlines()
    #         with open(new_rst_file, 'w') as f:
    #             for line in lines:
    #                 if line.startswith('.. literalinclude::') and 'examples/' in line:
    #                     f.write(line.replace(f'examples/{current_name}.yaml', f'examples/{new_name}.yaml'))
    #                 else:
    #                     f.write(line)
    #     if self.examples_rst:
    #         modify_toctree(self.examples_rst, action='rename', index=index, new_entry=new_name)

    # def set_author(self, index:int, author_name:str, author_email:str):
    #     """
    #     Set the author information for an example in the examples RST file.
        
    #     Parameters
    #     ----------
    #     index : int
    #         The 1-based index of the example to set the author for.
    #     author_name : str
    #         The name of the author.
    #     author_email : str
    #         The email of the author.
    #     """
    #     if self.examples_folder_path:
    #         current_name = get_name_from_toctree(self.examples_rst, index)
    #         current_rst_file = os.path.join(self.examples_folder_path, f'{current_name}.rst')
    #         if not os.path.exists(current_rst_file):
    #             raise FileNotFoundError(f'Current example file {current_rst_file} does not exist. Cannot set author for non-existing example files.')
    #         with open(current_rst_file, 'a') as f:
    #             f.write(f'\n\n')
    #             f.write(example_footer(author_name, author_email))