# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
Manages examples for ``pestifer run`` within the documentation.
This module provides the `SphinxExampleManager` class, which allows for the management of examples in the documentation.
It includes methods to add, insert, update, and delete examples, as well as to shift title indices in the examples RST file.
It also provides functionality to modify the TOC tree for the examples RST file.  The only way :class:`SphinxExampleManager` understands the **order** of the examples is the order they are referenced in the TOC tree of the examples RST file.
This is typically the file `docs/source/examples.rst` in the package.
"""
import os
import logging
logger = logging.getLogger(__name__)
from .toctree_util import modify_toctree,get_name_from_toctree
from ..core.examplemanager import Example
from ..core.stringthings import example_footer

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
        self.examples_folder_name = examples_folder_name
        self.examples_rst = os.path.join(docs_source_path, examples_rst_name) if docs_source_path else None
        assert os.path.isfile(self.examples_rst), f'Examples RST file {self.examples_rst} does not exist or is not a file'
        self.examples_source_path = os.path.join(docs_source_path, examples_folder_name) if docs_source_path else None
        assert os.path.isdir(self.examples_source_path), f'Examples source path {self.examples_source_path} does not exist or is not a directory'
    
    # def add_example(self, example:Example):
    #     """
    #     Add a new example to the documentation.
    #     This method creates a new example file in the examples source path and updates the examples RST file.

    #     Parameters
    #     ----------
    #     example : Example
    #         The example to add. The example must have a valid index (greater than 0).  It should have been set by the :class:`ExampleManager` prior to calling this method to reflect its position at the end of the current examples list.
    #     """
    #     if example.index==0:
    #         raise ValueError('Example index must be greater than 0. Use insert_example to insert an example at a specific index.')
    #     if self.examples_source_path:
    #         destination_file= os.path.join(self.examples_source_path, f'{example.name}.rst')
    #         if not os.path.exists(destination_file):
    #            with open(destination_file, 'w') as f:
    #                f.write(example.new_rst())
    #         else:
    #             raise FileExistsError(f'Example file {destination_file} already exists. Cannot overwrite existing example files.')
    #     if self.examples_rst:
    #         modify_toctree(self.examples_rst, action='add', new_entry=example.name)

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
            The path to the newly created RST file.
        """
        title= f"Example {example.index}: {example.description}"
        yaml_file_name=example.name if example.name.endswith('.yaml') else f"{example.name}.yaml"
        if example.pdbID.startswith('P') and len(example.pdbID) > 4:
            # this is likely an alphafold ID
            url= f"https://alphafold.com/api/prediction/{example.pdbID}"
        else:
            url= f"https://www.rcsb.org/structure/{example.pdbID}"
        basestring= f".. _example {example.name}:\n\n{title}\n" \
               f"{'-' * len(title)}\n\n" \
               f"`PDB ID {example.pdbID} <{url}>`_ is...\n\n" \
               f"This example demonstrates that ...\n\n" \
               f".. literalinclude:: ../../../pestifer/resources/examples/{yaml_file_name}\n" \
               f"    :language: yaml\n\n"
        if example.author_email and example.author_name:
            basestring += example_footer(example.author_name, example.author_email)
        return basestring
    
    def title_index_shift(self, start_index:int, shift:int):
        """
        Shift the index of titles in the examples RST file.
        
        Parameters
        ----------
        start_index : int
            The starting index (1-based) from which to shift titles.
        shift : int
            The amount by which to shift the indices.
        """
        existing_rst_files = [f for f in os.listdir(self.examples_source_path) if f.endswith('.rst')]
        for existing_file in existing_rst_files:
            with open(os.path.join(self.examples_source_path, existing_file), 'r') as f:
                lines = f.readlines()
            with open(os.path.join(self.examples_source_path, existing_file), 'w') as f:
                for line in lines:
                    if line.startswith('Example'):
                        try:
                            index = int(line.split(':')[0].split(' ')[1])
                            if index >= start_index:
                                new_index = index + shift
                                new_title = line.replace(f'Example {index}', f'Example {new_index}')
                                f.write(new_title)
                            else:
                                f.write(line)
                        except ValueError:
                            f.write(line)
                    else:
                        f.write(line)

    def insert_example(self, index:int, example:Example):
        """
        Insert an example at a specific index in the examples RST file.
        
        Parameters
        ----------
        index : int
            The 1-based index label of the inserted example.
        example : Example
            The example to insert.
        """
        destination_file = os.path.join(self.examples_source_path, f'{example.name}.rst')
        if not os.path.exists(destination_file):
            self.title_index_shift(index, 1)  # Shift existing titles to make space for the new example
            with open(destination_file, 'w') as f:
                f.write(self.new_rst(example))
        else:
            raise FileExistsError(f'Example file {destination_file} already exists. Cannot overwrite existing example files.')
        if self.examples_rst:
            modify_toctree(self.examples_rst, action='insert', new_entry=example.name, index=index, common_prefix=self.examples_folder_name)

    def update_example(self, index: int, example:Example):
        """
        Update an existing example's documentation with content in example
        
        Parameters
        ----------
        index : int
            The 1-based index of the example to update.
        example : Example
            The example whose content will be used to update the existing example
        """
        current_name = get_name_from_toctree(self.examples_rst, index)
        current_rst_file = os.path.join(self.examples_source_path, f'{current_name}.rst')
        if not os.path.exists(current_rst_file):
            raise FileNotFoundError(f'Current example file {current_rst_file} does not exist. Cannot update non-existing example files.')
        updated_name = example.name
        if current_name != updated_name:
            # rename the current example rst file to include the index
            updated_rst_file = os.path.join(self.examples_source_path, f'{updated_name}.rst')
            os.rename(current_rst_file, updated_rst_file)
            logger.debug(f'Renamed {current_rst_file} to {updated_rst_file}')
            if self.examples_rst:
                modify_toctree(self.examples_rst, action='update', target=current_name, new_entry=example.name)
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
                        line = ''
                elif 'PDB ID' in line:
                    # replace the PDB ID with the new example's PDB ID
                    if example.pdbID:
                        url = f"https://alphafold.com/api/prediction/{example.pdbID}" if example.pdbID.startswith('P') and len(example.pdbID) > 4 else f"https://www.rcsb.org/structure/{example.pdbID}"
                        line = f"`PDB ID {example.pdbID} <{url}>`_\n"
                elif 'Alphafold ID' in line:
                    # replace the Alphafold ID with the new example's PDB ID
                    if example.pdbID:
                        url = f"https://alphafold.com/api/prediction/{example.pdbID}" if example.pdbID.startswith('P') and len(example.pdbID) > 4 else f"https://www.rcsb.org/structure/{example.pdbID}"
                        line = f"`Alphafold ID {example.pdbID} <{url}>`_\n"
                f.write(line)

    def delete_example(self, example:Example):
        """
        Delete an example from the examples RST file.  This method first searches the TOC tree for the example's name, then deletes the corresponding file in the examples source path (typically docs/source/examples/name.rst). Since the examples are ordered by their position in the TOC tree, this method will also shift the indices of subsequent examples down by one, by editing the titles of their respective RST files.
        
        Parameters
        ----------
        example : Example
            The example to delete.
        """
        destination_file = os.path.join(self.examples_source_path, f'{example.name}.rst')
        if os.path.exists(destination_file):
            os.remove(destination_file)
            self.title_index_shift(example.index, -1)
        else:
            raise FileNotFoundError(f'Example file {destination_file} does not exist. Cannot delete non-existing example files.')
        if self.examples_rst:
            modify_toctree(self.examples_rst, action='delete', target=example.name)

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
    #     if self.examples_source_path:
    #         current_name = get_name_from_toctree(self.examples_rst, index)
    #         current_rst_file = os.path.join(self.examples_source_path, f'{current_name}.rst')
    #         if not os.path.exists(current_rst_file):
    #             raise FileNotFoundError(f'Current example file {current_rst_file} does not exist. Cannot rename non-existing example files.')
    #         new_rst_file = os.path.join(self.examples_source_path, f'{new_name}.rst')
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
    #     if self.examples_source_path:
    #         current_name = get_name_from_toctree(self.examples_rst, index)
    #         current_rst_file = os.path.join(self.examples_source_path, f'{current_name}.rst')
    #         if not os.path.exists(current_rst_file):
    #             raise FileNotFoundError(f'Current example file {current_rst_file} does not exist. Cannot set author for non-existing example files.')
    #         with open(current_rst_file, 'a') as f:
    #             f.write(f'\n\n')
    #             f.write(example_footer(author_name, author_email))