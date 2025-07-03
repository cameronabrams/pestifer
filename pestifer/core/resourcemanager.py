# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
Defines the ``ResourceManager`` class for managing access to :mod:`pestifer.resources`.
"""

import glob
import os
import shutil
import yaml
from .. import resources
from ..charmmff.charmmffcontent import CHARMMFFContent
from .labels import Labels

import logging
logger=logging.getLogger(__name__)

class ExampleManager:
    """
    A class for managing example input files in pestifer.
    This class provides methods to check out example YAML files, report the list of examples, and manage example resources.
    
    Parameters
    ----------
    example_path : str
        The path to the directory containing example input files. This directory should contain an ``info.yaml`` file
        that describes the examples available in that directory.
    """

    def __init__(self,example_path):
        if not os.path.isdir(example_path):
            raise FileNotFoundError(f'Example path {example_path} does not exist or is not a directory')
        self.path=example_path
        info_file=os.path.join(example_path,'info.yaml')
        if not os.path.isfile(info_file):
            raise FileNotFoundError(f'Example path {example_path} does not contain info.yaml')
        with open(info_file,'r') as f:
            self.info=yaml.safe_load(f)
        if 'input-files' not in self.info:
            raise KeyError(f'info.yaml in {example_path} does not contain input-files key')
        self.examples_list=self.info['input-files']
    
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
        example_yaml=self.examples_list[index-1]['name']
        example_yaml_path=os.path.join(self.path,example_yaml)
        if not os.path.isfile(example_yaml_path):
            raise FileNotFoundError(f'Example YAML file {example_yaml_path} does not exist')
        shutil.copy(example_yaml_path,os.getcwd())
        logger.info(f'Checked out example {index} from {self.path} to current working directory')
        return example_yaml # return the name of the copied file
    
    def report_examples_list(self,header=False,formatter=r'{:>7s}    {:<30s}    {}'):
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
            report_lines = [formatter.format('Index','Name','Description')+'\n']
        else:
            report_lines = []
        for i, example in enumerate(self.examples_list):
            name = example.get('name', 'Unnamed Example').split('.')[0]  # remove file extension
            description = example.get('description', 'No description provided')
            index_str= f'{(i+1):>2d}'
            report_lines.append(formatter.format(index_str, name, description))
        return ''.join(report_lines)

    def append_example(self,name:str,description:str,pdbID:str):
        """
        Append a new example to the examples list.

        Parameters
        ----------
        name : str
            The name of the example file to append. This should be a valid file path.
        description : str
            A description of the example.
        pdbID : str
            The PDB ID associated with the example. 

        Raises
        ------
        ValueError
            If the name, description, or pdbID is not provided. 
        FileNotFoundError
            If the example file does not exist at the specified path.
        """
        if not name or not description or not pdbID:
            raise ValueError('Name, description, and pdbID must be provided')
        new_example = {
            'name': name,
            'description': description,
            'pdbID': pdbID
        }
        if not os.path.isfile(name):
            raise FileNotFoundError(f'Example file {name} does not exist')
        new_example['name'] = os.path.basename(name)
        shutil.copy(name, self.path)
        self.examples_list.append(new_example)
        self.info['input-files'] = self.examples_list
        info_file=os.path.join(self.path,'info.yaml')
        with open(info_file,'w') as f:
            yaml.safe_dump(self.info,f)
        logger.info(f'Added new example: {name}')

    def delete_example(self,index:int):
        """
        Delete an example from the examples list by its index.
        
        Parameters
        ----------
        index : int
            The index of the example to delete (1-based).

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
        example=self.examples_list.pop(real_index)
        info_file=os.path.join(self.path,'info.yaml')
        with open(info_file,'w') as f:
            yaml.safe_dump(self.info,f)
        example_file=os.path.join(self.path,example['name'])
        if os.path.isfile(example_file):
            os.remove(example_file)
        logger.info(f'Deleted example {index}: {example["name"]}')

    def insert_example(self,index:int,name:str,description:str,pdbID:str):
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

        Raises
        ------
        IndexError
            If the index is out of range for the examples list.
        ValueError
            If the name, description, or pdbID is not provided.
        FileNotFoundError
            If the example file does not exist at the specified path.
        """
        real_index = index - 1  # convert to zero-based index
        if real_index < 0 or real_index > len(self.examples_list):
            raise IndexError(f'Index {index} is out of range for examples list of length {len(self.examples_list)}')
        if not name or not description or not pdbID:
            raise ValueError('Name, description, and pdbID must be provided')
        new_example = {
            'name': name,
            'description': description,
            'pdbID': pdbID
        }
        if not os.path.isfile(name):
            raise FileNotFoundError(f'Example file {name} does not exist')
        new_example['name'] = os.path.basename(name)
        shutil.copy(name, self.path)
        self.examples_list.insert(real_index, new_example)
        self.info['input-files'] = self.examples_list
        info_file=os.path.join(self.path,'info.yaml')
        with open(info_file,'w') as f:
            yaml.safe_dump(self.info,f)
        logger.info(f'Inserted new example at index {index}: {name}')   

class ResourceManager:
    """
    A class for managing pestifer's built-in resources.
    This class provides access to various resources such as CHARMM force fields, example input files, and TcL scripts.
    It also manages the ``ycleptic`` configuration file used for user-specific configurations.
    The resources are organized into directories, and the class provides methods to access these resources.
    """
    
    base_resources=['charmmff','examples','tcl','ycleptic']
    """
    A list of base resources that are managed by the ResourceManager.
    These resources include CHARMM force fields, example input files, TcL scripts, and the ycleptic configuration file.
    """
    
    ignored_resources=['__pycache__','_archive','bash']
    """
    A list of resource directories that are ignored by the ResourceManager.
    These directories are not considered part of the base resources and are excluded from resource management.
    """

    def __init__(self):
        self.resources_path=os.path.dirname(resources.__file__)
        self.resource_dirs=[x for x in glob.glob(os.path.join(self.resources_path,'*')) if os.path.isdir(x) 
                            and not os.path.basename(x) in ResourceManager.ignored_resources]
        assert all([x in [os.path.basename(_) for _ in self.resource_dirs] for x in ResourceManager.base_resources]),f'some resources seem to be missing'
        self.ycleptic_configdir=os.path.join(self.resources_path,'ycleptic')
        ycleptic_files=glob.glob(os.path.join(self.ycleptic_configdir,'*'))
        assert len(ycleptic_files)==1,f'Too many config files in {self.ycleptic_configdir}: {ycleptic_files}'
        self.ycleptic_config=ycleptic_files[0]
        self.resource_path={}
        for r in ResourceManager.base_resources:
            self.resource_path[r]=os.path.join(self.resources_path,r)
            if not os.path.isdir(self.resource_path[r]):
                raise FileNotFoundError(f'Resource {r} not found at {self.resource_path[r]} -- your installation is likely incomplete')
        self.charmmff_content=CHARMMFFContent(self.resource_path['charmmff'])
        self.pdbrepository=self.charmmff_content.pdbrepository
        self.example_manager=ExampleManager(self.resource_path['examples'])
        self.labels=Labels

    def __str__(self):
        cp=os.path.commonpath(list(self.resource_path.values()))
        retstr=f'Pestifer resources are found under\n    {cp}\n'
        for r,p in self.resource_path.items():
            retstr+=f'        {p.replace(cp+os.sep,"")+os.sep}\n'
        return retstr

    def show(self,out_stream=print,components={},fullnames=False,missing_fullnames={}):
        """
        Show a summary of the specified components.
        
        Parameters
        ----------
        out_stream : callable, optional
            A callable that takes a string and outputs it. Default is `print`.
        components : dict, optional
            A dictionary specifying the components to show. The keys are the component names
            and the values are the specifications for each component.
        fullnames : bool, optional
            If True, show full names for the components. Default is False.
        missing_fullnames : dict, optional
            A dictionary of missing full names for the components. Default is an empty dictionary.
        """
        for c,spec in components.items():
            if not c in self.base_resources:
                logger.warning(f'{c} is not a base resource; expected one of {", ".join(self.base_resources)}')
            path=self.get_resource_path(c)
            if c=='examples':
                out_stream(f'\nExamples:\n\n{self.example_manager.report_examples_list(header=True)}')
            elif c=='charmmff':
                if 'toppar' in spec:
                    out_stream(f'{self.charmmff_content.tarfilename}')
                if 'pdb' in spec:
                    self.pdbrepository.show(out_stream,fullnames=fullnames,missing_fullnames=missing_fullnames)
                if 'custom' in spec:
                    path=self.get_charmmff_customdir()
                    with open(os.path.join(path,'00PESTIFER-README.txt'),'r') as f:
                        msg=f.read()
                    out_stream(msg)
            elif c=='tcl':
                path=self.get_tcldir()
                with open(os.path.join(path,'00PESTIFER-README.txt'),'r') as f:
                    msg=f.read()
                out_stream(msg)

    def get_ycleptic_config(self):
        """
        Get the path to the ``ycleptic`` configuration file.
        This file is used for defining the formats of user-specific configurations and is located in the resources directory.
        
        Returns
        -------
        str
            The path to the ``ycleptic`` configuration file.
        """
        return self.ycleptic_config
    
    def get_resource_path(self,r):
        """
        Get the path to a specific resource.

        Parameters
        ----------
        r : str
            The name of the resource to get the path for. This should be one of the base resources defined in ``ResourceManager.base_resources``.
        """
        return self.resource_path.get(r,None)
            
    def get_charmmff_customdir(self):
        """
        Get the path to the custom CHARMM force field directory.
        This directory is used for storing custom CHARMM force field files that are not part of the standard distribution.

        Returns
        -------
        str
            The path to the custom CHARMM force field directory.
        """
        return os.path.join(self.resource_path['charmmff'],'custom')

    def get_tcldir(self):
        """
        Get the path to the TcL scripts and packages directory.
        
        Returns
        -------
        str
            The path to the TcL scripts and packages directory.
        """
        return self.resource_path['tcl']
    
    def update_atomselect_macros(self):
        """
        Update the atomselect macros in the macros.tcl file based on the content of labels.py.
        This method updates atomselect macros for residues that aren't already part of
        atomselect macros in a base VMD installation.
        """
        p=os.path.join(self.resource_path['tcl'],'macros.tcl')
        with open(p,'w') as f:
            f.write('## Updates atomselect macros to account for residues not in the macros in a base VMD installation.\n')
            f.write('## This file is automatically generated by ``pestifer.core.resourcemanager.update_atomselect_macros().``\n')
            f.write('## Do not edit this file directly; changes will be overwritten.\n')
            f.write('## Instead, edit `pestifer.core.labels.py` and run ``pestifer modify-package --update-atomselect-macros.``\n')
            f.write('## The Tcl proc ``update_atomselect_macro`` is defined in the PestiferUtil Tcl package.\n\n')
            self.labels.update_atomselect_macros(f)

    def get_tcl_pkgdir(self):
        """
        Get the path to the TcL package directory.
        This directory contains TcL packages that are used by pestifer for various tasks, such as loading additional functionality or libraries.

        Returns
        -------
        str
            The path to the TcL package directory.
        """
        return os.path.join(self.resource_path['tcl'],'pkg')
    
    def get_tcl_scriptsdir(self):
        """
        Get the path to the TcL scripts directory.
        This directory contains TcL scripts that are used by pestifer for various tasks, such as setting up simulations or processing data.

        Returns
        -------
        str
            The path to the TcL scripts directory.
        """
        return os.path.join(self.resource_path['tcl'],'scripts')
    
    def update_pdbrepository(self,user_pdbrepository=[]):
        """
        Update the PDB repository with user-defined paths.

        Parameters
        ----------
        user_pdbrepository : list, optional
            A list of paths to user-defined PDB collections. These paths will be added to the PDB repository.
            If not provided, an empty list is used, meaning no user-defined paths will be added.
        """
        for path in user_pdbrepository:
            logger.info(f'Adding user PDB collection: {path}')
            self.charmmff_content.pdbrepository.add_path(path)

    def update_charmmff(self,tarball='',user_custom_directory=None):
        """
        Update the CHARMM force field content with a tarball and/or a user-defined custom directory

        Parameters
        ----------
        tarball : str, optional
            The path to a tarball containing CHARMM force field files. If provided, the tarball will be loaded into the CHARMM force field content.
            If not provided, no tarball will be loaded beyond the default CHARMM force field files.
        user_custom_directory : str, optional
            The path to a user-defined custom directory containing CHARMM force field files. If provided, this directory will be added to the CHARMM force field content.
            If not provided, no custom directory will be added. 
        """
        if tarball:
            self.charmmff_content.load_charmmff(tarball)
        if user_custom_directory:
            self.charmmff_content.add_custom_directory(user_custom_directory)