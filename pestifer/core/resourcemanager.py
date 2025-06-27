# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" Defines the ResourceManager class for managing access to pestifer's built-in resources 
"""
import glob
import os
import shutil
import yaml
from .. import resources
from ..charmmff.charmmffcontent import CHARMMFFContent

import logging
logger=logging.getLogger(__name__)

class ExampleManager:
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
        # copies requested example YAML input file to CWD
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
        # return a string representation of the examples list
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
        # append a new example to the examples list
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
        # delete an example from the examples list
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
        # insert a new example at the specified index
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
    base_resources=['charmmff','examples','tcl','ycleptic']
    ignored_resources=['__pycache__','_archive','bash']
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

    def __str__(self):
        cp=os.path.commonpath(list(self.resource_path.values()))
        retstr=f'Pestifer resources are found under\n    {cp}\n'
        for r,p in self.resource_path.items():
            retstr+=f'        {p.replace(cp+os.sep,"")+os.sep}\n'
        return retstr

    def show(self,out_stream=print,components={},fullnames=False,missing_fullnames={}):
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
        return self.ycleptic_config
    
    def get_resource_path(self,r):
        return self.resource_path.get(r,None)
            
    def get_charmmff_customdir(self):
        return os.path.join(self.resource_path['charmmff'],'custom')

    def get_tcldir(self):
        return self.resource_path['tcl']
    
    def get_tcl_pkgdir(self):
        return os.path.join(self.resource_path['tcl'],'pkg')
    
    def get_tcl_scriptsdir(self):
        return os.path.join(self.resource_path['tcl'],'scripts')
    
    def update_pdbrepository(self,user_pdbrepository=[]):
        for path in user_pdbrepository:
            logger.info(f'Adding user PDB collection: {path}')
            self.charmmff_content.pdbrepository.add_path(path)

    def update_charmmff(self,tarball='',user_custom_directory=None):
        if tarball:
            self.charmmff_content.load_charmmff(tarball)
        if user_custom_directory:
            self.charmmff_content.add_custom_directory(user_custom_directory)