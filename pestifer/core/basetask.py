# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
This module defines the :class:`BaseTask` class, which serves as a base class for all task types in the pestifer framework.
The :class:`BaseTask` class provides a common interface and some basic functionality for tasks, including task initialization,
logging, state management, and file handling. It is intended to be subclassed by specific task types that implement their own
specific behavior and functionality.
The :class:`BaseTask` class includes methods for task execution, state variable management, file collection, and task name management.
It also provides a mechanism for inheriting state from prior tasks, saving the task state to a YAML file, and copying the state to files based on specified extensions.
The class also includes methods for converting coordinate files to PDB files and vice versa, as well as creating constraint PDB files based on task specifications.
The :class:`BaseTask` class is designed to be flexible and extensible, allowing for the creation of various task types that can perform different operations on molecular structures.

Available tasks that inherit from :class:`BaseTask` include all those in the :mod:`pestifer.tasks` subpackage.

"""
import logging
import os
import shutil
import yaml

from .baseobj import BaseObj
from .stringthings import FileCollector

logger=logging.getLogger(__name__)
logging.getLogger("matplotlib").setLevel(logging.WARNING)

class BaseTask(BaseObj):
    """ 
    A base class for Tasks.
    This class is intended to be subclassed by specific task types.
    It provides a common interface and some basic functionality for tasks.

    Parameters
    ----------
    config_specs : dict, optional
        A dictionary of configuration specifications. This can include task-specific settings.
    controller_specs : dict, optional
        A dictionary of specifications provided by the controller that created this task.
        It can include attributes like `prior`, `index`, `scripters`, `taskname`, and `controller_index`.

    """
    
    req_attr=BaseObj.req_attr+['specs','config','index','prior','scripters','taskname','controller_index']
    """ 
    List of required attributes as strings. These attributes must have values assigned at instance creation time, otherwise an AssertionError is raised.
    The required attributes are:
    
        * ``controller_index``: index of the controller that created this task
        * ``specs``: dictionary of specifications; under ``directives`` in yaml
        * ``scripters``: dictionary of ``FileWriters``
        * ``prior``: identifier of prior task in sequence
        * ``index``: unique integer index of task in run
        * ``config``: access to the run config
        * ``taskname``: caller-supplied name of task
    """

    yaml_header='generic_task'
    """
    The yaml header for this task type.  This is used to identify the task type in yaml files.
    """
    
    init_msg_options=['INITIATED','STARTED','BEGUN','SET IN MOTION','KICKED OFF','LIT','SPANKED ON THE BOTTOM']
    """
    A list of message options that are used to log the initiation of the task.
    These messages are used to log the start of the task in a consistent manner.
    The messages are used to log the start of the task in a consistent manner.
    They are used to log the start of the task in a consistent manner.
    """

    def __init__(self,config_specs={},controller_specs={}):
        """
        Constructor for the BaseTask class.
        """
        specs=config_specs.copy()
        prior=controller_specs.get('prior',None)
        index=controller_specs.get('index',0)
        scripters=controller_specs.get('scripters',{})
        taskname=controller_specs.get('taskname','generic_task')
        config=controller_specs.get('config',{})
        controller_index=controller_specs.get('controller_index',0)
        if prior:
            index=prior.index+1
        logger.debug(f'Creating task {taskname} with index {index}')
        input_dict = {
            'index':index,
            'scripters':scripters,
            'prior':prior,
            'specs':specs,
            'config':config,
            'taskname':taskname,
            'controller_index':controller_index
        }
        super().__init__(input_dict)
        self.subtaskcount=0
        self.statevars={}
        self.FC=FileCollector()
        self.result=0

    def override_taskname(self,taskname):
        """
        Override the task name with a new name.  Certain situations may require a Controller to override the task name.
        This method is used to change the task name after the task has been created.

        Parameters
        ----------
        taskname : str
            The new task name to set.
        """
        logger.debug(f'Overriding task name {self.taskname} to {taskname}')
        self.taskname=taskname

    def do(self):
        """
        This a stub method that should be overridden by subclasses.
        It is intended to be the main method that performs the task's operations.

        Returns
        -------
        int
            The result of the task execution. By default, it returns 0, indicating success.
        """
        return self.result

    def __str__(self):
        return f'{self.controller_index} - {self.index} - {self.taskname} [has prior {self.prior!=None}]'

    def log_message(self,message,**kwargs):
        """ 
        Logs a message with the task's index and name.

        Parameters
        ----------
        message : str
            The message to log.
        **kwargs : dict
            Additional keyword arguments that can be used to add extra information to the log message.
            These are typically used to provide additional context or details about the task's execution.
        """
        extra=''
        for k,v in kwargs.items():
            if v:
                extra+=f' ({k}: {v})'
        mtoks=[x.strip() for x in [x.upper() for x in message.split()]]
        if not any([x in self.init_msg_options for x in mtoks]):
            extra+=f' (result: {self.result})'
        logger.info(f'Controller {self.controller_index:02} Task {self.index:02} \'{self.taskname}\' {message} {extra}')

    def get_keepfiles(self):
        """ 
        Returns a list of files that should be kept after the task is done.
        """
        if hasattr(self,'keepfiles'):
            return self.keepfiles
        else:   
            return []
        
    def next_basename(self,*obj):
        """
        Generates a new basename for the task based on the controller index, task index, subtask count, and task name.
        If the ``specs`` dictionary contains a 'basename' key, it overrides the default basename.
        The basename is constructed in the format: ``controller_index-task_index-subtask_count_taskname-label``.
        If ``obj`` is provided, it appends a label to the basename based on the first element of ``obj``.

        Parameters
        ----------
        *obj : tuple
            Optional tuple of objects. If provided, the first element is used to create a label for the basename.
        """

        label=''
        if len(obj)==1 and len(obj[0])>0:
            label=f'-{obj[0]}'
        default_basename=f'{self.controller_index:02d}-{self.index:02d}-{self.subtaskcount:02d}_{self.taskname}{label}'
        overwrite_basename=self.specs.get('basename',None)
        if overwrite_basename:
            logger.debug(f'Overriding basename {default_basename} with {overwrite_basename}')
            self.basename=overwrite_basename
        else:
            self.basename=default_basename
        self.subtaskcount+=1
    
    def coor_to_pdb(self):
        """
        Converts the coordinate file to a PDB file using VMD.
        This method creates a new VMD script to perform the conversion.
        It uses the ``namdbin2pdb`` Tcl proc to convert the coordinate file to a PDB file.
        The resulting PDB file is named based on the task's basename.
        """
        vm=self.scripters['vmd']
        vm.newscript(f'{self.basename}-coor2pdb')
        psf=self.statevars['psf']
        coor=self.statevars['coor']
        vm.addline(f'namdbin2pdb {psf} {coor} {self.basename}.pdb')
        vm.writescript()
        vm.runscript()

    def pdb_to_coor(self):
        """
        Converts the PDB file to a coordinate file using VMD.
        This method creates a new VMD script to perform the conversion.
        It uses the ``pdb2namdbin`` Tcl proc to convert the PDB file to a coordinate file.
        The resulting coordinate file is named based on the task's basename.
        """
        vm=self.scripters['vmd']
        vm.newscript(f'{self.basename}-pdb2coor')
        pdb=self.statevars['pdb']
        vm.addline(f'pdb2namdbin {pdb} {self.basename}.coor')
        vm.writescript()
        vm.runscript()

    def make_constraint_pdb(self,specs,statekey='consref'):
        """
        Creates a PDB file with constraints based on the specifications provided.
        This method generates a VMD script that selects atoms based on the specifications and writes a PDB file with the specified constraints.
        The resulting PDB file is named based on the task's basename and the state key.
        This allows the task to keep track of the constraints PDB file in its state.
        The method uses the VMD script writer to create and run the script that generates the constraints PDB file.
        The script selects all atoms, sets their occupancy to 0.0, selects the constrained atoms based on the specifications, and sets their occupancy to the specified force constant.
        Finally, it writes the PDB file with the specified name.
        If the 'consref' key is not provided in the specifications, it defaults to ``f'{self.basename}-constraints.pdb'``.
        The resulting PDB file is then stored in the task's state variables under the specified statekey.
        If the specified atoms are not found in the PDB, an error will be raised.

        Parameters
        ----------
        specs : dict
            A dictionary containing the specifications for the constraints. It should include:

                - ``k``: The force constant for the constraints (default is taken from the config).
                - ``atoms``: A string defining the atoms to be constrained (default is ``all``).
                - ``consref``: The name of the output PDB file for constraints (default is ``f'{self.basename}-constraints.pdb'``).

        statekey : str, optional
            The key under which the resulting PDB file will be stored in the state variables. Default is ``consref``.
        """
        vm=self.scripters['vmd']
        pdb=self.statevars['pdb']
        force_constant=specs.get('k',self.config['user']['namd']['harmonic']['spring_constant'])
        constrained_atoms_def=specs.get('atoms','all')
        logger.debug(f'constraint spec: {specs["atoms"]}')
        c_pdb=specs.get('consref','')
        if not c_pdb:
            c_pdb=f'{self.basename}-constraints.pdb'
        vm.newscript(f'{self.basename}-make-constraint-pdb')
        vm.addline(f'mol new {pdb}')
        vm.addline(f'set a [atomselect top all]')
        vm.addline(f'$a set occupancy 0.0')
        vm.addline(f'set c [atomselect top "{constrained_atoms_def}"]')
        vm.addline(f'$c set occupancy {force_constant}')
        vm.addline(f'$a writepdb {c_pdb}')
        vm.writescript()
        vm.runscript()
        self.update_statevars(statekey,c_pdb,mode='file')
    
    def inherit_state(self):
        """
        Inherits the state variables from a prior task if it exists.
        This method checks if there is a prior task and, if so, copies its state variables
        to the current task's state variables. This allows the current task to use the state from
        the prior task, which can be useful for tasks that depend on the results of previous tasks. 
        If there is no prior task, the state variables remain empty.
        """
        if self.prior:
            logger.debug(f'Inheriting state from prior task {self.prior.taskname}')
            self.statevars=self.prior.statevars.copy()

    def save_state(self,exts=[]):
        """
        Saves the current state of the task to files based on the specified extensions.
        This method updates the state variables with the file paths for each specified extension.
        It checks if the files exist and updates the state variables accordingly.
        If the file does not exist, it raises a FileNotFoundError in strict mode or
        logs a debug message in permissive mode.

        If ``pdb`` is in the list of extensions and ``vel`` is also in the state variables,
        it assumes that the binary velocity file is stale and removes it from the state variables.
        If ``coor`` is in the list of extensions and ``pdb`` is not, it converts the coordinate file
        to a PDB file and updates the state variables accordingly.
        If ``pdb`` is in the list of extensions and ``coor`` is not, it converts the PDB file to a coordinate file
        and updates the state variables accordingly.
        
        Parameters
        ----------
        exts : list
            A list of file extensions to save the state for. Each extension corresponds to a specific file type.
            The extensions should be provided without the leading dot (e.g., ``psf``, ``pdb``, ``coor``, ``vel``, ``xsc``).
            The method will update the state variables with the file paths for each extension.
            If the extension is ``xsc``, it uses ``permissive`` mode, allowing the file to be optional.
            For other extensions, it uses ``strict`` mode, meaning the file must exist.

        """
        for ext in exts:
            mode='strict'
            if ext=='xsc':
                mode='permissive'
            self.update_statevars(ext,f'{self.basename}.{ext}',vtype='file',mode=mode)
        if 'pdb' in exts and 'vel' in self.statevars: # if pdb has changed, assume binary velocity file is stale
            del self.statevars['vel']
        if 'coor' in exts and 'pdb' not in exts:
            self.coor_to_pdb()
            self.update_statevars('pdb',f'{self.basename}.pdb',vtype='file')
        elif 'pdb' in exts and 'coor' not in exts:
            self.pdb_to_coor()
            self.update_statevars('coor',f'{self.basename}.coor',vtype='file')

    def update_statevars(self,key,value,vtype='',mode='strict'):
        """ 
        Updates the state variables with a new key-value pair.
        This method checks if the value is a file and whether it exists based on the specified mode.
        If the value is a file and does not exist, it raises a FileNotFoundError in strict mode or logs a debug message in permissive mode.
        If the value is not a file or exists, it updates the state variables with the new key-value pair.   
        
        Parameters
        ----------
        key : str
            The key under which the value will be stored in the state variables.
        value : str
            The value to be stored in the state variables. This can be a file path or any other string.
        vtype : str, optional
            The type of the value. If 'file', it indicates that the value should be treated as a file path.
            If not specified, it defaults to an empty string, meaning the value is not treated as a file.
        mode : str, optional
            The mode for checking the existence of the file. It can be ``strict`` or ``permissive``.
            In ``strict`` mode, if the file does not exist, a FileNotFoundError is raised.
            In ``permissive`` mode, if the file does not exist, a debug message is logged instead.    
        """
        if vtype=='file':
            if not os.path.exists(value):
                if mode=='strict':
                    raise FileNotFoundError(f'Expected file "{value}" not found.')
                else:
                    logger.debug(f'Attempt to update state {key} with "{value}" failed since {value} is not a file.')
            else:
                self.statevars[key]=value
        else:
            self.statevars[key]=value

    def write_statefile(self):
        """
        Writes the current state variables to a YAML file.
        This method creates a YAML file with the current state variables of the task.
        The filename is determined by the ``statefile`` key in the task's specifications.
        If the ``statefile`` key is not present, it defaults to ``{self.basename}-state.yaml``.
        The state variables are written in YAML format, allowing for easy serialization and deserialization
        of the task's state.
        """
        statefilename=self.specs.get('statefile',f'{self.basename}-state.yaml')
        with open(statefilename,'w') as f:
            yaml.dump(self.statevars,f)

    def copy_state(self,exts=[]):
        """
        Copies the state variables to files based on the specified extensions.
        This method iterates over the provided extensions and checks if each extension is present in the state variables.
        If the extension is found, it copies the corresponding file to a new file named ``{self.basename}.{ext}``.
        The copied file is then updated in the state variables with the new file path.
        This allows the task to inherit the state from previous tasks or to save the current state to files for later use. 

        Parameters
        ----------  
        exts : list
            A list of file extensions to copy the state for. Each extension corresponds to a specific file type.
            The extensions should be provided without the leading dot (e.g., ``psf``, ``pdb``, ``coor``, ``vel``, ``xsc``).
            The method will copy the files corresponding to each extension and update the state variables accordingly.
        """
        copies=[]
        for ext in exts:
            if ext in self.statevars:
                ffile=f'{self.basename}.{ext}'
                shutil.copy(self.statevars[ext],ffile)
                self.update_statevars(ext,ffile,vtype='file')
                copies.append(ffile)
        logger.debug(f'Copy state inherited {", ".join(copies)}')