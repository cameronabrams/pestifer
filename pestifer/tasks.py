# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" Task definitions

    This module defines several types of tasks that pestifer can 
    perform sequentially as part of a run.

    All tasks inherit BaseTask, and any task that is specified in the 
    config must have a yaml_header member.

"""
from .basemod import BaseMod
import logging
logger=logging.getLogger(__name__)
import shutil
import os
# from argparse import Namespace
import yaml
from .util import is_periodic
from .molecule import Molecule
from .chainids import ChainIDManager
from .colvars import *
from .stringthings import FileCollector
from .modcontainer import ModContainer
from .command import Command

class BaseTask(BaseMod):
    """ A base class for Tasks.
    
    Attributes
    ----------
    req_attr: list
        * specs: dictionary of specifications; under 'directives' in yaml
        * writers: dictionary of FileWriters
        * prior: name of prior task in sequence
        * index: unique integer index of task in run
        * config: access to the run config
        * taskname: caller-supplied name of task
    
    subtaskcount: int
        a count of subtasks in the task
    
    basename: string
        a basename used in file-naming conventions for the task

    statevars: dict
        dictionary that captures the current state of the run

    FC: FileCollector
        a structure for tracking files created by the task

    Methods
    -------

    next_basename(*obj):
        determines the next basename based on the task index, subtask
        label (obj[0]) and a subtask count

    coor_to_pdb():
        Generates and executes a VMD script to generate a new PDB
        file from an existing psf and coor file.
    
    pdb_to_coor(): 
        Generates and executes a VMD script to generate a new coor
        file from an existing psf and PDB file.
    
    make_constraint_pdb(dict): 
        Generates and executes a VMD script (based on the built-in
        make_constraint_pdb script) that generates a PDB file with 
        the appropriate attribute-tagging based on the directives in 
        the input dict.

    inherit_state(): 
        copies the statevar dict from the previous task onto this
        task's statevar dict

    save_state(list): 
        Based on the current basename, stores the current names of 
        all files whose extension types are in the list
    
    update_statevars(key,value): 
        Updates this task's statevars dict subject to some controls.

    """
    req_attr=BaseMod.req_attr+['specs','config','index','prior','writers','taskname']
    yaml_header='generic_task'
    # exts=['.psf','.pdb','.coor','.xsc'] # extensions of files that can be transferred from one task to the next
    _taskcount=0

    def __init__(self,input_dict,taskname,config,writers,prior):
        specs=input_dict.copy()
        input_dict = {
            'index':BaseTask._taskcount,
            'writers': writers,
            'prior':prior,
            'specs':specs,
            'config':config,
            'taskname':taskname
        }
        # for k,v in specs.items():
        #     if not v or (type(v)==str and v.lower()=='none'):
        #         specs[k]={}
        super().__init__(input_dict)
        BaseTask._taskcount+=1
        self.subtaskcount=0
        self.statevars={}
        self.FC=FileCollector()

    def next_basename(self,*obj):
        label=''
        if len(obj)==1 and len(obj[0])>0:
            label=f'-{obj[0]}'
        default_basename=f'{self.index:02d}-{self.subtaskcount:02d}-{self.taskname}{label}'
        overwrite_basename=self.specs.get('basename',None)
        if overwrite_basename:
            self.basename=overwrite_basename
        else:
            self.basename=default_basename
        self.subtaskcount+=1
    
    def coor_to_pdb(self):
        vm=self.writers['vmd']
        vm.newscript(f'{self.basename}-coor2pdb')
        vm.usescript('namdbin2pdb')
        vm.writescript()
        psf=self.statevars['psf']
        coor=self.statevars['coor']
        vm.runscript(psf=psf,coor=coor,pdb=f'{self.basename}.pdb')

    def pdb_to_coor(self):
        vm=self.writers['vmd']
        vm.newscript(f'{self.basename}-pdb2coor')
        vm.usescript('pdb2namdbin')
        vm.writescript()
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        vm.runscript(psf=psf,pdb=pdb,coor=f'{self.basename}.coor')

    def make_constraint_pdb(self,specs,statekey='consref'):
        vm=self.writers['vmd']
        pdb=self.statevars['pdb']
        force_constant=specs.get('k',self.config['user']['namd2']['harmonic']['spring_constant'])
        vm.newscript(f'{self.basename}-make-constraint-pdb')
        vm.usescript('make_constraint_pdb')
        vm.writescript()
        logger.debug(f'constraint spec: {specs["atoms"]}')
        c_pdb=specs.get('consref','')
        if not c_pdb:
            c_pdb=f'{self.basename}-constraints.pdb'
        vm.runscript(
            pdb=pdb,
            refpdb=c_pdb,
            constrained_atoms_def=','.join(specs['atoms'].split()),
            force_constant=force_constant)
        self.update_statevars(statekey,c_pdb,mode='file')
    
    def inherit_state(self):
        if self.prior:
            logger.debug(f'Inheriting state from prior task {self.prior.taskname}')
            self.statevars=self.prior.statevars.copy()

    def save_state(self,exts=[]):
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
        statefilename=self.specs.get('statefile',f'{self.basename}-state.yaml')
        with open(statefilename,'w') as f:
            yaml.dump(self.statevars,f)

    def copy_state(self,exts=[]):
        copies=[]
        for ext in exts:
            if ext in self.statevars:
                ffile=f'{self.basename}.{ext}'
                shutil.copy(self.statevars[ext],ffile)
                self.update_statevars(ext,ffile,vtype='file')
                copies.append(ffile)
        logger.debug(f'Copy state inherited {", ".join(copies)}')
        

class MDTask(BaseTask):
    """ A class for handling all NAMD2 runs
    
    Attributes
    ----------
    yaml_header(str):
        the name used to declare this task in an input yaml file

    Methods
    -------
    
    do():
        Inherits state and performs namd2 run based on specs; updates run state

    namd2run(): 
        Generates the NAMD2 config file based on specs and then executes NAMD2

    """
    yaml_header='md'

    def do(self):
        logger.info(f'Task {self.taskname} ({self.specs["ensemble"]}) {self.index:02d} initiated')
        self.inherit_state()            
        self.namd2run()
        self.save_state(exts=['coor','vel','xsc'])
        logger.info(f'Task {self.taskname} {self.index:02d} complete')

    def copy_charmmpar_local(self):
        local_names=[]
        na=self.writers['namd2']
        namd2_params_abs=na.standard_charmmff_parfiles+na.custom_charmmff_parfiles
        for nf in namd2_params_abs:
            d,n=os.path.split(nf)
            shutil.copy(nf,n)
            local_names.append(n)
        return local_names

    def namd2run(self,baselabel='',absolute_paths=True,extras={},script_only=False):
        specs=self.specs
        logger.debug(f'md task specs {specs}')
        ensemble=specs['ensemble']
        if not baselabel:
            self.next_basename(ensemble)
        else:
            self.next_basename(baselabel)
        params={}
        namd_global_params=self.config['user']['namd2']
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        coor=self.statevars.get('coor',None)
        vel=self.statevars.get('vel',None)
        xsc=self.statevars.get('xsc',None)
        cell=self.statevars.get('cell',None)
        firsttimestep=self.statevars.get('firsttimestep',0)
        self.statevars['periodic']=is_periodic(cell,xsc)

        temperature=specs['temperature']
        pressure=specs['pressure']
        params['tcl']=[]
        params['tcl'].append(f'set temperature {temperature}')

        nsteps=specs['nsteps']
        dcdfreq=specs['dcdfreq']
        xstfreq=specs['xstfreq']

        constraints=specs.get('constraints',{})
        other_params=specs.get('other_parameters',{})
        params.update(namd_global_params['generic'])
        params['structure']=psf
        params['coordinates']=pdb
        params['temperature']='$temperature'

        if coor:
            params['bincoordinates']=coor
        if vel:
            params['binvelocities']=vel
            del params['temperature']
        
        na=self.writers['namd2']
        if absolute_paths:
            params['parameters']=na.standard_charmmff_parfiles+na.custom_charmmff_parfiles
        else:
            params['parameters']=[]
            namd2_params_abs=na.standard_charmmff_parfiles+na.custom_charmmff_parfiles
            for nf in namd2_params_abs:
                d,n=os.path.split(nf)
                params['parameters'].append(n)
        
        if xsc or cell:
            if xsc:
                params['extendedSystem']=xsc
            else:
                params['tcl'].append(f'source {cell}')
            if ensemble=='NPT' and xstfreq:
                params['xstfreq']=xstfreq
        
        if self.statevars['periodic']:
            params.update(namd_global_params['solvated'])

        else:
            params.update(namd_global_params['vacuum'])
        
        if ensemble in ['NPT','NVT']:
            params.update(namd_global_params['thermostat'])
            if ensemble=='NPT':
                if not self.statevars['periodic']:
                    raise Exception(f'Cannot use barostat on a system without PBCs')
                params['tcl'].append(f'set pressure {pressure}')
                params.update(namd_global_params['barostat'])
        params['outputName']=f'{self.basename}'

        if dcdfreq:
            params['dcdfreq']=dcdfreq
        if 'restartfreq' in specs:
            params['restartfreq']=specs['restartfreq']
        if constraints:
            self.make_constraint_pdb(constraints)
            params['constraints']='on'
            params['consref']=self.statevars['consref']
            params['conskfile']=self.statevars['consref']
            params['conskcol']='O'
        params.update(other_params)
        params.update(extras)
        params['firsttimestep']=firsttimestep
        if ensemble=='minimize':
            assert specs['minimize']>0,f'Error: you must specify how many minimization cycles'
            params['minimize']=specs['minimize']
        elif ensemble in ['NVT','NPT']:
            assert nsteps>0,f'Error: you must specify how many time steps to run'
            params['run']=nsteps
        
        na.newscript(self.basename)
        na.writescript(params)
        if not script_only:
            na.runscript()
            if ensemble!='minimize':
                self.update_statevars('firsttimestep',firsttimestep+nsteps)
            else:
                self.update_statevars('firsttimestep',firsttimestep+specs['minimize'])
            if cell: # this is a use-once statevar
                del self.statevars['cell']
        return params

class PsfgenTask(BaseTask):
    """ A class for handling invocations of psfgen
    
    Attributes
    ----------
    yaml_header(str) 

    Methods
    -------
    do(): 
        Based on specs, reads in input PDB/mmCIF file, generates parsed Molecule instances, generates
        the psfgen script, and executes VMD to perform the psfgen run.

    """
    yaml_header='psfgen'
    def __init__(self,input_dict,taskname,config,writers,prior):
        super().__init__(input_dict,taskname,config,writers,prior)
        self.chainIDmanager=ChainIDManager(format=self.specs['source']['file_format'])

    def do(self):
        logger.info(f'Task {self.taskname} {self.index:02d} initiated')
        self.inherit_state()
        logger.debug('Injesting molecule(s)')
        self.injest_molecules()
        self.statevars['base_molecule']=self.base_molecule
        logger.debug('Running first psfgen')
        self.psfgen()
        # we now have a full coordinate set, so we can do coormods
        self.coormods()
        min_loop_length=self.specs['source']['sequence']['loops']['min_loop_length']
        nloops=self.base_molecule.has_loops(min_loop_length=min_loop_length)*self.base_molecule.num_images()
        if nloops>0:
            logger.debug(f'Declashing {nloops} loops')
            self.declash_loops(self.specs['source']['sequence']['loops'])
        logger.info(f'Task {self.taskname} {self.index:02d} complete')

    def coormods(self):
        if sum([len(x) for x in self.mods.coormods.__dict__.values()])>0:
            ba=self.base_molecule.active_biological_assembly
            logger.debug(f'Performing coormods')
            self.next_basename('coormods')
            vm=self.writers['vmd']
            vm.newscript(self.basename)
            psf=self.statevars['psf']
            pdb=self.statevars['pdb']
            vm.load_psf_pdb(psf,pdb,new_molid_varname='mCM')
            for transform in ba.transforms:
                self.mods.coormods.crotations.write_TcL(vm,chainIDmap=transform.chainIDmap)
            vm.write_pdb(self.basename,'mCM')
            vm.endscript()
            vm.writescript()
            vm.runscript()
            self.save_state(exts=['pdb'])

    def psfgen(self):
        self.next_basename('build')
        pg=self.writers['psfgen']
        pg.newscript(self.basename)
        pg.topo_aliases()
        pg.set_molecule(self.base_molecule,altcoords=self.specs['source'].get('altcoords',None))
        pg.describe_molecule(self.base_molecule)
        pg.complete(self.basename)
        pg.endscript()
        pg.writescript()
        pg.runscript()
        self.save_state(exts=['psf','pdb'])
        pg.cleanup(cleanup=self.specs['cleanup'])

    def declash_loops(self,specs):
        mol=self.base_molecule
        cycles=specs['declash']['maxcycles']
        self.update_statevars('min_loop_length',specs['min_loop_length'])
        if not mol.has_loops() or not cycles:
            return
        self.next_basename('declash')
        vt=self.writers['vmd']
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        vt.newscript(self.basename)
        vt.load_psf_pdb(psf,pdb,new_molid_varname='mLL')
        mol.write_loop_lines(vt,cycles=cycles,min_length=specs['min_loop_length'])
        vt.write_pdb(self.basename,'mLL')
        vt.endscript()
        vt.writescript()
        vt.runscript()
        self.save_state(exts=['pdb'])

    def injest_molecules(self):
        specs=self.specs
        self.molecules={}
        self.source_specs=specs['source']
        logger.debug(f'User-input modspecs {self.specs["mods"]}')
        self.mods=ModContainer(self.specs['mods'])
        # self.pdbs=self.mods.report_pdbs()
        # self.usermod_specs=specs['mods']
        # logger.debug(f'user mods at injest_molecules {self.mods.__dict__}')
        self.molecules[self.source_specs['id']]=Molecule(source=self.source_specs,usermods=self.mods,chainIDmanager=self.chainIDmanager).activate_biological_assembly(self.source_specs['biological_assembly'])
        self.base_molecule=self.molecules[self.source_specs['id']]
        # self.statevars['min_loop_length']=self.source_specs['sequence']['loops']['min_loop_length']
        # for p in self.pdbs:
        #     self.molecules[p]=Molecule(source=p)

class DomainSwapTask(MDTask):
    yaml_header='domainswap'

    # def __init__(self,input_dict,taskname,config,writers,prior):
    #     super().__init__(input_dict,taskname,config,writers,prior)

    def do(self):
        logger.info(f'Task {self.taskname} {self.index:02d} initiated')
        self.inherit_state()
        logger.debug(f'Generating inputs for domain swap')
        self.make_inputs()
        logger.debug(f'Running NAMD to execute domain swap')
        self.namd2run(baselabel='domainswap-run',extras={'colvars':'on','colvarsconfig':self.statevars['cv']})
        self.save_state(exts=['vel','coor'])

    def make_inputs(self):
        specs=self.specs
        self.next_basename('domainswap-prep')
        vm=self.writers['vmd']
        vm.newscript(self.basename)
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        vm.usescript('domainswap')
        vm.writescript()
        vm.runscript(
            psf=psf,
            pdb=pdb,
            swap_domain_def=','.join(specs['swap_domain_def'].split()),
            anchor_domain_def=','.join(specs['anchor_domain_def'].split()),
            chain_swap_pairs=':'.join([','.join(x) for x in specs['chain_directional_swaps']]),
            force_constant=specs['force_constant'],
            target_numsteps=specs['target_numsteps'],
            cv=f'{self.basename}-cv.inp',
            refpdb=f'{self.basename}-ref.pdb')
        self.update_statevars('cv',f'{self.basename}-cv.inp',vtype='file')
        
class LigateTask(MDTask):
    yaml_header='ligate'
    # statevars={}
    # def __init__(self,input_dict,taskname,config,writers,prior):
    #     super().__init__(input_dict,taskname,config,writers,prior)

    def do(self):
        logger.info(f'Task {self.taskname} {self.index:02d} initiated')
        self.inherit_state()
        self.base_molecule=self.statevars['base_molecule']
        if not self.base_molecule.has_loops(min_loop_length=self.statevars['min_loop_length']):
            logger.info('No loops. Ligation bypassed.')
            return
        logger.debug('Storing sequence gaps.')
        self.write_gaps()
        logger.debug('Measuring gap distances.')
        self.measure_distances(self.specs['steer'])
        logger.debug('Steering loop C-termini toward their partner N-termini')
        self.do_steered_md(self.specs['steer'])
        self.save_state(exts=['coor','vel'])
        logger.debug('Connecting loop C-termini to their partner N-termini')
        self.connect()
        self.save_state(exts=['psf','pdb'])
        logger.info(f'Task {self.taskname} {self.index:02d} complete')

    def write_gaps(self):
        self.next_basename('gaps')
        mol=self.base_molecule
        datafile=f'{self.basename}.inp'
        writer=self.writers['data']
        writer.newfile(datafile)
        mol.write_gaps(writer)
        writer.writefile()
        self.update_statevars('data',datafile,vtype='file')

    def measure_distances(self,specs):
        receiver_flexible_zone_radius=specs['receiver_flexible_zone_radius']
        comment_chars='#!$'
        self.next_basename('measure')
        resultsfile=f'{self.basename}.dat'
        vm=self.writers['vmd']
        vm.newscript(self.basename)
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        vm.usescript('measure_bonds')
        vm.writescript()
        datafile=self.statevars['data']
        vm.runscript(psf=psf,pdb=pdb,i=datafile,opdb=f'{self.basename}.pdb',o=resultsfile,rfzr=receiver_flexible_zone_radius)
        self.update_statevars('fixedref',f'{self.basename}.pdb',vtype='file')
        self.update_statevars('results',resultsfile,vtype='file')
        with open(resultsfile,'r') as f:
            datalines=f.read().split('\n')
        self.gaps=[]
        for line in datalines:
            if len(line)>0 and not line[0] in comment_chars:
                data=line.split()
                thisgap={
                    'segname':data[0],
                    'serial_i':int(data[1]),
                    'serial_j':int(data[2]),
                    'distance':float(data[3])
                }
                self.gaps.append(thisgap)

    def do_steered_md(self,specs):
        self.next_basename('steer')
        writer=self.writers['data']
        writer.newfile(f'{self.basename}-cv.inp')
        for i,g in enumerate(self.gaps):
            g['name']=f'GAP{i:02d}'
            declare_distance_cv_atoms(g,writer)
        for i,g in enumerate(self.gaps):
            g['k']=specs['force_constant']
            g['targ_distance']=specs['target_distance']
            g['targ_numsteps']=specs['nsteps']
            declare_harmonic_distance_bias(g,writer)
        writer.writefile()
        savespecs=self.specs
        self.specs=specs
        self.namd2run(extras={        
            'fixedatoms':'on',
            'fixedatomsfile':self.statevars['fixedref'],
            'fixedatomscol': 'O',
            'colvars': 'on',
            'colvarsconfig': f'{self.basename}-cv.inp'
        })
        self.specs=savespecs

    def connect(self):
        self.write_connect_patches()
        self.connect_gaps()

    def write_connect_patches(self):
        self.next_basename('gap_patches')
        mol=self.base_molecule
        datafile=f'{self.basename}.inp'
        writer=self.writers['data']
        writer.newfile(datafile)
        mol.write_connect_patches(writer)
        writer.writefile()
        self.update_statevars('data',datafile,vtype='file')

    def connect_gaps(self):
        self.next_basename('heal')
        pg=self.writers['psfgen']
        pg.newscript(self.basename)
        pg.topo_aliases()
        topfile=os.path.join(self.config.charmmff_custom_path,'mylink.top')
        pg.addline(f'topology {topfile}')
        pg.usescript('loop_closure')
        pg.writescript()
        patchfile=self.statevars['data']
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        pg.runscript(psf=psf,pdb=pdb,p=patchfile,o=self.basename)
        self.save_state(exts=['psf','pdb'])

class SolvateTask(BaseTask):
    yaml_header='solvate'
    # opt_attr=Task.opt_attr+[yaml_header]

    def do(self):
        self.statevars=self.prior.statevars.copy()
        if 'xsc' in self.statevars:
            del self.statevars['xsc']
        self.next_basename()
        logger.info(f'Task {self.taskname} {self.index:02d} initiated')
        vt=self.writers['vmd']
        vt.newscript(self.basename)
        vt.usescript('solv')
        vt.writescript()
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        vt.runscript(o=self.basename,pdb=pdb,psf=psf,pad=self.specs['pad'])
        self.save_state(exts=['psf','pdb'])
        self.update_statevars('cell',f'{self.basename}_cell.tcl',vtype='file')
        logger.info(f'Task {self.taskname} {self.index:02d} complete')

# class RelaxTask(MDTask):
#     yaml_header='relax'
#     # opt_attr=Task.opt_attr+[yaml_header]
#     def do(self):
#         logger.info(f'Task {self.taskname} {self.index:02d} initiated')
#         self.inherit_state()
#         self.next_basename('relax')
#         self.namd2script(basename,self.namd2prep(basename,self.specs))
#         self.relax(basename)
#         logger.info(f'Task {self.taskname} {self.index:02d} complete')

# class MinimizeTask(Task):
#     yaml_header='minimize'
#     opt_attr=Task.opt_attr+[yaml_header]
#     def do(self):
#         logger.info(f'Task {self.taskname} {self.index:02d} initiated')
#         self.statevars=self.prior.statevars.copy()
#         basename=self.next_basename('minimize')
#         self.minimize(basename,self.specs)
#         logger.info(f'Task {self.taskname} {self.index:02d} complete')

class ManipulateTask(BaseTask):
    yaml_header='manipulate'
    # opt_attr=Task.opt_attr+[yaml_header]
    def do(self):
        logger.info(f'Task {self.taskname} {self.index:02d} initiated')
        if self.prior:
            logger.debug(f'Task {self.taskname} prior {self.prior.taskname}')
            self.statevars=self.prior.statevars.copy()
        self.coormods(self.specs['mods'])
        # self.minimize(self.specs['minimize'])

    def coormods(self,specs):
        self.mods=ModContainer(specs)
        if self.mods.coormods:
            logger.debug(f'performing coormods')
            self.next_basename('coormods')
            vm=self.writers['vmd']
            vm.newscript(self.basename)
            psf=self.statevars['psf']
            pdb=self.statevars['pdb']
            vm.load_psf_pdb(psf,pdb,new_molid_varname='mCM')
            self.mods.coormods.crotations.write_TcL(vm)
            vm.write_pdb(self.basename,'mCM')
            vm.endscript()
            vm.writescript()
            vm.runscript()
            self.save_state(exts=['pdb'])

class TerminateTask(MDTask):
    yaml_header='terminate'
    # opt_attr=Task.opt_attr+[yaml_header]
    def do(self):
        logger.info(f'Task {self.taskname} {self.index:02d} initiated')
        self.inherit_state()
        self.next_basename()
        self.copy_state(exts=['psf','pdb','coor','xsc','vel'])
        self.write_chainmaps()
        self.write_statefile()
        self.make_package()
        logger.info(f'Task {self.taskname} {self.index:02d} complete')

    def write_chainmaps(self):
        bm=self.statevars.get('base_molecule',None)
        if bm:
            maps=bm.get_chainmaps()
            with open(self.specs['chainmapfile'],'w') as f:
                yaml.dump(maps,f)
            del self.statevars['base_molecule']

    def make_package(self):
        specs=self.specs.get('package',{})
        # logger.debug(f'make_package specs {specs}')
        if not specs:
            return
        self.inherit_state()
        self.FC.clear()  # populate a file collector to make the tarball
        logger.debug(f'Packaging for namd2 using basename {self.basename}')
        savespecs=self.specs
        self.specs=specs
        params=self.namd2run(script_only=True,absolute_paths=False)
        self.specs=savespecs
        self.FC.append(f'{self.basename}.namd')
        constraints=specs.get('constraints',{})
        if constraints:
            self.make_constraint_pdb(constraints)
            self.FC.append(self.statevars['consref'])
        local_params=self.copy_charmmpar_local()
        for n in local_params:
            self.FC.append(n)
        for ext in ['psf','pdb','coor','xsc','vel']:
            self.FC.append(self.statevars[ext])

        if specs["topogromacs"]:
            logger.debug(f'running topogromacs')
            with open(f'{self.basename}_par.inp','w') as f:
                for pf in params['parameters']:
                    f.write(f'parameters {pf}\n')
            vt=self.writers['vmd']
            vt.newscript(f'{self.basename}_tg')
            vt.usescript('tg')
            vt.writescript()
            psf=self.statevars['psf']
            pdb=self.statevars['pdb']
            inputname=os.path.splitext(self.statevars['coor'])[0]
            vt.runscript(o=self.basename,pdb=pdb,psf=psf,i=inputname,parinp=f'{self.basename}_par.inp',ospf=f'{self.basename}_tg.psf',opdb=f'{self.basename}_tg.pdb',top=f'{self.basename}_topogromacs.top',cellfile=f'{self.basename}_cell.inp')
            with open(f'{self.basename}_cell.inp','r') as f:
                box=f.read().split()
            boxstr=' '.join(box)
            c=Command(f'gmx editconf -f {self.basename}_tg.pdb -o {self.basename}_topogromacs.pdb -box {boxstr}')
            c.run()
            self.FC.append(f'{self.basename}_topogromacs.pdb')
            self.FC.append(f'{self.basename}_topogromacs.top')
        self.FC.tarball(specs["basename"])
