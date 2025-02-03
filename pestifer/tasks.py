# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" Task definitions

    This module defines several types of tasks that pestifer can 
    perform sequentially as part of a run.

    All tasks inherit BaseTask, and any task that is specified in the 
    config must have a yaml_header member.

"""
import logging
import os
import shutil
import yaml
import matplotlib.pyplot as plt 
import networkx as nx
import numpy as np
import pandas as pd

from copy import deepcopy
from scipy.constants import physical_constants

from .basemod import BaseMod
from .chainidmanager import ChainIDManager
from .colvars import *
from .command import Command
from .molecule import Molecule
from .modmanager import ModManager
from .mods import CleavageSite, CleavageSiteList
from .progress import PestiferProgress
from .psf import PSFContents
from .ring import ring_check
from .stringthings import FileCollector
from .util import is_periodic, cell_from_xsc, cell_to_xsc

logger=logging.getLogger(__name__)
logging.getLogger("matplotlib").setLevel(logging.WARNING)

g_per_amu=physical_constants['atomic mass constant'][0]*1000
A_per_cm=1.e8
A3_per_cm3=A_per_cm**3

class BaseTask(BaseMod):
    """ A base class for Tasks.
    
    Attributes
    ----------
    req_attr: list
        * `specs`: dictionary of specifications; under '`directives`' in yaml
        * `writers`: dictionary of `FileWriters`
        * `prior`: identifier of prior task in sequence
        * `index`: unique integer index of task in run
        * `config`: access to the run config
        * `taskname`: caller-supplied name of task
    
    `subtaskcount`: int
        a count of subtasks in the task
    
    `basename`: string
        a basename used in file-naming conventions for the task

    `statevars`: dict
        dictionary that captures the current state of the run

    `FC`: `FileCollector`
        a structure for tracking files created by the task

    Methods
    -------

    `next_basename(*obj)`:
        determines the next basename based on the task index, subtask
        label (obj[0]) and a subtask count

    `coor_to_pdb()`:
        Generates and executes a VMD script to generate a new PDB
        file from an existing psf and coor file.
    
    `pdb_to_coor()`: 
        Generates and executes a VMD script to generate a new coor
        file from an existing psf and PDB file.
    
    `make_constraint_pdb(dict)`: 
        Generates and executes a VMD script that generates a PDB file 
        with the appropriate attribute-tagging based on the directives 
        in the input dict.

    `inherit_state()`: 
        copies the statevar dict from the previous task onto this
        task's statevar dict

    `save_state(list)`: 
        Based on the current basename, stores the current names of 
        all files whose extension types are in the list
    
    `update_statevars(key,value)`: 
        Updates this task's statevars dict subject to some controls.

    """
    req_attr=BaseMod.req_attr+['specs','config','index','prior','writers','taskname']
    yaml_header='generic_task'
    _taskcount=0
    init_msg_options=['INITIATED','STARTED','BEGUN','SET IN MOTION','KICKED OFF','LIT','SPANKED ON THE BOTTOM']

    def __init__(self,input_dict,taskname,config,writers,prior):
        specs=input_dict.copy()
        BaseTask._taskcount=specs.get('index',BaseTask._taskcount)
        logger.debug(f'Creating task {taskname} with index {BaseTask._taskcount}')
        input_dict = {
            'index':BaseTask._taskcount,
            'writers': writers,
            'prior':prior,
            'specs':specs,
            'config':config,
            'taskname':taskname
        }
        super().__init__(input_dict)
        BaseTask._taskcount+=1
        self.subtaskcount=0
        self.statevars={}
        self.FC=FileCollector()
        self.result=0

    def do(self):
        return self.result

    def log_message(self,message,**kwargs):
        extra=''
        for k,v in kwargs.items():
            if v:
                extra+=f' ({k}: {v})'
        mtoks=[x.strip() for x in [x.upper() for x in message.split()]]
        if not any([x in self.init_msg_options for x in mtoks]):
            extra+=f' (result: {self.result})'
        logger.info(f'Task {self.index:02d} \'{self.taskname}\'{extra} {message}')

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
        psf=self.statevars['psf']
        coor=self.statevars['coor']
        vm.addline(f'namdbin2pdb {psf} {coor} {self.basename}.pdb')
        vm.writescript()
        vm.runscript()

    def pdb_to_coor(self):
        vm=self.writers['vmd']
        vm.newscript(f'{self.basename}-pdb2coor')
        pdb=self.statevars['pdb']
        vm.addline(f'pdb2namdbin {pdb} {self.basename}.coor')
        vm.writescript()
        vm.runscript()

    def make_constraint_pdb(self,specs,statekey='consref'):
        vm=self.writers['vmd']
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
        
class RestartTask(BaseTask):
    """ This task only resets the workflow chain to named values of the psf, pdb, xsc, and coor files """
    yaml_header='restart'
    def do(self):
        self.log_message('initiated')
        self.next_basename()
        exts_actual=[]
        for ext in ['psf','coor','pdb','xsc']:
            fname=self.specs.get(ext,'')
            if fname:
                shutil.copy(fname,f'{self.basename}.{ext}')
                exts_actual.append(ext)
        self.save_state(exts=exts_actual)
        self.log_message('complete')
        self.result=0
        return super().do()

class MDTask(BaseTask):
    """ A class for handling all NAMD runs
    
    Attributes
    ----------
    `yaml_header(str)`:
        the name used to declare this task in an input yaml file

    Methods
    -------
    
    `do()`:
        Inherits state and performs namd run based on specs; updates run state

    `copy_charmmpar_local()`:abbr:
        Copies all charmm parameter files from pestifer's resource library to
        the current working directory; returns the list of file names that
        were copied

    `namdrun()`: 
        Generates the NAMD config file based on specs and then executes NAMD

    """
    yaml_header='md'

    def do(self):
        self.log_message('initiated',ensemble=self.specs.get('ensemble',None))
        self.inherit_state()            
        self.result=self.namdrun()
        if self.result==0: self.save_state(exts=['coor','vel','xsc'])
        self.log_message('complete',ensemble=self.specs.get('ensemble',None))
        return super().do()

    def namdrun(self,baselabel='',extras={},script_only=False,**kwargs):
        specs=self.specs
        logger.debug(f'md task specs {specs}')
        ensemble=specs['ensemble']
        if not baselabel:
            self.next_basename(ensemble)
        else:
            self.next_basename(baselabel)
        
        params={}
        namd_global_params=self.config['user']['namd']
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        coor=self.statevars.get('coor',None)
        vel=self.statevars.get('vel',None)
        xsc=self.statevars.get('xsc',None)
        # cell=self.statevars.get('cell',None)
        firsttimestep=self.statevars.get('firsttimestep',0)
        self.statevars['periodic']=is_periodic(xsc)

        temperature=specs['temperature']
        if ensemble=='NPT':
            pressure=specs['pressure']
        params['tcl']=[]
        params['tcl'].append(f'set temperature {temperature}')

        nsteps=specs['nsteps']
        dcdfreq=specs['dcdfreq']
        xstfreq=specs['xstfreq']

        constraints=specs.get('constraints',{})
        other_params=specs.get('other_parameters',{})
        colvars=specs.get('colvar_specs',{})
        params.update(namd_global_params['generic'])
        params['structure']=psf
        params['coordinates']=pdb
        params['temperature']='$temperature'

        if coor:
            params['bincoordinates']=coor
        if vel:
            params['binvelocities']=vel
            del params['temperature']
        
        na=self.writers['namd']
        na.update_par()
        self.local_parameter_files=na.copy_charmm_par()
        params['parameters']=self.local_parameter_files
        
        if xsc:
            params['extendedSystem']=xsc
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
        if colvars:
            writer=self.writers['data']
            writer.newfile(f'{self.basename}-cv.inp')
            colvar_writer(colvars,writer,pdb=pdb)
            writer.writefile()
            params['colvars']='on'
            params['colvarsconfig']=f'{self.basename}-cv.inp'

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
            local_execution_only=not self.statevars['periodic']
            single_gpu_only=kwargs.get('single_gpu_only',False)
            result=na.runscript(single_molecule=(not self.statevars['periodic']),local_execution_only=local_execution_only,single_gpu_only=single_gpu_only)
            if result!=0:
                return -1
            inherited_etitles=[]
            if self.prior and self.prior.taskname=='md' and hasattr(self.prior,'mdlog'):
                inherited_etitles=self.prior.mdlog.etitles
                logger.debug(f'Offering these etitles: {inherited_etitles}')
            self.mdlog=na.getlog(inherited_etitles=inherited_etitles)

            if os.path.exists(f'{self.basename}.xst'):
                self.xstlog=na.getxst()
                logger.debug(f'{self.xstlog.df.shape[0]} xst entries stored')

            if ensemble!='minimize':
                self.update_statevars('firsttimestep',firsttimestep+nsteps)
            else:
                self.update_statevars('firsttimestep',firsttimestep+specs['minimize'])
            if specs.get('remove_parfiles',False):
                na.remove_charmm_par(params['parameters'])
        return 0

class MDPlotTask(BaseTask):
    """ A class for making plots of energy-like quantities from a series of one or more NAMD 
        runs """
    yaml_header='mdplot'
    def do(self):
        self.log_message('initiated')
        self.inherit_state()            
        datasources=[]
        xstsources=[]
        root=self.prior
        while root.prior!=None and root.prior.taskname=='md':
            if hasattr(root.prior,'mdlog'):
                datasources.append(root.prior.mdlog.edata)
            if hasattr(root.prior,'xstlog'):
                xstsources.append(root.prior.xstlog.df)
            root=root.prior
        logger.debug(f'concatenating energy-like data from {len(datasources)} sequential logs')
        edata=pd.concat(datasources[::-1])
        savedata=self.specs.get('savedata',None)
        xstdata=pd.concat(xstsources[::-1])
        if savedata:
            logger.debug(f'Saving energy-like data to {savedata}.')
            try:
                edata.to_csv(savedata,header=True,index=False)
            except:
                logger.debug(f'For some reason, could not write a dataframe to csv')
                logger.debug(edata.iloc[:3,:].to_string())
            logger.debug(f'Saving cell data to xst-{savedata}.')
            try:
                xstdata.to_csv(f'xst-{savedata}',header=True,index=False)
            except:
                logger.debug(f'For some reason, could not write a dataframe to csv')
                logger.debug(xstdata.iloc[:3,:].to_string())
        
        traces=self.specs.get('traces',[])
        legend=self.specs.get('legend',False)
        basename=self.specs.get('basename','myplot')
        for trace in traces:
            unitspecs=[]
            figsize=self.specs.get('figsize',(9,6))
            fig,ax=plt.subplots(1,1,figsize=figsize)
            if type(trace)!=list:
                tracelist=[trace]
            else:
                tracelist=trace
            for t_i in tracelist:
                unitspec=self.specs.get('units',{}).get(t_i,'*')
                if unitspec=='*':
                    units=1.0
                else:
                    if unitspec=='g_per_cc':
                        units=g_per_amu*A3_per_cm3
                    else:
                        logger.debug(f'Unitspec "{unitspec}" not recognized.')
                        units=1.0
                unitspecs.append(unitspec)
                if t_i.upper() in edata:
                    key=t_i.upper()
                    ax.plot(edata['TS'],edata[key]*units,label=key.title())
                elif t_i in xstdata:
                    ax.plot(xstdata['step'],xstdata[t_i]*units,label=t_i)
            ax.set_xlabel('time step')
            tracename=','.join(tracelist)
            ax.set_ylabel(tracename+' ('+','.join([_ for _ in unitspecs if _!='*'])+')')
            if legend:
                plt.legend()
            tracename='-'.join(tracelist)
            plt.savefig(f'{basename}-{tracename}.png',bbox_inches='tight')
            plt.clf()
        self.log_message('complete')
        self.result=0
        return super().do()

class PsfgenTask(BaseTask):
    """ A class for handling invocations of psfgen which create a molecule from a base PDB/mmCIF file
    or from a PSF file generated previously by psfgen
    
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
        self.molecules={}

    def do(self):
        self.log_message('initiated')
        self.inherit_state()
        logger.debug('Injesting molecule(s)')
        self.injest_molecules()
        self.statevars['base_molecule']=self.base_molecule
        logger.debug(f'base mol num images {self.base_molecule.num_images()}')
        logger.debug('Running first psfgen')
        self.result=self.psfgen()
        if self.result!=0:
            return super().do()
        # we now have a full coordinate set, so we can do coormods
        self.coormods()
        # min_loop_length=0
        min_loop_length=self.specs['source'].get('sequence',{}).get('loops',{}).get('min_loop_length',0)
        self.update_statevars('min_loop_length',min_loop_length)
        nloops=self.base_molecule.has_loops(min_loop_length=min_loop_length)*self.base_molecule.num_images()
        if nloops>0 and self.specs['source']['sequence']['loops']['declash']['maxcycles']>0:
            logger.debug(f'Declashing {nloops} loops')
            self.declash_loops(self.specs['source']['sequence']['loops'])
        nglycans=self.base_molecule.nglycans()*self.base_molecule.num_images()
        if nglycans>0 and self.specs['source']['sequence']['glycans']['declash']['maxcycles']>0:
            logger.debug(f'Declashing {nglycans} glycan segments')
            self.declash_glycans(self.specs['source']['sequence']['glycans'])
        self.log_message('complete')
        return super().do()

    def coormods(self):
        coormods=self.modmanager.get('coormods',{})
        logger.debug(f'psfgen task has {len(coormods)} coormods')
        logger.debug(f'{coormods}')
        ba=self.base_molecule.active_biological_assembly
        if coormods:
            logger.debug(f'performing coormods')
            for modtype,modlist in coormods.items():
                if len(modlist)>0:
                    self.next_basename(modtype)
                    vm=self.writers['vmd']
                    vm.newscript(self.basename,packages=['PestiferCRot'])
                    psf=self.statevars['psf']
                    pdb=self.statevars['pdb']
                    vm.load_psf_pdb(psf,pdb,new_molid_varname='mCM')
                    for transform in ba.transforms:
                        modlist.write_TcL(vm,chainIDmap=transform.chainIDmap)
                    vm.write_pdb(self.basename,'mCM')
                    vm.writescript()
                    vm.runscript()
                    self.save_state(exts=['pdb'])

    def psfgen(self):
        self.next_basename('build')
        pg=self.writers['psfgen']
        pg.newscript(self.basename,packages=['PestiferCRot'])
        pg.set_molecule(self.base_molecule,altcoords=self.specs.get('source',{}).get('altcoords',None))
        pg.describe_molecule(self.base_molecule)
        pg.writescript(self.basename)
        result=pg.runscript()
        if result!=0:
            return result
        if self.specs.get('remove_charmmtop',False):
            pg.remove_charmm_top(pg.topologies)
        self.save_state(exts=['psf','pdb'])
        self.strip_remarks()
        return 0
        
    def strip_remarks(self):
        pdb=self.statevars['pdb']
        c=Command(f'grep -v ^REMARK {pdb} > tmp').run()
        shutil.move('tmp',pdb)

    def declash_loops(self,specs):
        mol=self.base_molecule
        cycles=specs['declash']['maxcycles']
        if not mol.has_loops() or not cycles:
            logger.debug(f'Loop declashing is intentionally not done.')
            return
        self.next_basename('declash-loops')
        vt=self.writers['vmd']
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        vt.newscript(self.basename,packages=['PestiferDeclash'])
        vt.load_psf_pdb(psf,pdb,new_molid_varname='mLL')
        mol.write_loop_lines(vt,cycles=cycles,min_length=specs['min_loop_length'],include_c_termini=specs['declash']['include_C_termini'])
        vt.write_pdb(self.basename,'mLL')
        vt.writescript()
        vt.runscript()
        self.save_state(exts=['pdb'])

    def declash_glycans(self,specs):
        mol=self.base_molecule
        cycles=specs['declash']['maxcycles']
        clashdist=specs['declash']['clashdist']
        if not mol.nglycans() or not cycles:
            logger.debug(f'Glycan declashing is intentionally not done.')
            return
        self.next_basename('declash-glycans')
        outpdb=f'{self.basename}.pdb'
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        vt=self.writers['vmd']
        vt.newscript(self.basename,packages=['PestiferDeclash'])
        vt.addline(f'mol new {psf}')
        vt.addline(f'mol addfile {pdb} waitfor all')
        vt.addline(f'set a [atomselect top all]')
        vt.addline(f'set molid [molinfo top get id]')
        nglycan=self._write_glycans(vt)
        vt.addline(f'vmdcon -info "Declashing $nglycans glycans; clashdist {clashdist}; maxcycles {cycles}"')
        vt.addline(r'for {set i 0} {$i<$nglycans} {incr i} {')
        vt.addline(f'   declash_pendant $molid $glycan_idx($i) $rbonds($i) $movers($i) {cycles} {clashdist}')
        vt.addline(r'}')
        vt.addline(f'$a writepdb {outpdb}')
        vt.writescript()
        logger.debug(f'Declashing {nglycan} glycans')
        vt.runscript(progress_title='declash-glycans')
        self.save_state(exts=['pdb'])

    def _write_glycans(self,fw):
        psf=self.statevars['psf']
        logger.debug(f'Injesting {psf}')
        struct=PSFContents(psf,parse_topology=['bonds'])
        logger.debug(f'Making graph structure of glycan atoms...')
        glycanatoms=struct.atoms.get(segtype='glycan')
        logger.debug(f'{len(glycanatoms)} total glycan atoms')
        glycangraph=glycanatoms.graph()
        G=[glycangraph.subgraph(c).copy() for c in nx.connected_components(glycangraph)]
        logger.debug(f'Preparing declash input for {len(G)} glycans')
        fw.addline(f'set nglycans {len(G)}')
        for i,g in enumerate(G):
            logger.debug(f'Glycan {i} has {len(g)} atoms')
            serials=[x.serial for x in g]
            for at in g:
                lig_ser=[x.serial for x in at.ligands]
                for k,ls in enumerate(lig_ser):
                    if not ls in serials:
                        at.is_root=True
                        rp=at.ligands[k]
                        logger.debug(f'-> Atom {str(at)} is the root, bound to atom {str(rp)}')
            indices=' '.join([str(x.serial-1) for x in g])
            fw.comment(f'Glycan {i}:')
            fw.addline(f'set glycan_idx({i}) [list {indices}]')
            fw.addline(f'set rbonds({i}) [list]')
            fw.addline(f'set movers({i}) [list]')
            for bond in nx.bridges(g):
                ai,aj=bond
                if not (ai.isH() or aj.isH()) and not ai.is_pep(aj):
                    g.remove_edge(ai,aj)
                    S=[g.subgraph(c).copy() for c in nx.connected_components(g)]
                    assert len(S)==2,f'Bond {ai.serial-1}-{aj.serial-1} when cut makes more than 2 components'
                    for sg in S:
                        is_root=any([hasattr(x,'is_root') for x in sg])
                        if not is_root:
                            if ai in sg:
                                sg.remove_node(ai)
                            if aj in sg:
                                sg.remove_node(aj)
                            if len(sg)>1 or (len(sg)==1 and not [x for x in sg.nodes][0].isH()):
                                mover_serials=[x.serial for x in sg]
                                mover_indices=" ".join([str(x-1) for x in mover_serials])
                                logger.debug(f'{str(ai)}--{str(aj)} is a rotatable bridging bond')
                                fw.addline(f'lappend rbonds({i}) [list {ai.serial-1} {aj.serial-1}]')
                                logger.debug(f'  -> movers: {" ".join([str(x) for x in sg])}')
                                fw.addline(f'lappend movers({i}) [list {mover_indices}]')
                    g.add_edge(ai,aj)
        return len(G)

    def injest_molecules(self):
        specs=self.specs
        self.source_specs=specs['source']
        logger.debug(f'User-input modspecs {self.specs["mods"]}')
        self.modmanager=ModManager(self.specs['mods'])
        seqmods=self.modmanager.get('seqmods',{})
        logger.debug(f'Injesting seqmods {seqmods}')
        if 'grafts' in seqmods:
            logger.debug(f'looking for graft sources to injest')
            Grafts=seqmods['grafts']
            for g in Grafts:
                if not g.source_pdbid in self.molecules:
                    logger.debug(f'Injesting graft source {g.source_pdbid}')
                    this_source={
                        'id':g.source_pdbid,
                        'file_format':'PDB'
                    }
                    self.molecules[g.source_pdbid]=Molecule(source=this_source)
                g.activate(deepcopy(self.molecules[g.source_pdbid]))
        self.chainIDmanager=ChainIDManager(
            format=self.source_specs['file_format'],
            transform_reserves=self.source_specs.get('transform_reserves',{}),
            remap=self.source_specs.get('remap_chainIDs',{}))
        self.base_molecule=Molecule(source=self.source_specs,
                                    modmanager=self.modmanager,
                                    chainIDmanager=self.chainIDmanager).activate_biological_assembly(self.source_specs['biological_assembly'])
        if self.source_specs.get('id',{}):
            key=self.source_specs['id']
        elif self.source_specs.get('prebuilt',{}):
            key=f'{self.source_specs["prebuilt"]["psf"]}-{self.source_specs["prebuilt"]["pdb"]}'
            xsc=self.source_specs.get('xsc','')
            if xsc:
                self.update_statevars('xsc',xsc)
        elif self.source_specs.get('alphafold',{}):
            key=f'{self.source_specs["alphafold"]}'
        else:
            raise Exception(f'The "source" directive of "psfgen" must have "id" , "prebuilt", or "alphafold"')
        self.molecules[key]=self.base_molecule
        for molid,molecule in self.molecules.items():
            logger.debug(f'Molecule "{molid}": {molecule.num_atoms()} atoms in {molecule.num_residues()} residues; {molecule.num_segments()} segments.')

    def update_molecule(self):
        """Updates all segments of the base molecule based on the 
           current coordinate file.  All ssbonds and links are 
           carried forward. No biological assembly beyond the apparent
           asymmetric unit is assumed. This should permit generation
           of a new psfgen script based on this new molecule to 
           recreate it or modify it.
        """
        # get the key of the base_molecule
        logger.debug(f'{self.taskname} has {len(self.molecules)} entries in its molecules dict')
        base_key='base'
        for k,v in self.molecules:
            if v==self.base_molecule:
                base_key=k
        # assert base_key!='UNSET',f'Cannot update a non-existent base molecule'
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        xsc=self.statevars.get('xsc','')
        source={
            'prebuilt': {
                'psf':psf,
                'pdb':pdb,
                'xsc':xsc
            }
        }
        if hasattr(self,'chainIDManager') and hasattr(self,'modmanager'):
            updated_molecule=Molecule(source=source,chainIDmanager=self.chainIDmanager,modmanager=self.modmanager).activate_biological_assembly(0)
        else:
            updated_molecule=Molecule(source=source).activate_biological_assembly(0)

        self.molecules[base_key]=updated_molecule
        self.base_molecule=updated_molecule

class CleaveTask(PsfgenTask):
    yaml_header='cleave'
    def do(self):
        self.log_message('initiated')
        self.inherit_state()
        cleavage_sites=CleavageSiteList([CleavageSite(x) for x in self.specs['sites']])
        # update base molecule to the point an inferential psfgen call could reproduce it, up to ssbonds and links
        self.base_molecule=self.statevars['base_molecule']
        self.update_molecule()
        self.base_molecule.cleave_chains(cleavage_sites)
        self.result=self.psfgen()
        # self.save_state(exts=['psf','pdb']) # already done in psfgen()
        self.log_message('complete')
        return self.result

class DomainSwapTask(MDTask):
    yaml_header='domainswap'

    def do(self):
        self.log_message('initiated')
        self.inherit_state()
        logger.debug(f'Generating inputs for domain swap')
        self.make_inputs()
        logger.debug(f'Running NAMD to execute domain swap')
        self.result=self.namdrun(baselabel='domainswap-run',extras={'colvars':'on','colvarsconfig':self.statevars['cv']},single_gpu_only=True)
        if self.result!=0:
            return self.result
        self.save_state(exts=['vel','coor'])
        self.log_message('complete')
        return self.result

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

    def do(self):
        self.log_message('initiated')
        self.inherit_state()
        self.base_molecule=self.statevars['base_molecule']
        if not self.base_molecule.has_loops(min_loop_length=self.statevars['min_loop_length']):
            self.log_message('bypassed')
            return
        logger.debug('Storing sequence gaps.')
        self.write_gaps()
        logger.debug('Measuring gap distances.')
        self.measure_distances(self.specs['steer'])
        logger.debug('Steering loop C-termini toward their partner N-termini')
        self.result=self.do_steered_md(self.specs['steer'])
        if self.result!=0:
            return self.result
        self.save_state(exts=['coor','vel'])
        logger.debug('Connecting loop C-termini to their partner N-termini')
        self.result=self.connect()
        if self.result!=0:
            return self.result
        self.save_state(exts=['psf','pdb'])
        self.log_message('complete')
        return self.result
    
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
        comment_chars='#!$'
        self.next_basename('measure')
        vm=self.writers['vmd']
        vm.newscript(self.basename)
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        datafile=self.statevars['data']
        opdb=f'{self.basename}.pdb'
        receiver_flexible_zone_radius=specs.get('receiver_flexible_zone_radius',0.0)
        resultsfile=f'{self.basename}.dat'
        vm.addline(f'measure_bonds {psf} {pdb} {datafile} {opdb} {resultsfile} {receiver_flexible_zone_radius} ')
        vm.writescript()
        vm.runscript()
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
            g['colvars']=f'GAP{i:02d}'
            declare_distance_cv_atoms(g,writer)
        for i,g in enumerate(self.gaps):
            g['forceConstant']=specs['force_constant']
            g['targ_distance']=specs['target_distance']
            g['targ_numsteps']=specs['nsteps']
            declare_single_harmonic_distance_bias(g,writer)
        writer.writefile()
        savespecs=self.specs
        self.specs=specs
        result=self.namdrun(extras={        
            'fixedatoms':'on',
            'fixedatomsfile':self.statevars['fixedref'],
            'fixedatomscol': 'O',
            'colvars': 'on',
            'colvarsconfig': f'{self.basename}-cv.inp'
            },single_gpu_only=True)
        self.specs=savespecs
        return result

    def connect(self):
        self.write_connect_patches()
        result=self.connect_gaps()
        return result

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
        # pg.topo_aliases()
        topfile=os.path.join(self.config.charmmff_custom_path,'mylink.top')
        pg.addline(f'topology {topfile}')
        # pg.usescript('loop_closure')
        patchfile=self.statevars['data']
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        pg.load_project(psf,pdb)
        pg.addline(f'source {patchfile}')
        pg.writescript(self.basename,guesscoord=True,regenerate=True)
        result=pg.runscript()
        if result==0:
            self.save_state(exts=['psf','pdb'])
        return result
    
class SolvateTask(BaseTask):
    yaml_header='solvate'

    def do(self):
        self.log_message('initiated')
        self.inherit_state()
        self.next_basename()
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        xsc=self.statevars.get('xsc',None)
        use_minmax=True
        if xsc is not None:
            box,origin=cell_from_xsc(xsc)
            if box is not None and origin is not None:
                use_minmax=False
                basisvec=box.diagonal()
                LL=origin-0.5*basisvec
                UR=origin+0.5*basisvec
        if use_minmax:
            p=PDBParser(PDBcode=os.path.splitext(pdb)[0]).parse()
            x=np.array([a.x for a in p.parsed['ATOM']])
            y=np.array([a.y for a in p.parsed['ATOM']])
            z=np.array([a.z for a in p.parsed['ATOM']])
            minmax=np.array([[x.min(),y.min(),z.min()],[x.max(),y.max(),z.max()]])
            spans=minmax[1]-minmax[0]
            maxspan=spans.max()
            cubic=self.specs.get('cubic',False)
            sympad=np.array([0.,0.,0.])
            if cubic:
                sympad=0.5*(maxspan-spans)
            pad=self.specs.get('pad',10.0)
            LL=minmax[0]-pad-sympad
            UR=minmax[1]+pad+sympad
            basisvec=UR-LL
            origin=0.5*(UR+LL)
            cell_to_xsc(np.diag(basisvec),origin,f'{self.basename}.xsc')

        ll_tcl=r'{ '+' '.join([str(_) for _ in LL.tolist()])+r' }'
        ur_tcl=r'{ '+' '.join([str(_) for _ in UR.tolist()])+r' }'
        box_tcl=r'{ '+ll_tcl+' '+ur_tcl+r' }'
        vt=self.writers['vmd']
        vt.newscript(self.basename)
        vt.addline( 'package require solvate')
        vt.addline( 'package require autoionize')
        vt.addline( 'psfcontext mixedcase')
        vt.addline(f'mol new {psf}')
        vt.addline(f'mol addfile {pdb} waitfor all')
        vt.addline(f'solvate {psf} {pdb} -minmax {box_tcl} -o {self.basename}_solv')
        ai_args=[]
        sc=self.specs.get('salt_con',None)
        if sc is None:
            ai_args.append('-neutralize')
        else:
            ai_args.append(f'-sc {sc}')
        cation=self.specs.get('cation',None)
        if cation is not None:
            ai_args.append(f'-cation {cation}')
        anion=self.specs.get('anion',None)
        if anion is not None:
            ai_args.append(f'-cation {anion}')
        vt.addline(f'autoionize -psf {self.basename}_solv.psf -pdb {self.basename}_solv.pdb {" ".join(ai_args)} -o {self.basename}')
        vt.writescript()
        self.result=vt.runscript()
        if self.result!=0:
            return super().do()
        self.save_state(exts=['psf','pdb','xsc'])
        self.log_message('complete')
        return super().do()

class DesolvateTask(BaseTask):
    yaml_header='desolvate'
    def do(self):
        self.catdcd=self.config.shell_commands['catdcd']
        self.do_idx_psf_gen()
        self.do_dcd_prune()

    def do_idx_psf_gen(self):
        self.log_message('initiated')
        # self.inherit_state()
        self.next_basename()
        psf=self.specs['psf']
        keepatselstr=self.specs['keepatselstr']
        idx_outfile=self.specs['idx_outfile']
        psf_outfile=self.specs['psf_outfile']
        vt=self.writers['vmd']
        vt.newscript(self.basename)
        vt.addline( 'package require psfgen')
        vt.addline(f'mol new {psf}')
        vt.addline(f'set keepsel [atomselect top "{keepatselstr}"]')
        vt.addline( 'set keepsegid [lsort -unique [$keepsel get segid]]')
        vt.addline(f'set dumpsel [atomselect top "not ({keepatselstr})"]')
        vt.addline(f'set dumpsegid [lsort -unique [$dumpsel get segid]]')
        vt.addline(f'vmdcon -info "Writing {idx_outfile}"')
        vt.addline(f'set fp [open "{idx_outfile}" "w"]')
        vt.addline( 'puts $fp "[$keepsel get index]"')
        vt.addline( 'close $fp')
        vt.addline( 'vmdcon -info "Keeping segids $keepsegid"')
        vt.addline( 'vmdcon -info "Dumping segids $dumpsegid"')
        vt.addline(f'readpsf {psf}')
        vt.addline(r'if {[IsIntersectionEmpty $keepsegid $dumpsegid]} {')
        vt.addline(r'   foreach badsegid $dumpsegid {')
        vt.addline( '       delatom $badsegid')
        vt.addline(r'   }')
        vt.addline(r'}')
        vt.addline(f'vmdcon -info "Writing {psf_outfile}"')
        vt.addline(f'writepsf {psf_outfile}')
        vt.writescript()
        self.result=vt.runscript(progress_title='psfidx')
        
    def do_dcd_prune(self):
        idx_outfile=self.specs['idx_outfile']
        dcd_outfile=self.specs['dcd_outfile']
        dcd_infiles=self.specs['dcd_infiles']
        dcd_stride=self.specs['dcd_stride']
        progress_struct=PestiferProgress(name='catdcd',track_stdout=False)
        c=Command(f'{self.catdcd} -i {idx_outfile} -stride {dcd_stride} -o {dcd_outfile} {" ".join(dcd_infiles)}')
        c.run(progress=progress_struct)

class ManipulateTask(BaseTask):
    yaml_header='manipulate'
    def do(self):
        self.log_message('initiated')
        if self.prior:
            logger.debug(f'Task {self.taskname} prior {self.prior.taskname}')
            self.inherit_state()
        logger.debug(f'manipulate {self.specs["mods"]}')
        self.modmanager=ModManager(self.specs['mods'])
        self.result=self.coormods()
        self.log_message('complete')
        return super().do()

    def coormods(self):
        coormods=self.modmanager.get('coormods',{})
        if coormods:
            logger.debug(f'performing coormods')
            for modtype,modlist in coormods.items():
                self.next_basename(modtype)
                vm=self.writers['vmd']
                vm.newscript(self.basename,packages=['Orient'])
                psf=self.statevars['psf']
                pdb=self.statevars['pdb']
                vm.load_psf_pdb(psf,pdb,new_molid_varname='mCM')
                modlist.write_TcL(vm)
                vm.write_pdb(self.basename,'mCM')
                vm.writescript()
                result=vm.runscript()
                if result!=0:
                    return result
                self.save_state(exts=['pdb'])
        return 0

class TerminateTask(MDTask):  #need to inherit for namdrun() method
    yaml_header='terminate'
    
    def do(self):
        self.log_message('initiated')
        self.inherit_state()
        self.next_basename()
        self.copy_state(exts=['psf','pdb','coor','xsc','vel'])
        self.write_chainmaps()
        self.write_statefile()
        self.result=self.make_package()
        self.log_message('complete')
        return self.result

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
            return 0
        self.inherit_state()
        self.FC.clear()  # populate a file collector to make the tarball
        logger.debug(f'Packaging for namd using basename {self.basename}')
        savespecs=self.specs
        self.specs=specs
        params={}
        result=self.namdrun(script_only=True)
        self.specs=savespecs
        self.FC.append(f'{self.basename}.namd')
        constraints=specs.get('constraints',{})
        if constraints:
            self.make_constraint_pdb(constraints)
            self.FC.append(self.statevars['consref'])
        local_params=self.local_parameter_files
        for n in local_params:
            self.FC.append(n)
        for ext in ['psf','pdb','coor','xsc','vel']:
            if ext in self.statevars:
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
        return 0

class RingCheckTask(BaseTask):
    """ A class for checking for pierced rings
    
    Attributes
    ----------
    yaml_header(str) 

    Methods
    -------
    do(): 
        Based on specs, executes ring piercing check and deletes offending residues

    """
    yaml_header='ring_check'

    def do(self):
        self.log_message('initiated')
        self.inherit_state()
        psf=self.statevars.get('psf',None)
        pdb=self.statevars.get('pdb',None)
        xsc=self.statevars.get('xsc',None)
        cutoff=self.specs.get('cutoff',3.5)
        segtypes=self.specs.get('segtypes',['lipid'])
        delete_these=self.specs.get('delete','piercee')
        npiercings=ring_check(psf,pdb,xsc,cutoff=cutoff,segtypes=segtypes)
        if npiercings:
            ess='s' if len(npiercings)>1 else ''
            if delete_these=="none":
                logger.debug(f'No action taken regarding {len(npiercings)} pierced-ring configuration{ess}')
                for r in npiercings:
                    logger.debug(f'  Piercing of {r["piercee"]["segname"]}-{r["piercee"]["resid"]} by {r["piercer"]["segname"]}-{r["piercer"]["resid"]}')
            else:
                self.next_basename('ring_check')
                pg=self.writers['psfgen']
                pg.newscript(self.basename)
                pg.load_project(psf,pdb)
                logger.debug(f'Deleting all {delete_these}s from {len(npiercings)} pierced-ring configuration{ess}')
                for r in npiercings:
                    logger.debug(f'   Deleting segname {r[delete_these]["segname"]} residue {r[delete_these]["resid"]}')
                    pg.addline(f'delatom {r[delete_these]["segname"]} {r[delete_these]["resid"]}')
                pg.writescript(self.basename)
                pg.runscript()
                self.save_state(exts=['psf','pdb'])
        self.log_message('complete')
        self.result=0
        return super().do()




