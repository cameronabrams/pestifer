"""

.. module:: tasks
   :synopsis: defines all tasks
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
from .basemod import BaseMod
import logging
logger=logging.getLogger(__name__)
import shutil
import os
import yaml
import numpy as np
from .util import inspect_classes, is_periodic
from .molecule import Molecule
from .chainids import ChainIDManager
from .colvars import *
from .stringthings import FileCollector

class Task(BaseMod):
    req_attr=BaseMod.req_attr+['specs','index','prior','writers','config','taskname']
    yaml_header='generic_task'
    default_specs={}
    exts=['.psf','.pdb','.coor','.xsc']
    _taskcount=0
    def __init__(self,input_dict,taskname,config,writers,prior):
        specs=input_dict.copy()
        for k,v in self.__class__.default_specs.items():
            if not k in specs:
                specs[k]=v
        input_dict = {
            'index':Task._taskcount,
            'config':config,
            'writers': writers,
            'prior':prior,
            'specs':specs,
            'taskname':taskname
        }
        for k,v in specs.items():
            if not v or (type(v)==str and v.lower()=='none'):
                specs[k]={}
        super().__init__(input_dict)
        Task._taskcount+=1
        self.subtaskcount=0
        self.statevars={}
        self.FC=FileCollector()

    def update_statefile(self,key,filename,mode='strict'):
        if os.path.exists(filename):
            self.statevars[key]=filename
        else:
            if mode=='strict':
                raise FileNotFoundError(f'{filename} not found')

    def next_basename(self,*obj):
        label=''
        if len(obj)==1 and len(obj[0])>0:
            label=f'-{obj[0]}'
        basename=f'{self.index:02d}-{self.subtaskcount:02d}-{self.taskname}{label}'
        self.subtaskcount+=1
        return basename
            
    def coor_to_pdb(self,basename):
        vm=self.writers['vmd']
        vm.newscript(f'{basename}-coor2pdb')
        vm.usescript('namdbin2pdb')
        vm.writescript()
        psf=self.statevars['psf']
        coor=self.statevars['coor']
        vm.runscript(psf=psf,coor=coor,pdb=f'{basename}.pdb')
        self.update_statefile('pdb',f'{basename}.pdb')

    def pdb_to_coor(self,basename):
        vm=self.writers['vmd']
        vm.newscript(f'{basename}-pdb2coor')
        vm.usescript('pdb2namdbin')
        vm.writescript()
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        vm.runscript(psf=psf,pdb=pdb,coor=f'{basename}.coor')
        self.update_statefile('coor',f'{basename}.coor')

    def minimize(self,specs):
        namd_params=self.config.namd_params
        temperature=specs.get('temperature',None)
        basename=self.next_basename('minimize')
        nminsteps=specs.get('nminsteps',1000)
        dcdfreq=specs.get('dcdfreq',0)
        na=self.writers['namd2']
        if not temperature:
            temperature=namd_params['generic'].get('temperature',311)
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        coor=self.statevars.get('coor',None)
        vel=self.statevars.get('vel',None)
        xsc=self.statevars.get('xsc',None)
        cell=self.statevars.get('cell',None)
        self.statevars['periodic']=is_periodic(cell,xsc)
        params={}
        params['tcl']=[f'set temperature {temperature}']
        params['structure']=psf
        params['coordinates']=pdb
        if coor:
            params['bincoordinates']=coor
        if vel:
            params['binvelocities']=vel
        else:
            params['temperature']='$temperature'
        if xsc:
            params['extendedSystem']=xsc
        params['parameters']=na.standard_charmmparfiles+na.custom_charmmparfiles
        params.update(namd_params['generic'])

        if self.statevars['periodic']:
            params.update(namd_params['solvated'])
            if cell:
                params['tcl'].append(f'source {self.statevars["cell"]}')
        else:
            params.update(namd_params['vacuum'])
        params.update(namd_params['thermostat'])           
        params['outputName']=f'{basename}'
        if dcdfreq:
            params['dcdfreq']=dcdfreq
        params['firsttimestep']=0
        if nminsteps:
            params['minimize']=nminsteps
        na.newscript(basename)
        na.writescript(params)
        na.runscript()
        self.update_statefile('coor',f'{basename}.coor')
        self.update_statefile('xsc',f'{basename}.xsc',mode='permissive')
        if cell:
            del self.statevars['cell']
        self.coor_to_pdb(basename)

    def relax(self,specs,label=None):
        namd_params=self.config.namd_params
        temperature=specs.get('temperature',None)
        pressure=specs.get('pressure',None)
        if not temperature:
            temperature=namd_params['generic'].get('temperature',311)
        if label==None:
            basename=self.next_basename('relax')
        else:
            basename=self.next_basename(label)
        nminsteps=specs.get('nminsteps',0)
        nsteps=specs.get('nsteps',1000)
        dcdfreq=specs.get('dcdfreq',0)
        xstfreq=specs.get('xstfreq',0)
        na=self.writers['namd2']
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        params={}
        params['structure']=psf
        params['coordinates']=pdb
        coor=self.statevars.get('coor',None)
        vel=self.statevars.get('vel',None)
        xsc=self.statevars.get('xsc',None)
        self.statevars['periodic']=is_periodic(None,xsc)
        if coor:
            params['bincoordinates']=coor
        if vel:
            params['binvelocities']=vel
        params.update({'tcl':[f'set temperature {temperature}']})
        params['temperature']='$temperature'
        params['parameters']=na.standard_charmmparfiles+na.custom_charmmparfiles
        params.update(namd_params['generic'])
        if xsc:
            params['extendedSystem']=xsc
            if xstfreq:
                params['xstfreq']=xstfreq
        if self.statevars['periodic']:
            params.update(namd_params['solvated'])
        else:
            params.update(namd_params['vacuum'])
        params.update(namd_params['thermostat'])
        if pressure:
            if not self.statevars['periodic']:
                raise Exception(f'Cannot use barostat on a system without PBCs')
            params['tcl'].append(f'set pressure {pressure}')
            params.update(namd_params['barostat'])
        params['outputName']=f'{basename}'
        if dcdfreq:
            params['dcdfreq']=dcdfreq
        params['firsttimestep']=0
        if nminsteps:
            params['minimize']=nminsteps
        params['run']=nsteps
        na.newscript(basename)
        na.writescript(params)
        na.runscript()
        self.update_statefile('coor',f'{basename}.coor')
        self.update_statefile('vel',f'{basename}.vel')
        self.update_statefile('xsc',f'{basename}.xsc',mode='permissive')
        self.coor_to_pdb(basename)

class PsfgenTask(Task):
    yaml_header='psfgen'
    default_specs={'cleanup':False,'mods':{},'layloops':{},'minimize':{'nsteps':1000}}
    def __init__(self,input_dict,taskname,config,writers,prior):
        super().__init__(input_dict,taskname,config,writers,prior)
        self.chainIDmanager=ChainIDManager(format=config['rcsb_file_format'])
        self.mods=self.specs['mods']

    def do(self):
        logger.info(f'Task {self.taskname} {self.index:02d} initiated')
        if self.prior:
            logger.debug(f'... prior {self.prior.taskname}')
            self.statevars=self.prior.statevars.copy()
        logger.debug('Parsing modifications')
        self.modparse()
        logger.debug('Injesting molecule(s)')
        self.injest_molecules(self.specs['source'])
        self.statevars['base_molecule']=self.base_molecule
        logger.debug('Running first psfgen')
        self.psfgen()
        if self.base_molecule.has_loops():
            logger.debug('Declashing loops')
            self.declash_loops(self.specs.get('declash',{}))
        logger.debug('Minimizing')
        self.minimize(self.specs['minimize'])
        logger.info(f'Task {self.taskname} {self.index:02d} complete')

    def psfgen(self):
        basename=self.next_basename('build')
        pg=self.writers['psfgen']
        pg.newscript(basename)
        pg.topo_aliases()
        pg.set_molecule(self.base_molecule)
        pg.describe_molecule(self.base_molecule,self.mods)
        pg.complete(basename)
        pg.endscript()
        pg.writescript()
        pg.runscript()
        self.update_statefile('pdb',f'{basename}.pdb')
        self.update_statefile('psf',f'{basename}.psf')
        pg.cleanup(cleanup=self.specs['cleanup'])

    def declash_loops(self,specs):
        basename=self.next_basename('declash')
        mol=self.base_molecule
        if not mol.has_loops():
            return
        vt=self.writers['vmd']
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        vt.newscript(basename)
        vt.load_psf_pdb(psf,pdb,new_molid_varname='mLL')
        cycles=specs.get('maxcycles',100)
        mol.write_loop_lines(vt,cycles=cycles,min_length=specs.get('min_loop_length',4))
        vt.write_pdb(basename,'mLL')
        vt.endscript()
        vt.writescript()
        vt.runscript()
        self.update_statefile('pdb',f'{basename}.pdb')

    def modparse(self):
        mod_classes,modlist_classes=inspect_classes('pestifer.mods','List')
        retdict={}
        input_dict=self.specs.get('mods',{})
        for hdr,entries in input_dict.items():
            class_name=[name for name,cls in mod_classes.items() if cls.yaml_header==hdr][0]
            cls=mod_classes[class_name]
            LCls=modlist_classes.get(f'{class_name}List',list)
            for entry in entries:
                assert type(entry) in [dict,str]
                newmod=cls(entry)
                newmod.source='USER'
                if not hdr in retdict:
                    retdict[hdr]=LCls([])
                retdict[hdr].append(newmod)
        self.specs['mods']=retdict
        # TODO: gather names of all aux pdb files for
        self.pdbs=[]

    def injest_molecules(self,specs):
        self.molecules={}
        psf_exists=False
        self.basename=specs.get('rcsb',None)
        bioassemb=specs.get('biological_assembly',0)
        excludes=specs.get('exclude',{})
        self.molecules[self.basename]=Molecule(config=self.config,source=self.basename,chainIDmanager=self.chainIDmanager,excludes=excludes,use_psf=psf_exists).activate_biological_assembly(bioassemb)
        self.base_molecule=self.molecules[self.basename]
        for p in self.pdbs:
            self.molecules[p]=Molecule(config=self.config,source=p)

class LigateTask(Task):
    yaml_header='ligate'
    default_specs={'steer':{},'connect':{},'minimize':{}}
    statevars={}
    def __init__(self,input_dict,taskname,config,writers,prior):
        super().__init__(input_dict,taskname,config,writers,prior)

    def do(self):
        logger.info(f'Task {self.taskname} {self.index:02d} initiated')
        if self.prior:
            logger.debug(f'Task {self.taskname} prior {self.prior.taskname}')
            self.statevars=self.prior.statevars.copy()
        self.base_molecule=self.statevars['base_molecule']
        if not self.base_molecule.has_loops(min_length=self.specs.get('min_loop_length',4)):
            logger.info('No loops. Ligation bypassed.')
            return
        logger.debug('Steering loop ends')
        self.steerends(self.specs['steer'])
        logger.debug('Connecting loops')
        self.connect(self.specs['connect'])
        logger.debug('Minimizing')
        self.minimize(self.specs['minimize'])
        logger.info(f'Task {self.taskname} {self.index:02d} complete')

    def steerends(self,specs):
        logger.debug('...Writing gaps')
        self.write_gaps()
        logger.debug('...Measuring gap distances')
        self.measure_distances()
        logger.debug('...Doing steered MD')
        self.do_steered_md(specs)
    
    def write_gaps(self):
        basename=self.next_basename('gaps')
        mol=self.base_molecule
        datafile=f'{basename}.inp'
        writer=self.writers['data']
        writer.newfile(datafile)
        mol.write_gaps(writer)
        writer.writefile()
        self.update_statefile('data',datafile)

    def measure_distances(self):
        comment_chars='#!$'
        basename=self.next_basename('measure')
        resultsfile=f'{basename}.dat'
        vm=self.writers['vmd']
        vm.newscript(basename)
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        vm.usescript('measure_bonds')
        vm.writescript()
        datafile=self.statevars['data']
        vm.runscript(psf=psf,pdb=pdb,i=datafile,opdb=f'{basename}.pdb',o=resultsfile)
        self.update_statefile('fixedref',f'{basename}.pdb')
        self.update_statefile('results',resultsfile)
        with open(resultsfile,'r') as f:
            datalines=f.read().split('\n')
        gaps=[]
        for line in datalines:
            if len(line)>0 and not line[0] in comment_chars:
                data=line.split()
                thisgap={
                    'segname':data[0],
                    'serial_i':int(data[1]),
                    'serial_j':int(data[2]),
                    'distance':float(data[3])
                }
                gaps.append(thisgap)
        self.specs['gaps']=gaps

    def do_steered_md(self,specs):
        basename=self.next_basename('steer')
        nsteps=specs.get('nsteps',1000)
        dcdfreq=specs.get('dcdfreq',0)
        temperature=specs.get('temperature',310)
        writer=self.writers['data']
        writer.newfile(f'{basename}-cv.inp')
        for i,g in enumerate(self.specs['gaps']):
            g['name']=f'GAP{i:02d}'
            declare_distance_cv_atoms(g,writer)
        for i,g in enumerate(self.specs['gaps']):
            g['k']=specs.get('force_constant',20.0)
            g['targ_distance']=specs.get('target_distance',2.0)
            g['targ_numsteps']=nsteps
            declare_harmonic_distance_bias(g,writer)
        writer.writefile()
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        na=self.writers['namd2']
        params={'structure':psf,'coordinates':pdb}
        params.update({'tcl':[f'set temperature {temperature}']})
        params['temperature']='$temperature'
        params['parameters']=na.standard_charmmparfiles+na.custom_charmmparfiles
        namd_params=self.config.namd_params
        params.update(namd_params['generic'])
        params.update(namd_params['vacuum'])
        params.update(namd_params['thermostat'])

        params['outputName']=f'{basename}'
        if dcdfreq:
            params['dcdfreq']=dcdfreq
        params['firsttimestep']=0
        extras={
            'fixedatoms':'on',
            'fixedatomsfile':self.statevars['fixedref'],
            'fixedatomscol': 'O',
            'colvars': 'on',
            'colvarsconfig': f'{basename}-cv.inp'
        }
        params.update(extras)
        params['run']=nsteps
        na.newscript(basename)
        na.writescript(params)
        na.runscript()
        self.update_statefile('coor',f'{basename}.coor')
        self.update_statefile('vel',f'{basename}.vel')
        self.coor_to_pdb(basename)

    def connect(self,specs):
        self.write_connect_patches(specs)
        self.connect_gaps(specs)

    def write_connect_patches(self,specs):
        basename=self.next_basename('gap_patches')
        mol=self.base_molecule
        datafile=f'{basename}.inp'
        writer=self.writers['data']
        writer.newfile(datafile)
        mol.write_connect_patches(writer)
        writer.writefile()
        self.update_statefile('data',datafile)

    def connect_gaps(self,specs):
        basename=self.next_basename('heal')
        pg=self.writers['psfgen']
        pg.newscript(basename)
        pg.topo_aliases()
        topfile=os.path.join(self.config.charmm_custom_path,'mylink.top')
        pg.addline(f'topology {topfile}')
        pg.usescript('loop_closure')
        pg.writescript()
        patchfile=self.statevars['data']
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        pg.runscript(psf=psf,pdb=pdb,p=patchfile,o=basename)
        self.update_statefile('psf',f'{basename}.psf')
        self.update_statefile('pdb',f'{basename}.pdb')
        if 'vel' in self.statevars:
            del self.statevars['vel']
        self.pdb_to_coor(basename)

class SolvateTask(Task):
    yaml_header='solvate'
    opt_attr=Task.opt_attr+[yaml_header]
    default_specs={'solvate':{},'autoionize':{},'minimize':{'nminsteps':100}}
    def do(self):
        self.statevars=self.prior.statevars.copy()
        basename=self.next_basename()
        logger.info(f'Task {self.taskname} {self.index:02d} initiated')
        vt=self.writers['vmd']
        vt.newscript(basename)
        vt.usescript('solv')
        vt.writescript()
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        vt.runscript(o=basename,pdb=pdb,psf=psf)
        self.update_statefile('cell',f'{basename}_cell.tcl')
        self.update_statefile('psf',f'{basename}.psf')
        self.update_statefile('pdb',f'{basename}.pdb')
        for oldext in ['coor','vel','xsc']:
            if oldext in self.statevars:
                del self.statevars[oldext]
        self.minimize(self.specs['minimize'])
        logger.info(f'Task {self.taskname} {self.index:02d} complete')

class RelaxTask(Task):
    yaml_header='relax'
    opt_attr=Task.opt_attr+[yaml_header]
    def do(self):
        logger.info(f'Task {self.taskname} {self.index:02d} initiated')
        self.statevars=self.prior.statevars.copy()
        self.relax(self.specs,label='')
        logger.info(f'Task {self.taskname} {self.index:02d} complete')

class TerminateTask(Task):
    yaml_header='terminate'
    opt_attr=Task.opt_attr+[yaml_header]
    default_specs={'basename':'final','chainmapfile':'chainmaps.yaml','statefile':'states.yaml'}
    def do(self):
        logger.info(f'Task {self.taskname} {self.index:02d} initiated')
        self.statevars=self.prior.statevars.copy()
        for ext in self.exts+['.vel']:
            aext=ext[1:]
            if aext in self.statevars:
                ffile=f'{self.specs["basename"]}{ext}'
                shutil.copy(self.statevars[aext],ffile)
                self.update_statefile(aext,ffile)
        bm=self.statevars.get('base_molecule',None)
        if bm:
            maps=bm.get_chainmaps()
            with open(self.specs['chainmapfile'],'w') as f:
                yaml.dump(maps,f)
            del self.statevars['base_molecule']
        with open(self.specs["statefile"],'w') as f:
            yaml.dump(self.statevars,f)
        logger.info(f'Task {self.taskname} {self.index:02d} complete')
