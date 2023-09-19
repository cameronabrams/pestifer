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
from argparse import Namespace
import yaml
from .util import inspect_classes, is_periodic
from .molecule import Molecule
from .chainids import ChainIDManager
from .colvars import *
from .stringthings import FileCollector
from .modcontainer import ModContainer

class Task(BaseMod):
    req_attr=BaseMod.req_attr+['specs','config','index','prior','writers','taskname']
    yaml_header='generic_task'
    exts=['.psf','.pdb','.coor','.xsc'] # extensions of files that can be transferred from one task to the next
    _taskcount=0
    def __init__(self,input_dict,taskname,config,writers,prior):
        specs=input_dict.copy()
        input_dict = {
            'index':Task._taskcount,
            'writers': writers,
            'prior':prior,
            'specs':specs,
            'config':config,
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
        namd_params=self.config['user']['namd2']
        basename=self.next_basename('minimize')
        nminsteps=specs['nminsteps']
        dcdfreq=specs['dcdfreq']
        na=self.writers['namd2']
        temperature=namd_params['generic']['temperature']
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
        params['parameters']=na.standard_charmmff_parfiles+na.custom_charmmff_parfiles
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

    def namd2prep(self,basename,specs,absolute_paths=True):
        namd_params=self.config['user']['namd2']
        ensemble=specs['ensemble']
        temperature=specs['temperature']
        pressure=specs['pressure']

        nminsteps=specs['nminsteps']
        nsteps=specs['nsteps']
        dcdfreq=specs['dcdfreq']
        xstfreq=specs['xstfreq']
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
        if absolute_paths:
            params['parameters']=na.standard_charmmff_parfiles+na.custom_charmmff_parfiles
        else:
            params['parameters']=[]
            namd2_params_abs=na.standard_charmmff_parfiles+na.custom_charmmff_parfiles
            for nf in namd2_params_abs:
                d,n=os.path.split(nf)
                params['parameters'].append(n)
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
        if ensemble=='NPT':
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
        return params

    def namd2script(self,basename,params):
        na=self.writers['namd2']
        na.newscript(basename)
        na.writescript(params)

    def relax(self,basename):
        na=self.writers['namd2']
        na.runscript()
        self.update_statefile('coor',f'{basename}.coor')
        self.update_statefile('vel',f'{basename}.vel')
        self.update_statefile('xsc',f'{basename}.xsc',mode='permissive')
        self.coor_to_pdb(basename)

class PsfgenTask(Task):
    yaml_header='psfgen'
    def __init__(self,input_dict,taskname,config,writers,prior):
        super().__init__(input_dict,taskname,config,writers,prior)
        self.chainIDmanager=ChainIDManager(format=self.specs['source']['file_format'])

    def do(self):
        logger.info(f'Task {self.taskname} {self.index:02d} initiated')
        if self.prior:
            logger.debug(f'... prior {self.prior.taskname}')
            self.statevars=self.prior.statevars.copy()
        logger.debug('Injesting molecule(s)')
        self.injest_molecules()
        self.statevars['base_molecule']=self.base_molecule
        logger.debug('Running first psfgen')
        self.psfgen()
        nloops=self.base_molecule.has_loops(min_loop_length=self.statevars['min_loop_length'])*self.base_molecule.num_images()
        if nloops>0:
            logger.debug(f'Declashing {nloops} loops')
            self.declash_loops(self.specs['source']['sequence']['loops'])
        logger.debug('Minimizing')
        self.minimize(self.specs['minimize'])
        logger.info(f'Task {self.taskname} {self.index:02d} complete')

    def psfgen(self):
        basename=self.next_basename('build')
        pg=self.writers['psfgen']
        pg.newscript(basename)
        pg.topo_aliases()
        pg.set_molecule(self.base_molecule)
        pg.describe_molecule(self.base_molecule)
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
        cycles=specs['declash']['maxcycles']
        mol.write_loop_lines(vt,cycles=cycles,min_length=specs['min_loop_length'])
        vt.write_pdb(basename,'mLL')
        vt.endscript()
        vt.writescript()
        vt.runscript()
        self.update_statefile('pdb',f'{basename}.pdb')

    def injest_molecules(self):
        specs=self.specs
        self.molecules={}
        self.source_specs=specs['source']
        logger.debug(f'user-input modspecs {self.specs["mods"]}')
        self.mods=ModContainer(self.specs['mods'])
        # self.pdbs=self.mods.report_pdbs()
        # self.usermod_specs=specs['mods']
        logger.debug(f'user mods at injest_molecules {self.mods.__dict__}')
        self.molecules[self.source_specs['id']]=Molecule(source=self.source_specs,usermods=self.mods,chainIDmanager=self.chainIDmanager).activate_biological_assembly(self.source_specs['biological_assembly'])
        self.base_molecule=self.molecules[self.source_specs['id']]
        self.statevars['min_loop_length']=self.source_specs['sequence']['loops']['min_loop_length']
        # for p in self.pdbs:
        #     self.molecules[p]=Molecule(source=p)

class LigateTask(Task):
    yaml_header='ligate'
    statevars={}
    def __init__(self,input_dict,taskname,config,writers,prior):
        super().__init__(input_dict,taskname,config,writers,prior)

    def do(self):
        logger.info(f'Task {self.taskname} {self.index:02d} initiated')
        if self.prior:
            logger.debug(f'Task {self.taskname} prior {self.prior.taskname}')
            self.statevars=self.prior.statevars.copy()
        self.base_molecule=self.statevars['base_molecule']
        if not self.base_molecule.has_loops(min_loop_length=self.statevars['min_loop_length']):
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
        self.measure_distances(specs)
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

    def measure_distances(self,specs):
        receiver_flexible_zone_radius=specs['receiver_flexible_zone_radius']
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
        vm.runscript(psf=psf,pdb=pdb,i=datafile,opdb=f'{basename}.pdb',o=resultsfile,rfzr=receiver_flexible_zone_radius)
        self.update_statefile('fixedref',f'{basename}.pdb')
        self.update_statefile('results',resultsfile)
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
        basename=self.next_basename('steer')
        nsteps=specs['nsteps']
        dcdfreq=specs['dcdfreq']
        temperature=specs['temperature']
        writer=self.writers['data']
        writer.newfile(f'{basename}-cv.inp')
        for i,g in enumerate(self.gaps):
            g['name']=f'GAP{i:02d}'
            declare_distance_cv_atoms(g,writer)
        for i,g in enumerate(self.gaps):
            g['k']=specs['force_constant']
            g['targ_distance']=specs['target_distance']
            g['targ_numsteps']=nsteps
            declare_harmonic_distance_bias(g,writer)
        writer.writefile()
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        na=self.writers['namd2']
        params={'structure':psf,'coordinates':pdb}
        params.update({'tcl':[f'set temperature {temperature}']})
        params['temperature']='$temperature'
        params['parameters']=na.standard_charmmff_parfiles+na.custom_charmmff_parfiles
        logger.debug(f'Parameter files: {params["parameters"]}')
        namd_params=self.config['user']['namd2']
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
        topfile=os.path.join(self.config.charmmff_custom_path,'mylink.top')
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
        vt.runscript(o=basename,pdb=pdb,psf=psf,pad=self.specs['pad'])
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
        basename=self.next_basename('relax')
        self.namd2script(basename,self.namd2prep(basename,self.specs))
        self.relax(basename)
        logger.info(f'Task {self.taskname} {self.index:02d} complete')

class TerminateTask(Task):
    yaml_header='terminate'
    opt_attr=Task.opt_attr+[yaml_header]
    def do(self):
        logger.info(f'Task {self.taskname} {self.index:02d} initiated')
        self.statevars=self.prior.statevars.copy()
        self.writeresults()
        if 'package' in self.specs:
            self.package(self.specs['package'])

    def writeresults(self):
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

    def package(self,specs):
        if not specs:
            return
        self.statevars=self.prior.statevars.copy()
        self.FC.clear()
        basename=specs["basename"]
        params=self.namd2prep(basename,specs)
        local_params=[]
        for nf in params["parameters"]:
            d,n=os.path.split(nf)
            shutil.copy(nf,n)
            local_params.append(n)
            self.FC.append(n)
        for ext in self.exts+['.vel']:
            aext=ext[1:]
            self.FC.append(self.statevars[aext])
        if os.path.exists(self.statevars['vel']):
            del params['temperature']
        params['parameters']=local_params
        self.namd2script(basename,params)
        self.FC.append(f'{basename}.namd')
        self.FC.tarball(specs["basename"])

        