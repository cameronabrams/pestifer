from .basemod import BaseMod
import logging
logger=logging.getLogger(__name__)
import shutil
import os
from .util import inspect_classes
from .molecule import Molecule
from .chainids import ChainIDManager
from .colvars import *

class Task(BaseMod):
    req_attr=BaseMod.req_attr+['index','prior','specs','taskname','writers','config']
    yaml_header='generic_task'
    default_specs={}
    exts=['.psf','.pdb','.coor','.xst']
    _taskcount=0
    def __init__(self,input_dict,config,writers,prior):
        assert len(input_dict)==1
        taskname=list(input_dict.keys())[0]
        specval=list(input_dict.values())[0]
        specs={} if not specval else specval
        for k,v in self.__class__.default_specs.items():
            if not k in specs:
                specs[k]=v
        input_dict.update({
            'index':Task._taskcount,
            'config':config,
            'writer': writers,
            'prior':prior,
            'specs':specs,
            'taskname':taskname
        })
        super().__init__(input_dict)
        Task._taskcount+=1

    def prior_statefiles(self):
        if not self.prior: return {}
        return self.prior.statefiles
    
    def pass_thru(self):
        for fn in self.priors().values():
            if os.path.exists(fn):
                bn,ext=os.path.splitext(fn)
                shutil.copyfile(fn,f'{self.basename}{ext}')
            else:
                logger.debug(f'Task {self.index} ({self.taskname}): prior task\'s {fn} not found')
        else:
            logger.debug(f'Task {self.index} ({self.taskname}) empty passthru')
    
    def relax(self,basename):
        specs=self.specs
        temperature=specs.get('temperature',300)
        outputbasename=f'{basename}-relaxed'
        nsteps=specs.get('nsteps',1000)
        psf=self.statefiles['psf']
        pdb=self.statefiles['pdb']
        coor=self.statefile.get('coor',None)
        params={'structure':psf}
        if coor:
            params['coordinates']=coor
        else:
            params['coordinates']=pdb
        params.update({'tcl':[f'set temperature {temperature}']})
        charmmpar=self.config['StdCharmmParam']
        params['parameters']=charmmpar
        namd_params=self.config.namd_params
        params.update(namd_params['generic'])
        params.update(namd_params['vacuum'])
        params.update(namd_params['thermostat'])           
        params['outputName']=f'{outputbasename}'
        params['firsttimestep']=0
        params['run']=nsteps
        na=self.writers['namd2']
        na.newscript(basename)
        na.writescript(params)
        na.runscript()
        self.statefiles['coor']=f'{outputbasename}.coor'
        vm=self.writers['vmd']
        vm.newscript(f'{outputbasename}-coor2pdb')
        vm.usescript('namdbin2pdb')
        psf=self.statefiles['psf']
        vm.runscript(psf,f'{outputbasename}.coor',f'{outputbasename}.pdb')
        self.statefiles['pdb']=f'{outputbasename}.pdb'

class PsfgenTask(Task):
    yaml_header='psfgen'
    req_attr=Task.req_attr+['source']
    opt_attr=Task.opt_attr+['mods','ligation','cleanup']
    default_specs={'cleanup':False,'mods':None,'ligation':None}
    statefiles={}
    def __init__(self):
        self.chainIDmanager=ChainIDManager()
        super().__init__()

    def do(self):
        self.modparse()
        self.injest_molecules()
        basename=f'{self.index:02d}-{self.taskname}'
        self.psfgen(basename)
        self.ligate(f'{basename}-ligate')
        self.relax(f'{basename}-relax')

    def psfgen(self,basename):
        pg=self.writers['psfgen']
        pg.newscript(basename)
        pg.topo_aliases()
        pg.set_molecule(self.base_molecule)
        pg.describe_molecule(self.base_molecule,self.mod_dict)
        self.statefiles.update(pg.complete())
        pg.endscript()
        pg.writescript()
        pg.runscript()
        pg.cleanup(cleanup=self.specs['cleanup'])

    def ligate(self,basename):
        if not 'ligation' in self.specs:
            return
        specs=self.specs['ligation']
        self.layloops(specs,f'{basename}-lay')
        self.steerends(specs,f'{basename}-steer')
        self.connect(specs,f'{basename}-connect')

    def layloops(self,specs,basename):
        mol=self.base_molecule
        vt=self.writers['vmd']
        psf=self.statefiles['psf']
        pdb=self.statefiles['pdb']
        vt.newscript(basename)
        vt.load_psf_pdb(psf,pdb,new_molid_varname='mLL')
        cycles=specs.get('cycles',100)
        sac_n=specs.get('min_loop_length',4)
        mol.write_loop_lines(vt,cycles=cycles,sac_n=sac_n)
        self.statefiles['pdb']=f'{basename}-layed.pdb'
        vt.write_pdb(self.statefiles['pdb'],'mLL')
        vt.endscript()
        vt.writescript()
        vt.runscript()

    def steerends(self,specs,basename):
        self.write_gaps(specs,f'{basename}-gaps')
        self.measure_distances(specs,f'{basename}-measure')
        self.do_steered_md(specs,f'{basename}-smd')
    
    def write_gaps(self,basename):
        mol=self.base_molecule
        datafile=f'{basename}-gaps.inp'
        writer=self.writers['data']
        writer.newfile(datafile)
        mol.write_gaps(writer)
        writer.writefile()
        self.statefiles['datafile']=datafile

    def measure_distances(self,specs,basename):
        vm=self.writers['vmd']
        vm.newscript(basename)
        psf=self.statefiles['psf']
        pdb=self.statefiles['pdb']
        vm.usescript('measure_bonds')
        vm.writescript()
        datafile=self.statefiles['data']
        self.statefiles['fixedref']=f'{basename}-fixedref.pdb'
        vm.runscript(psf,pdb,datafile,self.statefiles['fixedref'])
        with open(datafile,'r') as f:
            datalines=f.read().split('\n')
        gaps=[]
        for line in datalines:
            data=line.split()
            thisgap={
                'chainID':data[0],
                'serial_i':int(data[1]),
                'serial_j':int(data[2]),
                'distance':float(data[3])
            }
            gaps.append(thisgap)
        self.specs['gaps']=gaps

    def do_steered_md(self,specs,basename):
        smdspecs=specs.get('smdclose',{})
        nsteps=smdspecs.get('nsteps',1000)
        temperature=smdspecs.get('temperature',310)
        writer=self.writers['data']
        writer.newfile('cv.inp')
        for i,g in enumerate(self.specs['gaps']):
            g['name']=f'GAP{i:02d}'
            declare_distance_cv_atoms(g,writer)
        for i,g in enumerate(self.specs['gaps']):
            g['k']=smdspecs.get('force_constant',20.0)
            g['targ_distance']=smdspecs.get('target_distance',2.0)
            g['targ_numsteps']=nsteps
            declare_harmonic_distance_bias(g,writer)
        writer.writefile()
        psf=self.statefiles['psf']
        pdb=self.statefiles['pdb']

        params={'structure':psf,'coordinates':pdb}
        params.update({'tcl':[f'set temperature {temperature}']})
        charmmpar=self.config['StdCharmmParam']
        params['parameters']=charmmpar
        namd_params=self.config.namd_params
        params.update(namd_params['generic'])
        params.update(namd_params['vacuum'])
        params.update(namd_params['thermostat'])
        extras={
            'fixedatoms':'on',
            'fixedatomsfile':self.statefiles['fixedref'],
            'fixedatomscol': 'O',
            'colvars': 'on',
            'colvarsconfig': 'cv.inp'
        }
        params.update(extras)
        params['outputName']=f'{basename}-steered'
        params['firsttimestep']=0
        params['run']=nsteps
        na=self.writers['namd2']
        na.newscript(basename)
        na.writescript(params)
        na.runscript()
        self.statefiles['coor']=f'{basename}-steered.coor'
        vm=self.writers['vmd']
        vm.newscript(f'{basename}-coor2pdb')
        vm.usescript('namdbin2pdb')
        psf=self.statefiles['psf']
        vm.runscript(psf,f'{basename}-steered.coor',f'{basename}-steered.pdb')
        self.statefiles['pdb']=f'{basename}-steered.pdb'

    def connect(self,specs,basename):
        self.write_heal_patches(specs,f'{basename}-heal-patches')
        self.heal_gaps(specs,f'{basename}-heal-gaps')

    def write_heal_patches(self,specs,basename):
        mol=self.base_molecule
        datafile=f'{basename}.inp'
        writer=self.writers['data']
        writer.newfile(datafile)
        mol.write_heal_patches(writer)
        writer.writefile()
        self.statefiles['datafile']=datafile

    def heal_gaps(self,specs,basename):
        vm=self.writers['psfgen']
        vm.newscript(f'{basename}-loop-closure')
        topfile=os.path.join(self.config.charmm_toppar_path,'mylink.top')
        vm.addline(f'topology {topfile}')
        vm.usescript('loop_closure')
        vm.writescript()
        patchfile=self.statefiles['datafile']
        psf=self.statefiles['psf']
        pdb=self.statefiles['pdb']
        outputbasename=f'{basename}-healed'
        vm.runscript(psf,pdb,patchfile,f'{outputbasename}.psf',f'{outputbasename}.pdb')
        self.statefiles['psf']=f'{outputbasename}.psf'
        self.statefiles['pdb']=f'{outputbasename}.pdb'

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

    def injest_molecules(self):
        self.molecules={}
        psf_exists=False
        if type(self.source)==dict:
            self.basename=self.source.get('rcsb',None)
            bioassemb=self.source.get('biological_assembly',None)
        elif type(self.source)==str:
            self.basename=self.source
            assert os.path.exists(f'{self.basename}.pdb')
            self.pdb=f'{self.basename}.pdb'
            psf_exists=os.path.exists(f'{self.basename}.psf')
            if psf_exists:
                self.psf=f'{self.basename}.psf'
            bioassemb=0
        self.molecules[self.basename]=Molecule(config=self.config,source=self.basename,use_psf=psf_exists).activate_biological_assembly(bioassemb,self.chainIDmanager)
        self.base_molecule=self.molecules[self.basename]
        for p in self.pdbs:
            self.molecules[p]=Molecule(config=self.config,source=p)

class SolvateTask(Task):
    yaml_header='solvate'
    opt_attr=Task.opt_attr+[yaml_header]
    def do(self):
        logger.info(f'Solvate task initiated')
        self.statefiles=self.prior_statefiles()
        basename=f'{self.index:02d}-{self.taskname}'
        vt=self.writers['vmd']
        vt.newscript(basename)
        vt.usescript('solv')
        vt.writescript()
        psf=self.statefiles['psf']
        pdb=self.statefiles['pdb']
        logger.debug(f'solvate: inputs {psf} {pdb} to {vt.basename}.tcl')
        vt.runscript(o=basename,pdb=pdb,psf=psf)
        self.statefiles['pdb']=f'{basename}.pdb'
        self.statefiles['psf']=f'{basename}.psf'

class RelaxTask(Task):
    yaml_header='relax'
    opt_attr=Task.opt_attr+[yaml_header]
    def do(self):
        self.statefiles=self.prior_statefiles()
        basename=f'{self.index:02d}-{self.taskname}'
        self.relax(basename)

