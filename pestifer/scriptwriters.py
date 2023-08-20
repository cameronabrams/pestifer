"""

.. module:: psfgen
   :synopsis: Manages script generation and execution of psfgen under vmd
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
logger=logging.getLogger(__name__)
from .command import Command
import os
from .molecule import Molecule
from .stringthings import ByteCollector, FileCollector, my_logger
import datetime
import shutil

class Filewriter:
    def __init__(self):
        self.B=ByteCollector()

    def newfile(self,filename):
        self.filename=filename
        self.B.reset()

    def injest_file(self,filename):
        self.B.injest_file(filename)

    def addline(self,data):
        self.B.addline(data)

    def writefile(self):
        with open(self.filename,'w') as f:
            f.write(str(self.B))

class Scriptwriter(Filewriter):
    def __init__(self):
        self.F=FileCollector()
        self.default_ext='.tcl'
        self.default_script=f'pestifer-script{self.default_ext}'
        self.scriptname=self.default_script
        super().__init__()
    
    def newscript(self,basename=None):
        timestampstr=datetime.datetime.today().ctime()
        if basename:
            self.basename=basename
        else:
            self.basename=os.path.splitext(self.default_script)[0]
        self.scriptname=f'{self.basename}{self.default_ext}'
        self.newfile(self.scriptname)
        self.B.banner(f'pestifer: {self.basename}{self.default_ext}')
        self.B.banner(f'Created {timestampstr}')

    def writescript(self):
        self.writefile()

class VMD(Scriptwriter):
    def __init__(self,config):
        super().__init__()
        self.config=config
        self.tcl_path=config.tcl_path
        self.vmd_startup=config.vmd_startup_script
        self.default_script=config['vmd_scriptname']

    def usescript(self,scriptbasename):
        scriptname=os.path.join(self.tcl_path,f'scripts/{scriptbasename}.tcl')
        timestampstr=datetime.datetime.today().ctime()
        if not os.path.exists(scriptname):
            raise FileNotFoundError(f'Pestifer script {scriptbasename}.tcl is not found.')
        self.B.banner(f'Injested from {scriptbasename}, {timestampstr}')
        self.injest_file(scriptname)

    def loadmodule(self,modulename):
        procname=os.path.join(self.tcl_path,f'proc/{modulename}.tcl')
        self.B.addline(f'source {procname}')

    def endscript(self):
        self.B.addline('exit')
        self.B.banner('END PESTIFER VMD SCRIPT')
        self.B.banner('Thank you for using pestifer!')

    def set_molecule(self,mol):
        mol.molid_varname=f'm{mol.molid}'
        self.B.addline(f'mol new {mol.source}.pdb waitfor all')
        self.B.addline(f'set {mol.molid_varname} [molinfo top get id]')

    def load_psf_pdb(self,*objs,new_molid_varname='mX'):
        if len(objs)==1:
            basename=objs[0]
            pdb=f'{basename}.pdb'
            psf=f'{basename}.psf'
        else:
            psf,pdb=objs
        self.B.addline(f'mol new {psf}')
        self.B.addline(f'mol addfile {pdb} waitfor all')
        self.B.addline(f'set {new_molid_varname} [molinfo top get id]')

    def write_pdb(self,basename,molid_varname):
        self.B.addline(f'set X [atomselect ${molid_varname} all]')
        self.B.addline(f'$X writepdb {basename}.pdb')

    def center_molecule(self,mol:Molecule):
        self.B.banner('Centering molecule')
        self.B.addline(f'set a [atomselect ${mol.molid_varname} all]')
        self.B.addline(f'set or [measure center $a weight mass]')
        self.B.addline(f'$a moveby [vescale -i $or]')
        self.B.addline(f'$a delete')

    def reset_molecule_orientation(self,mol:Molecule,specs):
        selspec=specs.get('selspec',{})
        if not selspec:
            return
        group1=selspec.get('group1','')
        group2=selspec.get('group2','')
        if not group1 or not group2:
            return
        self.B.banner('Resetting molecule orientation')
        self.B.addline(f'set g1 [measure center [atomselect ${mol.molid_varname} "protein and {group1}] weight mass]')
        self.B.addline(f'set g2 [measure center [atomselect ${mol.molid_varname} "protein and {group2}] weight mass]')
        self.B.addline('set pi 3.1415928')
        self.B.addline('set dv [vecsub $g1 $g2]')
        self.B.addline('set d [veclength $dv]')
        self.B.addline('set cp [expr [lindex $dv 0]/$d]')
        self.B.addline('set sp [expr [lindex $dv 1]/$d]')
        self.B.addline('set p [expr acos($cp)]')
        self.B.addline('if { [expr $sp < 0.0] } {')
        self.B.addline('    set p [expr 2*$pi-$p]')
        self.B.addline('}')
        self.B.addline('set ct [exr [lindex $dv 2]/$d]')
        self.B.addline('set t [expr acos($ct)]')
        self.B.addline(f'set a [atomselect ${mol.molid_varname} all]')
        self.B.addline('$a move [transaxis z [expr -1 * $p] rad]')
        self.B.addline('$a move [transaxis y [expr -1 * $t] rad]')
        self.B.addline('$a delete')
        self.B.banner('Done resetting molecule orientation')

    def runscript(self,*args,**options):
        assert hasattr(self,'scriptname'),f'No scriptname set.'
        c=Command(f'{self.config["vmd"]} -dispdev text -startup {self.vmd_startup} -e {self.scriptname} -args -respath {self.tcl_path}',**options)
        c.run()
        self.logname=f'{self.basename}.log'
        with open(self.logname,'w') as f:
            my_logger(f'STDOUT from "{c.command}"',f.write)
            f.write(c.stdout+'\n')
            my_logger(f'STDERR from "{c.command}"',f.write)
            f.write(c.stderr+'\n')
            my_logger(f'END OF LOG',f.write)
    
    def cleanup(self,cleanup=False):
        if cleanup:
            logger.info(f'Post-execution clean-up: {len(self.F)} files removed.')
            for file in self.F:
                if os.path.exists(file):
                    os.remove(file)

class Psfgen(VMD):
    def __init__(self,config):
        super().__init__(config)
        self.pestifer_charmmpath=config.charmm_toppar_path
        self.user_charmm_topparpath=config.user_charmm_toppar_path
        self.default_script=config['psfgen_scriptname']

    def newscript(self,basename=None):
        super().newscript(basename=basename)
        self.B.addline('package require psfgen')
        self.B.addline('psfcontext mixedcase')

    def topo_aliases(self):
        for t in self.config['StdCharmmTopo']:
            ft=os.path.join(self.user_charmm_topparpath,t)
            self.B.addline(f'topology {ft}')
        for lt in self.config['LocalCharmmTopo']:
            ft=os.path.join(self.pestifer_charmmpath,lt)
            self.B.addline(f'topology {ft}')
        for pdba in self.config['PDBAliases']:
            self.B.addline(f'pdbalias {pdba}')
        for k,v in self.config['PDB_to_CHARMM_Resnames'].items():
            self.B.addline(f'set RESDICT({k}) {v}')
        for k,v in self.config['PDB_to_CHARMM_Atomnames'].items():
            self.B.addline(f'set ANAMEDICT({k}) {v}')
        self.B.banner('END HEADER')

    def load_project(self,*objs):
        if len(objs)==1:
            basename=objs[0]
            self.B.addline(f'readpsf {basename}.psf pdb {basename}.pdb')
        else:
            psf,pdb=objs
            self.B.addline(f'readpsf {psf} pdb {pdb}')
        
    def describe_molecule(self,mol:Molecule,mods):
        molid_varname=f'm{mol.molid}'
        self.B.addline(f'mol top ${molid_varname}')
        mol.write_TcL(self.B,mods,file_collector=self.F)

    def complete(self,statename):
        self.B.addline('guesscoord')
        self.B.addline('regenerate angles dihedrals')
        psf=f'{statename}.psf'
        pdb=f'{statename}.pdb'
        self.B.addline(f'writepsf cmap {psf}')
        self.B.addline(f'writepdb {pdb}')
        return {'psf':psf,'pdb':pdb}

class NAMD2(Scriptwriter):
    def __init__(self,config):
        super().__init__()
        self.config=config
        self.templates_path=config.namd_template_path
        self.default_basename=config['namd2_configbasename']
        self.max_cpu_count=os.cpu_count()
        self.default_ext='.namd'
        self.default_script=f'{self.default_basename}{self.default_ext}'
        
    def newscript(self,basename=None):
        super().newscript(basename)
        self.scriptname=f'{basename}{self.default_ext}'
        self.addline('NAMD2 script')

    def writescript(self,params):
        for k,v in params.items():
            if type(v)==list:
                for val in v:
                    self.addline(f'{k} {val}')
            else:
                self.addline(f'{k} {v}')
        super().writescript()

    def runscript(self):
        assert hasattr(self,'scriptname'),f'No scriptname set.'
        c=Command(f'{self.config["charmrun"]} +p {self.max_cpu_count} {self.config["namd2"]} {self.scriptname}')
        c.run()
        self.logname=f'{self.basename}.log'
        with open(self.logname,'w') as f:
            my_logger(f'STDOUT from "{c.command}"',f.write)
            f.write(c.stdout+'\n')
            my_logger(f'STDERR from "{c.command}"',f.write)
            f.write(c.stderr+'\n')
            my_logger(f'END OF LOG',f.write)