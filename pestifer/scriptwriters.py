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

class Scriptwriter:
    def __init__(self):
        self.B=ByteCollector()
        self.F=FileCollector()
        self.default_script='pestifer-script.tcl'

    def addline(self,data):
        self.B.addline(data)

    def newscript(self,basename=None):
        timestampstr=datetime.datetime.today().ctime()
        if basename:
            self.basename=basename
        else:
            self.basename=os.path.splitext(self.default_script)[0]
        self.B.reset()
        self.B.banner(f'PESTIFER SCRIPT {self.basename}.tcl')
        self.B.banner(f'Created {timestampstr}')

class VMD(Scriptwriter):
    def __init__(self,config):
        self.config=config
        self.tcl_path=config.tcl_path
        self.vmd_startup=config.vmd_startup_script
        self.default_script=config['vmd_scriptname']
        super().__init__()

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

    def writescript(self):
        self.scriptname=f'{self.basename}.tcl'
        with open(self.scriptname,'w') as f:
            f.write(str(self.B))

    def runscript(self):
        assert hasattr(self,'scriptname'),f'No scriptname set.'
        c=Command(f'{self.config["vmd"]} -dispdev text -startup {self.vmd_startup} -e {self.scriptname} -args -respath {self.tcl_path}')
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

    def topo_aliases(self):
        self.B.addline('psfcontext mixedcase')
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

    def transform_postmods(self):
        self.B.addline('guesscoord')
        self.B.addline('regenerate angles dihedrals')
        psf=f'{self.basename}.psf'
        pdb=f'{self.basename}.pdb'
        self.B.addline(f'writepsf cmap {psf}')
        self.B.addline(f'writepdb {pdb}')

class NAMD2(Scriptwriter):
    def __init__(self,config):
        self.config=config
        self.templates_path=config.namd_template_path
        self.default_script=config['namd2_configname']
        super().__init__()