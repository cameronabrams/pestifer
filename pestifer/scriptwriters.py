# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Scriptwriters
"""
import logging
logger=logging.getLogger(__name__)
from .command import Command
import os
from .util import reduce_intlist
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

    def banner(self,data):
        self.B.banner(data)

    def comment(self,data):
        self.B.comment(data)

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
        self.banner(f'pestifer: {self.basename}{self.default_ext}')
        self.banner(f'Created {timestampstr}')

    def writescript(self):
        self.writefile()

    def addfile(self,filename):
        self.F.append(filename)

class VMD(Scriptwriter):
    def __init__(self,config):
        super().__init__()
        self.config=config
        self.vmd=config.vmd
        self.tcl_path=config.tcl_path
        self.tcl_proc_path=config.tcl_proc_path
        self.tcl_script_path=config.tcl_script_path
        self.vmd_startup=config.vmd_startup_script

    def usescript(self,scriptbasename):
        scriptname=os.path.join(self.tcl_script_path,f'{scriptbasename}.tcl')
        timestampstr=datetime.datetime.today().ctime()
        if not os.path.exists(scriptname):
            raise FileNotFoundError(f'Pestifer script {scriptbasename}.tcl is not found.')
        self.banner(f'Begin {scriptbasename}, {timestampstr}')
        self.injest_file(scriptname)
        self.banner(f'End {scriptbasename}')

    def endscript(self):
        self.addline('exit')
        self.banner('END PESTIFER VMD SCRIPT')
        self.banner('Thank you for using pestifer!')

    def set_molecule(self,mol,altcoords=None):
        mol.molid_varname=f'm{mol.molid}'
        ext='.pdb' if mol.rcsb_file_format=='PDB' else '.cif'
        self.addline(f'mol new {mol.sourcespecs["id"]}{ext} waitfor all')
        self.addline(f'set {mol.molid_varname} [molinfo top get id]')
        if altcoords:
            self.addline(f'mol addfile {altcoords}')
            self.addline(f'animate delete beg 0 end 0 {mol.molid_varname}')
        if mol.rcsb_file_format=='mmCIF':
            # VMD appends a "1" to any two-letter chain ID from a CIF file,
            # so let's undo that
            # also, CIF HETATM records have a "." for residue number, which
            # VMD interprets as 0, so we'll replace that with the residue numbers
            # assigned by pestifer when it reads in the cif file
            au=mol.asymmetric_unit
            residues=au.residues
            uCIDs=residues.uniqattrs(['chainID'])['chainID']
            self.comment('Resetting chains and resids for this CIF-source molecule')
            for c in uCIDs:
                chain=residues.filter(chainID=c)
                resids=[]
                for x in chain:
                    resids.extend([str(y.resseqnum) for y in x.atoms])
                residlist=' '.join(resids)
                serials=chain.atom_serials(as_type=int)
                vmd_red_list=reduce_intlist(serials)
                self.addline(f'set a [atomselect ${mol.molid_varname} "serial {vmd_red_list}"]')
                self.addline(f'$a set chain {c}')
                self.addline(f'$a set resid [ list {residlist} ]')
                self.addline(f'$a delete')
            self.comment('Done.')

    def backup_selection(self,selname,dataholder='data',attributes=['chain','x','y','z','resid','resname','name']):
        self.addline(f'set {dataholder} [ backup ${selname} [ list {" ".join(attributes)} ] ]')
    
    def restore_selection(self,selname,dataholder='data',attributes=['chain','x','y','z','resid','resname','name']):
        self.addline(f'restore ${selname} [ list {" ".join(attributes)} ]  ${dataholder}')

    def load_psf_pdb(self,*objs,new_molid_varname='mX'):
        if len(objs)==1:
            basename=objs[0]
            pdb=f'{basename}.pdb'
            psf=f'{basename}.psf'
        else:
            psf,pdb=objs
        self.addline(f'mol new {psf}')
        self.addline(f'mol addfile {pdb} waitfor all')
        self.addline(f'set {new_molid_varname} [molinfo top get id]')
        self.molid_varname=new_molid_varname

    def write_pdb(self,basename,molid_varname):
        self.addline(f'set X [atomselect ${molid_varname} all]')
        self.addline(f'$X writepdb {basename}.pdb')

    def center_molecule(self,mol):
        self.banner('Centering molecule')
        self.addline(f'set a [atomselect ${mol.molid_varname} all]')
        self.addline(f'set or [measure center $a weight mass]')
        self.addline(f'$a moveby [vescale -i $or]')
        self.addline(f'$a delete')

    # def reset_molecule_orientation(self,mol,specs):
    #     selspec=specs.get('selspec',{})
    #     if not selspec:
    #         return
    #     group1=selspec.get('group1','')
    #     group2=selspec.get('group2','')
    #     if not group1 or not group2:
    #         return
    #     self.banner('Resetting molecule orientation')
    #     self.addline(f'set g1 [measure center [atomselect ${mol.molid_varname} "protein and {group1}] weight mass]')
    #     self.addline(f'set g2 [measure center [atomselect ${mol.molid_varname} "protein and {group2}] weight mass]')
    #     self.addline('set pi 3.1415928')
    #     self.addline('set dv [vecsub $g1 $g2]')
    #     self.addline('set d [veclength $dv]')
    #     self.addline('set cp [expr [lindex $dv 0]/$d]')
    #     self.addline('set sp [expr [lindex $dv 1]/$d]')
    #     self.addline('set p [expr acos($cp)]')
    #     self.addline('if { [expr $sp < 0.0] } {')
    #     self.addline('    set p [expr 2*$pi-$p]')
    #     self.addline('}')
    #     self.addline('set ct [exr [lindex $dv 2]/$d]')
    #     self.addline('set t [expr acos($ct)]')
    #     self.addline(f'set a [atomselect ${mol.molid_varname} all]')
    #     self.addline('$a move [transaxis z [expr -1 * $p] rad]')
    #     self.addline('$a move [transaxis y [expr -1 * $t] rad]')
    #     self.addline('$a delete')
    #     self.banner('Done resetting molecule orientation')

    def runscript(self,*args,**options):
        assert hasattr(self,'scriptname'),f'No scriptname set.'
        c=Command(f'{self.vmd} -dispdev text -startup {self.vmd_startup} -e {self.scriptname} -args -respath {self.tcl_proc_path}',**options)
        c.run()
        self.logname=f'{self.basename}.log'
        with open(self.logname,'w') as f:
            my_logger(f'STDOUT from "{c.c}"',f.write)
            f.write(c.stdout+'\n')
            my_logger(f'STDERR from "{c.c}"',f.write)
            f.write(c.stderr+'\n')
            my_logger(f'END OF LOG',f.write)
    
    def cleanup(self,cleanup=False):
        if cleanup:
            nremoved=len(self.F)
            self.F.flush()
            logger.info(f'Post-execution clean-up: {nremoved} files removed.')

class Psfgen(VMD):
    def __init__(self,config):
        super().__init__(config)
        self.charmmff_config=config['user']['charmmff']
        self.psfgen_config=config['user']['psfgen']
        self.pestifer_charmmff_toppar_path=config.charmmff_toppar_path
        self.pestifer_charmmff_custom_path=config.charmmff_custom_path
        self.user_charmmff_toppar_path=config.user_charmmff_toppar_path
        # self.default_script=config['psfgen_scriptname']

    def newscript(self,basename=None):
        super().newscript(basename=basename)
        self.addline('package require psfgen')
        self.addline('psfcontext mixedcase')

    def topo_aliases(self):
        for t in self.charmmff_config['standard']['topologies']:
            ft=os.path.join(self.pestifer_charmmff_toppar_path,t)
            self.addline(f'topology {ft}')
        for lt in self.charmmff_config['custom']['topologies']:
            ft=os.path.join(self.pestifer_charmmff_custom_path,lt)
            self.addline(f'topology {ft}')
        for pdba in self.psfgen_config['aliases']:
            self.addline(f'pdbalias {pdba}')
        self.banner('END HEADER')

    def load_project(self,*objs):
        if len(objs)==1:
            basename=objs[0]
            self.addline(f'readpsf {basename}.psf pdb {basename}.pdb')
        else:
            psf,pdb=objs
            self.addline(f'readpsf {psf} pdb {pdb}')
        
    def describe_molecule(self,mol):
        molid_varname=f'm{mol.molid}'
        self.addline(f'mol top ${molid_varname}')
        mol.write_TcL(self)

    def complete(self,statename):
        self.addline('guesscoord')
        self.addline('regenerate angles dihedrals')
        psf=f'{statename}.psf'
        pdb=f'{statename}.pdb'
        self.addline(f'writepsf cmap {psf}')
        self.addline(f'writepdb {pdb}')

class NAMD2(Scriptwriter):
    def __init__(self,config):
        super().__init__()
        self.charmmff_config=config['user']['charmmff']
        self.charmrun=config.charmrun
        self.namd2=config.namd2
        self.namd2_config=config['user']['namd2']
        self.max_cpu_count=os.cpu_count()
        self.default_ext='.namd'
        if config.user_charmmff_toppar_path:
            self.standard_charmmff_parfiles=[os.path.join(config.user_charmmff_toppar_path,x) for x in self.charmmff_config['standard']['parameters']]
        else:
            self.standard_charmmff_parfiles=[os.path.join(config.charmmff_toppar_path,x) for x in self.charmmff_config['standard']['parameters']]
        self.custom_charmmff_parfiles=[os.path.join(config.charmmff_custom_path,x) for x in self.charmmff_config['custom']['parameters']]

    def write_parcommands(self,filename,fetch=False):
        with open(filename,'w') as f:
            for pf in self.standard_charmmff_parfiles+self.custom_charmmff_parfiles:
                d,bn=os.path.split(pf)
                if fetch:
                    shutil.copy(pf,bn)
                f.write(f'parameters {bn}\n')

    def newscript(self,basename=None):
        super().newscript(basename)
        self.scriptname=f'{basename}{self.default_ext}'
        self.banner('NAMD2 script')

    def writescript(self,params):
        for k,v in params.items():
            if type(v)==list:
                for val in v:
                    if k=='tcl':
                        self.addline(val)
                    else:
                        self.addline(f'{k} {val}')
            else:
                if k=='tcl':
                    self.addline(v)
                else:
                    self.addline(f'{k} {v}')
        super().writescript()

    def runscript(self):
        assert hasattr(self,'scriptname'),f'No scriptname set.'
        c=Command(f'{self.charmrun} +p {self.max_cpu_count} {self.namd2} {self.scriptname}')
        c.run()
        self.logname=f'{self.basename}.log'
        with open(self.logname,'w') as f:
            my_logger(f'STDOUT from "{c.command}"',f.write)
            f.write(c.stdout+'\n')
            my_logger(f'STDERR from "{c.command}"',f.write)
            f.write(c.stderr+'\n')
            my_logger(f'END OF LOG',f.write)