# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" Defines the Filewriter class and several descendents specialized for writing
    scripts for VMD, psfgen, and NAMD
"""
import os
import datetime
import logging

from .command import Command
from .config import Config
from .stringthings import ByteCollector, FileCollector
from .util.colors import *
from .util.logparsers import VMDLog, NAMDLog, PackmolLog, PsfgenLog
from .util.progress import NAMDProgress, PsfgenProgress, PestiferProgress, PackmolProgress
from .util.util import reduce_intlist

logger=logging.getLogger(__name__)

class Filewriter:
    def __init__(self,comment_char='#'):
        self.B=ByteCollector(comment_char=comment_char)
        self.is_written=False
        self.indent=''

    def newfile(self,filename):
        self.filename=filename
        self.B.reset()
        self.is_written=False
        return self

    def injest_file(self,filename):
        self.B.injest_file(filename)

    def addline(self,data,indents=0):
        self.B.addline(indents*self.indent+data)

    def banner(self,data):
        self.B.banner(data)

    def comment(self,data):
        self.B.comment(data)

    def has_statement(self,statement):
        return self.B.has_statement(statement)

    def writefile(self,force=False):
        if not self.is_written or force:
            with open(self.filename,'w') as f:
                f.write(str(self.B))
            self.is_written=True
        else:
            logger.debug(f'{self.filename} has already been written')

class TcLScriptwriter(Filewriter):
    def __init__(self,config):
        self.config=config
        self.progress=self.config.progress
        self.F=FileCollector()
        self.default_ext='.tcl'
        self.default_script=f'pestifer-script{self.default_ext}'
        self.scriptname=self.default_script
        super().__init__(comment_char='#')
    
    def newscript(self,basename=None):
        timestampstr=datetime.datetime.today().ctime()
        if basename:
            self.basename=basename
        else:
            self.basename=os.path.splitext(self.default_script)[0]
        self.scriptname=f'{self.basename}{self.default_ext}'
        self.newfile(self.scriptname)
        msg=f'{__package__}: {self.basename}{self.default_ext}'
        self.comment(msg)
        self.banner(f'Created {timestampstr}')

    def writescript(self):
        self.writefile()

    def addfile(self,filename):
        self.F.append(filename)

class VMD(TcLScriptwriter):
    def __init__(self,config:Config):
        super().__init__(config)
        self.vmd=config.shell_commands['vmd']
        self.tcl_root=config.tcl_root
        self.tcl_pkg_path=config.tcl_pkg_path
        self.tcl_script_path=config.tcl_script_path
        self.vmd_startup=config.vmd_startup_script

    def newscript(self,basename=None,packages=[]):
        super().newscript(basename=basename)
        self.packages=packages
        for p in self.packages:
            self.addline(f'package require {p}')
            self.addline(f'namespace import {p}::*')

    def usescript(self,scriptbasename,local=False,add_banners=False):
        if not local:
            scriptname=os.path.join(self.tcl_script_path,f'{scriptbasename}.tcl')
        else:
            scriptname=os.path.join('.',f'{scriptbasename}.tcl')
        timestampstr=datetime.datetime.today().ctime()
        if not os.path.exists(scriptname):
            raise FileNotFoundError(f'Pestifer script {scriptbasename}.tcl is not found.')
        if add_banners:
            self.banner(f'Begin {scriptbasename}, {timestampstr}')
        self.injest_file(scriptname)
        if add_banners:
            self.banner(f'End {scriptbasename}')

    def writescript(self,force_exit=False):
        if not self.has_statement('exit') or force_exit:
            self.addline('exit')
        self.banner(f'END {__package__.upper()} VMD SCRIPT')
        self.banner(f'Thank you for using {__package__}!')
        super().writescript()

    def set_molecule(self,mol,altcoords=None):
        mol.molid_varname=f'm{mol.molid}'
        ext='.pdb' if mol.rcsb_file_format=='PDB' else '.cif'
        if mol.sourcespecs.get('id',{}):
            self.addline(f'mol new {mol.sourcespecs["id"]}{ext} waitfor all')
        elif mol.sourcespecs.get('prebuilt',{}):
            pdb=mol.sourcespecs['prebuilt']['pdb']
            self.addline(f'mol new {pdb} waitfor all autobonds off')
        elif mol.sourcespecs.get('alphafold',{}):
            pdb=f'{mol.sourcespecs["alphafold"]}.pdb'
            self.addline(f'mol new {pdb} waitfor all')
        else:
            raise Exception(f'Molecule specs has none of "id", "prebuilt", or "alphafold"')
        self.addline(f'set {mol.molid_varname} [molinfo top get id]')
        self.addline(f'set nf [molinfo ${mol.molid_varname} get numframes]')
        self.addline(r'if { $nf > 1 } { animate delete beg 0 end [expr $nf - 2] $'+mol.molid_varname+r' }')
        if altcoords:
            mol.set_coords(altcoords)
            self.addline(f'mol addfile {altcoords} waitfor all')
            self.addline(f'set nf [molinfo ${mol.molid_varname} get numframes]')
            self.addline(r'animate delete beg 0 end [expr $nf - 2] $'+mol.molid_varname+r' }')
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
                self.addline(f'set TMP [atomselect ${mol.molid_varname} "serial {vmd_red_list}"]')
                self.addline(f'$TMP set chain {c}')
                self.addline(f'$TMP set resid [ list {residlist} ]')
                self.addline(f'$TMP delete')
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
        self.addline(f'set TMP [atomselect ${molid_varname} all]')
        self.addline(f'$TMP writepdb {basename}.pdb')
        self.addline(f'$TMP delete')

    def center_molecule(self,mol):
        self.banner('Centering molecule')
        self.addline(f'set TMP [atomselect ${mol.molid_varname} all]')
        self.addline(f'set or [measure center $TMP weight mass]')
        self.addline(f'$TMP moveby [vescale -1 $or]')
        self.addline(f'$TMP delete')

    def shift_coords(self,factors):
        self.banner('Shifting')
        
    def runscript(self,*args,**options):
        assert hasattr(self,'scriptname'),f'No scriptname set.'
        self.logname=f'{self.basename}.log'
        self.logparser=VMDLog(basename=self.basename)
        logger.debug(f'Log file: {self.logname}')
        c=Command(f'{self.vmd} -dispdev text -startup {self.vmd_startup} -e {self.scriptname} -args --tcl-root {self.tcl_root}',**options)
        progress_struct=None
        progress_title=options.get('progress_title','')
        if self.progress and progress_title!='':
            progress_struct=PestiferProgress(name=progress_title)
            self.logparser.enable_progress_bar(progress_struct)
        return c.run(logfile=self.logname,logparser=self.logparser)
    
    def cleanup(self,cleanup=False):
        if cleanup:
            nremoved=len(self.F)
            self.F.flush()
            logger.info(f'Post-execution clean-up: {nremoved} files removed.')

class Psfgen(VMD):
    def __init__(self,config):
        super().__init__(config)
        self.charmmff=config.RM.charmmff_content
        self.charmmff_config=config['user']['charmmff']
        self.psfgen_config=config['user']['psfgen']

    def fetch_standard_charmm_topologies(self):
        topology_local=[]
        for t in self.charmmff_config['standard']['topologies']+self.charmmff_config['custom']['topologies']:
            self.charmmff.copy_charmmfile_local(t)
            topology_local.append(t)
        # the order of the topologies is imporant, and the top_all35_ethers.rtf should be last
        if 'top_all35_ethers.rtf' in topology_local:
            topology_local.remove('top_all35_ethers.rtf')
            topology_local.append('top_all35_ethers.rtf')
        logger.debug(f'local topologies: {topology_local}')
        return topology_local

    def newscript(self,basename=None,packages=[],additional_topologies=[]):
        super().newscript(basename=basename,packages=packages)
        self.addline('package require psfgen')
        self.addline('psfcontext mixedcase')
        self.topologies=self.fetch_standard_charmm_topologies()
        for at in additional_topologies:
            assert os.sep not in at,f'Topology file {at} must not contain a path.'
            self.charmmff.copy_charmmfile_local(at)
            self.topologies.append(at)
        self.topologies=list(set(self.topologies))  # remove duplicates
        if 'top_all35_ethers.rtf' in self.topologies:
            self.topologies.remove('top_all35_ethers.rtf')
            self.topologies.append('top_all35_ethers.rtf')
        for t in self.topologies:
            self.addline(f'topology {t}')
        for pdba in self.psfgen_config['aliases']:
            alias_tokens=pdba.split()
            if alias_tokens[1]=='*': # wild-card for all protein residues
                for protein_resname in self.psfgen_config['segtypes']['protein']['resnames']:
                    this_pdba=f'{alias_tokens[0]} {protein_resname} {alias_tokens[2]} {alias_tokens[3]}'
                    self.addline(f'pdbalias {this_pdba}')
            else:
                self.addline(f'pdbalias {pdba}')

    def atomselect_macros(self):
        """Converts all resnames in the base config by segtype into updated atomselect macros
        This is usually only done one time, or whenever developer updates base config
        The command 'inittcl' is the only way this method is ever called
        """
        for segtypename,segtyperec in self.psfgen_config['segtypes'].items():
            logger.debug(f'Searching base record of segtype {segtypename} for resname atomselect macro')
            logger.debug(', '.join(segtyperec.keys()))
            if segtyperec.get('macro',False)==True:
                new_macro='resname '+' '.join(segtyperec['resnames'])
                self.addline(f'update_atomselect_macro {segtypename} "{new_macro}" 0')

    def load_project(self,*objs):
        logger.debug(f'load_project called with {len(objs)} arguments')
        for _ in objs:
            logger.debug(f'   {_}')
        if len(objs)==1:
            basename=objs[0]
            self.addline(f'readpsf {basename}.psf pdb {basename}.pdb')
        else:
            psf,pdb=objs
            self.addline(f'readpsf {psf} pdb {pdb}')
        
    def describe_molecule(self,mol):
        self.molid_varname=f'm{mol.molid}'
        self.addline(f'mol top ${self.molid_varname}')
        mol.write_TcL(self)

    def writescript(self,statename,guesscoord=True,regenerate=True,writepsf=True,writepdb=True,force_exit=False):
        if guesscoord:
            self.addline('guesscoord')
        if regenerate:
            self.addline('regenerate angles dihedrals')
        if writepsf:
            self.addline(f'writepsf cmap {statename}.psf')
        if writepdb:
            self.addline(f'writepdb {statename}.pdb')
        super().writescript(force_exit=force_exit)

    def runscript(self,*args,**options):
        assert hasattr(self,'scriptname'),f'No scriptname set.'
        self.logname=f'{self.basename}.log'
        self.logparser=PsfgenLog(basename=self.basename)
        logger.debug(f'Log file: {self.logname}')
        clean_options=options.copy()
        for k,v in options.items():
            if v=='':
                clean_options.pop(k)
        c=Command(f'{self.vmd} -dispdev text -startup {self.vmd_startup} -e {self.scriptname} -args --tcl-root {self.tcl_root}',**clean_options)
        progress_struct=None
        if self.progress:
            progress_struct=PsfgenProgress()
            self.logparser.enable_progress_bar(progress_struct)
        return c.run(logfile=self.logname,logparser=self.logparser)
    
class NAMD(TcLScriptwriter):
    def __init__(self,config:Config):
        super().__init__(config)
        self.charmmff=config.RM.charmmff_content
        self.charmmff_config=config['user']['charmmff']
        self.charmrun=config.shell_commands['charmrun']
        self.namd=config.shell_commands['namd3']
        self.namdgpu=config.shell_commands['namd3gpu']
        self.namd_type=config.namd_type
        self.namd_config=config['user']['namd']
        self.namd_version=int(self.namd_config['namd-version'])
        logger.debug(f'Using NAMD v {self.namd_version}')
        self.local_ncpus=config.local_ncpus
        self.ncpus=config.ncpus
        self.gpu_devices=config.gpu_devices
        self.slurmvars=config.slurmvars
        self.ngpus=config.ngpus
        if self.namd_version==2:
            self.namd_deprecates={}
        else:
            self.namd_deprecates=config.namd_deprecates
        logger.debug(f'{self.ncpus} cpus are available for namd')
        self.default_ext='.namd'
    
    def fetch_standard_charmm_parameters(self):
        parameters_local=[]
        for t in list(set(self.charmmff_config['standard']['parameters']+self.charmmff_config['custom']['parameters'])):
            self.charmmff.copy_charmmfile_local(t)
            parameters_local.append(t)
        logger.debug(f'local parameters: {parameters_local}')
        return parameters_local

    def newscript(self,basename=None,addl_paramfiles=[]):
        super().newscript(basename)
        self.scriptname=f'{basename}{self.default_ext}'
        self.banner('NAMD script')
        self.parameters=self.fetch_standard_charmm_parameters()
        for at in addl_paramfiles:
            if not at in self.parameters:
                self.charmmff.copy_charmmfile_local(at)
                self.parameters.append(at)
        for p in self.parameters:
            assert os.sep not in p
            self.addline(f'parameters {p}')

    def writescript(self,params,cpu_override=False):
        logger.debug(f'params: {params}')
        tailers=['minimize','run','numsteps']
        for k,v in params.items():
            if k not in tailers:
                if type(v)==list:
                    for val in v:
                        if k=='tcl':
                            self.addline(val)
                        else:
                            self.addline(f'{self.namd_deprecates.get(k,k)} {val}')
                else:
                    if k=='tcl':
                        self.addline(v)
                    else:
                        self.addline(f'{self.namd_deprecates.get(k,k)} {v}')
        if self.namd_type=='gpu' and not cpu_override:
            for k,v in self.namd_config['gpu-resident'].items():
                self.addline(f'{k} {v}')
        for t in tailers:
            if t in params:
                self.addline(f'{self.namd_deprecates.get(t,t)} {params[t]}')
        super().writescript()

    def runscript(self,**kwargs):
        assert hasattr(self,'scriptname'),f'No scriptname set.'
        if kwargs.get('local_execution_only',False):
            use_cpu_count=self.local_ncpus
        else:
            use_cpu_count=self.ncpus
        if kwargs.get('single_gpu_only',False):
            use_gpu_count=1
            use_gpu_devices='0'
        else:
            use_gpu_count=self.ngpus
            use_gpu_devices=self.gpu_devices
        if self.namd_type=='cpu' or kwargs.get('cpu_override',False):
            c=Command(f'{self.charmrun} +p {use_cpu_count} {self.namd} {self.scriptname}')
        elif self.namd_type=='gpu':
            if len(self.slurmvars)>0:
                use_cpu_count=8 if use_gpu_count==1 else (use_gpu_count-1)*8 + 8-(use_gpu_count-1)
            else:
                use_cpu_count=self.local_ncpus
            c=Command(f'{self.namdgpu} +p{use_cpu_count} +setcpuaffinity +devices {use_gpu_devices} {self.scriptname}')
        self.logname=f'{self.basename}.log'
        self.logparser=NAMDLog(basename=self.basename)
        progress_struct=None
        if self.progress:
            logger.debug(f'NAMD runscript using progress')
            progress_struct=NAMDProgress()
            self.logparser.enable_progress_bar(progress_struct)
        else:
            logger.debug(f'NAMD runscript NOT using progress')
        return c.run(logfile=self.logname,logparser=self.logparser)

class PackmolInputWriter(Filewriter):
    def __init__(self,config):
        super().__init__(comment_char='#')
        self.indent=4*' '
        self.config=config
        self.progress=self.config.progress
        self.F=FileCollector()
        self.default_ext='.inp'
        self.default_script=f'packmol{self.default_ext}'
        self.scriptname=self.default_script

    def newscript(self,basename=None):
        timestampstr=datetime.datetime.today().ctime()
        if basename:
            self.basename=basename
        else:
            self.basename=os.path.splitext(self.default_script)[0]
        self.scriptname=f'{self.basename}{self.default_ext}'
        self.newfile(self.scriptname)
        self.banner(f'{__package__}: {self.basename}{self.default_ext}')
        self.banner(f'Created {timestampstr}')

    def writescript(self):
        self.writefile()

    def runscript(self,*args,**options):
        assert hasattr(self,'scriptname'),f'No scriptname set.'
        self.logname=f'{self.basename}.log'
        self.logparser=PackmolLog(basename=self.basename)
        logger.debug(f'Log file: {self.logname}')
        cmd=Command(f'{self.config.shell_commands["packmol"]} < {self.scriptname}')
        progress_struct=None
        if self.progress:
            progress_struct=PackmolProgress()
            self.logparser.enable_progress_bar(progress_struct)
        return cmd.run(ignore_codes=[173],logfile=self.logname,logparser=self.logparser)
    