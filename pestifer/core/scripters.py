# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
Defines the :class:`Filewriter` class and several descendents specialized for writing
and running scripts for VMD, psfgen, NAMD, and packmol.  These include:

- :class:`TcLScripter`: A base class for creating and running Tcl scripts.
- :class:`VMDScripter`: A class for creating and running VMD Tcl scripts.
- :class:`PsfgenScripter`: A class for creating and running psfgen Tcl scripts.
- :class:`NAMDScripter`: A class for creating and running NAMD Tcl scripts.
- :class:`PackmolScripter`: A class for creating and running Packmol scripts.
"""
import os
import datetime
import logging

from .command import Command
from .config import Config
from .stringthings import ByteCollector, FileCollector
from ..util.colors import *
from ..util.logparsers import VMDLog, NAMDLog, PackmolLog, PsfgenLog
from ..util.progress import NAMDProgress, PsfgenProgress, PestiferProgress, PackmolProgress
from ..util.util import reduce_intlist

logger=logging.getLogger(__name__)

class Filewriter:
    """
    A class for writing files with specific formats.
    
    Parameters
    ----------
    comment_char : str, optional
        The character used for comments in the script. Default is '#'.
    """
    def __init__(self,comment_char='#'):
        self.B=ByteCollector(comment_char=comment_char)
        self.is_written=False
        self.indent=''

    def newfile(self,filename):
        """
        Initializes a new file for writing.
        
        Parameters
        ----------
        filename : str
            The name of the file to be created or overwritten.
        """
        self.filename=filename
        self.B.reset()
        self.is_written=False
        return self

    def ingest_file(self,filename):
        """
        Ingest a file into the script writer.

        Parameters
        ----------
        filename : str
            The name of the file to be ingested. The file's contents will be added to the script.

        """
        self.B.ingest_file(filename)

    def addline(self,data,indents=0):
        """
        Add a line of code to the script.
        
        Parameters
        ----------
        data : str
            The line of code to be added.
        indents : int, optional
            The number of indentation levels to apply. Default is 0.
        """
        self.B.addline(indents*self.indent+data)

    def banner(self,data):
        """
        Add a banner comment to the script.

        Parameters
        ----------
        data : str
            The banner comment to be added.
        """
        self.B.banner(data)

    def comment(self,data):
        """
        Add a comment to the script.

        Parameters
        ----------
        data : str
            The comment to be added.
        """
        self.B.comment(data)

    def has_statement(self,statement):
        """
        Check if a specific statement is present in the script.
        
        Parameters
        ----------
        statement : str
            The statement to check for.
        """
        return self.B.has_statement(statement)

    def writefile(self,force=False):
        """
        Write the contents of the script to a file.
        
        Parameters
        ----------
        force : bool, optional
            If True, overwrite the file even if it has already been written.
        """
        if not self.is_written or force:
            with open(self.filename,'w') as f:
                f.write(str(self.B))
            self.is_written=True
        else:
            logger.debug(f'{self.filename} has already been written')

class TcLScripter(Filewriter):
    """
    This class extends the Filewriter class to provide functionality for creating and managing Tcl scripts.

    Parameters
    ----------
    config : Config
        The configuration object containing settings for the script.
    """
    def __init__(self,config):
        self.config=config
        self.progress=self.config.progress
        self.F=FileCollector()
        self.default_ext='.tcl'
        self.default_script=f'pestifer-script{self.default_ext}'
        self.scriptname=self.default_script
        super().__init__(comment_char='#')
    
    def newscript(self,basename=None):
        """
        Initialize a new Tcl script with a specified basename.
        If no basename is provided, a default script name is used.

        Parameters
        ----------
        basename : str, optional
            The base name for the script file. If not provided, a default name is used.
        """
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
        """
        Writes the Tcl script to the file specified by the ``scriptname`` attribute.
        """
        self.writefile()

    def addfile(self,filename):
        """
        Add a file to the list of files associated with the script; these might be inputs or outputs (mostly outputs).
        This method appends the specified filename to the internal ``FileCollector``.

        Parameters
        ----------
        filename : str
            The name of the file to be added.
        """
        self.F.append(filename)

class VMDScripter(TcLScripter):
    """
    This class extends the TcLScripter class to provide functionality for creating and managing VMD scripts.

    Parameters
    ----------
    config : Config
        The configuration object containing settings for the script.
    """ 

    def __init__(self,config:Config):
        super().__init__(config)
        self.vmd=config.shell_commands['vmd']
        self.tcl_root=config.tcl_root
        self.tcl_pkg_path=config.tcl_pkg_path
        self.tcl_script_path=config.tcl_script_path
        self.vmd_startup=config.vmd_startup_script
        self.indent=' '*4

    def newscript(self,basename=None,packages=[]):
        """
        Initialize a new VMD script with a specified basename and packages.
        If no basename is provided, a default script name is used.
        
        Parameters
        ----------
        basename : str, optional
            The base name for the script file. If not provided, a default name is used.
        packages : list, optional
            A list of Tcl packages to be imported in the script.
        """
        super().newscript(basename=basename)
        self.packages=packages
        for p in self.packages:
            self.addline(f'package require {p}')
            self.addline(f'namespace import {p}::*')

    def usescript(self,scriptbasename,local=False,add_banners=False):
        """
        Ingest a VMD Tcl script file into the current script.
        This method reads the specified script file and adds its contents to the current script.

        Parameters
        ----------
        scriptbasename : str
            The base name of the script file to be ingested (without extension).
        local : bool, optional
            If True, the script is searched for in the current directory; otherwise, it is searched for in the Tcl script path.
        add_banners : bool, optional
            If True, banners are added to the script indicating the start and end of the ingested script.
        """
        if not local:
            scriptname=os.path.join(self.tcl_script_path,f'{scriptbasename}.tcl')
        else:
            scriptname=os.path.join('.',f'{scriptbasename}.tcl')
        timestampstr=datetime.datetime.today().ctime()
        if not os.path.exists(scriptname):
            raise FileNotFoundError(f'Pestifer script {scriptbasename}.tcl is not found.')
        if add_banners:
            self.banner(f'Begin {scriptbasename}, {timestampstr}')
        self.ingest_file(scriptname)
        if add_banners:
            self.banner(f'End {scriptbasename}')

    def writescript(self,force_exit=False):
        """
        Finalize the VMD script by adding an exit statement and banners.
        This method checks if an exit statement is already present in the script; if not, it adds one.

        Parameters
        ----------
        force_exit : bool, optional
            If True, the exit statement is added even if one is already present.
        """
        if not self.has_statement('exit') or force_exit:
            self.addline('exit')
        self.banner(f'END {__package__.upper()} VMD SCRIPT')
        self.banner(f'Thank you for using {__package__}!')
        super().writescript()

    def set_molecule(self,mol,altcoords=None):
        """
        Set up a VMD molecule from the given molecule object.
    
        Parameters
        ----------
        mol : Molecule
            The molecule object to be set up in VMD.
        altcoords : str, optional
            Alternative coordinates file to be loaded for the molecule.
        """
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
        """
        Backup a selection of atoms in VMD to a data holder variable.
        This method creates a backup of the specified selection by storing its attributes in a data holder variable.
        
        Parameters
        ----------
        selname : str
            The name of the selection to be backed up.
        dataholder : str, optional
            The name of the variable where the backup data will be stored. Default is ``data``.
        attributes : list, optional
            A list of attributes to be included in the backup. Default includes ``chain``, ``x``, ``y``, ``z``, ``resid``, ``resname``, and ``name``.
        """

        self.addline(f'set {dataholder} [ backup ${selname} [ list {" ".join(attributes)} ] ]')
    
    def restore_selection(self,selname,dataholder='data',attributes=['chain','x','y','z','resid','resname','name']):
        """
        Restore a selection of atoms in VMD from a data holder variable.
        This method restores the specified selection by retrieving its attributes from the data holder variable.

        Parameters
        ----------
        selname : str
            The name of the selection to be restored.
        dataholder : str, optional
            The name of the variable from which the backup data will be retrieved. Default is ``data``.
        attributes : list, optional
            A list of attributes to be restored. Default includes ``chain``, ``x``, ``y``, ``z``, ``resid``, ``resname``, and ``name``.
        """
        self.addline(f'restore ${selname} [ list {" ".join(attributes)} ]  ${dataholder}')

    def load_psf_pdb(self,*objs,new_molid_varname='mX'):
        """
        Load a PSF and PDB file into VMD.
        This method takes one or two objects (PSF and PDB files) and loads them into VMD.
        If only one object is provided, it is assumed to be the basename of both the PSF and PDB files.
        If two objects are provided, the first is treated as the PSF file and the second as the PDB file.

        Parameters
        ----------
        objs : tuple
            A tuple containing either one or two objects:

            - If one object, it should be the basename of both the PSF and PDB files.
            - If two objects, the first should be the PSF file and the second should be the PDB file.
        
        new_molid_varname : str, optional
            The variable name to be used for the new molecule ID in VMD. Default is 'mX'.
        """
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
        """
        Write the current VMD molecule to a PDB file.
        This method exports the current molecule in VMD to a PDB file with the specified basename.

        Parameters
        ----------
        basename : str
            The base name for the PDB file to be created.
        molid_varname : str
            The variable name of the molecule ID in VMD. This is used to select the molecule
            for writing to the PDB file.
        """
        self.addline(f'set TMP [atomselect ${molid_varname} all]')
        self.addline(f'$TMP writepdb {basename}.pdb')
        self.addline(f'$TMP delete')
        
    def runscript(self,*args,**options):
        """
        Run the VMD script using the VMD command line interface.
        This method constructs a command to execute VMD with the specified script and options.
        
        Parameters
        ----------
        args : tuple
            Additional arguments to be passed to the VMD command.
        options : dict
            A dictionary of options to be passed to the VMD command.
        
        Returns
        -------
        Command
            The command object that was constructed to run the VMD script, after it has been executed.
        """
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
        """
        Perform post-execution clean-up by flushing the file collector.
        If ``cleanup`` is True, it will remove all files collected during the script execution.

        Parameters
        ----------
        cleanup : bool, optional
            If True, the method will remove all files collected during the script execution. Default is False
        """
        if cleanup:
            nremoved=len(self.F)
            self.F.flush()
            logger.info(f'Post-execution clean-up: {nremoved} files removed.')

class PsfgenScripter(VMDScripter):
    """
    This class extends the VMDScripter class to provide functionality for creating and managing psfgen scripts.

    Parameters
    ----------
    config : Config
        The configuration object containing settings for the script.
    """
    def __init__(self,config):
        super().__init__(config)
        self.charmmff=config.RM.charmmff_content
        self.charmmff_config=config['user']['charmmff']
        self.psfgen_config=config['user']['psfgen']
        self.postregencommands=ByteCollector()

    def addpostregenerateline(self,line):
        """
        Add a line to the post-regenerate commands section of the psfgen script.

        Parameters
        ----------
        line : str
            The line of Tcl code to be added to the post-regenerate commands.
        """
        self.postregencommands.addline(line)

    def fetch_standard_charmm_topologies(self):
        """
        Fetch the standard CHARMM topologies from the configuration and copy them to the local directory.
        This method retrieves the topologies defined in the CHARMM force field configuration and copies them to 
        the local directory for use in the psfgen script.

        Returns
        -------
        list
            A list of topology files that have been copied to the local directory.
        """
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
        """
        Initialize a new psfgen script with a specified basename and packages.
        If no basename is provided, a default script name is used.

        Parameters
        ----------
        basename : str, optional
            The base name for the script file. If not provided, a default name is used.
        packages : list, optional
            A list of Tcl packages to be imported in the script.
        additional_topologies : list, optional
            A list of additional topology files to be included in the script. These files should not contain
            a path (i.e., they should be in the current working directory).
        """
        super().newscript(basename=basename,packages=packages)
        self.addline('package require psfgen')
        self.addline('psfcontext mixedcase')
        self.topologies=self.fetch_standard_charmm_topologies()
        for at in sorted(additional_topologies):
            assert os.sep not in at,f'Topology file {at} must not contain a path.'
            if not at in self.topologies:
                self.charmmff.copy_charmmfile_local(at)
                self.topologies.append(at)
        # self.topologies=list(set(self.topologies))  # remove duplicates
        # if 'top_all35_ethers.rtf' in self.topologies:
        #     self.topologies.remove('top_all35_ethers.rtf')
        #     self.topologies.append('top_all35_ethers.rtf')
        for t in self.topologies:
            self.addline(f'topology {t}')
        # logger.debug(f'psfgen aliases: {self.psfgen_config["aliases"]}')
        for alias_type,alias_list in self.psfgen_config['aliases'].items():
            logger.debug(f'Adding {len(alias_list)} {alias_type} aliases to psfgen script')
            for pdba in alias_list:
                alias_tokens=pdba.split()
                if alias_tokens[1]=='*': # wild-card for all protein residues
                    for protein_resname in self.psfgen_config['segtypes']['protein']['resnames']:
                        this_pdba=f'{alias_tokens[0]} {protein_resname} {alias_tokens[2]} {alias_tokens[3]}'
                        self.addline(f'pdbalias {alias_type} {this_pdba}')
                else:
                    self.addline(f'pdbalias {alias_type} {pdba}')

    def atomselect_macros(self):
        """
        Converts all resnames in the base config by segtype into updated atomselect macros
        This is usually only done one time, or whenever developer updates base config
        The subcommand ``inittcl`` is the only way this method is ever called (development use only)
        """
        for segtypename,segtyperec in self.psfgen_config['segtypes'].items():
            logger.debug(f'Searching base record of segtype {segtypename} for resname atomselect macro')
            logger.debug(', '.join(segtyperec.keys()))
            if segtyperec.get('macro',False)==True:
                new_macro='resname '+' '.join(segtyperec['resnames'])
                self.addline(f'update_atomselect_macro {segtypename} "{new_macro}" 0')

    def load_project(self,*objs):
        """
        Load a project into psfgen by reading the PSF and PDB files.
        This method takes one or two objects (PSF and PDB files) and loads them into psfgen.
        If only one object is provided, it is assumed to be the basename of both the PSF and PDB files.
        If two objects are provided, the first is treated as the PSF file and the second as the PDB file.
        
        Parameters
        ----------
        objs : tuple
            A tuple containing either one or two objects:

            - If one object, it should be the basename of both the PSF and PDB files.
            - If two objects, the first should be the PSF file and the second should be the PDB file.
        
        """
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
        """
        Describe a molecule in the psfgen script.
        This method sets up the molecule in the psfgen script by writing the necessary Tcl commands
        to create a new molecule and load its topology and coordinates.
        
        Parameters
        ----------
        mol : Molecule
            The molecule object to be described in the psfgen script.
        """
        self.molid_varname=f'm{mol.molid}'
        self.addline(f'mol top ${self.molid_varname}')
        mol.write_TcL(self)
    
    def writescript(self,statename,guesscoord=True,regenerate=True,writepsf=True,writepdb=True,force_exit=False):
        """
        Finalize the psfgen script by adding commands to write the PSF and PDB
        files, and optionally regenerate angles and dihedrals, and guess coordinates.
        This method checks if an exit statement is already present in the script; if not, it adds one.

        Parameters
        ----------
        statename : str
            The name of the state to be used for writing the PSF and PDB files.
        guesscoord : bool, optional
            If True, the script will include a command to guess coordinates. Default is True.
        regenerate : bool, optional
            If True, the script will include a command to regenerate angles and dihedrals. Default is True.
        writepsf : bool, optional
            If True, the script will include a command to write the PSF file. Default is True.
        writepdb : bool, optional
            If True, the script will include a command to write the PDB file. Default is True.
        force_exit : bool, optional
            If True, the script will include an exit statement even if one is already present. Default is False.
        """
        if guesscoord:
            self.addline('guesscoord')
        if regenerate:
            self.addline('regenerate angles dihedrals')
        if len(self.postregencommands.byte_collector)>0:
            self.B.write(self.postregencommands.byte_collector)
        if writepsf:
            self.addline(f'writepsf cmap {statename}.psf')
        if writepdb:
            self.addline(f'writepdb {statename}.pdb')
        super().writescript(force_exit=force_exit)

    def runscript(self,*args,**options):
        """
        Run the psfgen script using the VMD command line interface.
        This method constructs a command to execute VMD with the specified script and options.
        
        Parameters
        ----------
        args : tuple
            Additional arguments to be passed to the VMD command.
        options : dict
            A dictionary of options to be passed to the VMD command.
        
        Returns
        -------
        Command
            The command object that was constructed to run the psfgen script, after it has been executed.
        """
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
        else:
            logger.debug('Progress bar is disabled for psfgen script')
        return c.run(logfile=self.logname,logparser=self.logparser)
    
class NAMDScripter(TcLScripter):
    """
    This class extends the TcLScripter class to provide functionality for creating and managing NAMD scripts (which NAMD refers to as "configurations" -- this is not to be confused with Pestifer's Config class).

    Parameters
    ----------
    config : Config
        The configuration object containing settings for the script.
    """
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
        """
        Fetch the standard CHARMM parameters from the configuration and copy them to the local directory.
        This method retrieves the parameters defined in the CHARMM force field configuration and copies them to
        the local directory for use in the NAMD script.
        
        Returns
        -------
        list
            A list of parameter files that have been copied to the local directory.
        """
        logger.debug('Fetching standard CHARMM parameters')
        parameters_local=[]
        for t in self.charmmff_config['standard']['parameters']:
            self.charmmff.copy_charmmfile_local(t)
            parameters_local.append(t)
        for t in self.charmmff_config['custom']['parameters']:
            if t not in parameters_local:
                self.charmmff.copy_charmmfile_local(t)
                parameters_local.append(t)
        logger.debug(f'local parameters: {parameters_local}')
        return parameters_local

    def newscript(self,basename=None,addl_paramfiles=[]):
        """
        Initialize a new NAMD script with a specified basename and additional parameter files.
        If no basename is provided, a default script name is used.

        Parameters
        ----------
        basename : str, optional
            The base name for the script file. If not provided, a default name is used.
        addl_paramfiles : list, optional
            A list of additional parameter files to be included in the script. These files should not contain
            a path (i.e., they should be in the current working directory).
        """
        super().newscript(basename)
        self.scriptname=f'{basename}{self.default_ext}'
        self.banner('NAMD script')
        self.parameters=self.fetch_standard_charmm_parameters()
        for at in sorted(addl_paramfiles):
            if not at in self.parameters:
                self.charmmff.copy_charmmfile_local(at)
                self.parameters.append(at)
        for p in self.parameters:
            assert os.sep not in p
            self.addline(f'parameters {p}')

    def writescript(self,params,cpu_override=False):
        """
        Write the NAMD script based on the provided parameters.
        This method constructs the NAMD script by adding lines based on the parameters provided.
        If the NAMD type is 'gpu' and `cpu_override` is False, it will also add GPU-specific configurations.
        
        Parameters
        ----------
        params : dict
            A dictionary of parameters to include in the NAMD script.
        cpu_override : bool, optional
            If True, the script will be written with CPU-specific configurations, even if the NAMD type is 'gpu'.
        """
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
        """
        Run the NAMD script using the NAMD command line interface.
        This method constructs a command to execute NAMD with the specified script and options.

        Parameters
        ----------
        kwargs : dict
            A dictionary of options to be passed to the NAMD command. This can include:
            - ``local_execution_only``: If True, use the local CPU count instead of the configured CPU count.
            - ``single_gpu_only``: If True, use only one GPU device.
            - ``cpu_override``: If True, force the use of CPU settings even if the NAMD type is ``gpu``.
        """
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

class PackmolScripter(Filewriter):
    """
    This class extends the Filewriter class to provide functionality for creating and managing Packmol scripts.
    
    Parameters
    ----------
    config : Config
        The configuration object containing settings for the script.
    """
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
        """
        Create a new Packmol input script.

        Parameters
        ----------
        basename : str, optional
            The base name for the script file (without extension). If not provided, a default name will be used.
        """
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
        """
        Write the Packmol input script to the file specified in the ``scriptname`` attribute.
        """
        self.writefile()

    def runscript(self,*args,**options):
        """
        Run the Packmol script using the specified shell command.
        This method constructs a command to execute Packmol with the specified script and options.
        """
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
    