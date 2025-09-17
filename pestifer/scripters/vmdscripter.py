# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import datetime
import logging
import os

from .tclscripter import TcLScripter

from ..core.command import Command

from ..logparsers import VMDLogParser

from ..molecule.molecule import Molecule
from ..molecule.residue import ResidueList

from ..objs.crot import Crot, CrotList
from ..objs.orient import Orient, OrientList
from ..objs.rottrans import RotTrans, RotTransList

from ..util.util import reduce_intlist
from ..util.progress import PestiferProgress

logger = logging.getLogger(__name__)

class VMDScripter(TcLScripter):
    """
    This class extends the TcLScripter class to provide functionality for creating and managing VMD scripts.

    Parameters
    ----------
    config : Config
        The configuration object containing settings for the script.
    """ 

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.vmd = kwargs.get('vmd')
        self.tcl_root = kwargs.get('tcl_root')
        self.tcl_pkg_path = kwargs.get('tcl_pkg_path')
        self.tcl_script_path = kwargs.get('tcl_script_path')
        self.vmd_startup = kwargs.get('vmd_startup_script')
        self.indent = ' ' * 4

    def newscript(self, basename=None, packages=[]):
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
        self.packages = packages
        for p in self.packages:
            self.addline(f'package require {p}')
            self.addline(f'namespace import {p}::*')

    def usescript(self, scriptbasename, local=False, add_banners=False):
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
            scriptname = os.path.join(self.tcl_script_path, f'{scriptbasename}.tcl')
        else:
            scriptname = os.path.join('.', f'{scriptbasename}.tcl')
        timestampstr = datetime.datetime.today().ctime()
        if not os.path.exists(scriptname):
            raise FileNotFoundError(f'Pestifer script {scriptbasename}.tcl is not found.')
        if add_banners:
            self.banner(f'Begin {scriptbasename}, {timestampstr}')
        self.ingest_file(scriptname)
        if add_banners:
            self.banner(f'End {scriptbasename}')

    def writescript(self, force_exit=False):
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

    def set_molecule(self, mol: Molecule, altcoords=None):
        """
        Set up a VMD molecule from the given molecule object.
    
        Parameters
        ----------
        mol : Molecule
            The molecule object to be set up in VMD.
        altcoords : str, optional
            Alternative coordinates file to be loaded for the molecule.
        """
        self.molid_varname = f'm{mol.molid}'
        ext = '.pdb' if mol.rcsb_file_format == 'PDB' else '.cif'
        if mol.sourcespecs.get('source_id', {}):
            self.addline(f'mol new {mol.sourcespecs["source_id"]}{ext} waitfor all')
        elif mol.sourcespecs.get('prebuilt', {}):
            pdb = mol.sourcespecs['prebuilt']['pdb']
            self.addline(f'mol new {pdb} waitfor all autobonds off')
        else:
            raise ValueError(f'Molecule specs has neither "source_id" nor "prebuilt" specified')
        self.addline(f'set {self.molid_varname} [molinfo top get id]')
        self.addline(f'set nf [molinfo ${self.molid_varname} get numframes]')
        self.addline(r'if { $nf > 1 } { animate delete beg 0 end [expr $nf - 2] $'+self.molid_varname+r' }')
        if altcoords:
            mol.set_coords(altcoords)
            self.addline(f'mol addfile {altcoords} waitfor all')
            self.addline(f'set nf [molinfo ${self.molid_varname} get numframes]')
            self.addline(r'animate delete beg 0 end [expr $nf - 2] $'+self.molid_varname+r' }')
        if mol.rcsb_file_format == 'mmCIF':
            # VMD appends a "1" to any two-letter chain ID from a CIF file,
            # so let's undo that
            # also, CIF HETATM records have a "." for residue number, which
            # VMD interprets as 0, so we'll replace that with the residue numbers
            # assigned by pestifer when it reads in the cif file
            au = mol.asymmetric_unit
            residues = au.segments.collect_residues()
            uCIDs = residues.uniqattrs(['chainID'])['chainID']
            self.comment('Resetting chains and resids for this CIF-source molecule')
            for c in uCIDs:
                chain: ResidueList = residues.filter(lambda x: x.chainID == c)
                resids = []
                for x in chain.data:
                    resids.extend([str(y.resid.resid) for y in x.atoms])
                residlist = ' '.join(resids)
                serials = chain.atom_serials(as_type=int)
                vmd_red_list = reduce_intlist(serials)
                self.addline(f'set TMP [atomselect ${self.molid_varname} "serial {vmd_red_list}"]')
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
        if len(objs) == 1:
            basename = objs[0]
            pdb = f'{basename}.pdb'
            psf = f'{basename}.psf'
        else:
            psf, pdb = objs
        self.addline(f'mol new {psf}')
        self.addline(f'mol addfile {pdb} waitfor all')
        self.addline(f'set {new_molid_varname} [molinfo top get id]')
        self.molid_varname = new_molid_varname

    def write_pdb(self, basename, molid_varname):
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

    def runscript(self, *args, **options):
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
        assert hasattr(self, 'scriptname'), f'No scriptname set.'
        self.logname = f'{self.basename}.log'
        self.logparser = VMDLogParser(basename=self.basename)
        logger.debug(f'Log file: {self.logname}')
        c = Command(f'{self.vmd} -dispdev text -startup {self.vmd_startup} -e {self.scriptname} -args --tcl-root {self.tcl_root}', **options)
        progress_struct = None
        progress_title = options.get('progress_title', '')
        if self.progress and progress_title != '':
            progress_struct = PestiferProgress(name=progress_title, color=options.get('progress_color','fuchsia'))
            self.logparser.enable_progress_bar(progress_struct)
        return c.run(logfile=self.logname, logparser=self.logparser)

    def cleanup(self, cleanup=False):
        """
        Perform post-execution clean-up by flushing the file collector.
        If ``cleanup`` is True, it will remove all files collected during the script execution.

        Parameters
        ----------
        cleanup : bool, optional
            If True, the method will remove all files collected during the script execution. Default is False
        """
        if cleanup:
            logger.debug(f'Post-execution clean-up: {len(self.F)} files removed.')
            self.F.flush()

    def write_crots(self, crots: CrotList, chainIDmap: dict = {}):
        for crot in crots.data:
            self.write_crot(crot, chainIDmap=chainIDmap)

    def write_orients(self, orients: OrientList):
        for orient in orients.data:
            self.write_orient(orient)

    def write_rottranslist(self, rottranslist: RotTransList):
        for rt in rottranslist.data:
            self.write_rottrans(rt)

    def write_crot(self, crot: Crot, chainIDmap: dict = {}):
        """
        Write a CROT object to the VMD script.
        This method generates VMD commands to create a CROT object based on the provided CROT data.

        Parameters
        ----------
        crot : Crot
            The CROT object containing the data to be written to the script.
        chainIDmap : dict, optional
            A mapping of original chain IDs to new chain IDs. Default is an empty dictionary.
        """
        
        the_chainID = chainIDmap.get(crot.chainID, crot.chainID)
        molid_varname = self.molid_varname
        molid = f'${molid_varname}'
        # endIsCterm=kwargs.get('endIsCterm',True)
        if crot.angle in ['PHI', 'PSI', 'OMEGA']:
            self.addline('set r1 [[atomselect {} "chain {} and resid {} and name CA"] get residue]'.format(molid, the_chainID, crot.resid1.resid))
            self.addline('set r2 [[atomselect {} "chain {} and resid {} and name CA"] get residue]'.format(molid, the_chainID, crot.resid2.resid))
            if crot.resid1 <= crot.resid2:
                direction = 'C'
            else:
                direction = 'N'
            self.addline(f'brot {molid} $r1 $r2 {crot.angle.lower()} {direction} {crot.degrees}')
            # if endIsCterm:
            #     W.addline('Crot_{}_toCterm $r1 $r2 {} {} {}'.format(self.angle.lower(),the_chainID,molid,self.degrees))
            # else:
            #     W.addline('Crot_{} $r1 $r2 {} {} {}'.format(self.angle.lower(),the_chainID,molid,self.degrees))
        elif crot.angle in ['CHI1','CHI2']:  # this is a side-chain bond
            self.addline('set r1 [[atomselect {} "chain {} and resid {} and name CA"] get residue]'.format(molid,the_chainID,crot.resid1.resid))
            self.addline(f'brot {molid} $r1 -1 {crot.angle[:-1].lower()} {crot.angle[-1]} {crot.degrees}')
        elif crot.angle=='ANGLEIJK':
            self.addline('set rotsel [atomselect {} "segname {}"]'.format(molid,crot.segnamejk))
            self.addline('set ri [lindex [[atomselect {} "segname {} and resid {} and name {}"] get {{x y z}}] 0]'.format(molid,crot.segnamei,crot.residi.resid,crot.atomi))
            self.addline('set rj [lindex [[atomselect {} "segname {} and resid {} and name {}"] get {{x y z}}] 0]'.format(molid,crot.segnamejk,crot.residj.resid,crot.atomj))
            self.addline('set rk [lindex [[atomselect {} "segname {} and resid {} and name {}"] get {{x y z}}] 0]'.format(molid,crot.segnamejk,crot.residk.resid,crot.atomk))
            self.addline('set rij [vecsub $ri $rj]')
            self.addline('set rjk [vecsub $rj $rk]')
            self.addline('set cijk [veccross $rij $rjk]')
            self.addline('$rotsel move [trans center $rj origin $rj axis $cijk {} degrees]'.format(crot.degrees))
        elif crot.angle=='ALPHA':
            self.addline('set r1 [[atomselect {} "chain {} and resid {} and name CA"] get residue]'.format(molid,the_chainID,crot.resid1.resid))
            self.addline('set r2 [[atomselect {} "chain {} and resid {} and name CA"] get residue]'.format(molid,the_chainID,crot.resid2.resid))
            self.addline('set rterm [[atomselect {} "chain {} and resid {} and name CA"] get residue]'.format(molid,the_chainID,crot.resid3.resid))
            self.addline('fold_alpha $r1 $r2 $rterm {}'.format(molid))

    def write_orient(self, orient: Orient):
        """
        Write a VMD Orient object to the script.
        This method generates the Tcl commands to orient a coordinate set in VMD based on the provided Orient object.

        Parameters
        ----------
        orient : Orient
            The Orient object containing the orientation data.
        """
        self.addline('set a [atomselect top all]')
        self.addline('set I [draw principalaxes $a]')
        adict=dict(x=r'{1 0 0}',y=r'{0 1 0}',z=r'{0 0 1}')
        self.addline(f'set A [orient $a [lindex $I 2] {adict[orient.axis]}]')
        self.addline(r'$a move $A')
        if hasattr(orient,'refatom'):
            self.addline(r'set com [measure center $a]')
            self.addline(r'$a moveby [vecscale $com -1]')
            self.addline(f'set z [[atomselect top "name {orient.refatom}"] get z]')
            self.addline(r'if { $z < 0.0 } {')
            self.addline(r'   $a move [transaxis x 180 degrees]')
            self.addline(r'}')

    def write_rottrans(self, rottrans: RotTrans):
        molid_varname=self.molid_varname
        molid=f'${molid_varname}'
        self.addline(f'set mover [atomselect {molid} all]')
        if rottrans.movetype=='TRANS':
            self.addline(f'$mover moveby [list {rottrans.x} {rottrans.y} {rottrans.z}]')
        elif rottrans.movetype=='ROT':
            self.addline(f'set COM [measure center $mover weight mass]')
            self.addline(f'$mover move [trans origin $COM axis {rottrans.axis} {rottrans.angle}]')

    def write_protein_loop_lines(self, mol: Molecule , cycles=100, **options):
        """
        Write Tcl commands to declare loops in the asymmetric unit that are in the ``MISSING``
        state and have a length greater than or equal to the specified minimum length.
        This method iterates through the segments of the asymmetric unit and generates Tcl commands
        to declare the loops.
        
        Parameters
        ----------
        mol : Molecule
            An instance of the Molecule class that contains the structural information.
        cycles : int, optional
            The number of cycles for the loop declaration. Default is 100.
        options : dict, optional
            Additional options for loop declaration. It can include:

            - ``min_length``: The minimum length of loops to consider. Default is 4.
            - ``include_c_termini``: If True, includes C-terminal loops in the declaration. If False, skips C-terminal loops. Default is False.
        """
        ba = mol.active_biological_assembly
        au = mol.asymmetric_unit
        min_length = mol.min_loop_length
        include_c_termini = options.get('include_c_termini', False)
        for S in au.segments.data:
            # chainID=S.chainID
            if S.segtype in 'protein':
                asymm_segname = S.segname
                n_subsegs = len(S.subsegments)
                for b in S.subsegments.data:
                    is_c_terminus = (S.subsegments.index(b) == (n_subsegs - 1))
                    is_processible = b.state == 'MISSING' and b.num_items() >= min_length
                    if is_processible and (not include_c_termini) and is_c_terminus:
                        logger.debug(f'A.U. C-terminal loop {b.pstr()} declashing is skipped')
                        is_processible = False
                    if is_processible:
                        reslist = [f'{r.resid.resid}' for r in S.residues.data[b.bounds[0]:b.bounds[1] + 1]]
                        tcllist = '[list ' + ' '.join(reslist) + ']'
                        for transform in ba.transforms:
                            cm = transform.chainIDmap
                            act_segID = cm.get(asymm_segname, asymm_segname)
                            self.addline(f'declash_loop $mLL {act_segID} {tcllist} {cycles}')
