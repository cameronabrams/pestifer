# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import datetime
import logging
import os
import subprocess

from .tcl import TcLScripter

from ..core.command import Command

from ..logparsers import VMDLogParser

from ..molecule.molecule import Molecule
from ..molecule.residue import ResidueList

from ..objs.align import Align
from ..objs.transfer_coords import TransferCoords
from ..objs.crot import Crot, CrotList
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
        has_bare_exit = any(
            l.strip() == 'exit'
            for l in self.B.byte_collector.split('\n')
            if l.strip() and not l.strip().startswith('#')
        )
        if not has_bare_exit or force_exit:
            self.addline('exit')
        self.banner(f'END PESTIFER VMD SCRIPT')
        self.banner(f'Thank you for using Pestifer!')
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
        # VMD runs in pestifer's own session (not a new one): a detached session strips
        # the controlling terminal, and an rlwrap-wrapping VMD launcher then exits without
        # running the script.  Feeding /dev/null on stdin also keeps such a launcher from
        # choosing rlwrap at all (VMD -dispdev text -e … reads no stdin).
        return c.run(logfile=self.logname, logparser=self.logparser,
                     new_session=False, stdin=subprocess.DEVNULL)

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

    def write_rottranslist(self, rottranslist: RotTransList):
        for rt in rottranslist.data:
            self.write_rottrans(rt)

    def write_crot(self, crot: Crot, chainIDmap: dict = {}, molid: str = None):
        """
        Write a CROT object to the VMD script.
        This method generates VMD commands to create a CROT object based on the provided CROT data.

        Parameters
        ----------
        crot : Crot
            The CROT object containing the data to be written to the script.
        chainIDmap : dict, optional
            A mapping of original chain IDs to new chain IDs. Default is an empty dictionary.
        molid : str, optional
            The molid variable name to operate on; defaults to the scripter's current molecule.
        """

        the_chainID = chainIDmap.get(crot.chainID, crot.chainID)
        molid = f'${self.molid_varname}' if not molid else f'${molid}'
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
            logger.warning('ANGLEIJK is deprecated as an irotation (it is a rigid-body segment '
                           'rotation, not an internal-coordinate one); use the transrot '
                           'AXISANGLE movetype instead.')
            self.addline('set rotsel [atomselect {} "segname {}"]'.format(molid,crot.segnamejk))
            selI = f'segname {crot.segnamei} and resid {crot.residi.resid} and name {crot.atomi}'
            selJ = f'segname {crot.segnamejk} and resid {crot.residj.resid} and name {crot.atomj}'
            selK = f'segname {crot.segnamejk} and resid {crot.residk.resid} and name {crot.atomk}'
            self._write_axisangle_rotation('rotsel', selI, selJ, selK, crot.degrees, molid)
        elif crot.angle=='ALPHA':
            self.addline('set r1 [[atomselect {} "chain {} and resid {} and name CA"] get residue]'.format(molid,the_chainID,crot.resid1.resid))
            self.addline('set r2 [[atomselect {} "chain {} and resid {} and name CA"] get residue]'.format(molid,the_chainID,crot.resid2.resid))
            self.addline('set rterm [[atomselect {} "chain {} and resid {} and name CA"] get residue]'.format(molid,the_chainID,crot.resid3.resid))
            self.addline('fold_alpha $r1 $r2 $rterm {}'.format(molid))
        elif crot.angle=='GLYCAN_PENDANT':
            self.addline('glycan_pendant_rotate {} {} {} {} {} {} {}'.format(molid,the_chainID,crot.residi.resid,crot.residj.resid,crot.atomi,crot.atomj,crot.degrees))
        else:
            raise ValueError(f'Unknown CROT angle type: {crot.angle}')

    def write_rottrans(self, rottrans: RotTrans, molid: str = None):
        """Emit VMD commands for a rigid-body transformation of a (possibly partial) fragment.

        The fragment is named by ``rottrans.sel`` (a VMD atomselection, default ``all``).  A
        rigid-body transform is only meaningful when the fragment is fully disconnected from the
        rest of the system, so for any selection other than ``all`` a disconnection guard is
        emitted: if any selected atom is bonded to an atom outside the selection, the script prints
        a ``PESTIFER-ERROR`` and exits non-zero rather than silently stretching the crossing bond.
        Rotations are performed about the fragment's own center of mass.
        """
        molid = f'${self.molid_varname}' if not molid else f'${molid}'
        sel = rottrans.sel or 'all'
        self.addline(f'set mover [atomselect {molid} "{sel}"]')
        if sel != 'all':
            self._write_disconnection_guard(sel)
        if rottrans.movetype=='TRANS':
            self.addline(f'$mover moveby [list {rottrans.x} {rottrans.y} {rottrans.z}]')
        elif rottrans.movetype=='ROT':
            # rotate about the selection's center of mass.  VMD's `trans center {c}` rotates
            # about c; `trans origin {c}` instead *moves* c to the global origin (and rotates
            # there), which is not what we want.
            self.addline(f'set COM [measure center $mover weight mass]')
            self.addline(f'$mover move [trans center $COM axis {rottrans.axis} {rottrans.angle}]')
        elif rottrans.movetype=='ALIGN':
            self._write_align(rottrans, molid)
        elif rottrans.movetype=='AXISANGLE':
            selI, selJ, selK = rottrans.axis_atoms
            self._write_axisangle_rotation('mover', selI, selJ, selK, rottrans.angle, molid)

    def _write_axisangle_rotation(self, mover_var: str, selI: str, selJ: str, selK: str,
                                  degrees, molid: str):
        """Rotate the fragment held by the atomselect variable *mover_var* about the axis normal to
        the angle defined by three atom selections.

        The pivot is the (mass-unweighted) center of ``selJ`` and the axis is
        ``(rI - rJ) x (rJ - rK)`` -- i.e. the normal to the plane of the three points.  Each of
        ``selI``/``selJ``/``selK`` is a VMD atomselection (normally selecting a single atom).  This
        is the shared engine for the ``transrot`` ``AXISANGLE`` movetype and the deprecated
        ``ANGLEIJK`` irotation.
        """
        self.addline(f'set _ri [measure center [atomselect {molid} "{selI}"]]')
        self.addline(f'set _rj [measure center [atomselect {molid} "{selJ}"]]')
        self.addline(f'set _rk [measure center [atomselect {molid} "{selK}"]]')
        self.addline(f'set _axis [veccross [vecsub $_ri $_rj] [vecsub $_rj $_rk]]')
        self.addline(f'${mover_var} move [trans center $_rj origin $_rj axis $_axis {degrees} degrees]')

    def _write_align(self, rottrans: RotTrans, molid: str):
        """Emit the minimal (roll-free) rotation carrying ``source`` onto ``target``, about the
        fragment's center of mass.

        The rotation axis is ``source x target`` and the angle is ``atan2(|source x target|,
        source . target)`` -- the unique smallest rotation between the two directions.  Degenerate
        cases are handled: already-aligned vectors produce no rotation, and antiparallel vectors are
        flipped 180 degrees about an arbitrary perpendicular axis.
        """
        self._write_align_vector(rottrans.source, molid, '_src')
        self._write_align_vector(rottrans.target, molid, '_tgt')
        self.addline(f'set _src [vecnorm $_src]')
        self.addline(f'set _tgt [vecnorm $_tgt]')
        self.addline(f'set _cross [veccross $_src $_tgt]')
        self.addline(f'set _sin [veclength $_cross]')
        self.addline(f'set _cos [vecdot $_src $_tgt]')
        self.addline(f'set COM [measure center $mover weight mass]')
        self.addline(f'if {{ $_sin < 1.0e-6 }} {{')
        self.addline(f'  if {{ $_cos < 0.0 }} {{')
        self.addline(f'    set _perp [veccross $_src {{1.0 0.0 0.0}}]')
        self.addline(f'    if {{ [veclength $_perp] < 1.0e-6 }} {{ set _perp [veccross $_src {{0.0 1.0 0.0}}] }}')
        self.addline(f'    $mover move [trans center $COM axis [vecnorm $_perp] 180.0]')
        self.addline(f'  }}')
        self.addline(f'}} else {{')
        self.addline(f'  set _angle [expr {{atan2($_sin, $_cos) * 180.0 / 3.14159265358979}}]')
        self.addline(f'  $mover move [trans center $COM axis [vecnorm $_cross] $_angle]')
        self.addline(f'}}')

    def _write_align_vector(self, spec: list, molid: str, varname: str):
        """Set Tcl variable ``varname`` to an ALIGN vector.

        ``spec`` is either a literal ``[x, y, z]`` (three numbers) or a pair of atomselections
        ``[selA, selB]`` whose mass-weighted centers define the vector (from ``selA`` to ``selB``).
        """
        if len(spec) == 3 and all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in spec):
            x, y, z = spec
            self.addline(f'set {varname} [list {x} {y} {z}]')
        else:
            selA, selB = spec
            self.addline(f'set _va [atomselect {molid} "{selA}"]')
            self.addline(f'set _vb [atomselect {molid} "{selB}"]')
            self.addline(f'set {varname} [vecsub [measure center $_vb weight mass] [measure center $_va weight mass]]')
            self.addline(f'$_va delete')
            self.addline(f'$_vb delete')

    def _write_disconnection_guard(self, sel: str):
        """Emit a Tcl guard that hard-errors unless ``$mover`` is fully disconnected.

        Relies on ``$mover`` already being an atomselection.  ``getbonds`` returns one bonded-index
        list per selected atom (in selection order); if any bonded index is absent from the set of
        selected indices, a bond crosses the selection boundary and the transform would deform the
        molecule.
        """
        self.addline(f'set _in_sel [dict create]')
        self.addline(f'foreach _i [$mover get index] {{ dict set _in_sel $_i 1 }}')
        self.addline(f'foreach _bl [$mover getbonds] {{')
        self.addline(f'  foreach _b $_bl {{')
        self.addline(f'    if {{ ![dict exists $_in_sel $_b] }} {{')
        self.addline(f'      puts "PESTIFER-ERROR: transrot selection \\"{sel}\\" is not fully disconnected (a bond crosses the selection boundary); refusing to apply a rigid-body transform."')
        self.addline(f'      exit 1')
        self.addline(f'    }}')
        self.addline(f'  }}')
        self.addline(f'}}')

    def write_align(self, align: Align):
        """Write VMD commands to align the pipeline system to a reference coordinate file.

        Uses ``measure fit`` to compute the least-RMSD homogeneous transformation
        between ``mobile_sel`` atoms in the pipeline molecule and ``ref_sel`` atoms
        in the reference molecule, then applies that matrix to ``apply_to`` atoms.
        The reference molecule is deleted after the fit so it leaves no state.
        """
        molid_varname = self.molid_varname
        molid = f'${molid_varname}'
        ref_sel = align.ref_sel if align.ref_sel is not None else align.mobile_sel
        ref_pdb = align.effective_ref_pdb
        if align.ref_psf:
            self.addline(f'mol load psf {{{align.ref_psf}}} pdb {{{ref_pdb}}}')
        else:
            self.addline(f'mol load pdb {{{ref_pdb}}}')
        self.addline('set _ref_mol [molinfo top get id]')
        self.addline(f'set _mob_fit [atomselect {molid} "{align.mobile_sel}"]')
        self.addline(f'set _ref_fit [atomselect $_ref_mol "{ref_sel}"]')
        self.addline('set _n_mob [$_mob_fit num]')
        self.addline('set _n_ref [$_ref_fit num]')
        self.addline('if { $_n_mob != $_n_ref } {')
        self.addline(f'    puts "PESTIFER-ERROR: align: mobile_sel \\"{align.mobile_sel}\\" has $_n_mob atoms but ref_sel \\"{ref_sel}\\" has $_n_ref atoms"')
        self.addline('    exit 1')
        self.addline('}')
        self.addline('set _rmsd_before [measure rmsd $_mob_fit $_ref_fit]')
        self.addline(f'vmdcon -info "Align RMSD before: $_rmsd_before"')
        self.addline('set _M [measure fit $_mob_fit $_ref_fit]')
        self.addline(f'set _mover [atomselect {molid} "{align.apply_to}"]')
        self.addline('$_mover move $_M')
        self.addline('set _rmsd_after [measure rmsd $_mob_fit $_ref_fit]')
        self.addline(f'vmdcon -info "Align RMSD after: $_rmsd_after"')
        self.addline('mol delete $_ref_mol')
        self.addline('$_mob_fit delete')
        self.addline('$_ref_fit delete')
        self.addline('$_mover delete')

    def write_transfer_coords(self, tc: TransferCoords):
        """Write VMD commands to copy coordinates from a donor system onto the pipeline system.

        If ``tc.pre_align`` is true, the entire donor is first rigidly fitted to the
        pipeline system using ``align_donor_sel`` / ``align_mobile_sel`` before the
        coordinate transfer is performed.  Both the alignment and transfer selections
        are checked for congruency before use.
        """
        molid_varname = self.molid_varname
        molid = f'${molid_varname}'
        if tc.donor_psf:
            self.addline(f'mol load psf {{{tc.donor_psf}}} pdb {{{tc.donor_pdb}}}')
        else:
            self.addline(f'mol load pdb {{{tc.donor_pdb}}}')
        self.addline('set _donor_mol [molinfo top get id]')
        if tc.pre_align:
            self.addline(f'set _aln_donor [atomselect $_donor_mol "{tc.align_donor_sel}"]')
            self.addline(f'set _aln_mob [atomselect {molid} "{tc.align_mobile_sel}"]')
            self.addline('set _n_aln_donor [$_aln_donor num]')
            self.addline('set _n_aln_mob [$_aln_mob num]')
            self.addline('if { $_n_aln_donor != $_n_aln_mob } {')
            self.addline(f'    puts "PESTIFER-ERROR: transfer_coords align: align_donor_sel \\"{tc.align_donor_sel}\\" has $_n_aln_donor atoms but align_mobile_sel \\"{tc.align_mobile_sel}\\" has $_n_aln_mob atoms"')
            self.addline('    exit 1')
            self.addline('}')
            self.addline('set _aln_M [measure fit $_aln_donor $_aln_mob]')
            self.addline('set _aln_all [atomselect $_donor_mol all]')
            self.addline('$_aln_all move $_aln_M')
            self.addline('$_aln_donor delete')
            self.addline('$_aln_mob delete')
            self.addline('$_aln_all delete')
        self.addline(f'set _donor_sel [atomselect $_donor_mol "{tc.donor_sel}"]')
        self.addline(f'set _mob_sel [atomselect {molid} "{tc.mobile_sel}"]')
        self.addline('set _n_donor [$_donor_sel num]')
        self.addline('set _n_mob [$_mob_sel num]')
        self.addline('if { $_n_donor != $_n_mob } {')
        self.addline(f'    puts "PESTIFER-ERROR: transfer_coords: donor_sel \\"{tc.donor_sel}\\" has $_n_donor atoms but mobile_sel \\"{tc.mobile_sel}\\" has $_n_mob atoms"')
        self.addline('    exit 1')
        self.addline('}')
        self.addline('$_mob_sel set {x y z} [$_donor_sel get {x y z}]')
        self.addline('mol delete $_donor_mol')
        self.addline('$_donor_sel delete')
        self.addline('$_mob_sel delete')

