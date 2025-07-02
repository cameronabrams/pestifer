# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Desolvate task for the Pestifer framework.  This is normally used as part of the ``pestifer desolvate`` command to process existing DCD files.
This task is responsible for generating an index file and a PSF file from a given PSF and PDB file,
and then pruning a DCD file based on the generated index.
It inherits from the :class:`BaseTask <pestifer.core.basetask.BaseTask>` class and uses the VMD scripter to create the necessary scripts for processing.

Usage is described in the :ref:`subs_desolvate` documentation.
"""
from ..core.basetask import BaseTask
from ..core.command import Command
from ..util.progress import PestiferProgress
import os

class DesolvateTask(BaseTask):
    """
    DesolvateTask class for processing DCD files and generating index and PSF files.
    """
    yaml_header='desolvate'
    """
    YAML header for the DesolvateTask, used to identify the task in configuration files as
    part of a ``tasks`` list.  This is not the normal use case for Desolvate.
    """
    def do(self):
        self.catdcd=self.config.shell_commands['catdcd']
        self.do_idx_psf_gen()
        self.do_dcd_prune()

    def do_idx_psf_gen(self):
        """
        Generate an index file and a PSF file from the given PSF and PDB files.
        This method uses the VMD scripter to create a script that performs the necessary operations
        to generate the index file and the PSF file based on the specified selection criteria.
        The generated script is then executed to produce the output files."""
        self.log_message('initiated')
        # self.inherit_state()
        self.next_basename()
        psf=self.specs['psf']
        pdb=self.specs['pdb']
        keepatselstr=self.specs['keepatselstr']
        idx_outfile=self.specs['idx_outfile']
        psf_outfile=self.specs['psf_outfile']
        if pdb:
            pdb_outfile=os.path.splitext(psf_outfile)[0]+'.pdb'
        vt=self.scripters['vmd']
        vt.newscript(self.basename)
        vt.addline( 'package require psfgen')
        vt.addline(f'mol new {psf}')
        if pdb:
            vt.addline(f'mol addfile {pdb}')
        vt.addline(f'set keepsel [atomselect top "{keepatselstr}"]')
        vt.addline( 'set keepsegid [lsort -unique [$keepsel get segid]]')
        vt.addline(f'set dumpsel [atomselect top "not ({keepatselstr})"]')
        vt.addline(f'set dumpsegid [lsort -unique [$dumpsel get segid]]')
        vt.addline(f'vmdcon -info "Writing {idx_outfile}"')
        vt.addline(f'set fp [open "{idx_outfile}" "w"]')
        vt.addline( 'puts $fp "[$keepsel get index]"')
        vt.addline( 'close $fp')
        if pdb:
            vt.addline(f'$keepsel writepdb {pdb_outfile}')
        vt.addline( 'vmdcon -info "Keeping segids $keepsegid"')
        vt.addline( 'vmdcon -info "Dumping segids $dumpsegid"')
        vt.addline(f'readpsf {psf}')
        vt.addline(r'if {[IsIntersectionEmpty $keepsegid $dumpsegid]} {')
        vt.addline(r'   foreach badsegid $dumpsegid {')
        vt.addline( '       delatom $badsegid')
        vt.addline(r'   }')
        vt.addline(r'}')
        vt.addline(f'vmdcon -info "Writing {psf_outfile}"')
        vt.addline(f'writepsf {psf_outfile}')
        vt.writescript()
        self.result=vt.runscript(progress_title='psfidx')
        
    def do_dcd_prune(self):
        """
        Prune a DCD file based on the generated index file.
        This method uses the `catdcd` command to create a new DCD file that
        contains only the frames corresponding to the indices specified in the index file.
        The new DCD file is created with a specified stride, and the original DCD files
        are concatenated into the new DCD file.
        """
        idx_outfile=self.specs['idx_outfile']
        dcd_outfile=self.specs['dcd_outfile']
        dcd_infiles=self.specs['dcd_infiles']
        dcd_stride=self.specs['dcd_stride']
        progress_struct=PestiferProgress(name='catdcd',track_stdout=False)
        c=Command(f'{self.catdcd} -i {idx_outfile} -stride {dcd_stride} -o {dcd_outfile} {" ".join(dcd_infiles)}')
        c.run(progress=progress_struct)
