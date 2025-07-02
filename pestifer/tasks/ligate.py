# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`LigateTask` class for ligating loops in molecular dynamics simulations.
This class is a descendant of the :class:`MDTask <pestifer.tasks.md.MDTask>` class and is used to ligate loops
in a molecular structure using the NAMD molecular dynamics engine.
It measures the distances between loop termini, steers them toward each other, and connects them using a specified patch.
The resulting structure is then saved as a PSF/PDB files.

Usage is described in the :ref:`subs_runtasks_ligate` documentation.
"""
import logging
from .md import MDTask
from ..util.namdcolvars import declare_distance_cv_atoms, declare_single_harmonic_distance_bias

logger=logging.getLogger(__name__)

class LigateTask(MDTask):
    """
    LigateTask class for ligating loops in molecular dynamics simulations.
    """
    yaml_header='ligate'
    """
    YAML header for the LigateTask, used to identify the task in configuration files as part of a ``tasks`` list.
    """
    def do(self):
        """
        Execute the ligate task. This method checks if the base molecule has loops,
        measures the distances between loop termini, steers them toward each other,
        and connects them using a specified patch. The resulting structure is saved as a PSF/PDB fileset.
        It also writes the gaps to a data file and measures the distances between the termini of the loops.
        If the base molecule does not have loops, the task is bypassed.
        The method returns the result of the NAMD run, which is 0 on success or a non-zero error code on failure.
        If the task is successful, it saves the state of the simulation with the specified extensions.
        If the task is bypassed, it logs a message and returns without performing any operations.
        """
        self.log_message('initiated')
        self.inherit_state()
        self.base_molecule=self.statevars['base_molecule']
        if not self.base_molecule.has_loops(min_loop_length=self.statevars['min_loop_length']):
            self.log_message('bypassed')
            return
        logger.debug('Storing sequence gaps.')
        self.write_gaps()
        logger.debug('Measuring gap distances.')
        steering_specs=self.specs.get('steer',{})
        if not steering_specs:
            logger.debug(f'No steering specifications for ligate task; this is a bug; bypassing')
            return
        self.measure_distances(steering_specs)
        logger.debug('Steering loop C-termini toward their partner N-termini')
        self.result=self.do_steered_md(self.specs['steer'])
        if self.result!=0:
            return self.result
        self.save_state(exts=['coor','vel'])
        logger.debug('Connecting loop C-termini to their partner N-termini')
        connect_specs=self.specs.get('connect',{})
        self.result=self.connect(connect_specs)
        if self.result!=0:
            return self.result
        self.save_state(exts=['psf','pdb'])
        self.log_message('complete')
        return self.result
    
    def write_gaps(self):
        """
        Write the gaps in the base molecule to a data file.
        """
        self.next_basename('gaps')
        mol=self.base_molecule
        datafile=f'{self.basename}.inp'
        writer=self.scripters['data']
        writer.newfile(datafile)
        mol.write_gaps(writer)
        writer.writefile()
        self.update_statevars('data',datafile,vtype='file')

    def measure_distances(self,specs):
        """
        Measure the distances between loop termini.
        
        Parameters
        ----------
        specs : dict
            Specifications for the measurement, including the radius of the flexible zone around the receiver.
            This method uses the VMD scripter to create a script that measures the distances between
            the termini of the loops in the base molecule. It generates a data file containing the distances
            and saves the results in a specified output file. The method also updates the state variables
            with the results and the fixed reference structure.
        """
        comment_chars='#!$'
        self.next_basename('measure')
        vm=self.scripters['vmd']
        vm.newscript(self.basename)
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        datafile=self.statevars['data']
        opdb=f'{self.basename}.pdb'
        receiver_flexible_zone_radius=specs.get('receiver_flexible_zone_radius',0.0)
        resultsfile=f'{self.basename}.dat'
        vm.addline(f'measure_bonds {psf} {pdb} {datafile} {opdb} {resultsfile} {receiver_flexible_zone_radius} ')
        vm.writescript()
        vm.runscript()
        self.update_statevars('fixedref',f'{self.basename}.pdb',vtype='file')
        self.update_statevars('results',resultsfile,vtype='file')
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
        """
        Perform steered molecular dynamics to steer the loop termini toward each other.
        """
        self.next_basename('steer')
        writer=self.scripters['data']
        writer.newfile(f'{self.basename}-cv.inp')
        for i,g in enumerate(self.gaps):
            g['colvars']=f'GAP{i:02d}'
            declare_distance_cv_atoms(g,writer)
        for i,g in enumerate(self.gaps):
            g['forceConstant']=specs['force_constant']
            g['targ_distance']=specs['target_distance']
            g['targ_numsteps']=specs['nsteps']
            declare_single_harmonic_distance_bias(g,writer)
        writer.writefile()
        savespecs=self.specs
        self.specs=specs
        result=self.namdrun(extras={        
            'fixedatoms':'on',
            'fixedatomsfile':self.statevars['fixedref'],
            'fixedatomscol': 'O',
            'colvars': 'on',
            'colvarsconfig': f'{self.basename}-cv.inp'
            },single_gpu_only=True)
        self.specs=savespecs
        return result

    def connect(self,connect_specs):
        """
        Connect the loop termini using the specified ``LINK`` patch.
        """
        logger.debug(f'Connect specs: {connect_specs} (unused)')
        self.write_connect_patches()
        result=self.connect_gaps()
        return result

    def write_connect_patches(self):
        """
        Write the connect patches to a data file.
        This method generates a data file that contains the connection patches for the loop termini.
        It uses the base molecule to write the connection patches and updates the state variables with the new data file.
        The data file is named based on the current basename, which is generated by the ``next_basename`` method.
        """
        self.next_basename('gap_patches')
        mol=self.base_molecule
        datafile=f'{self.basename}.inp'
        writer=self.scripters['data']
        writer.newfile(datafile)
        mol.write_connect_patches(writer)
        writer.writefile()
        self.update_statevars('data',datafile,vtype='file')

    def connect_gaps(self):
        """
        Connect the gaps in the loop termini.
        This method uses the psfgen scripter to create a script that connects the gaps in the loop termini
        using the specified patch. It generates a new PSF file and PDB file based on the current state of the base molecule
        and the connection patches defined in the data file. The script is then executed, and if successful,
        the resulting PSF and PDB files are saved in the current state.
        
        Returns
        -------
        int
            The result of the psfgen script execution. A return value of 0 indicates success, while any other value indicates failure.
        """
        self.next_basename('heal')
        pg=self.scripters['psfgen']
        pg.newscript(self.basename)
        CC=self.config.RM.charmmff_content
        CC.copy_charmmfile_local('pestifer.top')
        pg.addline(f'topology pestifer.top')
        patchfile=self.statevars['data']
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        pg.load_project(psf,pdb)
        pg.addline(f'source {patchfile}')
        pg.writescript(self.basename,guesscoord=True,regenerate=True)
        result=pg.runscript()
        if result==0:
            self.save_state(exts=['psf','pdb'])
        return result
    