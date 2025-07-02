# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`CleaveTask` class for cleaving chains in a molecular structure.
This class is a descendant of the :class:`PsfgenTask <pestifer.psfgen.PsfgenTask>` class and is used
to cleave chains in a molecular structure based on specified cleavage sites.
It reads the cleavage sites from the task specifications, updates the base molecule to a state
that can be reproduced by an inferential psfgen call, and then performs the cleavage operation.
The resulting structure is then processed by the psfgen method, and the final result is saved
as a PSF file.

Usage is described in the :ref:`subs_runtasks_cleave` documentation.
"""
from .psfgen import PsfgenTask
from ..objs.cleavagesite import CleavageSite, CleavageSiteList

class CleaveTask(PsfgenTask):
    """
    CleaveTask class for cleaving chains in a molecular structure.
    This class inherits from the :class:`PsfgenTask <.psfgen.PsfgenTask>` class and is used to cleave chains
    in a molecular structure based on specified :class:`cleavage sites <pestifer.objs.cleavagesite.CleavageSite>`.
    Because cleaving chains is a topological operation, ``psfgen`` is required, and the system's PSF file is updated.
    
    Cleavage sites are listed in the ``sites`` entry of the task specifications.

    Example
    -------
    ::

        tasks:
          - (prior tasks...)
          - cleave:
              sites:
                - A:10-11   
                - B:5-6
          - (subsequent tasks...)
    """
    
    yaml_header='cleave'
    """
    YAML header for the CleaveTask, used to identify the task in configuration files as part of a ``tasks`` list.
    """

    def do(self):
        """
        Execute the cleave task.
        This method reads the cleavage sites from the task specifications, updates the base molecule to a state
        that can be reproduced by an inferential psfgen call, and then performs the cleavage operation(s).
        The resulting structure is then processed by the psfgen method, and the final result is saved as a PSF/PDB fileset.
        """
        self.log_message('initiated')
        self.inherit_state()
        cleavage_sites=CleavageSiteList([CleavageSite(x) for x in self.specs['sites']])
        # update base molecule to the point an inferential psfgen call could reproduce it, up to ssbonds and links
        self.base_molecule=self.statevars['base_molecule']
        self.update_molecule()
        self.base_molecule.cleave_chains(cleavage_sites)
        self.result=self.psfgen()
        # self.save_state(exts=['psf','pdb']) # already done in psfgen()
        self.log_message('complete')
        return self.result