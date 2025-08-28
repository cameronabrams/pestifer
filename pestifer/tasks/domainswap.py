# Author: Cameron F. Abrams, <cfa22@drexel.edu>

"""
Definition of the :class:`DomainSwapTask` class for performing domain swaps in molecular dynamics simulations.
This class is a descendant of the :class:`MDTask <pestifer.tasks.md.MDTask>` class and is used to perform domain swaps
using the NAMD molecular dynamics engine.
It generates the necessary input files for a domain swap operation, runs NAMD to execute the swap,
and saves the resulting state of the simulation.

It uses the :ref:`tcl-domainswap` Tcl script.  Usage is described in the :ref:`config_ref tasks domainswap` documentation.
"""

import logging

from .mdtask import MDTask

from ..core.artifacts import *
from ..scripters import VMDScripter

logger=logging.getLogger(__name__)

class DomainSwapTask(MDTask):
    """
    DomainSwapTask class for performing domain swaps in molecular dynamics simulations.
    """

    _yaml_header = 'domainswap'

    def do(self):
        self.make_inputs()
        self.result = self.namdrun(
            baselabel = 'domainswap-run', 
            extras = {'colvars': 'on', 'colvarsconfig': self.statevars['cv']}, 
            single_gpu_only = True)
        return self.result

    def make_inputs(self):
        specs = self.specs
        self.next_basename('domainswap-prep')
        vm: VMDScripter = self.get_scripter('vmd')
        vm.newscript(self.basename)
        state: StateArtifacts = self.get_current_artifact('state')
        vm.usescript('domainswap')
        vm.writescript()
        vm.runscript(
            psf = state.psf.name,
            pdb = state.pdb.name,
            swap_domain_def = ','.join(specs['swap_domain_def'].split()),
            anchor_domain_def = ','.join(specs['anchor_domain_def'].split()),
            chain_swap_pairs = ':'.join([','.join(x) for x in specs['chain_directional_swaps']]),
            force_constant = specs['force_constant'],
            target_numsteps = specs['target_numsteps'],
            cv = f'{self.basename}-cv.in',
            refpdb = f'{self.basename}-ref.pdb'
        )
        self.register(f'{self.basename}-cv', key='in', artifact_type=NAMDColvarsConfigArtifact)
        self.register(self.basename, key='tcl', artifact_type=VMDScriptArtifact)
        self.register(f'{self.basename}-ref', key='refpdb', artifact_type=PDBFileArtifact)
        self.register({self.basename}, key='log', artifact_type=LogFileArtifact)
