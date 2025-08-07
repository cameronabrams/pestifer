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

from pathlib import Path

from .md import MDTask

from ..core.artifacts import PDBFile, NAMDColvarsConfig, VMDScript, VMDLogFile
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
        psf: Path = self.get_current_artifact_path('psf')
        pdb: Path = self.get_current_artifact_path('pdb')
        vm.usescript('domainswap')
        vm.writescript()
        vm.runscript(
            psf = psf.name,
            pdb = pdb.name,
            swap_domain_def = ','.join(specs['swap_domain_def'].split()),
            anchor_domain_def = ','.join(specs['anchor_domain_def'].split()),
            chain_swap_pairs = ':'.join([','.join(x) for x in specs['chain_directional_swaps']]),
            force_constant = specs['force_constant'],
            target_numsteps = specs['target_numsteps'],
            cv = f'{self.basename}-cv.in',
            refpdb = f'{self.basename}-ref.pdb'
        )
        self.register_current_artifact(NAMDColvarsConfig(f'{self.basename}-cv'))
        self.register_current_artifact(VMDScript(self.basename))
        self.register_current_artifact(PDBFile(f'{self.basename}-ref'), key='refpdb')
        self.register_current_artifact(NAMDColvarsConfig(f'{self.basename}-cv'))
        self.register_current_artifact(VMDLogFile(f'{self.basename}'))
