# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`SolvateTask` class for solvating a molecular structure.
This class is a descendant of the :class:`BaseTask <pestifer.core.basetask.BaseTask>` class and is used to solvate a molecular structure
using the VMD solvate and autoionize packages.
It generates a solvated PDB and PSF file, optionally adding ions based on specified salt concentration and ion types.
The task can also handle cubic or rectangular boxes for solvation based on the provided specifications.
The task reads the input PDB and PSF files, calculates the bounding box for the solvation,
and generates the necessary Tcl commands to perform the solvation and ionization.
The resulting solvated structure is saved with the specified basename, and the state is updated with the new PSF, PDB, and XSC files.

Usage is described in the :ref:`subs_runtasks_solvate` documentation.
"""
import numpy as np

from pidibble.pdbparse import PDBParser

from .basetask import VMDTask
from ..core.artifacts import *
from ..util.util import cell_from_xsc, cell_to_xsc
from ..scripters import VMDScripter

class SolvateTask(VMDTask):
    """
    SolvateTask class for solvating a molecular structure.
    """
    _yaml_header = 'solvate'
    """
    YAML header for the SolvateTask, used to identify the task in configuration files as part of a ``tasks`` list.
    This header is used to declare SolvateTask objects in YAML task lists.
    """
    def do(self):
        """
        Execute the solvate task.
        This method initializes the task, inherits the state from prior tasks, and performs the solvation and ionization.
        """
        self.next_basename()
        state: StateArtifacts = self.get_current_artifact('state')
        # psf: Path = self.get_current_artifact_path('psf')
        # pdb: Path = self.get_current_artifact_path('pdb')
        # xsc: Path = self.get_current_artifact_path('xsc')
        # self.stash_current_artifact('vel')
        use_minmax = True
        if state.xsc is not None:
            box, origin = cell_from_xsc(state.xsc.name)
            if box is not None and origin is not None:
                use_minmax = False
                basisvec = box.diagonal()
                LL = origin - 0.5 * basisvec
                UR = origin + 0.5 * basisvec
            xsc = state.xsc
        if use_minmax:
            p = PDBParser(filepath=state.pdb.name).parse()
            x = np.array([a.x for a in p.parsed['ATOM']])
            y = np.array([a.y for a in p.parsed['ATOM']])
            z = np.array([a.z for a in p.parsed['ATOM']])
            minmax = np.array([[x.min(), y.min(), z.min()], [x.max(), y.max(), z.max()]])
            spans = minmax[1] - minmax[0]
            maxspan = spans.max()
            cubic = self.specs.get('cubic', False)
            sympad = np.array([0., 0., 0.])
            if cubic:
                sympad = 0.5 * (maxspan - spans)
            pad = self.specs.get('pad', 10.0)
            LL = minmax[0] - pad - sympad
            UR = minmax[1] + pad + sympad
            basisvec = UR - LL
            origin = 0.5 * (UR + LL)
            cell_to_xsc(np.diag(basisvec), origin, f'{self.basename}.xsc')
            xsc = NAMDXscFileArtifact(self.basename)

        ll_tcl = r'{ ' + ' '.join([str(_) for _ in LL.tolist()]) + r' }'
        ur_tcl = r'{ ' + ' '.join([str(_) for _ in UR.tolist()]) + r' }'
        box_tcl = r'{ ' + ll_tcl + ' ' + ur_tcl + r' }'
        vt: VMDScripter = self.scripters['vmd']
        vt.newscript(self.basename)
        vt.addline( 'package require solvate')
        vt.addline( 'package require autoionize')
        vt.addline( 'psfcontext mixedcase')
        vt.addline(f'mol new {state.psf.name}')
        vt.addline(f'mol addfile {state.pdb.name} waitfor all')
        vt.addline(f'solvate {state.psf.name} {state.pdb.name} -minmax {box_tcl} -o {self.basename}_solv')
        ai_args=[]
        sc=self.specs.get('salt_con',None)
        if sc is None:
            ai_args.append('-neutralize')
        else:
            ai_args.append(f'-sc {sc}')
        cation=self.specs.get('cation',None)
        if cation is not None:
            ai_args.append(f'-cation {cation}')
        anion=self.specs.get('anion',None)
        if anion is not None:
            ai_args.append(f'-anion {anion}')
        vt.addline(f'autoionize -psf {self.basename}_solv.psf -pdb {self.basename}_solv.pdb {" ".join(ai_args)} -o {self.basename}')
        vt.writescript()
        self.register(self.basename, key='tcl', artifact_type=VMDScriptArtifact)
        self.result = vt.runscript(progress_title='solvate', progress_color='fluorescent_blue')
        self.register(self.basename, key='log', artifact_type=VMDLogFileArtifact)
        self.register(self.basename+'_solv', key='log', artifact_type=VMDLogFileArtifact)
        if self.result != 0:
            return self.result
        psf1 = PSFFileArtifact(f'{self.basename}_solv.psf')
        pdb1 = PDBFileArtifact(f'{self.basename}_solv.pdb')
        self.register(dict(psf=psf1, pdb=pdb1, xsc=xsc), key='state', artifact_type=StateArtifacts)
        self.pdb_to_coor(f'{self.basename}.pdb')
        self.register(dict(pdb=PDBFileArtifact(f'{self.basename}', pytestable=True), coor=NAMDCoorFileArtifact(f'{self.basename}'), psf=PSFFileArtifact(f'{self.basename}', pytestable=True), xsc=xsc), key='state', artifact_type=StateArtifacts)
        return self.result