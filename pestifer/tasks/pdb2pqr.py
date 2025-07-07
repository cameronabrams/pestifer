# Author: Cameron F. Abrams, <cfa22@drexel.edu>

"""
Definition of the :class:`PDB2PQRTask` class for preparing PDB files for PQR generation using PDB2PQR.
This class is a descendant of the :class:`PsfgenTask <pestifer.tasks.psfgen.PsfgenTask>` class and is used to prepare PDB files,
generate PQR files with titration states, and apply necessary modifications to the molecular structure.

Usage is described in the :ref:`subs_runtasks_pdb2pqr` documentation.

The documentation of the pdb2pqr command is available at `readthedocs <https://pdb2pqr.readthedocs.io/en/latest/>`_.
"""

import logging
import os
from ..core.command import Command
from ..core.objmanager import ObjManager
from ..molecule.atom import Atom, AtomList
from ..molecule.chainidmanager import ChainIDManager
from .psfgen import PsfgenTask
from ..objs.patch import Patch
from ..objs.mutation import Mutation
from ..util.logparsers import PDB2PQRLog
from ..util.progress import PDB2PQRProgress
from pidibble.pdbparse import PDBParser

logger=logging.getLogger(__name__)

class PDB2PQRTask(PsfgenTask):
    """
    PDB2PQRTask class for preparing PDB files for PQR generation using PDB2PQR.
    """
    yaml_header='pdb2pqr'
    """
    YAML header for the PDB2PQRTask, used to identify the task in configuration
    files as part of a ``tasks`` list.
    """
    def do(self):
        """
        Execute the PDB2PQR task.
        """
        self.log_message('initiated')
        self.inherit_state()
        self.next_basename()
        self.base_molecule=self.statevars['base_molecule']
        pH=self.specs.get('pH',7.0)
        xsc=self.statevars.get('xsc',None)

        prep_result=self.prep_input()
        self.run_pdb2pqr(pH=pH)
        self.prep_mods()
        self.chainIDmanager=ChainIDManager()
        self.update_molecule()
        self.psfgen()
        self.save_state(exts=['psf','pdb'])
        self.log_message('complete')

    def prep_input(self):
        """
        Prepare the input PDB file for PDB2PQR processing.
        This method writes a hydrogen-free PDB file with histidine resnames reverted to HIS and
        water resnames changed to HOH. It uses the VMD scripter to create a script that performs these modifications.
        The resulting PDB file is saved with the basename followed by `_pprep.pdb`.
        It also updates the state variables with the new PDB file.
        The method returns the result of the VMD script execution.
        
        Returns
        -------
        int
            The result of the VMD script execution, typically 0 on success.
        """
        # writes a hydrogen-free PDB file with histidine resnames all reverted to HIS and water resnames to HOH
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        vt=self.scripters['vmd']
        vt.newscript(self.basename)
        vt.addline(f'mol new {psf}')
        vt.addline(f'mol addfile {pdb} waitfor all')
        vt.addline('set a [atomselect top noh]')
        vt.addline('set w [atomselect top water]')
        vt.addline('set h [atomselect top "resname HSD HSE HSP HIS"]')
        vt.addline('set z [atomselect top "resname ZN2"]')
        vt.addline('$w set resname HOH')
        vt.addline('$h set resname HIS')
        vt.addline('$z set resname ZN')
        vt.addline(f'$a writepdb {self.basename}_pprep.pdb')
        vt.writescript()
        result=vt.runscript()
        return result
    
    def run_pdb2pqr(self,pH=7.0):
        """
        Run PDB2PQR on the prepared PDB file to generate a PQR file
        with titration states. This method constructs a command to run PDB2PQR with the specified parameters,
        including the force field, output force field, pH, titration state method, and output file names.
        It uses the `Command` class to execute the PDB2PQR command and parse the results using the `PDB2PQRLog` parser.
        The method also updates the metadata with the protonation states of histidines based on the pKa values.
        
        Parameters
        ----------
        pH : float, optional
            The pH value to be used for titration state calculations. Default is 7.
        """
        # runs PDB2PQR on the prepped PDB file, generating a PQR file with titration states, and parses the results
        c=Command(f'pdb2pqr --ff CHARMM --ffout CHARMM --with-ph {pH} --titration-state-method propka --pdb-output {self.basename}_pqr.pdb {self.basename}_pprep.pdb {self.basename}.pqr')
        self.log_parser=PDB2PQRLog(basename=f'{self.basename}_run')
        PS=PDB2PQRProgress()
        self.log_parser.enable_progress_bar(PS)

        c.run(logfile=f'{self.basename}_run.log',log_stderr=True,logparser=self.log_parser)
        self.log_parser.metadata['pka_table']['protonated'] = self.log_parser.metadata['pka_table']['respka']>pH
        logger.debug(f'PDB2PQR run completed; pka_table:\n{self.log_parser.metadata['pka_table'].to_string()}')

        self.log_parser.metadata['histidines']={k:[] for k in ['HSD','HSE','HSP']}
        if os.path.exists(f'{self.basename}_pqr.pdb'):
            pqr=PDBParser(PDBcode=f'{self.basename}_pqr').parse().parsed
            AL=AtomList([Atom(x) for x in pqr['ATOM']])
            self.prot_deprot_dict={}
            for a in AL:
                if a.resname in ['HSD','HSE','HSP']:
                    if not a.resseqnum in self.log_parser.metadata['histidines'][a.resname]:
                        self.log_parser.metadata['histidines'][a.resname].append((a.chainID,a.resseqnum))

        logger.debug(f'Protonated/deprotonated histidines: {self.log_parser.metadata["histidines"]}')

        pka_table=self.log_parser.metadata['pka_table']
        allhismods=self.log_parser.metadata['histidines']
        # replace the resname field in the pka_table with the histidine resnames
        for hisname in ['HSD','HSE','HSP']:
            if hisname in allhismods:
                for chain_resnum in allhismods[hisname]:
                    chainID,resnum=chain_resnum
                    pka_table.loc[(pka_table['resnum']==resnum) & (pka_table['reschain']==chainID),'resname']=hisname

    def prep_mods(self):
        """
        Prepare modifications based on the PKA table from the log parser.
        This method iterates through the PKA table and creates `Patch` objects for residues that require modifications.
        It handles different residue types such as histidines, acidic residues (ASP, GLU), basic residues (TYR, SER, LYS, ARG), and terminal residues (N+, C-).
        The patches are created based on the protonation state of the residues and their resnames.
        The patches are stored in the `objmanager` for later use.
        """
        self.objmanager=ObjManager()
        # columns
        # resname  resnum reschain  respka  resmodelpka resatomtype  protonated
        for i,row in self.log_parser.metadata['pka_table'].iterrows():
            patchname=None
            mut_resname=None
            if row['resname'] in ['HSE','HSP']: # HSD is the default histidine!
                # use a mutation to change HSD to HSE or HSP
                if row['resname']=='HSE':
                    mut_resname='HSE'
                elif row['resname']=='HSP':
                    mut_resname='HSP'
            elif row['resname'] in ['ASP','GLU']: # Acidic residues
                if row['protonated']:
                    if row['resname']=='ASP':
                        patchname='ASPP'
                    elif row['resname']=='GLU':
                        patchname='GLUP'
                    else:
                        patchname=None
            elif row['resname'] in ['TYR','SER','LYS','ARG']: # basic residues
                if not row['protonated']: # deprotonate!
                    if row['resname']=='TYR':
                        patchname='TYRO' # provided in pestifer.top
                    elif row['resname']=='SER':
                        patchname='SERD'
                    elif row['resname']=='LYS':
                        patchname='LSN'
                    elif row['resname']=='ARG':
                        patchname='RN2'
            elif row['resname'] == 'N+':
                if not row['protonated']:
                    patchname='NNEU'
            elif row['resname'] == 'C-':
                if row['protonated']:
                    patchname='CNEU'

            if not patchname and not mut_resname:
                logger.debug(f'No patch or mutation will be applied for {row["resname"]} at {row["reschain"]}:{row["resnum"]}')
                continue
            if patchname:
                shortcode=f"{patchname}:{row['reschain']}:{row['resnum']}"
                pat=Patch(shortcode)
                self.objmanager.ingest(pat)
            if mut_resname:
                shortcode=f"{row['reschain']}:{row['resname']},{row['resnum']},{mut_resname}"
                mut=Mutation(shortcode)
                self.objmanager.ingest(mut)