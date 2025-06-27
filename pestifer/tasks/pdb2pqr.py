# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import logging
import os
from ..core.command import Command
from ..core.objmanager import ObjManager
from ..molecule.atom import Atom, AtomList
from ..molecule.chainidmanager import ChainIDManager
from .psfgen import PsfgenTask
from ..objs.patch import Patch
from ..util.logparsers import PDB2PQRLog
from ..util.progress import PDB2PQRProgress
from pidibble.pdbparse import PDBParser

logger=logging.getLogger(__name__)

class PDB2PQRTask(PsfgenTask):
    yaml_header='pdb2pqr'

    def do(self):
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
        # writes a hydrogen-free PDB file with histidine resnames all reverted to HIS and water resnames to HOH
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        vt=self.writers['vmd']
        vt.newscript(self.basename)
        vt.addline(f'mol new {psf}')
        vt.addline(f'mol addfile {pdb} waitfor all')
        vt.addline('set a [atomselect top noh]')
        vt.addline('set w [atomselect top water]')
        vt.addline('set h [atomselect top "resname HSD HSE HSP HIS"]')
        vt.addline('$w set resname HOH')
        vt.addline('$h set resname HIS')
        vt.addline(f'$a writepdb {self.basename}_pprep.pdb')
        vt.writescript()
        result=vt.runscript()
        return result
    
    def run_pdb2pqr(self,pH=7.0):
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
                        self.log_parser.metadata['histidines'][a.resname].append(a.resseqnum)

        logger.debug(f'Protonated/deprotonated histidines: {self.log_parser.metadata["histidines"]}')

        pka_table=self.log_parser.metadata['pka_table']
        allhismods=self.log_parser.metadata['histidines']
        # replace the resname field in the pka_table with the histidine resnames
        for hisname in ['HSD','HSE','HSP']:
            if hisname in allhismods:
                for resnum in allhismods[hisname]:
                    pka_table.loc[(pka_table['resnum']==resnum) & (pka_table['reschain']=='A'),'resname']=hisname

    def prep_mods(self):
        self.objmanager=ObjManager()
        # columns
        # resname  resnum reschain  respka  resmodelpka resatomtype  protonated
        for i,row in self.log_parser.metadata['pka_table'].iterrows():
            patchname=None
            if row['resname'] in ['HSE','HSP']: # HSD is the default histidine!
                # make the Patch shortcode
                if row['resname']=='HSE':
                    patchname='HS2'
                elif row['resname']=='HSP':
                    patchname='HSPP' # provided in pestifer.top
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

            if not patchname:
                logger.debug(f'No patch will be applied for {row["resname"]} at {row["reschain"]}:{row["resnum"]}')
                continue
            shortcode=f"{patchname}:{row['reschain']}:{row['resnum']}"
            pat=Patch(shortcode)
            self.objmanager.injest(pat)
        