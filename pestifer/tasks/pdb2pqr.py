# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import logging
import os
from ..command import Command
from ..basetask import BaseTask
from ..util.logparsers import PDB2PQRLog
from ..util.progress import PDB2PQRProgress
from ..atom import Atom, AtomList
from ..residue import Residue, ResidueList
from pidibble.pdbparse import PDBParser

logger=logging.getLogger(__name__)

class PDB2PQRTask(BaseTask):
    yaml_header='pdb2pqr'

    def do(self):
        self.log_message('initiated')
        self.inherit_state()
        self.next_basename()
        pH=self.specs.get('pH',7.0)
        xsc=self.statevars.get('xsc',None)

        prep_result=self.prep_input()
        self.run_pdb2pqr(pH=pH)

        if self.result!=0:
            return super().do()
        
        # self.save_state(exts=['psf','pdb','xsc'])
        self.log_message('complete')
        return super().do()

    def prep_input(self):
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

        