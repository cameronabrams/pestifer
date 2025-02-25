# Author: Cameron F. Abrams, <cfa22@drexel.edu>

from .psfgen import PsfgenTask
from ..objs.cleavagesite import CleavageSite, CleavageSiteList

class CleaveTask(PsfgenTask):
    yaml_header='cleave'
    def do(self):
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