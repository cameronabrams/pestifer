"""

.. module:: psfgen
   :synopsis: Manages script generation and execution of psfgen under vmd
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
logger=logging.getLogger(__name__)
from .command import Command
import os
from .config import ConfigGetParam
from .molecule import Molecule

class Psfgen:
    def __init__(self,resources):
        self.resources=resources
        tcl_path=resources.ResourcePaths['tcl']
        package_charmmdir=resources.ResourcePaths['charmm']
        system_charmm_toppardir=resources.system_charmm_toppardir
        self.script_name=ConfigGetParam('psfgen_script_name')
        self.script='#### BEGIN PESTIFER PSFGEN ####\n'
        self.script+='#### BEGIN HEADER ####\n'
        self.script+=f'source {tcl_path}/modules/src/loopmc.tcl\n'
        self.script+=f'source {tcl_path}/vmdrc.tcl\n'
        self.script+='package require psfgen\n'
        self.script+='psfcontxt mixedcase\n'
        for t in ConfigGetParam('StdCharmmTopo'):
            ft=os.path.join(system_charmm_toppardir,t)
            self.script+=f'topology {ft}\n'
        for lt in ConfigGetParam('LocalCharmmTopo'):
            ft=os.path.join(package_charmmdir,lt)
            self.script+=f'topology {ft}\n'
        for pdba in ConfigGetParam('PDBAliases'):
            self.script+=f'{pdba}\n'
        for k,v in ConfigGetParam('PDB_to_CHARMM_Resnames').items():
            self.script+=f'set RESDICT({k}) {v}\n'
        for k,v in ConfigGetParam('PDB_to_CHARMM_Atomnames').items():
            self.script+=f'set ANAMEDICT({k}) {v}\n'
        self.script+=f'#### END HEADER ####\n'

    def write(self):
        with open(self.script_name,'w') as f:
            f.write(self.script)
            f.write('exit\n')
            f.write('#### END PESTIFER PSFGEN ####\n')

    def cleanfiles(self):
        if os.path.exists(self.script_name):
            os.remove(self.script_name)

    def set_molecule(self,mol):
        molid_varname=f'm{mol.molid}'
        self.script+=f'mol new {mol.pdb_code} waitfor all\n'
        self.script+=f'set {molid_varname} [molinfo top get id]\n'
        
    def describe_molecule(self,mol:Molecule,mods):
        # if mol.cif:
        #     # need to renumber and rechain to user specs
        #     self.script+=f'set ciftmp [atomselect ${molid_varname} all]\n'
        #     self.script+='$ciftmp set chain [list {}]\n'.format(" ".join([_.chainID for _ in mol.Atoms]))
        #     self.script+'$ciftmp set resid [list {}]\n'.format(" ".join([str(_.resseqnum) for _ in self.Atoms]))
        molid_varname=f'm{mol.molid}'
        self.script+=f'mol top ${molid_varname}\n'
        self.script+=mol.write_TcL()

    def transform_postmods(self,outputs):
        self.script+='guesscoord\n'
        self.script+='regenerate angles dihedrals\n'
        psf=outputs['psf']
        pdb=outputs['pdb']
        self.script+=f'writepsf cmap {psf}\n'
        self.script+=f'writepdb {pdb}\n'

    def global_postmods(self):
        pass

    def run(self):
        c=Command(f'{self.resources.vmd} -dispdev text -e {self.script_name}')
        o,e=c.run()
        return o,e