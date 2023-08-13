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
from .stringthings import ByteCollector

class Psfgen:
    def __init__(self,resources,post_build_directives={}):
        self.resources=resources
        self.tcl_path=resources.ResourcePaths['tcl']
        self.package_charmmdir=resources.ResourcePaths['charmm']
        self.system_charmm_toppardir=resources.system_charmm_toppardir
        self.post_build_directives=post_build_directives
        self.script_name=ConfigGetParam('psfgen_script_name')
        self.B=ByteCollector()

    def beginscript(self):
        self.B.banner('BEGIN PESTIFER PSFGEN')
        self.B.banner('BEGIN HEADER')
        self.B.addline(f'source {self.tcl_path}/modules/src/loopmc.tcl')
        self.B.addline(f'source {self.tcl_path}/vmdrc.tcl')
        self.B.addline('package require psfgen')
        self.B.addline('psfcontext mixedcase')
        for t in ConfigGetParam('StdCharmmTopo'):
            ft=os.path.join(self.system_charmm_toppardir,t)
            self.B.addline(f'topology {ft}')
        for lt in ConfigGetParam('LocalCharmmTopo'):
            ft=os.path.join(self.package_charmmdir,lt)
            self.B.addline(f'topology {ft}')
        for pdba in ConfigGetParam('PDBAliases'):
            self.B.addline(f'pdbalias {pdba}')
        for k,v in ConfigGetParam('PDB_to_CHARMM_Resnames').items():
            self.B.addline(f'set RESDICT({k}) {v}')
        for k,v in ConfigGetParam('PDB_to_CHARMM_Atomnames').items():
            self.B.addline(f'set ANAMEDICT({k}) {v}')
        self.B.banner('END HEADER')

    def endscript(self):
        with open(self.script_name,'w') as f:
            self.B.addline('exit')
            self.B.banner('END PESTIFER PSFGEN')
            self.B.banner('Thank you for using pestifer!')
            f.write(self.B.byte_collector)

    def cleanfiles(self):
        if os.path.exists(self.script_name):
            os.remove(self.script_name)

    def set_molecule(self,mol):
        molid_varname=f'm{mol.molid}'
        self.B.addline(f'mol new {mol.pdb_code} waitfor all')
        self.B.addline(f'set {molid_varname} [molinfo top get id]')
        
    def describe_molecule(self,mol:Molecule,mods):
        # if mol.cif:
        #     # need to renumber and rechain to user specs
        #     self.B.addline(f'set ciftmp [atomselect ${molid_varname} all]')
        #     self.B.addline('$ciftmp set chain [list {}]').format(" ".join([_.chainID for _ in mol.Atoms]))
        #     self.script+'$ciftmp set resid [list {}]').format(" ".join([str(_.resseqnum) for _ in self.Atoms]))
        molid_varname=f'm{mol.molid}'
        self.B.addline(f'mol top ${molid_varname}')
        self.B.addline(mol.write_TcL(mods))

    def transform_postmods(self,outputs):
        self.B.addline('guesscoord')
        self.B.addline('regenerate angles dihedrals')
        psf=outputs['psf']
        pdb=outputs['pdb']
        self.B.addline(f'writepsf cmap {psf}')
        self.B.addline(f'writepdb {pdb}')

    def global_postmods(self):
        pass

    def run(self):
        c=Command(f'{self.resources.vmd} -dispdev text -e {self.script_name}')
        o,e=c.run()
        return o,e