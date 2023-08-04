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
from .resourcemanager import ResourcesGetPath as RGP
from .resourcemanager import ResourcesGet as RG
from .molecule import Molecule

class Psfgen:
    def __init__(self):
        self.script_name=ConfigGetParam('psfgen_script_name')
        self.script='#### BEGIN PESTIFER PSFGEN ####\n'
        self.script+='#### BEGIN HEADER ####\n'
        self.script+=f'source {RGP("tcl")}/modules/src/loopmc.tcl\n'
        self.script+=f'source {RGP("tcl")}/vmdrc.tcl\n'
        self.script+='package require psfgen\n'
        self.script+='psfcontxt mixedcase\n'
        for t in ConfigGetParam('StdCharmmTopo'):
            ft=os.path.join(RG('system_charmm_toppardir'),t)
            self.script+=f'topology {ft}\n'
        for lt in ConfigGetParam('LocalCharmmTopo'):
            ft=os.path.join(RGP('charmm'),lt)
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

    def describe_molecule(self,mol:Molecule):
        molid_varname=f'm{mol.molid}'
        self.script+=f'mol new {mol.pdb_code} waitfor all\n'
        self.script+=f'set {molid_varname} [molinfo top get id]\n'
        # if mol.cif:
        #     # need to renumber and rechain to user specs
        #     self.script+=f'set ciftmp [atomselect ${molid_varname} all]\n'
        #     self.script+='$ciftmp set chain [list {}]\n'.format(" ".join([_.chainID for _ in mol.Atoms]))
        #     self.script+'$ciftmp set resid [list {}]\n'.format(" ".join([str(_.resseqnum) for _ in self.Atoms]))
        self.script+=f'mol top ${molid_varname}\n'
        self.script+=f'#### BEGIN SEGMENTS ####\n'
        # Loops=[]
        # """ For each BIOMT transformation in the requested biological assembly """
        # for t in mol.activeBiologicalAssembly.biomt:
        #     # Identify and construct segments
        #     Segments=t.md.MakeSegments()
        #     for s in Segments:
        #         # Generate the psfgen script stanza for this segment
        #         stanza,loops=s.psfgen_stanza(includeTerminalLoops=self.md.userMods['includeTerminalLoops'],tmat=t)
        #         Loops.extend(loops) # psfgen postprocessing needs loop info
        #         ''' CONSTRUCTION IN PROGRESS '''
        #         self.script.extend(stanza)
        self.script+=f'#### END SEGMENTS ####\n'

    def run(self):
        c=Command(f'{RG("VMD")} -dispdev text -e {self.script_name}')
        o,e=c.run()
        return o,e