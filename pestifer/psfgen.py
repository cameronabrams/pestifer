"""

.. module:: psfgen
   :synopsis: Manages script generation and execution of psfgen under vmd
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
logger=logging.getLogger(__name__)
from .command import Command
import os
from .config import Config
from .resourcemanager import ResourcesGetPath as RGP
from .resourcemanager import ResourcesGet as RG
from .molecule import Molecule

class Psfgen:
    def __init__(self,**options):
        self.script_name=options.get('psfgen_script_name','my_mkpsf.tcl')
        self.script=['#### BEGIN HEADER ####']
        self.script.append(f'source {RGP("tcl")}/modules/src/loopmc.tcl')
        self.script.append(f'source {RGP("tcl")}/vmdrc.tcl')
        self.script.append('package require psfgen')
        self.script.append('psfcontxt mixedcase')
        for t in options['StdCharmmTopo']:
            ft=os.path.join(RG('system_charmm_toppardir'),t)
            self.script.append(f'topology {ft}')
        for lt in options['LocalCharmmTopo']:
            ft=os.path.join(RG('charmm'),lt)
            self.script.append(f'topology {ft}')
        for pdba in options['PDBAliases']:
            self.script.append(pdba)
        for k,v in _ResNameDict_PDB_to_CHARMM_.items():
            self.script.append(f'set RESDICT({k}) {v}')
        for k,v in _PDBAtomNameDict_.items():
            self.script.append(f'set ANAMEDICT({k}) {v}')
        self.script.append('#### END HEADER ####')

    def describe_molecule(self,mol:Molecule):
        molid_varname=f'm{mol.molid}'
        self.script.append(f'mol new {mol.pdb} waitfor all')
        self.script.append(f'set {molid_varname} [molinfo top get id]')
        if mol.cif:
            # need to renumber and rechain to user specs
            self.script.append(f'set ciftmp [atomselect ${molid_varname} all]')
            self.script.append('$ciftmp set chain [list {}]'.format(" ".join([_.chainID for _ in mol.Atoms])))
            self.script.append('$ciftmp set resid [list {}]'.format(" ".join([str(_.resseqnum) for _ in self.Atoms])))
        self.script.append(f'mol top ${molid_varname}')
        self.script.append('#### BEGIN SEGMENTS ####')
        Loops=[]
        """ For each BIOMT transformation in the requested biological assembly """
        for t in mol.activeBiologicalAssembly.biomt:
            # Identify and construct segments
            Segments=t.md.MakeSegments()
            for s in Segments:
                # Generate the psfgen script stanza for this segment
                stanza,loops=s.psfgen_stanza(includeTerminalLoops=self.md.userMods['includeTerminalLoops'],tmat=t)
                Loops.extend(loops) # psfgen postprocessing needs loop info
                ''' CONSTRUCTION IN PROGRESS '''
                self.script.extend(stanza)

    def run(self):
        if self.script[-1]!='exit':
            self.script.append('exit')
        with open(self.script_name,'w') as f:
            f.write('\n'.join(self.script))
        c=Command(f'{self.RM.vmd_location} -dispdev text -e {self.script_name}')
        c.run()