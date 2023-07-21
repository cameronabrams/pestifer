"""

.. module:: psfgen
   :synopsis: Manages script generation and execution of psfgen under vmd
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
logger=logging.getLogger(__name__)
from .command import Command
import os
from .resources import ResourceManager
from .atom import _PDBAtomNameDict_
from .residue import _ResNameDict_PDB_to_CHARMM_
from .molecule import Molecule

_molidcounter_=0
class Psfgen:
    def __init__(self,resman:ResourceManager,**options):
        if not resman:
            logger.fatal(f'Cannot use psfgen without a resource manager.')
        self.resman=resman
        self.script_name=options.get('psfgen_script_name','my_mkpsf.tcl')
        self.script=['#### BEGIN HEADER ####']
        self.script.append(f'source {resman.resource_paths["tcl"]}/modules/src/loopmc.tcl')
        self.script.append(f'source {resman.resource_paths["tcl"]}/vmdrc.tcl')
        self.script.append('package require psfgen')
        self.script.append('psfcontxt mixedcase')
        for t in options['StdCharmmTopo']:
            ft=os.path.join(resman.system_charmm_toppardir,t)
            self.script.append(f'topology {ft}')
        for lt in options['LocalCharmmTopo']:
            ft=os.path.join(resman.resource_paths['charmm'],lt)
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
        for t in mol.activeBiologicalAssembly.biomt:
            Segments=t.md.MakeSegments()
            for s in Segments:
                stanza,loops=s.psfgen_stanza(includeTerminalLoops=self.md.userMods['includeTerminalLoops'],tmat=t)
                Loops.extend(loops) # psfgen postprocessing needs loop info
                ''' CONSTRUCTION IN PROGRESS '''
                self.script.extend(stanza)
        pass

    # def addline(self,line):
    #     self.script.append(line)

    def run(self):
        if self.script[-1]!='exit':
            self.script.append('exit')
        with open(self.script_name,'w') as f:
            f.write('\n'.join(self.script))
        c=Command(f'{self.resman.vmd_location} -dispdev text -e {self.script_name}')
        c.run()