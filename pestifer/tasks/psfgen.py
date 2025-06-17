# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import logging
import networkx as nx
import shutil

from copy import deepcopy

from ..basetask import BaseTask
from ..chainidmanager import ChainIDManager
from ..command import Command
from ..molecule import Molecule
from ..objmanager import ObjManager
from ..psfutil.psfcontents import PSFContents

logger=logging.getLogger(__name__)

class PsfgenTask(BaseTask):
    """ A class for handling invocations of psfgen which create a molecule from a base PDB/mmCIF file
    or from a PSF file generated previously by psfgen
    
    Attributes
    ----------
    yaml_header(str) 

    Methods
    -------
    do(): 
        Based on specs, reads in input PDB/mmCIF file, generates parsed Molecule instances, generates
        the psfgen script, and executes VMD to perform the psfgen run.

    """
    yaml_header='psfgen'
    def __init__(self,config_specs={},controller_specs={}):
        super().__init__(config_specs,controller_specs)
        self.molecules={}
        self.keepfiles=[]
        if self.specs.get('source',{}).get('prebuilt',{}):
            self.keepfiles=[self.specs['source']["prebuilt"]["psf"],self.specs['source']["prebuilt"]["pdb"]]
            xsc=self.specs['source'].get('xsc','')
            if xsc:
                self.keepfiles.append(xsc)

    def do(self):
        self.log_message('initiated')
        self.inherit_state()
        logger.debug('Injesting molecule(s)')
        self.injest_molecules()
        self.statevars['base_molecule']=self.base_molecule
        logger.debug(f'base mol num images {self.base_molecule.num_images()}')
        logger.debug('Running first psfgen')
        self.result=self.psfgen()
        if self.result!=0:
            return super().do()
        # we now have a full coordinate set, so we can do coormods
        self.coormods()
        # min_loop_length=0
        min_loop_length=self.specs['source'].get('sequence',{}).get('loops',{}).get('min_loop_length',0)
        self.update_statevars('min_loop_length',min_loop_length)
        nloops=self.base_molecule.has_loops(min_loop_length=min_loop_length)*self.base_molecule.num_images()
        if nloops>0 and self.specs['source']['sequence']['loops']['declash']['maxcycles']>0:
            logger.debug(f'Declashing {nloops} loops')
            self.declash_loops(self.specs['source']['sequence']['loops'])
        nglycans=self.base_molecule.nglycans()*self.base_molecule.num_images()
        if nglycans>0 and self.specs['source']['sequence']['glycans']['declash']['maxcycles']>0:
            logger.debug(f'Declashing {nglycans} glycan segments')
            self.declash_glycans(self.specs['source']['sequence']['glycans'])
        self.log_message('complete')
        return super().do()

    def coormods(self):
        coormods=self.objmanager.get('coord',{})
        logger.debug(f'psfgen task has {len(coormods)} coormods:')
        logger.debug(';'.join([str(_) for _ in coormods]))
        ba=self.base_molecule.active_biological_assembly
        if coormods:
            logger.debug(f'performing coormods')
            for objtype,objlist in coormods.items():
                if len(objlist)>0:
                    self.next_basename(objtype)
                    vm=self.writers['vmd']
                    packages=[]
                    if objtype=='crotations':
                        packages.append('PestiferCRot')
                    vm.newscript(self.basename,packages=packages)
                    psf=self.statevars['psf']
                    pdb=self.statevars['pdb']
                    vm.load_psf_pdb(psf,pdb,new_molid_varname='mCM')
                    for transform in ba.transforms:
                        objlist.write_TcL(vm,chainIDmap=transform.chainIDmap)
                    vm.write_pdb(self.basename,'mCM')
                    vm.writescript()
                    vm.runscript()
                    self.save_state(exts=['pdb'])

    def patch_topologies(self):
        objmanager=self.base_molecule.objmanager
        seqmods=objmanager.get('seq',{})
        patches=seqmods.get('patches',[])
        CC=self.config.RM.charmmff_content
        new_topfiles=set()
        for patch in patches:
            new_topfiles.add(CC.get_topfile_of_patchname(patch.patchname))
        return list(new_topfiles)

    def psfgen(self):
        self.next_basename('build')
        pg=self.writers['psfgen']
        pg.newscript(self.basename,packages=['PestiferCRot'],additional_topologies=self.patch_topologies())
        pg.set_molecule(self.base_molecule,altcoords=self.specs.get('source',{}).get('altcoords',None))
        pg.describe_molecule(self.base_molecule)
        pg.writescript(self.basename)
        result=pg.runscript()
        if result!=0:
            return result
        if self.specs.get('remove_charmmtop',False):
            pg.remove_charmm_top(pg.topologies)
        self.save_state(exts=['psf','pdb'])
        self.strip_remarks()
        return 0
        
    def strip_remarks(self):
        pdb=self.statevars['pdb']
        c=Command(f'grep -v ^REMARK {pdb} > tmp').run()
        shutil.move('tmp',pdb)

    def declash_loops(self,specs):
        mol=self.base_molecule
        cycles=specs['declash']['maxcycles']
        if not mol.has_loops() or not cycles:
            logger.debug(f'Loop declashing is intentionally not done.')
            return
        self.next_basename('declash-loops')
        vt=self.writers['vmd']
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        vt.newscript(self.basename,packages=['PestiferDeclash'])
        vt.load_psf_pdb(psf,pdb,new_molid_varname='mLL')
        mol.write_loop_lines(vt,cycles=cycles,min_length=specs['min_loop_length'],include_c_termini=specs['declash']['include_C_termini'])
        vt.write_pdb(self.basename,'mLL')
        vt.writescript()
        vt.runscript()
        self.save_state(exts=['pdb'])

    def declash_glycans(self,specs):
        mol=self.base_molecule
        cycles=specs['declash']['maxcycles']
        clashdist=specs['declash']['clashdist']
        if not mol.nglycans() or not cycles:
            logger.debug(f'Glycan declashing is intentionally not done.')
            return
        self.next_basename('declash-glycans')
        outpdb=f'{self.basename}.pdb'
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        vt=self.writers['vmd']
        vt.newscript(self.basename,packages=['PestiferDeclash'])
        vt.addline(f'mol new {psf}')
        vt.addline(f'mol addfile {pdb} waitfor all')
        vt.addline(f'set a [atomselect top all]')
        vt.addline(f'set molid [molinfo top get id]')
        nglycan=self._write_glycans(vt)
        vt.addline(f'vmdcon -info "Declashing $nglycans glycans; clashdist {clashdist}; maxcycles {cycles}"')
        vt.addline(r'for {set i 0} {$i<$nglycans} {incr i} {')
        vt.addline(f'   declash_pendant $molid $glycan_idx($i) $rbonds($i) $movers($i) {cycles} {clashdist}')
        vt.addline(r'}')
        vt.addline(f'$a writepdb {outpdb}')
        vt.writescript()
        logger.debug(f'Declashing {nglycan} glycans')
        vt.runscript(progress_title='declash-glycans')
        self.save_state(exts=['pdb'])

    def _write_glycans(self,fw):
        psf=self.statevars['psf']
        logger.debug(f'Injesting {psf}')
        struct=PSFContents(psf,parse_topology=['bonds'])
        logger.debug(f'Making graph structure of glycan atoms...')
        glycanatoms=struct.atoms.get(segtype='glycan')
        logger.debug(f'{len(glycanatoms)} total glycan atoms')
        glycangraph=glycanatoms.graph()
        G=[glycangraph.subgraph(c).copy() for c in nx.connected_components(glycangraph)]
        logger.debug(f'Preparing declash input for {len(G)} glycans')
        fw.addline(f'set nglycans {len(G)}')
        for i,g in enumerate(G):
            logger.debug(f'Glycan {i} has {len(g)} atoms')
            serials=[x.serial for x in g]
            for at in g:
                lig_ser=[x.serial for x in at.ligands]
                for k,ls in enumerate(lig_ser):
                    if not ls in serials:
                        at.is_root=True
                        rp=at.ligands[k]
                        logger.debug(f'-> Atom {str(at)} is the root, bound to atom {str(rp)}')
            indices=' '.join([str(x.serial-1) for x in g])
            fw.comment(f'Glycan {i}:')
            fw.addline(f'set glycan_idx({i}) [list {indices}]')
            fw.addline(f'set rbonds({i}) [list]')
            fw.addline(f'set movers({i}) [list]')
            for bond in nx.bridges(g):
                ai,aj=bond
                if not (ai.isH() or aj.isH()) and not ai.is_pep(aj):
                    g.remove_edge(ai,aj)
                    S=[g.subgraph(c).copy() for c in nx.connected_components(g)]
                    assert len(S)==2,f'Bond {ai.serial-1}-{aj.serial-1} when cut makes more than 2 components'
                    for sg in S:
                        is_root=any([hasattr(x,'is_root') for x in sg])
                        if not is_root:
                            if ai in sg:
                                sg.remove_node(ai)
                            if aj in sg:
                                sg.remove_node(aj)
                            if len(sg)>1 or (len(sg)==1 and not [x for x in sg.nodes][0].isH()):
                                mover_serials=[x.serial for x in sg]
                                mover_indices=" ".join([str(x-1) for x in mover_serials])
                                logger.debug(f'{str(ai)}--{str(aj)} is a rotatable bridging bond')
                                fw.addline(f'lappend rbonds({i}) [list {ai.serial-1} {aj.serial-1}]')
                                logger.debug(f'  -> movers: {" ".join([str(x) for x in sg])}')
                                fw.addline(f'lappend movers({i}) [list {mover_indices}]')
                    g.add_edge(ai,aj)
        return len(G)

    def injest_molecules(self):
        specs=self.specs
        self.source_specs=specs['source']
        logger.debug(f'User-input modspecs {self.specs["mods"]}')
        self.objmanager=ObjManager(self.specs['mods'])
        seqmods=self.objmanager.get('seq',{})
        logger.debug(f'Injesting seqmods {seqmods}')
        if 'grafts' in seqmods:
            logger.debug(f'looking for graft sources to injest')
            Grafts=seqmods['grafts']
            for g in Grafts:
                if not g.source_pdbid in self.molecules:
                    logger.debug(f'Injesting graft source {g.source_pdbid}')
                    this_source={
                        'id':g.source_pdbid,
                        'file_format':'PDB'
                    }
                    self.molecules[g.source_pdbid]=Molecule(source=this_source)
                g.activate(deepcopy(self.molecules[g.source_pdbid]))
        self.chainIDmanager=ChainIDManager(
            format=self.source_specs['file_format'],
            transform_reserves=self.source_specs.get('transform_reserves',{}),
            remap=self.source_specs.get('remap_chainIDs',{}))
        self.base_molecule=Molecule(source=self.source_specs,
                                    objmanager=self.objmanager,
                                    chainIDmanager=self.chainIDmanager).activate_biological_assembly(self.source_specs['biological_assembly'])
        if self.source_specs.get('id',{}):
            key=self.source_specs['id']
        elif self.source_specs.get('prebuilt',{}):
            key=f'{self.source_specs["prebuilt"]["psf"]}-{self.source_specs["prebuilt"]["pdb"]}'
            xsc=self.source_specs.get('xsc','')
            if xsc:
                self.update_statevars('xsc',xsc)
        elif self.source_specs.get('alphafold',{}):
            key=f'{self.source_specs["alphafold"]}'
        else:
            raise Exception(f'The "source" directive of "psfgen" must have "id" , "prebuilt", or "alphafold"')
        self.molecules[key]=self.base_molecule
        for molid,molecule in self.molecules.items():
            logger.debug(f'Molecule "{molid}": {molecule.num_atoms()} atoms in {molecule.num_residues()} residues; {molecule.num_segments()} segments.')

    def update_molecule(self):
        """Updates all segments of the base molecule based on the 
           current coordinate file.  All ssbonds and links are 
           carried forward. No biological assembly beyond the apparent
           asymmetric unit is assumed. This should permit generation
           of a new psfgen script based on this new molecule to 
           recreate it or modify it.
        """
        # get the key of the base_molecule
        logger.debug(f'{self.taskname} has {len(self.molecules)} entries in its molecules dict')
        base_key='base'
        for k,v in self.molecules:
            if v==self.base_molecule:
                base_key=k
        # assert base_key!='UNSET',f'Cannot update a non-existent base molecule'
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        xsc=self.statevars.get('xsc','')
        source={
            'prebuilt': {
                'psf':psf,
                'pdb':pdb,
                'xsc':xsc
            }
        }
        if hasattr(self,'chainIDManager') and hasattr(self,'objmanager'):
            updated_molecule=Molecule(source=source,chainIDmanager=self.chainIDmanager,objmanager=self.objmanager).activate_biological_assembly(0)
        else:
            updated_molecule=Molecule(source=source).activate_biological_assembly(0)

        self.molecules[base_key]=updated_molecule
        self.base_molecule=updated_molecule
