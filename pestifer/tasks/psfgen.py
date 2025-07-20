# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`PsfgenTask` class for handling invocations of psfgen which create a molecule from a base PDB/mmCIF file.

Usage is described in the :ref:`subs_runtasks_psfgen` documentation.

"""
import logging
import networkx as nx
import shutil
import os

from copy import deepcopy

from ..core.basetask import BaseTask
from ..molecule.chainidmanager import ChainIDManager
from ..core.command import Command
from ..molecule.molecule import Molecule
from ..core.objmanager import ObjManager
from ..core.pipeline import PipelineContext, PSFFile, PDBFile, XSCFile
from ..psfutil.psfatom import PSFAtomList
from ..psfutil.psfcontents import PSFContents

logger=logging.getLogger(__name__)

class PsfgenTask(BaseTask):
    """ 
    A class for handling invocations of psfgen which create a molecule from a base PDB/mmCIF file
    or from a PSF file generated previously by psfgen
    This class is a descendant of the :class:`BaseTask <pestifer.core.basetask.BaseTask>` class, and as such it has two parameters.

    Parameters
    ----------
    config_specs : dict
        Configuration specifications for the task.
    controller_specs : dict
        Controller specifications for the task.
    """

    yaml_header='psfgen'
    """
    YAML header for the PsfgenTask, used to identify the task in configuration files as part of a ``tasks`` list.
    """
    def __init__(self,ctx:PipelineContext,config_specs={},controller_specs={}):
        super().__init__(ctx,config_specs,controller_specs)
        self.molecules={}
        # self.keepfiles=[]
        # if self.specs.get('source',{}).get('prebuilt',{}):
        #     self.keepfiles=[self.specs['source']["prebuilt"]["psf"],self.specs['source']["prebuilt"]["pdb"]]
        #     xsc=self.specs['source'].get('xsc','')
        #     if xsc:
        #         self.keepfiles.append(xsc)

    def do(self):
        """
        Execute the psfgen task.
        This method initializes the task, ingests the base molecule, and runs the psfgen process.
        It also handles any necessary coormods and declashing of loops and glycans based on the task specifications.
        The results of the psfgen process are saved as a PSF/PDB fileset, and the state is updated accordingly.
        """
        self.log_message('initiated')
        self.inherit_artifacts()
        logger.debug('ingesting molecule(s)')
        self.ingest_molecules()
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
        self.nloops=self.base_molecule.has_loops(min_loop_length=min_loop_length)
        for segtype in ['protein','nucleicacid']:
            self.nloops[segtype]*=self.base_molecule.num_images()
        if self.nloops['protein']>0 and self.specs['source']['sequence']['loops']['declash']['maxcycles']>0:
            logger.debug(f'Declashing {self.nloops["protein"]} protein loops')
            self.declash_protein_loops(self.specs['source']['sequence']['loops'])
        if self.nloops['nucleicacid']>0 and self.specs['source']['sequence']['loops']['declash']['maxcycles']>0:
            logger.debug(f'Declashing {self.nloops["nucleicacid"]} nucleic acid loops')
            self.declash_na_loops(self.specs['source']['sequence']['loops'])
        nglycans=self.base_molecule.nglycans()*self.base_molecule.num_images()
        if nglycans>0 and self.specs['source']['sequence']['glycans']['declash']['maxcycles']>0:
            logger.debug(f'Declashing {nglycans} glycan segments')
            self.declash_glycans(self.specs['source']['sequence']['glycans'])
        self.log_message('complete')
        return super().do()

    def coormods(self):
        """
        Perform coordinate modifications based on the specifications provided in the task.
        """
        coormods=self.objmanager.get('coord',{})
        logger.debug(f'psfgen task has {len(coormods)} coormods:')
        logger.debug(';'.join([str(_) for _ in coormods]))
        ba=self.base_molecule.active_biological_assembly
        if coormods:
            logger.debug(f'performing coormods')
            for objtype,objlist in coormods.items():
                if len(objlist)>0:
                    self.next_basename(objtype)
                    vm=self.scripters['vmd']
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

    def resi_topologies(self):
        """
        Collect the topology files that are needed for the residues in the base molecule.
        """
        resis=set([x.resname for x in self.base_molecule.asymmetric_unit.residues])
        CC=self.config.RM.charmmff_content
        new_topfiles=set()
        for resname in resis:
            topfile=CC.get_topfile_of_resname(resname)
            if topfile:
                new_topfiles.add(topfile)
        return list(new_topfiles)

    def patch_topologies(self):
        """
        Collect the topology files that are needed for the patches and links in the base molecule.
        This method retrieves the CHARMMFF topology files associated with the patches and links defined 
        in the base molecule's object manager.
        """
        objmanager=self.base_molecule.objmanager
        seqmods=objmanager.get('seq',{})
        patches=seqmods.get('patches',[])
        CC=self.config.RM.charmmff_content
        new_topfiles=set()
        # logger.debug(f'New topologies: {new_topfiles}')
        for patch in patches:
            topfile=CC.get_topfile_of_patchname(patch.patchname)
            if topfile:
                new_topfiles.add(topfile)
        # logger.debug(f'New topologies: {new_topfiles}')
        topomods=objmanager.get('topol',{})
        links=topomods.get('links',[])
        for link in links:
            topfile=CC.get_topfile_of_patchname(link.patchname)
            if topfile:
                new_topfiles.add(topfile)
        # logger.debug(f'New topologies: {new_topfiles}')
        return list(new_topfiles)

    def psfgen(self):
        """
        Run the psfgen process to generate a PSF file from the base molecule.
        """
        self.next_basename('build')
        pg=self.scripters['psfgen']
        patch_topologies=self.patch_topologies()
        resi_topologies=self.resi_topologies()
        addl_topologies=list(set(patch_topologies+resi_topologies))
        pg.newscript(self.basename,packages=['PestiferCRot'],additional_topologies=addl_topologies)
        pg.set_molecule(self.base_molecule,altcoords=self.specs.get('source',{}).get('altcoords',None))
        pg.describe_molecule(self.base_molecule)
        pg.writescript(self.basename)
        result=pg.runscript()
        if result!=0:
            return result
        for ptop in addl_topologies:
            if ptop.endswith('.str'):
                if not self.statevars.get('charmmff_paramfiles',[]):
                    self.statevars['charmmff_paramfiles']=[]
                if ptop not in self.statevars['charmmff_paramfiles']:
                    logger.debug(f'Adding {ptop} to charmmff_paramfiles')
                    self.statevars['charmmff_paramfiles'].append(ptop)
        # register PSF, PDB, tcl, log, and all charmmff files in the pipeline context
        exts=['psf','pdb','tcl','log']
        for ext in exts:
            fname=f'{self.basename}.{ext}'
            if os.path.exists(fname):
                objtype=PSFFile if ext=='psf' else PDBFile if ext=='pdb' else XSCFile if ext=='xsc' else None
                if objtype:
                    self.ctx.register(
                        key=ext,
                        value=objtype(path=fname),
                        value_type=objtype,
                        produced_by=self,
                        type='intermediate',
                        propagate=ext in ['psf','pdb']
                    )
        for topology in pg.topologies:
            
            if topology.endswith('.str'):
                self.ctx.register(
                    key=topology,
                    value=topology,
                    value_type=str,
                    produced_by=self,
                    type='intermediate',
                    propagate=True
                )
        self.strip_remarks()
        return 0
        
    def strip_remarks(self):
        """
        Strip REMARK lines from the PDB file generated by psfgen.
        This method removes any REMARK lines from the PDB file to ensure that it contains only the relevant atomic coordinates and structure information.
        """
        pdb=self.get_current_artifact('pdb')
        if not pdb:
            logger.warning('No PDB file found to strip remarks from.')
            return
        if not os.path.exists(pdb):
            logger.warning(f'PDB file {pdb} does not exist, cannot strip remarks.')
            return
        logger.debug(f'Stripping REMARK lines from {pdb}')
        with open(pdb, 'r') as infile:
            lines = infile.readlines()
        with open(pdb, 'w') as outfile:
            for line in lines:
                if not line.startswith('REMARK'):
                    outfile.write(line)

    def declash_protein_loops(self,specs):
        """
        Declash loops in the base molecule using the custom ``PestiferDeclash`` TcL package.
        This method generates a VMD script to identify and declash loops in the molecular structure.
        It uses the ``PestiferDeclash`` package to perform the declashing operation, which involves identifying
        and modifying the coordinates of atoms in loop regions.
        """
        mol=self.base_molecule
        cycles=specs['declash']['maxcycles']
        if self.nloops['protein']==0 or not cycles:
            logger.debug(f'Protein loop declashing is intentionally not done.')
            return
        self.next_basename('declash-loops')
        vt=self.scripters['vmd']
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        vt.newscript(self.basename,packages=['PestiferDeclash'])
        vt.load_psf_pdb(psf,pdb,new_molid_varname='mLL')
        mol.write_protein_loop_lines(vt,cycles=cycles,min_length=specs['min_loop_length'],include_c_termini=specs['declash']['include_C_termini'])
        vt.write_pdb(self.basename,'mLL')
        vt.writescript()
        vt.runscript()
        self.ctx.register(
            key='pdb',
            value=PDBFile(path=f'{self.basename}.pdb'),
            value_type=PDBFile,
            produced_by=self,
            type='intermediate',
            propagate=True
        )

    def declash_na_loops(self,specs):
        """
        Declash nucleic acid loops in the base molecule using the custom ``PestiferDeclash`` TcL package.
        This method generates a VMD script to identify and declash nucleic acid loops in the molecular structure.
        It uses the ``PestiferDeclash`` package to perform the declashing operation, which involves identifying
        and modifying the coordinates of atoms in nucleic acid loop regions.
        """
        mol=self.base_molecule
        cycles=specs['declash']['maxcycles']
        clashdist=specs['declash']['clashdist']
        minlooplength=specs['min_loop_length']

        if self.nloops['nucleicacid'] == 0 or not cycles:
            logger.debug(f'Nucleic acid loop declashing is intentionally not done.')
            return
        self.next_basename('declash-na-loops')
        vt=self.scripters['vmd']
        psf=self.statevars['psf']
        pdb=self.statevars['pdb']
        outpdb=f'{self.basename}.pdb'
        vt.newscript(self.basename,packages=['PestiferDeclash'])
        vt.addline(f'mol new {psf}')
        vt.addline(f'mol addfile {pdb} waitfor all')
        vt.addline(f'set a [atomselect top all]')
        vt.addline(f'set molid [molinfo top get id]')
        nna=self._write_na_loops(vt,minlooplength=minlooplength)
        vt.addline(f'set nna {nna}')
        vt.addline(f'vmdcon -info "Declashing $nna nucleic acid loops; clashdist {clashdist}; maxcycles {cycles}"')
        vt.addline(r'for {set i 0} {$i<$nna} {incr i} {')
        vt.addline(f'   declash_pendant $molid $na_idx($i) $rbonds($i) $movers($i) {cycles} {clashdist}')
        vt.addline(r'}')
        vt.addline(f'$a writepdb {outpdb}')
        vt.writescript()
        logger.debug(f'Declashing {nna} nucleic acid loops')
        vt.runscript(progress_title='declash-nucleic-acid-loops')
        self.save_state(exts=['pdb'])

    def _write_na_loops(self,vt,**options):
        mol=self.base_molecule
        au=mol.asymmetric_unit
        psf=self.statevars['psf']
        logger.debug(f'ingesting {psf}')
        struct=PSFContents(psf,parse_topology=['bonds'])
        na_atoms=struct.atoms.get(segtype='nucleicacid')
        my_rep=list(set([(x.chainID,x.resseqnum) for x in na_atoms]))
        my_rep.sort(key=lambda x: (x[0],x[1]))
        logger.debug(f'Getting loops from {len(na_atoms)} nucleic acid atoms in PSF file {psf}')
        logger.debug(f'{my_rep}')
        min_length=options.get('minlooplength',4)
        include_c_termini=options.get('include_c_termini',False)
        i=0
        SL=[S for S in au.segments if S.segtype=='nucleicacid']
        for S in SL:
            asymm_segname=S.segname
            n_subsegs=len(S.subsegments)
            for b in S.subsegments:
                lr_resseqnum=S.residues[b.bounds[0]].resseqnum
                rr_resseqnum=S.residues[b.bounds[1]].resseqnum
                logger.debug(f'Processing subsegment {b.pstr()} for segname {asymm_segname} with bounds {lr_resseqnum}-{rr_resseqnum}')
                is_c_terminus=(S.subsegments.index(b)==(n_subsegs-1))
                is_processible=b.state=='MISSING' and b.num_items()>=min_length
                if is_processible and (not include_c_termini) and is_c_terminus:
                    logger.debug(f'A.U. C-terminal loop {b.pstr()} declashing is skipped')
                    is_processible=False
                if is_processible:
                    logger.debug(f'Processing loop {b.pstr()} {b.bounds} for segname {asymm_segname}')
                    loop_atoms=PSFAtomList([x for x in na_atoms if x.chainID==asymm_segname and x.resseqnum>=lr_resseqnum and x.resseqnum<=rr_resseqnum])
                    logger.debug(f'Loop {b.pstr()} has {len(loop_atoms)} atoms from PSFAtomList')
                    na_graph=loop_atoms.graph()
                    logger.debug(f'{na_graph}')
                    G=[na_graph.subgraph(c).copy() for c in nx.connected_components(na_graph)]
                    assert len(G)==1,f'NA loop {b.pstr()} has more than one connected component'
                    logger.debug(f'Loop {b.pstr()} has {len(loop_atoms)} atoms')
                    g=G[0]
                    serials=[x.serial for x in g]
                    for at in g:
                        lig_ser=[x.serial for x in at.ligands]
                        for k,ls in enumerate(lig_ser):
                            if not ls in serials:
                                at.is_root=True
                                rp=at.ligands[k]
                                logger.debug(f'-> Atom {str(at)} is the root, bound to atom {str(rp)}')
                    indices=' '.join([str(x.serial-1) for x in g])
                    vt.addline(f'set na_idx({i}) [list {indices}]')
                    vt.addline(f'set rbonds({i}) [list]')
                    vt.addline(f'set movers({i}) [list]')
                    for bond in nx.bridges(g):
                        ai,aj=bond
                        if not (ai.isH() or aj.isH()) and not ai.is_pep(aj):
                            g.remove_edge(ai,aj)
                            CC=[g.subgraph(c).copy() for c in nx.connected_components(g)]
                            assert len(CC)==2,f'Bond {ai.serial-1}-{aj.serial-1} when cut makes more than 2 components'
                            for sg in CC:
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
                                        vt.addline(f'lappend rbonds({i}) [list {ai.serial-1} {aj.serial-1}]')
                                        logger.debug(f'  -> movers: {" ".join([str(x) for x in sg])}')
                                        vt.addline(f'lappend movers({i}) [list {mover_indices}]')
                            g.add_edge(ai,aj)
                    i+=1
        return i

    def declash_glycans(self,specs):
        """
        Declash glycans in the base molecule using the custom ``PestiferDeclash`` TcL package.
        This method generates a VMD script to identify and declash glycans in the molecular structure.
        It uses the ``PestiferDeclash`` package to perform the declashing operation,
        which involves identifying and modifying the coordinates of atoms in glycan regions.
        """
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
        vt=self.scripters['vmd']
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
        logger.debug(f'ingesting {psf}')
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

    def ingest_molecules(self):
        """
        Ingests the base molecule from the specifications provided in the task.
        This method initializes the base molecule based on the source specifications,
        which can be a PDB file, a prebuilt PSF/PDB pair, or an AlphaFold model.
        It also handles any graft sources specified in the sequence modifications and
        activates the biological assembly of the base molecule."""

        # the prior task registered one of
        # - base_pdb, or
        # - base_cif, or 
        # - continuation pdb, psf, xsc files.
        artifacts=[p for p in self.ctx.artifacts if p.produced_by==self and p.key in ['base_coordinates','continuation_pdb','continuation_psf','continuation_xsc','continuation_vel']]
        this_source={}
        if len(artifacts)>0:
            logger.debug(f'Found {len(artifacts)} artifacts to ingest')
            for a in artifacts:
                if a.key=='base_coordinates':
                    basename,ext=os.path.splitext(a.value.path)
                    this_source['id']=basename
                    if ext.lower()=='.pdb':
                        this_source['file_format']='PDB'
                    elif ext.lower()=='.cif':
                        this_source['file_format']='mmCIF'
                    else:
                        raise RuntimeError(f'Unknown file format {ext} for base_coordinates artifact {a.value.path}')
                    break
                elif 'continuation' in a.key:
                    if not 'prebuilt' in this_source:
                        this_source['prebuilt']={}
                    if a.key=='continuation_pdb':
                        this_source['prebuilt']['pdb']=str(a.value.path)
                    elif a.key=='continuation_psf':
                        this_source['prebuilt']['psf']=str(a.value.path)
                    elif a.key=='continuation_xsc':
                        this_source['prebuilt']['xsc']=str(a.value.path)
                    elif a.key=='continuation_vel':
                        this_source['prebuilt']['vel']=str(a.value.path)
        else:
            raise RuntimeError(f'No artifacts found to ingest for {self.__class__.__name__}. Expected base_coordinates, or continuation_pdb/psf/xsc/vel files.')
        specs=self.specs
        self.source_specs=specs['source']
        assert not 'id' in self.source_specs,f'Version 2.0+ of Pestifer does not support "id" in source specs.  Psfgen task must inherit "base_coordinates" artifact or continuation artifacts from a prior task (fetch or continuation).'
        self.source_specs.update(this_source)
        logger.debug(f'User-input modspecs {self.specs["mods"]}')
        self.objmanager=ObjManager(self.specs['mods'])
        seqmods=self.objmanager.get('seq',{})
        logger.debug(f'ingesting seqmods {seqmods}')
        if 'grafts' in seqmods:
            logger.debug(f'looking for graft sources to ingest')
            Grafts=seqmods['grafts']
            for g in Grafts:
                if not g.source_pdbid in self.molecules:
                    logger.debug(f'ingesting graft source {g.source_pdbid}')
                    this_source={
                        'id':g.source_pdbid,
                        'file_format':'PDB'
                    }
                    self.molecules[g.source_pdbid]=Molecule(source=this_source)
                    self.ctx.register(
                        key=f'graft-{g.source_pdbid}',
                        value=self.molecules[g.source_pdbid],
                        value_type=Molecule,
                        produced_by=self,
                        type='instance',
                        propagate=True
                    )
                g.activate(deepcopy(self.molecules[g.source_pdbid]))
        self.chainIDmanager=ChainIDManager(
            format=self.source_specs['file_format'],
            transform_reserves=self.source_specs.get('transform_reserves',{}),
            remap=self.source_specs.get('remap_chainIDs',{}))
        self.base_molecule=Molecule(source=self.source_specs,
                                    objmanager=self.objmanager,
                                    chainIDmanager=self.chainIDmanager).activate_biological_assembly(self.source_specs['biological_assembly'])
        # register self.base_molecule in the pipeline context
        self.ctx.register(
            key='base_molecule',
            value=self.base_molecule,
            value_type=Molecule,
            produced_by=self,
            type='instance',
            propagate=True
        )
        for molid,molecule in self.molecules.items():
            logger.debug(f'Molecule "{molid}": {molecule.num_atoms()} atoms in {molecule.num_residues()} residues; {molecule.num_segments()} segments.')

    def update_molecule(self):
        """
        Updates all segments of the base molecule based on the 
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
        # get the psf, pdb, xsc from the pipeline context
        my_artifacts={x.key:x for x in self.ctx.artifacts if x.produced_by==self and x.key in ['psf','pdb','xsc','vel']}
        psf=my_artifacts.get('psf','')
        pdb=my_artifacts.get('pdb','')
        xsc=my_artifacts.get('xsc','')
        source={
            'prebuilt': {
                'psf':psf,
                'pdb':pdb,
                'xsc':xsc
            }
        }
        if hasattr(self,'chainIDmanager') and hasattr(self,'objmanager'):
            updated_molecule=Molecule(source=source,chainIDmanager=self.chainIDmanager,objmanager=self.objmanager).activate_biological_assembly(0)
        else:
            updated_molecule=Molecule(source=source).activate_biological_assembly(0)

        self.molecules[base_key]=updated_molecule
        self.base_molecule=updated_molecule
        self.ctx.register(
            key='base_molecule',
            value=self.base_molecule,
            value_type=Molecule,
            produced_by=self,
            type='instance',
            propagate=True
        )
