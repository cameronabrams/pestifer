# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import logging
import os

from .vmdscripter import VMDScripter

from ..core.command import Command
from ..core.labels import Labels
from ..core.objmanager import ObjManager
from ..util.stringthings import ByteCollector

from ..molecule.asymmetricunit import AsymmetricUnit
from ..molecule.atom import Atom, AtomList
from ..molecule.bioassemb import BioAssemb
from ..molecule.molecule import Molecule
from ..molecule.residue import Residue, ResidueList
from ..molecule.segment import Segment, SegmentList
from ..molecule.transform import Transform

from ..objs.cfusion import Cfusion, CfusionList
from ..objs.graft import Graft, GraftList
from ..objs.link import Link, LinkList
from ..objs.mutation import Mutation, MutationList
from ..objs.patch import Patch, PatchList
from ..objs.resid import ResID
from ..objs.ssbond import SSBond, SSBondList

from ..logparsers import PsfgenLogParser
from ..util.progress import PsfgenProgress
from ..util.stringthings import my_logger
from ..util.util import reduce_intlist

logger = logging.getLogger(__name__)

class PsfgenScripter(VMDScripter):
    """
    This class extends the VMDScripter class to provide functionality for creating and managing psfgen scripts.

    Parameters
    ----------
    config : Config
        The configuration object containing settings for the script.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.charmmff = kwargs.get('charmmff_content')
        self.charmmff_config = kwargs.get('charmmff_config')
        self.psfgen_config = kwargs.get('psfgen_config')
        self.postregencommands = ByteCollector()

    def addpostregenerateline(self,line):
        """
        Add a line to the post-regenerate commands section of the psfgen script.

        Parameters
        ----------
        line : str
            The line of Tcl code to be added to the post-regenerate commands.
        """
        self.postregencommands.addline(line)

    def fetch_standard_charmm_topologies(self):
        """
        Fetch the standard CHARMM topologies from the configuration and copy them to the local directory.
        This method retrieves the topologies defined in the CHARMM force field configuration and copies them to 
        the local directory for use in the psfgen script.

        Returns
        -------
        list
            A list of topology files that have been copied to the local directory.
        """
        topology_local = []
        # topologies are assumed correctly ordered in config file
        for t in self.charmmff_config['standard']['rtf'] + self.charmmff_config['custom']['rtf'] + self.charmmff_config['standard']['str'] + self.charmmff_config['custom']['str']:
            if not t in topology_local:
                self.charmmff.copy_charmmfile_local(t)
                topology_local.append(t)
        # the order of the topologies is imporant, and the top_all35_ethers.rtf should be last
        if 'top_all35_ethers.rtf' in topology_local:
            topology_local.remove('top_all35_ethers.rtf')
            topology_local.append('top_all35_ethers.rtf')
        logger.debug(f'Local CHARMMFF topologies:')
        my_logger(topology_local, logger.debug)
        return topology_local

    def newscript(self, basename: str = '', packages: list[str] = [], additional_topologies: list[str] = []):
        """
        Initialize a new psfgen script with a specified basename and packages.
        If no basename is provided, a default script name is used.

        Parameters
        ----------
        basename : str, optional
            The base name for the script file. If not provided, a default name is used.
        packages : list, optional
            A list of Tcl packages to be imported in the script.
        additional_topologies : list, optional
            A list of additional topology files to be included in the script. These files should not contain
            a path (i.e., they should be in the current working directory).  These can be either 'rtf' 
            or 'str' extension files.
        """
        super().newscript(basename=basename, packages=packages)
        self.addline('package require psfgen')
        self.addline('psfcontext mixedcase')
        self.topologies = self.fetch_standard_charmm_topologies()
        for at in sorted(additional_topologies):
            assert os.sep not in at, f'Topology file {at} must not contain a path.'
            if not at in self.topologies:
                self.charmmff.copy_charmmfile_local(at)
                self.topologies.append(at)
        # self.topologies=list(set(self.topologies))  # remove duplicates
        # if 'top_all35_ethers.rtf' in self.topologies:
        #     self.topologies.remove('top_all35_ethers.rtf')
        #     self.topologies.append('top_all35_ethers.rtf')
        for t in self.topologies:
            self.addline(f'topology {t}')
            # self.addfile(t) # appends this file to the scripters FileCollector for later cleanup
        # logger.debug(f'psfgen aliases: {self.psfgen_config["aliases"]}')
        for alias_type,alias_list in self.psfgen_config['aliases'].items():
            logger.debug(f'Adding {len(alias_list)} {alias_type} aliases to psfgen script')
            for pdba in alias_list:
                alias_tokens = pdba.split()
                if alias_tokens[1] == '*':  # wild-card for all protein residues
                    for protein_resname in self.psfgen_config['segtypes']['protein']['resnames']:
                        this_pdba = f'{alias_tokens[0]} {protein_resname} {alias_tokens[2]} {alias_tokens[3]}'
                        self.addline(f'pdbalias {alias_type} {this_pdba}')
                else:
                    self.addline(f'pdbalias {alias_type} {pdba}')

    def load_project(self,*objs):
        """
        Load a project into psfgen by reading the PSF and PDB files.
        This method takes one or two objects (PSF and PDB files) and loads them into psfgen.
        If only one object is provided, it is assumed to be the basename of both the PSF and PDB files.
        If two objects are provided, the first is treated as the PSF file and the second as the PDB file.
        
        Parameters
        ----------
        objs : tuple
            A tuple containing either one or two objects:

            - If one object, it should be the basename of both the PSF and PDB files.
            - If two objects, the first should be the PSF file and the second should be the PDB file.
        
        """
        logger.debug(f'load_project called with {len(objs)} arguments')
        for _ in objs:
            logger.debug(f'   {_}')
        if len(objs) == 1:
            basename = objs[0]
            self.addline(f'readpsf {basename}.psf pdb {basename}.pdb')
        else:
            psf, pdb = objs
            self.addline(f'readpsf {psf} pdb {pdb}')

    def describe_molecule(self, mol: Molecule):
        """
        Describe a molecule in the psfgen script.
        This method sets up the molecule in the psfgen script by writing the necessary Tcl commands
        to create a new molecule and load its topology and coordinates.
        
        Parameters
        ----------
        mol : Molecule
            The molecule object to be described in the psfgen script.
        """
        self.molid_varname = f'm{mol.molid}'
        self.addline(f'mol top ${self.molid_varname}')
        self.write_mol_tcl(mol)

    def write_mol_tcl(self, mol: Molecule):
        au: AsymmetricUnit = mol.asymmetric_unit
        segments: SegmentList = au.segments
        topomods = au.objmanager.get('topol', {})
        ssbonds: SSBondList = topomods.get('ssbonds', [])
        links: LinkList = topomods.get('links', [])
        ba: BioAssemb = mol.active_biological_assembly
        for transform in ba.transforms.data:
            self.banner(f'Transform {transform.index} begins')
            self.banner('The following mappings of A.U. asym ids is used:')
            for k,v in transform.chainIDmap.items():
                self.comment(f'A.U. chain {k}: Image chain {v}')
            self.banner('Segments follow')
            self.write_segments(segments, transform)
            self.banner('DISU patches follow')
            if ssbonds:
                self.write_ssbonds(ssbonds, transform)
            self.banner('LINK patches follow')
            if links:
                self.write_links(links, transform)
            self.banner(f'Transform {transform.index} ends')

    def write_segments(self, segments: SegmentList, transform: Transform = None):
        for segment in segments.data:
            self.comment(f'Segment {segment.segname} begins')
            self.write_segment(segment, transform)
            self.comment(f'Segment {segment.segname} ends')

    def write_segment(self, segment: Segment, transform: Transform = None):
        """
        Write a segment to the psfgen script.
        This method generates the Tcl commands to define a segment in the psfgen script,
        including its name, type, and any associated patches.

        Parameters
        ----------
        segment : Segment
            The segment object to be written to the psfgen script.
        transform : Transform, optional
            A transform object that may be applied to the segment. Default is None.
        """
        
        if segment.segtype == 'protein':
            self.write_polymer_stanza(segment, transform)
        elif segment.segtype == 'nucleicacid':
            self.write_nucleicacid_stanza(segment, transform)
        elif segment.segtype == 'glycan':
            self.write_glycan_stanza(segment, transform)
        else:
            self.write_generic_stanza(segment, transform)

    def write_polymer_stanza(self, segment: Segment, transform: Transform = None):
        """
        Write a polymer segment stanza to the psfgen script.
        This method generates the Tcl commands to define a polymer segment in the psfgen script,
        including its name, type, and any associated patches.

        Parameters
        ----------
        segment : Segment
            The polymer segment object to be written to the psfgen script.
        transform : Transform, optional
            A transform object that may be applied to the segment. Default is None.
        """
        
        chainIDmap: dict = transform.chainIDmap
        seglabel: str = segment.segname
        image_seglabel: str = chainIDmap.get(seglabel, seglabel)
        is_image: bool = image_seglabel != seglabel
        segtype: str = segment.segtype
        objmanager: ObjManager = segment.objmanager
        seqmods: dict = objmanager.get('seq', {})
        logger.debug(f'polymer_stanza {segtype} {seglabel}->{image_seglabel}')

        transform.register_mapping(segtype, image_seglabel, seglabel)

        loopspecs: dict = segment.specs.get('loops', {})
        sac_rn: str = loopspecs.get('sac_res_name', 'NOSAC')
        min_loop_length: int = loopspecs.get('min_loop_length', 0)
        build_all_terminal_loops: bool = segment.specs.get('include_terminal_loops', False)
        build_N_terminal_loop: bool = seglabel in segment.specs.get('build_zero_occupancy_N_termini', [])
        build_C_terminal_loop: bool = seglabel in segment.specs.get('build_zero_occupancy_C_termini', [])

        seg_mutations: MutationList = seqmods.get('mutations', MutationList([]))
        logger.debug(f'polymer_stanza for {segtype} segname {seglabel}; init mutations:')
        for m in seg_mutations:
            logger.debug(str(m))
        seg_Cfusions: CfusionList = seqmods.get('Cfusions', CfusionList([]))
        seg_patches: PatchList = seqmods.get('patches', PatchList([]))
        self.banner(f'Segment {image_seglabel} begins')
        logger.debug(f'Segment {image_seglabel} begins')
        for sf in seg_Cfusions.data:
            self.write_cfusion_presegment(sf)
        for i, b in enumerate(segment.subsegments.data):
            if b.state == 'RESOLVED':
                """ for a resolved subsegment, generate its pdb file """
                b.selname = f'{image_seglabel}{i:02d}'
                run = ResidueList(segment.residues[b.bounds[0]:b.bounds[1] + 1])
                b.pdb = f'segtype_polymer_{image_seglabel}_{run.data[0].resid.resid}_to_{run.data[-1].resid.resid}.pdb'
                logger.debug(f'Writing resolved subsegment {repr(b)} to {b.pdb}')
                self.addfile(b.pdb)
                serial_list = run.atom_serials(as_type=int)
                logger.debug(f'Last atom has serial {serial_list[-1]} ({run.data[-1].resname}{run.data[-1].resid.resid}):')
                at: Atom = segment.parent_molecule.asymmetric_unit.atoms.get(lambda x: x.serial == serial_list[-1])
                if hasattr(at, '__len__'):
                    for a in at:
                        logger.debug(f'whoops: {a.chainID} {a.resid.resid} {a.name} is {a.serial}')
                    raise Exception(f'More than one atom with serial {serial_list[-1]}??')
                logger.debug(f'    {at.serial} {at.resname} {at.name} in chain {at.chainID} residue {at.resname}{at.resid.resid}')
                assert at.resid == run.data[-1].resid
                vmd_red_list = reduce_intlist(serial_list)
                self.addline(f'set {b.selname} [atomselect $m{segment.parent_molecule.molid} "serial {vmd_red_list}"]')
                self.addline(f'${b.selname} set segname {image_seglabel}')
                if len(at.ORIGINAL_ATTRIBUTES) > 0:
                    if at.ORIGINAL_ATTRIBUTES["serial"] != at.serial:
                        self.banner(f'Atom with serial {at.ORIGINAL_ATTRIBUTES["serial"]} in PDB needs serial {at.serial} for VMD')
                """ Relabel chain ID and request coordinate transformation """
                if is_image:
                    self.backup_selection(b.selname, dataholder=f'{b.selname}_data')
                    self.addline(f'${b.selname} set chain {image_seglabel}')                 
                    self.addline(f'${b.selname} move {transform.write_TcL()}')
                self.addline(f'${b.selname} writepdb {b.pdb}')
            elif b.state == 'MISSING':
                if i == 0:
                    if build_N_terminal_loop or build_all_terminal_loops:
                        b.declare_buildable()
                elif i == (len(segment.subsegments) - 1):
                    if build_C_terminal_loop or build_all_terminal_loops:
                        b.declare_buildable()
                else:
                    b.declare_buildable()
        self.addline(f'segment {image_seglabel} '+'{')
        if segment.subsegments.data[0].state == 'MISSING' and not segment.subsegments.data[0].build:
            Nterminal_missing_subsegment = segment.subsegments.pop(0)
            logger.debug(f'Since terminal loops are not included, ignoring {str(Nterminal_missing_subsegment)}')
        if segment.subsegments.data[-1].state == 'MISSING' and not segment.subsegments.data[-1].build:
            Cterminal_missing_subsegment = segment.subsegments.pop(-1)
            logger.debug(f'Since terminal loops are not included, ignoring {str(Cterminal_missing_subsegment)}')
        for b in segment.subsegments.data:
            if b.state == 'RESOLVED':
                self.addline(f'pdb {b.pdb}', indents=1)
            elif b.state == 'MISSING' and b.build:
                for r in segment.residues[b.bounds[0]:b.bounds[1]+1]:
                    rname = Labels.charmm_resname_of_pdb_resname.get(r.resname, r.resname)
                    self.addline(f'residue {r.resid.resid} {rname} {image_seglabel}', indents=1)
                if b.num_items() >= min_loop_length and not b in [segment.subsegments[0], segment.subsegments[-1]]:
                    lrr = segment.residues[b.bounds[1]]
                    sac_resid = lrr.resid.increment()
                    b.sacres = Residue({'resname': sac_rn, 'resid': sac_resid, 'chainID': seglabel, 'segtype': segtype, 'segname': seglabel, 'resolved': False, 'atoms': AtomList([])})
                    logger.debug(f'Adding sac residue {str(b.sacres)}')
                    self.addline(f'residue {sac_resid.resid} {sac_rn} {image_seglabel}', indents=1)
        for cf in seg_Cfusions.data:
            self.write_cfusion_insegment(cf)
        self.write_mutations(seg_mutations)
        for p in seg_patches.data:
            if p.use_in_segment == 'first':
                self.addline(f'first {p.patchname}', indents=1)
            elif p.use_in_segment == 'last':
                self.addline(f'last {p.patchname}', indents=1)
        self.addline('}')
        self.banner(f'End segment {image_seglabel}')
        self.banner('Coordinate-specification commands')
        for b in segment.subsegments.data:
            if b.state == 'RESOLVED':
                self.comment(f'Subsegment {[segment.subsegments.index(b)]} is a resolved run')
                self.addline(f'coordpdb {b.pdb} {image_seglabel}')
        for b in segment.subsegments.data:
            if b.state == 'MISSING' and b.build:
                if segment.subsegments.index(b) > 0:  # only seed orientation for a loop that is not at the N-terminus
                    self.comment(f'Subsegment {[segment.subsegments.index(b)]}/{len(segment.subsegments)} is a missing loop')
                    this_run = ResidueList(segment.residues[b.bounds[0]:b.bounds[1] + 1])
                    prior_b = segment.subsegments[segment.subsegments.index(b) - 1]
                    self.comment(f'...attached to subsegment {segment.subsegments.index(prior_b)}')
                    prior_run = ResidueList(segment.residues[prior_b.bounds[0]:prior_b.bounds[1] + 1])
                    if segtype =='protein':
                        self.comment(f'Seeding orientation of model-built loop starting at {str(this_run[0])} from {str(prior_run[-1])}')
                        self.addline(f'{this_run.caco_str(prior_run, image_seglabel, transform.tmat)}')
        if len(seg_patches) > 0:
            self.banner(f'Patches for segment {image_seglabel}')
        for patch in seg_patches.data:
            if patch.use_in_segment == '' and not patch.use_after_regenerate:
                self.addline(f'patch {patch.patchname} {image_seglabel}:{patch.resid.resid}')
            if patch.use_after_regenerate:
                self.addpostregenerateline(f'patch {patch.patchname} {image_seglabel}:{patch.resid.resid}')
        for sc in seg_Cfusions:
            self.write_cfusion_postsegment(sc)
        self.banner('Intra-segmental terminal patches')
        if segtype == 'protein':
            for i, b in enumerate(segment.subsegments):
                # only non-terminal loops get the terminal patches
                if b.state == 'MISSING' and 0 < i < (len(segment.subsegments) - 1) and b.sacres is not None:
                    Cterm = segment.residues[b.bounds[1]]
                    self.addline(f'patch CTER {image_seglabel}:{Cterm.resid.resid}')
                    nextb = segment.subsegments.data[i + 1]
                    Nterm = segment.residues.data[nextb.bounds[0]]
                    patchname = 'NTER'
                    if Nterm.resname == 'PRO':
                        patchname = 'PROP'
                    elif Nterm.resname == 'GLY':
                        patchname = 'GLYP'
                    self.addline(f'patch {patchname} {image_seglabel}:{Nterm.resid.resid}')
                    logger.debug(f'deleting sacrificial residue {str(b.sacres)}')
                    self.addline(f'delatom {image_seglabel} {b.sacres.resid.resid}')
        self.banner('Restoring A.U. state for all resolved subsegments')
        for b in segment.subsegments.data:
            if b.state == 'RESOLVED':
                if is_image:
                    self.restore_selection(b.selname, dataholder=f'{b.selname}_data')
        self.banner(f'Segment {image_seglabel} ends')

    def write_glycan_stanza(self, segment: Segment, transform: Transform = None):
        """
        Write a glycan segment stanza to the psfgen script.
        This method generates the Tcl commands to define a glycan segment in the psfgen script,
        including its name, type, and any associated patches.

        Parameters
        ----------
        segment : Segment
            The glycan segment object to be written to the psfgen script.
        transform : Transform, optional
            A transform object that may be applied to the segment. Default is None.
        """
        self.write_generic_stanza(segment, transform)
    
    def write_nucleicacid_stanza(self, segment: Segment, transform: Transform = None):
        """
        Write a nucleic acid segment stanza to the psfgen script.
        This method generates the Tcl commands to define a nucleic acid segment in the psfgen script,
        including its name, type, and any associated patches.

        Parameters
        ----------
        segment : Segment
            The nucleic acid segment object to be written to the psfgen script.
        transform : Transform, optional
            A transform object that may be applied to the segment. Default is None.
        """
        self.write_polymer_stanza(segment, transform)

    def write_generic_stanza(self, segment: Segment, transform: Transform = None):
        """
        Write a generic segment stanza to the psfgen script.
        This method generates the Tcl commands to define a generic segment in the psfgen script,
        including its name, type, and any associated patches.

        Parameters
        ----------
        segment : Segment
            The generic segment object to be written to the psfgen script.
        transform : Transform, optional
            A transform object that may be applied to the segment. Default is None.
        """
        chainIDmap: dict = transform.chainIDmap if transform else {}
        seglabel: str = segment.segname
        image_seglabel: str = chainIDmap.get(seglabel, seglabel)
        is_image: bool = image_seglabel != seglabel
        segtype: str = segment.segtype
        objmanager: ObjManager = segment.objmanager
        seqmods: dict = objmanager.get('seq', {})
        topomods: dict = objmanager.get('topol', {})
        seg_grafts: GraftList = seqmods.get('grafts', GraftList([]))

        transform.register_mapping(segment.segtype, image_seglabel, seglabel)

        self.banner(f'Segment {image_seglabel} begins as image of {seglabel}')
        for g in seg_grafts:
            # g.write_pre_segment(W)
            self.write_graft_presegment(g)
        serial_list = segment.residues.atom_serials(as_type=int)
        resid_list = segment.residues.atom_resids(as_type=ResID)
        vmd_red_list = reduce_intlist(serial_list)
        pdb = f'segtype_generic_{image_seglabel}.pdb'
        selname = image_seglabel
        self.addfile(pdb)  # appends this file to the scripters FileCollector for later cleanup
        self.addline(f'set {selname} [atomselect $m{segment.parent_molecule.molid} "serial {vmd_red_list}"]')
        self.addline(f'${selname} set segname {image_seglabel}')
        if is_image:
            self.backup_selection(selname, dataholder=f'{selname}_data')
            self.addline(f'${selname} set chain {image_seglabel}')
            self.addline(f'${selname} move {transform.write_TcL()}')
        self.addline(f'${selname} set resid [list {" ".join([str(x) for x in resid_list])}]')
        self.addline(f'${selname} writepdb {pdb}')
        self.addline(f'segment {image_seglabel} '+'{')
        self.addline(f'first none', indents=1)
        self.addline(f'last none', indents=1)
        self.addline(f'pdb {pdb}', indents=1)
        for g in seg_grafts.data:
            # g.write_in_segment(self)
            self.addline (f'    pdb {g.segfile}')
        self.addline('}')
        self.addline(f'coordpdb {pdb} {image_seglabel}')
        for g in seg_grafts.data:
            # g.write_post_segment(self)
            self.addline(f'coordpdb {g.segfile} {g.chainID}')
        if is_image:
            self.banner(f'Restoring A.U. state for {seglabel}')
            self.restore_selection(selname, dataholder=f'{selname}_data')
        self.banner(f'Segment {image_seglabel} ends')

    def write_graft_presegment(self, G: Graft):
        """
        Writes the Tcl commands to create a graft segment in the Psfgen script.

        Parameters
        ----------
        G : Graft
            The Graft object representing the graft segment to be created.
            This method generates the Tcl commands to create a graft segment in the Psfgen script,
            including the selection of atoms and writing the segment to a PDB file.
        """

        self.comment(f'{str(G)}')
        self.addline(f'set topid [molinfo ${self.molid_varname} get id]')
        if os.path.exists(f'{G.source_pdbid}.pdb'):
            self.addline(f'mol new {G.source_pdbid}.pdb')
        else:
            self.addline(f'mol new {G.source_pdbid}')
        self.addline(f'set graftid [molinfo top get id]')
        self.addline(f'mol top $topid')

        alignment_target_resid_logic = f'(resid {G.target_root.resid}'
        if G.target_partner is not None:
            alignment_target_resid_logic += f' or resid {G.target_partner.resid}'
        alignment_target_resid_logic += ') and noh'

        self.addline(f'set target_sel [atomselect $topid "chain {G.residues[0].asym_chainID} and {alignment_target_resid_logic}"]')

        alignment_source_resid_logic = f'(resid {G.source_root.resid}'
        if G.source_partner is not None:
            alignment_source_resid_logic += f' or resid {G.source_partner.resid}'
        # if G.source_end is not None:
        #     alignment_source_resid_logic += f' or resid {G.source_end.resid}'
        alignment_source_resid_logic += ') and noh'

        self.addline(f'set source_sel [atomselect $graftid "chain {G.source_chainID} and {alignment_source_resid_logic}"]')
        self.addline(f'vmdcon -info "[$source_sel num] atoms in source, [$target_sel num] atoms in target"')

        # Now we need to select the residues that will be moved
        # The mover residues are those that are in the source segment *excluding* those used in the alignment

        movers_source_resid_logic = f'(not (resid {G.source_root.resid}'
        if G.source_partner is not None:
            movers_source_resid_logic += f' or resid {G.source_partner.resid}'
        movers_source_resid_logic += ')'
        if G.source_end is not None:
            movers_source_resid_logic += f' and resid <= {G.source_end.resid}'
        movers_source_resid_logic += ')'

        self.addline(f'set mover_sel [atomselect $graftid "chain {G.source_chainID} and {movers_source_resid_logic}"]')
        self.addline(f'vmdcon -info "[$mover_sel num] atoms will be moved"')

        # Now we need to align the mover residues to the target residues
        # We will use the transidentity command to get the transformation matrix
        # and then apply it to the mover residues

        # self.addline(f'set TT [transidentity]')
        self.addline(f'set TT [measure fit $source_sel $target_sel]')
        self.addline(f'vmdcon -info "Homog. trans. matrix: $TT"')
        self.addline(f'$mover_sel move $TT')
        G.segfile=f'graft{G.obj_id}.pdb'
        new_residlist=[]
        for y in G.donor_residues:
            new_residlist.extend([f'{y.resid.resid}' for x in y.atoms])  # r
        self.addline(f'$mover_sel set resid [list {" ".join(new_residlist)}]')
        self.addline(f'$mover_sel set chain {G.donor_residues[0].chainID}')
        self.addline(f'$mover_sel writepdb {G.segfile}')
        self.addfile(G.segfile)
        self.addline(f'$mover_sel delete')

    def write_cfusion_presegment(self, C: Cfusion):
        """
        Writes the Tcl commands to create a fusion segment in the Psfgen script.

        Parameters
        ----------
        C : Cfusion
            The Cfusion object representing the fusion segment to be created.
            This method generates the Tcl commands to create a fusion segment in the Psfgen script,
            including the selection of atoms and writing the segment to a PDB file.
        """

        self.addline(f'set topid [molinfo top get id]')
        self.addline(f'mol new {C.sourcefile}')
        self.addline(f'set cfusid [molinfo top get id]')
        self.addline(f'mol top $topid')
        self.addline(f'set fusres [atomselect $cfusid "protein and chain {C.sourceseg} and resid {C.resid1.resid} to {C.resid2.resid}"]')
        C.segfile=f'Cfusion{C.obj_id}_{C.sourceseg}_{C.resid1.resid}_to_{C.resid2.resid}.pdb'
        self.addline(f'$fusres writepdb {C.segfile}')
        self.addline(f'delete $cfusid')

    def write_cfusion_insegment(self, C: Cfusion):
        """
        Writes the Tcl commands to add the fusion into the active segment of the base molecule.

        Parameters
        ----------
        C : Cfusion
            The Cfusion object representing the fusion segment to be added.
        """
        self.addline(f'    pdb {C.segfile}')

    def write_cfusion_postsegment(self, C: Cfusion):
        """
        Writes the Tcl commands to finalize the fusion segment in the Psfgen script.

        Parameters
        ----------
        C : Cfusion
            The Cfusion object representing the fusion segment to be finalized.
        """
        self.addline(f'coordpdb {C.segfile} {C.chainID}')

    def write_mutations(self, mutations: MutationList):
        for mutation in mutations:
            self.write_mutation(mutation)

    def write_mutation(self, mutation: Mutation):
        self.addline(f'mutate {mutation.resid.resid} {mutation.newresname}', indents=1)

    def write_ssbonds(self, ssbonds: SSBondList, transform: Transform = None):
        for S in ssbonds:
            self.write_ssbond(S, transform)

    def write_ssbond(self, S: SSBond, transform: Transform = None):
        """
        Writes the Tcl commands to create a disulfide bond in the Psfgen script.

        Parameters
        ----------
        S : SSBond
            The SSBond object representing the disulfide bond to be created.
            This method generates the Tcl commands to create a disulfide bond in the Psfgen script,
            including the selection of atoms and writing the bond to a PDB file.
        """
        chainIDmap=transform.chainIDmap 
        # ok since these are only going to reference protein segments; protein segment names are the chain IDs
        logger.debug(f'writing patch for {str(S)}')
        c1 = chainIDmap.get(S.chainID1, S.chainID1)
        c2 = chainIDmap.get(S.chainID2, S.chainID2)
        r1 = S.resid1.resid
        r2 = S.resid2.resid
        self.addline(f'patch DISU {c1}:{r1} {c2}:{r2}')

    def write_links(self, links: LinkList, transform: Transform = None):
        for L in links:
            self.write_link(L, transform)

    def write_link(self, L: Link, transform: Transform = None):
        """
        Writes the Tcl commands to create a link in the Psfgen script.

        Parameters
        ----------
        L : Link
            The Link object representing the link to be created.
            This method generates the Tcl commands to create a link in the Psfgen script,
            including the selection of atoms and writing the link to a PDB file.
        """
        chainIDmap = transform.chainIDmap
        seg1 = L.residue1.chainID
        seg1 = chainIDmap.get(seg1, seg1)
        seg2 = L.residue2.chainID
        seg2 = chainIDmap.get(seg2, seg2)
        resid1 = L.residue1.resid
        resid2 = L.residue2.resid
        logger.debug(f'Link: {L.residue1.chainID}->{seg1}:{resid1} {L.residue2.chainID}->{seg2}:{resid2}')
        if not L.patchname == 'UNFOUND':
            write_post_regenerate = L.patchname in Patch._after_regenerate_patches
            if L.patchhead == 1:
                if write_post_regenerate:
                    self.addpostregenerateline(f'patch {L.patchname} {seg1}:{resid1} {seg2}:{resid2}')
                else:
                    self.addline(f'patch {L.patchname} {seg1}:{resid1} {seg2}:{resid2}')
            elif L.patchhead == 2:
                if write_post_regenerate:
                    self.addpostregenerateline(f'patch {L.patchname} {seg2}:{resid2} {seg1}:{resid1}')
                else:
                    self.addline(f'patch {L.patchname} {seg2}:{resid2} {seg1}:{resid1}')
        else:
            logger.warning(f'Could not identify patch for link: {str(L)}')
            self.comment(f'No patch found for {str(L)}')


    def writescript(self, statename, guesscoord=True, regenerate=True, writepsf=True, writepdb=True, force_exit=False):
        """
        Finalize the psfgen script by adding commands to write the PSF and PDB
        files, and optionally regenerate angles and dihedrals, and guess coordinates.
        This method checks if an exit statement is already present in the script; if not, it adds one.

        Parameters
        ----------
        statename : str
            The name of the state to be used for writing the PSF and PDB files.
        guesscoord : bool, optional
            If True, the script will include a command to guess coordinates. Default is True.
        regenerate : bool, optional
            If True, the script will include a command to regenerate angles and dihedrals. Default is True.
        writepsf : bool, optional
            If True, the script will include a command to write the PSF file. Default is True.
        writepdb : bool, optional
            If True, the script will include a command to write the PDB file. Default is True.
        force_exit : bool, optional
            If True, the script will include an exit statement even if one is already present. Default is False.
        """
        if guesscoord:
            self.addline('guesscoord')
        if regenerate:
            self.addline('regenerate angles dihedrals')
        if len(self.postregencommands.byte_collector)>0:
            self.B.write(self.postregencommands.byte_collector)
        if writepsf:
            self.addline(f'writepsf cmap {statename}.psf')
        if writepdb:
            self.addline(f'writepdb {statename}.pdb')
        super().writescript(force_exit=force_exit)

    def runscript(self, *args, **options):
        """
        Run the psfgen script using the VMD command line interface.
        This method constructs a command to execute VMD with the specified script and options.
        
        Parameters
        ----------
        args : tuple
            Additional arguments to be passed to the VMD command.
        options : dict
            A dictionary of options to be passed to the VMD command.
        
        Returns
        -------
        Command
            The command object that was constructed to run the psfgen script, after it has been executed.
        """
        assert hasattr(self, 'scriptname'), f'No scriptname set.'
        self.logname = f'{self.basename}.log'
        self.logparser = PsfgenLogParser(basename=self.basename)
        logger.debug(f'Log file: {self.logname}')
        clean_options = options.copy()
        for k, v in options.items():
            if v == '':
                clean_options.pop(k)
        c = Command(f'{self.vmd} -dispdev text -startup {self.vmd_startup} -e {self.scriptname} -args --tcl-root {self.tcl_root}', **clean_options)
        progress_struct = None
        if self.progress:
            progress_struct = PsfgenProgress()
            self.logparser.enable_progress_bar(progress_struct)
        else:
            logger.debug('Progress bar is disabled for psfgen script')
        result = c.run(logfile=self.logname, logparser=self.logparser)
        logger.debug(f'FileCollector:')
        my_logger(self.F, logger.debug)
        if not options.get('keep_tempfiles', False):
            self.F.flush()
        return result

