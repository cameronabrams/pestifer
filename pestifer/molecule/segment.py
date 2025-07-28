#Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
Defines the Segment class for generating and managing segments declared for psfgen
"""
from functools import singledispatchmethod
from typing import ClassVar, Optional, List, Any
from pydantic import Field, ConfigDict, PrivateAttr
import logging
logger=logging.getLogger(__name__)
from ..core.baseobj_new import BaseObj, BaseObjList
from ..objs.mutation import MutationList
from ..objs.cfusion import CfusionList
from ..objs.deletion import DeletionList
from ..objs.graft import GraftList
from ..objs.patch import PatchList
from ..core.labels import Labels
from ..core.objmanager import ObjManager
from .residue import Residue,ResidueList
from ..util.util import reduce_intlist
from ..core.scripters import PsfgenScripter
from .chainidmanager import ChainIDManager

class Segment(BaseObj):
    """
    A class for handling segments in a molecular structure.
    This class represents a segment defined by a list of residue indices and provides methods to check if a bond intersects the segment.
    It also provides methods to yield treadmilled versions of the segment and to check for equality with another segment.
    """
    model_config = ConfigDict(arbitrary_types_allowed=True)
    _required_fields = {'segtype','segname','chainID','residues','subsegments','parent_chain','specs'}
    _optional_fields = {'mutations','deletions','grafts','patches','attachments','psfgen_segname'}
    """
    Optional attributes for ``Segment``.
    """

    segtype: str = Field(..., description="The type of the segment (e.g., 'protein', 'nucleicacid', 'glycan').")
    segname: str = Field(..., description="The name of the segment as it is understood by psfgen and in a PSF file.")
    chainID: str = Field(..., description="The chain ID associated with the segment.")
    residues: ResidueList = Field(default_factory=ResidueList, description="A list of residues that make up the segment.")
    subsegments: list = Field(..., description="A list of subsegments within the segment, each with a state indicating whether it is resolved or missing.")
    parent_chain: str = Field(..., description="The chain ID of the parent chain to which this segment belongs.")
    specs: dict = Field(..., description="Specifications for the segment, including details like loops and terminal modifications.")

    mutations: Optional[MutationList] = Field(default=None, description="A list of mutations applied to the segment.")
    deletions: Optional[DeletionList] = Field(default=None, description="A list of deletions applied to the segment.")
    grafts: Optional[GraftList] = Field(default=None, description="A list of grafts applied to the segment.")
    patches: Optional[PatchList] = Field(default=None, description="A list of patches applied to the segment.")
    attachments: Optional[List] = Field(default=None, description="A list of attachments applied to the segment (Attachments are currently unimplemented).")
    psfgen_segname: Optional[str] = Field(default=None, description="The segment name as understood by psfgen.")

    _inheritable_objs: ClassVar[set] = {'mutations','patches','Cfusions','grafts'}
    """
    List of inheritable modifications for the segment.
    These modifications include mutations, patches, C-terminal fusions, and grafts.
    """

    _obj_manager: ObjManager = PrivateAttr(default_factory=ObjManager)
    """
    The object manager for the segment.
    This manager is used to track and manage the objects associated with the segment.
    """

    _parent_molecule: Any = PrivateAttr()

    class Adapter:
        def __init__(self, segtype: str = 'protein', segname: str = 'UNSET', chainID: str = 'UNSET', residues: ResidueList = None, subsegments: list = None, parent_chain: str = 'UNSET', specs: dict = None):
            self.segtype = segtype
            self.segname = segname
            self.chainID = chainID
            self.residues = residues if residues is not None else ResidueList([])
            self.subsegments = subsegments if subsegments is not None else []
            self.parent_chain = parent_chain
            self.specs = specs if specs is not None else {}
            
        @classmethod
        def from_dict(cls, data: dict):
            """
            Create an Adapter instance from a dictionary.
            
            Parameters
            ----------
            data : dict
                A dictionary containing the segment attributes.
                
            Returns
            -------
            Adapter
                An instance of the Adapter class initialized with the provided data.
            """
            return cls(**data)

        @classmethod
        def from_residue(cls, residue: Residue, segname: str = 'UNSET', specs: dict = None):
            """
            Create an Adapter instance from a Residue object.
            """
            apparent_chainID = residue.chainID
            apparent_segtype = residue.segtype
            apparent_segname = residue.chainID if segname == 'UNSET' else segname
            myRes = ResidueList([])
            myRes.append(residue)
            subsegments = myRes.state_bounds(lambda x: 'RESOLVED' if len(x.atoms) > 0 else 'MISSING')
            input_dict = {
                'specs': specs,
                'segtype': apparent_segtype,
                'segname': apparent_segname,
                'chainID': apparent_chainID,
                'residues': myRes,
                'subsegments': subsegments,
                'parent_chain': apparent_chainID
            }
            return cls(**input_dict)
        
        @classmethod
        def from_residue_list(cls, residues: ResidueList, segname: str = 'UNSET', specs: dict = None):
            """
            Create an Adapter instance from a ResidueList object.
            """
            apparent_chainID = residues[0].chainID
            apparent_segtype = residues[0].segtype
            apparent_segname = residues[0].chainID if segname == 'UNSET' else segname
            if apparent_segtype in ['protein', 'nucleicacid']:
                # a protein segment must have unique residue numbers
                assert residues.puniq(['resseqnum','insertion']),f'ChainID {apparent_chainID} has duplicate resseqnum-insertion!'
                # a protein segment may not have more than one protein chain
                assert all([x.chainID==residues[0].chainID for x in residues])
                residues.sort()
                subsegments=residues.state_bounds(lambda x: 'RESOLVED' if x.resolved else 'MISSING')
            else:
                logger.debug(f'Calling puniqify on residues of non-protein segment {segname}')
                residues.puniquify(fields=['resseqnum','insertion'],make_common=['chainID'])
                count=sum([1 for x in residues if hasattr(x,'_ORIGINAL_ATTRIBUTES')])
                if count>0:
                    logger.debug(f'{count} residue(s) were affected by puniquify:')
                    for x in residues:
                        if hasattr(x,'_ORIGINAL_ATTRIBUTES'):
                            logger.debug(f'    {x.chainID} {x.resname} {x.resseqnum}{x.insertion} was {x._ORIGINAL_ATTRIBUTES}')
                # this assumes residues are in a linear sequence?  not really..
                subsegments=residues.state_bounds(lambda x: 'RESOLVED' if len(x.atoms)>0 else 'MISSING')
            logger.debug(f'Segment {segname} has {len(residues)} residues across {len(subsegments)} subsegments')
            input_dict={
                'specs':specs,
                'segtype': apparent_segtype,
                'segname': apparent_segname,
                'chainID': apparent_chainID,
                'residues': residues,
                'subsegments': subsegments,
                'parent_chain': apparent_chainID
            }
            return cls(**input_dict)

        def to_dict(self):
            """
            Convert the Adapter instance to a dictionary.
            
            Returns
            -------
            dict
                A dictionary representation of the Adapter instance.
            """
            return {
                'segtype': self.segtype,
                'segname': self.segname,
                'chainID': self.chainID,
                'residues': self.residues,
                'subsegments': self.subsegments,
                'parent_chain': self.parent_chain,
                'specs': self.specs,
                'residues': self.residues
            }

        def __repr__(self):
            """
            Return a string representation of the Adapter instance.
            
            Returns
            -------
            str
                A string representation of the Adapter instance.
            """
            return f"Adapter(segtype={self.segtype}, segname={self.segname}, chainID={self.chainID}, residues={len(self.residues)}, subsegments={len(self.subsegments)}, parent_chain={self.parent_chain})"

    @BaseObj.from_input.register(Adapter)
    @classmethod
    def _from_adapter(cls, adapter: Adapter):
        input_dict = adapter.to_dict()
        return cls(**input_dict)

    @singledispatchmethod
    @classmethod
    def new(cls, input_obj: Any):
        raise NotImplementedError(f'Cannot create Segment from {type(input_obj)}')
    
    @new.register(dict)
    @classmethod
    def _new_from_dict(cls, input_dict: dict):
        """
        Create a new Segment instance from a dictionary.
        
        Parameters
        ----------
        input_dict : dict
            A dictionary containing the segment attributes.
            
        Returns
        -------
        Segment
            A new Segment instance initialized with the provided attributes.
        """
        adapter = cls.Adapter.from_dict(input_dict)
        return cls._from_adapter(adapter)
    
    @new.register(Residue)
    @classmethod
    def _new_from_residue(cls, res: Residue, segname: str = 'UNSET', specs: dict = {}):
        """
        Create a new Segment instance from a Residue object.
        
        Parameters
        ----------
        res : Residue
            The residue from which to create the segment.
            
        Returns
        -------
        Segment
            A new Segment instance initialized with the residue's attributes.
        """
        adapter = cls.Adapter.from_residue(res, segname=segname, specs=specs)
        logger.debug(f'Creating Segment from residue via adapter {adapter}')
        return cls._from_adapter(adapter)

    @new.register(ResidueList)
    @classmethod
    def _new_from_residue_list(cls, residues: ResidueList = ResidueList([]), segname: str = 'UNSET', specs: dict = {}):
        """
        Create a new Segment instance from a ResidueList object.
        
        Parameters
        ----------
        residues : ResidueList
            The list of residues from which to create the segment.
            
        Returns
        -------
        Segment
            A new Segment instance initialized with the residue list's attributes.
        """
        adapter = cls.Adapter.from_residue_list(residues, segname=segname, specs=specs)
        return cls._from_adapter(adapter)

    def describe(self) -> str:
        return f'<Segment {self.segname} of type {self.segtype} with {len(self.residues)} residues>'

    def __str__(self):
        return f'{self.segname}: type {self.segtype} chain {self.chainID} with {len(self.residues)} residues'

    def cleave(self,clv,daughter_chainID):
        """
        Cleave the segment at a specified cut location, creating a daughter segment.
        
        Parameters
        ----------
        clv : Cleavage
            The cut location specifying where to cleave the segment.
        daughter_chainID : str
            The chain ID to assign to the daughter segment.
        """
        assert clv.chainID==self.segname
        assert self.segtype=='protein'
        r2=self.residues.get(resseqnum=clv.resseqnum2,insertion=clv.insertion2)
        r2i=self.residues.index(r2)
        # These two slice operations create *copies* of lists of residues
        # The original list of residues contains elements that are referenced
        # by other structures, like links.   
        assert isinstance(self.residues,ResidueList)
        daughter_residues=self.residues[r2i:]
        assert isinstance(daughter_residues,ResidueList)
        parent_residues=self.residues[:r2i]
        assert daughter_residues[0] is self.residues[r2i]
        self.residues=parent_residues
        # this set_chainID call must act DEEPLY on structures in the list items
        # that have chainID attributes AND we have to revisit all links!
        daughter_residues.set(chainID=daughter_chainID)
        Dseg=Segment.new(daughter_residues,segname=daughter_chainID)
        assert isinstance(Dseg,Segment)
        Dseg._obj_manager=self._obj_manager.expel(daughter_residues)
        # Dseg.ancestor_obj=self.ancestor_obj
        return Dseg
    
    def write_TcL(self,W:PsfgenScripter,transform):
        """
        Write the Tcl commands to create a segment in the Psfgen script.
        This method generates the Tcl commands to create a new segment in the PsfgenScripter script
        based on the segment type and its specifications.
        
        Parameters
        ----------
        W : Psfgen
            The Psfgen script writer object to which the Tcl commands will be written.
        transform : Transform
            The transformation object containing mapping information.
        """

        if self.segtype=='protein':
            self.polymer_stanza(W,transform)
        elif self.segtype=='glycan':
            self.glycan_stanza(W,transform)
        elif self.segtype=='nucleicacid':
            self.nucleicacid_stanza(W,transform)
        else:
            self.generic_stanza(W,transform)
        
    def glycan_stanza(self,W:PsfgenScripter,transform):
        """
        Write the Tcl commands to create a glycan segment in the Psfgen script.

        Parameters
        ----------
        W : PsfgenScripter
            The PsfgenScripter script writer object to which the Tcl commands will be written.
        transform : Transform
            The transformation object containing mapping information.
        """
        self.generic_stanza(W,transform)
    
    def generic_stanza(self,W,transform):
        """
        Write the Tcl commands to create a generic segment in the Psfgen script.
        
        Parameters
        ----------
        W : PsfgenScripter
            The Psfgen script writer object to which the Tcl commands will be written.
        transform : Transform
            The transformation object containing mapping information.
        """
        assert len(self.subsegments)==1,'No missing atoms allowed in generic stanza'
        parent_molecule=self._parent_molecule
        chainIDmap=transform.chainIDmap
        seglabel=self.segname
        image_seglabel=chainIDmap.get(seglabel,seglabel)
        is_image=image_seglabel!=seglabel

        objmanager=self._obj_manager
        seqmods=objmanager.get('seq',{})
        seg_grafts=seqmods.get('grafts',GraftList([]))

        transform.register_mapping(self.segtype,image_seglabel,seglabel)

        W.banner(f'Segment {image_seglabel} begins as image of {seglabel}')
        for g in seg_grafts:
            g.write_pre_segment(W)
        serial_list=self.residues.atom_serials(as_type=int)
        resid_list=self.residues.atom_resseqnums(as_type=int)
        vmd_red_list=reduce_intlist(serial_list)
        pdb=f'segtype_generic_{image_seglabel}.pdb'
        selname=image_seglabel
        W.addfile(pdb) # appends this file to the scripters FileCollector for later cleanup
        W.addline(f'set {selname} [atomselect ${parent_molecule.molid_varname} "serial {vmd_red_list}"]')
        W.addline(f'${selname} set segname {image_seglabel}')
        if is_image:
            W.backup_selection(selname,dataholder=f'{selname}_data')
            W.addline(f'${selname} set chain {image_seglabel}')                 
            W.addline(f'${selname} move {transform.write_TcL()}')
        W.addline(f'${selname} set resid [list {" ".join([str(x) for x in resid_list])}]')
        W.addline(f'${selname} writepdb {pdb}')
        W.addline(f'segment {image_seglabel} '+'{')
        W.addline(f'first none',indents=1)
        W.addline(f'last none',indents=1)
        W.addline(f'pdb {pdb}',indents=1)
        for g in seg_grafts:
            g.write_in_segment(W)
        W.addline('}')
        W.addline(f'coordpdb {pdb} {image_seglabel}')
        for g in seg_grafts:
            g.write_post_segment(W)
        if is_image:
            W.banner(f'Restoring A.U. state for {seglabel}')
            W.restore_selection(selname,dataholder=f'{selname}_data')
        W.banner(f'Segment {image_seglabel} ends')

    def nucleicacid_stanza(self,W:PsfgenScripter,transform):
        """
        Write the Tcl commands to create a nucleic acid segment in the Psfgen script.  For now, this just redirects to a polymer stanza.

        Parameters
        ----------
        W : PsfgenScripter
            The Psfgen script writer object to which the Tcl commands will be written.
        transform : Transform
            The transformation object containing mapping information.
        """
        self.polymer_stanza(W,transform)

    def polymer_stanza(self,W:PsfgenScripter,transform):
        """
        Write the Tcl commands to create a polymer (protein or nucleic acid) 
        segment in the Psfgen script.

        Parameters
        ----------
        W : PsfgenScripter
            The Psfgen script writer object to which the Tcl commands will be written.
        transform : Transform
            The transformation object containing mapping information.
        """
        parent_molecule=self._parent_molecule
        chainIDmap=transform.chainIDmap
        seglabel=self.segname
        image_seglabel=chainIDmap.get(seglabel,seglabel)
        is_image=image_seglabel!=seglabel
        segtype=self.segtype
        objmanager=self._obj_manager
        seqmods=objmanager.get('seq',{})
        logger.debug(f'polymer_stanza {segtype} {seglabel}->{image_seglabel} seqmods: {seqmods}')

        transform.register_mapping(segtype,image_seglabel,seglabel)

        loopspecs=self.specs.get('loops',{})
        sac_rn=loopspecs.get('sac_res_name','NOSAC')
        min_loop_length=loopspecs.get('min_loop_length',0)
        build_all_terminal_loops=self.specs.get('include_terminal_loops',False)
        build_N_terminal_loop=seglabel in self.specs.get('build_zero_occupancy_N_termini',[])
        build_C_terminal_loop=seglabel in self.specs.get('build_zero_occupancy_C_termini',[])

        seg_mutations=seqmods.get('mutations',MutationList([]))
        logger.debug(f'polymer_stanza for {segtype} segname {seglabel}; init mutations:')
        for m in seg_mutations:
            logger.debug(str(m))
        seg_Cfusions=seqmods.get('Cfusions',CfusionList([]))
        seg_patches=seqmods.get('patches',PatchList([]))
        W.banner(f'Segment {image_seglabel} begins')
        logger.debug(f'Segment {image_seglabel} begins')
        for sf in seg_Cfusions:
            sf.write_pre_segment(W)
        for i,b in enumerate(self.subsegments):
            if b.state=='RESOLVED':
                """ for a resolved subsegment, generate its pdb file """
                b.selname=f'{image_seglabel}{i:02d}'
                run=ResidueList(self.residues[b.bounds[0]:b.bounds[1]+1])
                b.pdb=f'segtype_polymer_{image_seglabel}_{run[0].resseqnum}{run[0].insertion}_to_{run[-1].resseqnum}{run[-1].insertion}.pdb'
                W.addfile(b.pdb)
                serial_list=run.atom_serials(as_type=int)
                logger.debug(f'Last atom has serial {serial_list[-1]}')
                at=parent_molecule.asymmetric_unit.atoms.get(serial=serial_list[-1])
                if hasattr(at,'__len__'):
                    for a in at:
                        logger.debug(f'whoops: {a.chainID} {a.resseqnum} {a.name} is {a.serial}')
                    raise Exception(f'More than one atom with serial {serial_list[-1]}??')
                logger.debug(f'here it is {at.serial} {at.resname} {at.name}')
                assert at.resseqnum==run[-1].resseqnum
                vmd_red_list=reduce_intlist(serial_list)
                W.addline(f'set {b.selname} [atomselect ${parent_molecule.molid_varname} "serial {vmd_red_list}"]')
                W.addline(f'${b.selname} set segname {image_seglabel}')
                if hasattr(at,'_ORIGINAL_'):
                    if at._ORIGINAL_["serial"]!=at.serial:
                        W.banner(f'Atom with serial {at._ORIGINAL_["serial"]} in PDB needs serial {at.serial} for VMD')
                """ Relabel chain ID and request coordinate transformation """
                if is_image:
                    W.backup_selection(b.selname,dataholder=f'{b.selname}_data')
                    W.addline(f'${b.selname} set chain {image_seglabel}')                 
                    W.addline(f'${b.selname} move {transform.write_TcL()}')
                W.addline(f'${b.selname} writepdb {b.pdb}')
            elif b.state=='MISSING':
                if i==0:
                    if build_N_terminal_loop or build_all_terminal_loops:
                        b.declare_buildable()
                elif i==(len(self.subsegments)-1):
                    if build_C_terminal_loop or build_all_terminal_loops:
                        b.declare_buildable()
                else:
                    b.declare_buildable()
        W.addline(f'segment {image_seglabel} '+'{')
        if self.subsegments[0].state=='MISSING' and not self.subsegments[0].build:
            Nterminal_missing_subsegment=self.subsegments.pop(0)
            # logger.info(f'Since terminal loops are not included, ignoring {str(Nterminal_missing_subsegment)}')
        if self.subsegments[-1].state=='MISSING' and not self.subsegments[-1].build:
            Cterminal_missing_subsegment=self.subsegments.pop(-1)
            # logger.info(f'Since terminal loops are not included, ignoring {str(Cterminal_missing_subsegment)}')
        for b in self.subsegments:
            if b.state=='RESOLVED':
                W.addline(f'pdb {b.pdb}',indents=1)
            elif b.state=='MISSING' and b.build:
                for r in self.residues[b.bounds[0]:b.bounds[1]+1]:
                    rname=Labels.charmm_resname_of_pdb_resname.get(r.resname,r.resname)
                    W.addline(f'residue {r.resseqnum}{r.insertion} {rname} {image_seglabel}',indents=1)
                if b.num_items()>=min_loop_length and not b in [self.subsegments[0],self.subsegments[-1]]:
                    lrr=self.residues[b.bounds[1]]
                    sac_resseqnum=lrr.resseqnum
                    sac_insertion='A' if lrr.insertion in [' ',''] else chr(ord(lrr.insertion)+1)
                    assert sac_insertion<='Z',f'Residue {lrr.resseqnum} of chain {seglabel} already has too many insertion instances (last: {lrr.insertion}) to permit insertion of a sacrificial {sac_rn}'
                    b.sacres=Residue({'resname':sac_rn,'resseqnum':sac_resseqnum,'insertion':sac_insertion,'chainID':seglabel,'segtype':segtype,'resolved':False,'atoms':[]})
                    W.addline(f'residue {sac_resseqnum}{sac_insertion} {sac_rn} {image_seglabel}',indents=1)
        for cf in seg_Cfusions:
            cf.write_in_segment(W)
        for m in seg_mutations:
                W.comment(f'mutation source: {m.typekey}')
                W.addline(m.write_TcL())
        for p in seg_patches:
            if p.use_in_segment=='first':
                W.addline(f'first {p.patchname}',indents=1)
            elif p.use_in_segment=='last':
                W.addline(f'last {p.patchname}',indents=1)
        W.addline('}')
        W.banner(f'End segment {image_seglabel}')
        W.banner('Coordinate-specification commands')
        for b in self.subsegments:
            if b.state=='RESOLVED':
                W.comment(f'Subsegment {[self.subsegments.index(b)]} is a resolved run')
                W.addline(f'coordpdb {b.pdb} {image_seglabel}')
        for b in self.subsegments:
            if b.state=='MISSING' and b.build:
                if self.subsegments.index(b)>0: # only seed orientation for a loop that is not at the N-terminus
                    W.comment(f'Subsegment {[self.subsegments.index(b)]}/{len(self.subsegments)} is a missing loop')
                    this_run=ResidueList(self.residues[b.bounds[0]:b.bounds[1]+1])
                    prior_b=self.subsegments[self.subsegments.index(b)-1]
                    W.comment(f'...attached to subsegment {self.subsegments.index(prior_b)}')
                    prior_run=ResidueList(self.residues[prior_b.bounds[0]:prior_b.bounds[1]+1])
                    if segtype=='protein':
                        W.comment(f'Seeding orientation of model-built loop starting at {str(this_run[0])} from {str(prior_run[-1])}')
                        W.addline(f'{this_run.caco_str(prior_run,image_seglabel,parent_molecule.molid_varname,transform.tmat)}')
        if len(seg_patches)>0:
            W.banner(f'Patches for segment {image_seglabel}')
        for patch in seg_patches:
            if patch.use_in_segment=='' and not patch.use_after_regenerate:
                W.addline(f'patch {patch.patchname} {image_seglabel}:{patch.resseqnum}{patch.insertion}')
            if patch.use_after_regenerate:
                W.addpostregenerateline(f'patch {patch.patchname} {image_seglabel}:{patch.resseqnum}{patch.insertion}')
        for sc in seg_Cfusions:
            sc.write_post_segment(W)
        W.banner('Intra-segmental terminal patches')
        if segtype=='protein':
            for i,b in enumerate(self.subsegments):
                # only non-terminal loops get the terminal patches
                if b.state=='MISSING' and 0<i<(len(self.subsegments)-1) and hasattr(b,'sacres'):
                    Cterm=self.residues[b.bounds[1]]
                    W.addline(f'patch CTER {image_seglabel}:{Cterm.resseqnum}{Cterm.insertion}')
                    nextb=self.subsegments[i+1]
                    Nterm=self.residues[nextb.bounds[0]]
                    patchname='NTER'
                    if Nterm.resname=='PRO':
                        patchname='PROP'
                    elif Nterm.resname=='GLY':
                        patchname='GLYP'
                    W.addline(f'patch {patchname} {image_seglabel}:{Nterm.resseqnum}{Nterm.insertion}')
                    W.addline(f'delatom {image_seglabel} {b.sacres.resseqnum}{b.sacres.insertion}')
        W.banner('Restoring A.U. state for all resolved subsegments')
        for b in self.subsegments:
            if b.state=='RESOLVED':
                if is_image:
                    W.restore_selection(b.selname,dataholder=f'{b.selname}_data')
        W.banner(f'Segment {image_seglabel} ends')
        # THIS IS BAD self.objmanager.retire('seq')

class SegmentList(BaseObjList[Segment]):
    """
    A list of segments in a molecular structure.
    This class represents a collection of segments and provides methods to manage and manipulate them.
    It inherits from `AncestorAwareObjList` to maintain the context of the parent molecule.
    """

    def describe(self) -> str:
        return f'<SegmentList with {len(self)} segments>'

    @classmethod
    def generate(cls, seq_spec: dict = {}, residues: ResidueList = [], chainIDmanager: ChainIDManager = None):
        """
        Generate a SegmentList from a sequence specification and a list of residues.
        
        Parameters
        ----------
        seq_spec : dict
            A dictionary containing sequence specifications.
        residues : ResidueList
            A list of residues to be processed into segments.
        chainIDmanager : ChainIDManager
            An object managing chain IDs to ensure uniqueness.
        
        Returns
        -------
        SegmentList
            A new SegmentList instance containing the generated segments.
        """
        self = cls([])
        self.counters_by_segtype={}
        self.segnames=[]
        self.daughters={}
        self.segtype_of_segname={}
        assert all([x.segtype!='UNSET' for x in residues]),f'There are residues with UNSET segtype: {[(x.resname,x.chainID,x.resseqnum) for x in residues if x.segtype=="UNSET"]}'
        self.segtypes_ordered=[]
        for r in residues:
            if not r.chainID in self.segtype_of_segname:
                self.segtype_of_segname[r.chainID]=r.segtype
            if not r.segtype in self.segtypes_ordered:
                self.segtypes_ordered.append(r.segtype)
        logger.debug(f'Generating segments from list of {len(residues)} residues. segtypes detected: {self.segtypes_ordered}')
        initial_chainIDs=list(self.segtype_of_segname.keys())
        logger.debug(f'ChainIDs detected: {initial_chainIDs}')
        chainIDmanager.sandbag(initial_chainIDs)
        for stype in self.segtypes_ordered:
            self.counters_by_segtype[stype]=0
            res=residues.filter(segtype=stype)
            orig_chainIDs=res.uniqattrs(['chainID'])['chainID']
            logger.debug(f'Processing {len(res)} residues of segtype {stype} in {len(orig_chainIDs)} unique chainIDs')
            orig_res_groups={chainID: res.filter(chainID=chainID) for chainID in orig_chainIDs}
            logger.debug(f'Found {len(orig_res_groups)} original chainIDs for segtype {stype}: {list(orig_res_groups.keys())} {list([len(x) for x in orig_res_groups.values()])}')
            for chainID,c_res in orig_res_groups.items():
                assert isinstance(c_res,ResidueList),f'ChainID {chainID} residues are not a ResidueList: {type(c_res)}'
                this_chainID=chainID
                # c_res=res.filter(chainID=this_chainID)
                logger.debug(f'-> original chainID {chainID} in {len(c_res)} residues')
                this_chainID=chainIDmanager.check(this_chainID)
                if this_chainID != chainID:
                    logger.debug(f'{len(c_res)} residues with original chainID {chainID} and segtype {stype} are assigned new daughter chainID {this_chainID}')
                    if not chainID in self.daughters:
                        self.daughters[chainID]=[]
                    self.daughters[chainID].append(this_chainID)
                    c_res.set_chainIDs(this_chainID)
                    self.segtype_of_segname[this_chainID]=stype
                    if stype=='protein':
                        if chainID in seq_spec['build_zero_occupancy_C_termini']:
                            logger.debug(f'-> appending new chainID {this_chainID} to build_zero_occupancy_C_termini due to member {chainID}')
                            seq_spec['build_zero_occupancy_C_termini'].insert(seq_spec['build_zero_occupancy_C_termini'].index(chainID),this_chainID)
                        if chainID in seq_spec['build_zero_occupancy_N_termini']:
                            logger.debug(f'-> appending new chainID {this_chainID} to build_zero_occupancy_N_termini due to member {chainID}')
                            seq_spec['build_zero_occupancy_N_termini'].insert(seq_spec['build_zero_occupancy_N_termini'].index(chainID),this_chainID)
                        # seq_spec['build_zero_occupancy_C_termini'].remove(chainID)
                num_mis=sum([1 for x in c_res if len(x.atoms)==0])
                thisSeg=Segment.new(c_res, segname=this_chainID, specs=seq_spec)
                logger.debug(f'Made segment: stype {stype} chainID {this_chainID} segname {thisSeg.segname} ({num_mis} missing) (seq_spec {seq_spec})')
                self.append(thisSeg)
                self.segnames.append(thisSeg.segname)
                self.counters_by_segtype[stype]+=1
        return self
    
    def collect_residues(self):
        """
        Collect all residues from the segments in the list.

        Returns
        -------
        ResidueList
            A list of all residues associated with the segments."""
        residues = ResidueList([])
        for seg in self:
            residues.extend(seg.residues)
        return residues

    def collect_working_files(self):
        """
        Collect all working files (PDB files) from the segments in the list.

        Returns
        -------
        list
            A list of PDB files associated with the segments."""
        working_files=[]
        for seg in self:
            for subseg in seg.subsegments:
                if hasattr(subseg,'pdb'):
                    logger.debug(f'wf: {subseg.pdb}')
                    working_files.append(subseg.pdb)
        return working_files
    
    def write_TcL(self,W:PsfgenScripter,transform):
        """
        Write the Tcl commands to create all segments in the Psfgen script.
        
        Parameters
        ----------
        W : PsfgenScripter
            The Psfgen script writer object to which the Tcl commands will be written.
        transform : Transform
            The transformation object containing mapping information."""
        for seg in self:
            W.comment(f'Writing seg {seg.segname}')
            seg.write_TcL(W,transform)
            W.comment(f'Done with seg {seg.segname}')
    
    def get_segment_of_residue(self,residue):
        """
        Get the segment that contains a specified residue.

        Parameters
        ----------
        residue : Residue
            The residue for which to find the containing segment.

        Returns
        -------
        Segment or None
            The segment containing the specified residue, or None if not found.
        """
        for S in self:
            if residue in S.residues:
                return S
        return None

    def remove(self,item):
        """
        Remove a segment from the list and update associated attributes.
        
        Parameters
        ----------
        item : Segment
            The segment to remove from the list.
        """
        self.segnames.remove(item.segname)
        self.counters_by_segtype[self.segtype_of_segname[item.segname]]-=1
        return super().remove(item)

    def inherit_objs(self, objmanager:ObjManager):
        """
        Inherit objects from the object manager for each segment in the list.
        This method updates the object manager for each segment with inherited objects
        such as mutations, grafts, and patches.
        
        Parameters
        ----------
        objmanager : ObjectManager
            The object manager from which to inherit objects.
        """
        for item in self:
            item._obj_manager=objmanager.filter_copy(objnames=item._inheritable_objs,chainID=item.segname)
            counts = item._obj_manager.counts_by_header()
            logger.debug(f'Inherit objs: seg {item.segname} has {counts}')
