# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A graft refers to a set of residues that are added to the base molecule and sourced from a separate pdb. 
The graft's "target" is a residue that is congruent to the "reference" residue in the graft.  The graft is positioned by aligning all of
its atoms using a triple on the reference and the same triplet on the target as an alignment generator. The target residue's atoms
are replaced by those in the reference and the rest of the reference is incorporated.
"""
import logging
logger=logging.getLogger(__name__)
import os

from functools import singledispatchmethod

from ..core.baseobj import AncestorAwareObj, AncestorAwareObjList
from ..core.scripters import PsfgenScripter
from ..core.stringthings import ri_range
from .link import LinkList

class Graft(AncestorAwareObj):
    """
    A class for handling grafts.
    """

    req_attr=AncestorAwareObj.req_attr+['id','orig_chainID','orig_resseqnum1','orig_insertion1','orig_resseqnum2','orig_insertion2','source_pdbid','source_chainID','source_resseqnum1','source_insertion1','source_resseqnum2','source_insertion2','source_resseqnum3','source_insertion3']
    """
    Required attributes for a Graft object.
    These attributes must be provided when creating a Graft object.
    
    - ``id``: Unique identifier for the graft.
    - ``orig_chainID``: Chain ID of the target segment in the base molecule.
    - ``orig_resseqnum1``: N-terminal residue sequence number of the target segment.
    - ``orig_insertion1``: Insertion code of the N-terminal residue in the target segment.
    - ``orig_resseqnum2``: C-terminal residue sequence number of the target segment.
    - ``orig_insertion2``: Insertion code of the C-terminal residue in the target segment.
    - ``source_pdbid``: Basename of the source PDB file or PDB ID from which the graft is sourced.
    - ``source_chainID``: Chain ID in the source PDB file.
    - ``source_resseqnum1``: N-terminal residue sequence number of the source segment.
    - ``source_insertion1``: Insertion code of the N-terminal residue in the source segment.
    - ``source_resseqnum2``: C-terminal residue sequence number of the source segment.
    - ``source_insertion2``: Insertion code of the C-terminal residue in the source segment.
    - ``source_resseqnum3``: Optional third residue sequence number in the source segment.
    - ``source_insertion3``: Optional third insertion code in the source segment.
    """
    
    yaml_header='grafts'
    """
    YAML header for Graft objects.
    This header is used to identify Graft objects in YAML files.
    """
    
    objcat='seq'
    """
    Category of the Graft object.
    This categorization is used to group Graft objects in the object manager.
    """
    
    _Graft_counter=0
    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)

    @__init__.register(str)
    def _from_shortcode(self,shortcode):
        # shortcode format: "target:source"
        # target format:    "C_RRR[-SSS]"
        #  - C chainID 
        #  - RRR resid+insertion of first residue in target alignment basis
        #  - SSS optional resid+insertion of last residue in target alignment basis
        #    (if not present, only RRR is used)
        # source format:     "pdbid,C_RRR[#SSS][-TTT]"
        #  - pdbid basename of pdb file or pdb id
        #  - C chainID in source pdb
        #  - RRR resid+insertion of first residue in source alignment basis
        #  - SSS optional resid+insertion of last residue in target alignment basis
        #    (if not present, only RRR is used)
        #  - TTT optional resid+insertion completing RRR-TTT range for entire graft source
        #    (if not present, only RRR is the entire graft source)
        tokens=shortcode.split(':')
        assert len(tokens)==2,f'Malformed graft shortcode {shortcode}'
        target,source=tokens
        target_chainID,target_resrange=target.split('_')
        trr=ri_range(target_resrange)
        resseqnum1,insertion1=trr[0]
        if len(trr)>1:
            resseqnum2,insertion2=trr[1]
        else:
            resseqnum2,insertion2=resseqnum1,insertion1
        source=tokens[1].split(',')
        assert len(source)==2,f'Malformed graft spec {shortcode}'
        source_pdbid,source_chainresrange_spec=source
        source_chainID,source_resrange=source_chainresrange_spec.split('_')
        srr=ri_range(source_resrange)
        source_resseqnum1,source_insertion1=srr[0]
        if len(srr)>1:
            if len(srr)>2:
                source_resseqnum2,source_insertion2=srr[1]
                source_resseqnum3,source_insertion3=srr[2]
            else:
                if '#' in source_resrange:
                    source_resseqnum2,source_insertion2=srr[1]
                    source_resseqnum3,source_insertion3=source_resseqnum2,source_insertion2
                else:
                    source_resseqnum3,source_insertion3=srr[1]
                    source_resseqnum2,source_insertion2=source_resseqnum1,source_insertion1
        else:
            source_resseqnum2,source_insertion2=source_resseqnum1,source_insertion1
            source_resseqnum3,source_insertion3=source_resseqnum1,source_insertion1
        input_dict={
            'orig_chainID':target_chainID,
            'orig_resseqnum1':resseqnum1,
            'orig_insertion1':insertion1,
            'orig_resseqnum2':resseqnum2,
            'orig_insertion2':insertion2,
            'source_pdbid':source_pdbid,
            'source_chainID':source_chainID,
            'source_resseqnum1':source_resseqnum1,
            'source_insertion1':source_insertion1,
            'source_resseqnum2':source_resseqnum2,
            'source_insertion2':source_insertion2,
            'source_resseqnum3':source_resseqnum3,
            'source_insertion3':source_insertion3,
            'id':Graft._Graft_counter
        }
        logger.debug(f'Initializing graft {input_dict}')
        Graft._Graft_counter+=1
        self.residues=None
        super().__init__(input_dict)

    def activate(self,mol):
        """
        Activate the graft by linking it to a source molecule.

        Parameters
        ----------
        mol: Molecule
            The source molecule from which the graft will be sourced.
        """
        self.source_molecule=mol
        g_topomods=self.source_molecule.objmanager.get('topol',{})
        g_links=g_topomods.get('links',LinkList([]))
        self.source_seg=self.source_molecule.asymmetric_unit.segments.get(segname=self.source_chainID)
        self.mover_residues=type(self.source_seg.residues)([])
        self.index_residues=type(self.source_seg.residues)([])
        for residue in self.source_seg.residues:
            if residue>=f'{self.source_resseqnum1}{self.source_insertion1}' and residue<=f'{self.source_resseqnum2}{self.source_insertion2}':
                self.index_residues.append(residue)
            elif residue>f'{self.source_resseqnum2}{self.source_insertion2}' and residue<=f'{self.source_resseqnum3}{self.source_insertion3}':
                self.mover_residues.append(residue)
        # only ingest links that are internal to this set of residues or that link to index residues
        self.my_links=type(g_links)([])
        for l in g_links:
            if l.residue1 in self.mover_residues and l.residue2 in self.mover_residues:
                self.my_links.append(l)
            elif l.residue1 in self.index_residues and l.residue2 in self.mover_residues:
                self.my_links.append(l)
            elif l.residue2 in self.index_residues and l.residue1 in self.mover_residues:
                self.my_links.append(l)
        logger.debug(f'Activated graft {self.id}')

    def set_links(self,next_resseqnum):
        """
        Set the links for the graft based on the resolved receiver residues.

        Parameters
        ----------
        next_resseqnum: int
            The next available residue sequence number to be assigned to the graft residues.

        Returns
        -------
        next_resseqnum: int
            The updated next available residue sequence number after setting the links.
        """
        assert self.residues!=None
        logger.debug(f'Setting links for graft {self.id}')
        logger.debug(f'-> resolved receiver residues {[str(x) for x in self.residues]}')
        self.index_residues.set(chainID=self.residues[0].chainID)
        self.mover_residues.set(chainID=self.residues[0].chainID)
        for r in self.mover_residues:
            r.set(resseqnum=next_resseqnum)
            next_resseqnum+=1
        for l in self.my_links:
            if l.residue1 in self.index_residues and l.residue2 in self.mover_residues:
                l.residue1=self.residues[self.index_residues.index(l.residue1)]
            elif l.residue2 in self.index_residues and l.residue1 in self.mover_residues:
                l.residue2=self.residues[self.index_residues.index(l.residue2)]
        return next_resseqnum
    
    def __str__(self):
        res=f'Graft {self.id}: index residues {[str(x) for x in self.index_residues]} delivers {len(self.mover_residues)} residues to \n'
        res+=f'  target {[str(x) for x in self.residues]} along with {len(self.my_links)} internal links:\n'
        for l in self.my_links:
            res+=f'  -> link {str(l)}\n'
        return res

    def write_pre_segment(self,W:PsfgenScripter):
        """
        Write the pre-segment Tcl commands for the graft operation.
        This method generates the Tcl commands to prepare the graft operation in VMD.
        It sets up the molecule information, selects the appropriate atoms for the graft,
        and applies a transformation to align the graft with the target segment in the base molecule.
        It also writes the grafted segment to a PDB file.

        Parameters
        ----------
        W : PsfgenScripter
            The Psfgen script writer object to which the Tcl commands will be written.
        """
        W.comment(f'{str(self)}')
        W.addline(f'set topid [molinfo ${W.molid_varname} get id]')
        if os.path.exists(f'{self.source_pdbid}.pdb'):
            W.addline(f'mol new {self.source_pdbid}.pdb')
        else:
            W.addline(f'mol new {self.source_pdbid}')
        W.addline(f'set graftid [molinfo top get id]')
        W.addline(f'mol top $topid')
        W.addline(f'set mover_sel [atomselect $graftid "chain {self.source_chainID} and (resid {self.source_resseqnum1}{self.source_insertion1} to {self.source_resseqnum3}{self.source_insertion3}) and (not resid {self.source_resseqnum1}{self.source_insertion1} to {self.source_resseqnum2}{self.source_insertion2})"]')
        W.addline(f'set to_sel [atomselect ${W.molid_varname} "chain {self.orig_chainID} and resid {self.orig_resseqnum1}{self.orig_insertion1} to {self.orig_resseqnum2}{self.orig_insertion2}"]')
        W.addline(f'set from_sel [atomselect $graftid "chain {self.source_chainID} and resid {self.source_resseqnum1}{self.source_insertion1} to {self.source_resseqnum2}{self.source_insertion2}"]')
        W.addline(f'vmdcon -info "[$mover_sel num] atoms will move"')
        W.addline(f'vmdcon -info "    by comparing [$from_sel num] atoms to [$to_sel num] atoms"')
        W.addline(f'set TT [transidentity]')
        W.addline(f'set TT [measure fit $from_sel $to_sel]')
        W.addline(f'vmdcon -info "Homog. trans. matrix: $TT"')
        W.addline(f'$mover_sel move $TT')
        self.segfile=f'graft{self.id}.pdb'
        new_residlist=[]
        for y in self.mover_residues:
            new_residlist.extend([f'{y.resseqnum}' for x in y.atoms])
        W.addline(f'$mover_sel set resid [list {" ".join(new_residlist)}]')
        W.addline(f'$mover_sel set chain {self.mover_residues[0].chainID}')
        W.addline(f'$mover_sel writepdb {self.segfile}')
        W.addline(f'$mover_sel delete')

    def write_in_segment(self,W:PsfgenScripter):
        """
        Write the Tcl commands to add the graft into the active segment of the base molecule.
        This method generates the Tcl commands to add the graft segment to the active segment
        of the base molecule in the Psfgen script.
        
        Parameters
        ----------
        W : PsfgenScripter
            The Psfgen script writer object to which the Tcl commands will be written."""
        W.addline (f'    pdb {self.segfile}')

    def write_post_segment(self,W:PsfgenScripter):
        """
        Write the Tcl commands to finalize the graft segment in the Psfgen script.
        This method generates the Tcl commands to finalize the graft segment in the Psfgen script.

        Parameters
        ----------
        W : PsfgenScripter
            The Psfgen script writer object to which the Tcl commands will be written.
        """
        W.addline(f'coordpdb {self.segfile} {self.chainID}')

    def assign_residues(self,Residues):
        """
        Assign the residues for the graft based on the original chain ID and residue sequence numbers.
        This method scans the residues in the original chain and assigns those that fall within the specified range
        of residue sequence numbers and insertions to the graft's residues attribute.

        Parameters
        ----------
        Residues: ResidueList
            The list of residues from which the graft will be assigned.
        """
        assert self.residues==None
        logger.debug(f'Assigning receiver residues to graft {self.id} ({self.orig_resseqnum1}{self.orig_insertion1} to {self.orig_resseqnum2}{self.orig_insertion2})...')
        this_chain=Residues.filter(chainID=self.orig_chainID)
        logger.debug(f'...scanning {len(this_chain)} residues for ones to include in graft {self.id} [{self.orig_resseqnum1}{self.orig_insertion1}-{self.orig_resseqnum2}{self.orig_insertion2}]')
        self.residues=type(Residues)([])
        for r in this_chain:
            if r>=f'{self.orig_resseqnum1}{self.orig_insertion1}' and r<=f'{self.orig_resseqnum2}{self.orig_insertion2}':
                logger.debug(f'{r.resseqnum}{r.insertion} belongs')
                self.residues.append(r)
        logger.debug(f'...assigned {len(self.residues)} residues')

class GraftList(AncestorAwareObjList):
    """
    A class for handling lists of grafts.
    This class inherits from AncestorAwareObjList and provides methods to manage a list of Graft objects.
    It allows for the assignment of residue objects to each graft and the removal of grafts that do not have any assigned residues.
    It also handles the destruction of down-links from terminal residues in the grafts.
    """
    def assign_residues(self,Residues,Links):
        """
        Assign residue objects to each graft in the list.
        This method iterates through each graft in the list, assigns residues from the provided Residues
        object, and handles the removal of down-links from terminal residues.

        Parameters
        ----------
        Residues: ResidueList
            The list of residues from which the grafts will be assigned.
        Links: LinkList
            The list of links that may need to be modified based on the assigned residues.
        """
        logger.debug(f'Assigning residue objects to {len(self)} grafts')
        delete_us=[]
        down_group=[]
        down_links=[]
        for s in self:
            logger.debug(f' -> graft {s.id}')
            s.assign_residues(Residues)
            if s.residues:
                logger.debug(f'   -> {len(s.residues)} assigned')
                term_res=s.residues[-1]
                if len(term_res.down)>0:
                    logger.debug(f'-> terminal residue {term_res.chainID}_{term_res.resname}{term_res.resseqnum}{term_res.insertion}')
                    logger.debug(f'   links down to {", ".join([str(x) for x in term_res.down])}')
                    down_group,down_links=term_res.get_down_group()
                    logger.debug(f'Target residue down-links to a group of {len(down_group)} residues:')
                    if len(down_group)>0:
                        for dg in down_group:
                            logger.debug(f'  removing residue {str(dg)}')
                            Residues.remove(dg)
                    if len(down_links)>0:
                        for dl in down_links:
                            logger.debug(f'  removing link {str(dl)}')
                            Links.remove(dl)
                    if len(down_group)>1:
                        term_res.unlink(down_group[0])
                else:
                    logger.debug(f'-> no down-links to destroy.')

        delete_us=self.__class__([s for s in self if s.residues==None])
        for s in delete_us:
            self.remove(s)
        return delete_us,down_group,down_links

