# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import logging
logger=logging.getLogger(__name__)
import os

from functools import singledispatchmethod

from ..baseobj import AncestorAwareObj, AncestorAwareObjList
from ..scriptwriters import Psfgen
from ..stringthings import ri_range
from .link import LinkList

class Graft(AncestorAwareObj):
    """A class for handling grafts.  A graft refers to a set of residues that are added to the base molecule and sourced from a separate pdb. 
    The graft's "target" is a residue that is congruent to the "reference" residue in the graft.  The graft is positioned by aligning all of
    its atoms using a triple on the reference and the same triplet on the target as an alignment generator. The target residue's atoms
    are replaced by those in the reference and the rest of the reference is incorporated."""
    req_attr=AncestorAwareObj.req_attr+['id','orig_chainID','orig_resseqnum1','orig_insertion1','orig_resseqnum2','orig_insertion2','source_pdbid','source_chainID','source_resseqnum1','source_insertion1','source_resseqnum2','source_insertion2','source_resseqnum3','source_insertion3']
    yaml_header='grafts'
    objcat='seq'
    
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
        # only injest links that are internal to this set of residues or that link to index residues
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

    def write_pre_segment(self,W:Psfgen):
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
    def write_in_segment(self,W:Psfgen):
        W.addline (f'    pdb {self.segfile}')
    def write_post_segment(self,W:Psfgen):
        W.addline(f'coordpdb {self.segfile} {self.chainID}')
    # def execute_graft(self,target_residue,next_available_resid):
    #     graft_chainID=self.chainID 
    #     g_topomods=self.source_molecule.modmanager.get('topomods',{})
    #     g_links=g_topomods.get('links',LinkList([]))
    #     source_seg=self.source_molecule.asymmetric_unit.segments.get(segname=self.source_chainID)
    #     # build the list of new residues this graft contributes; the index residue is not included!
    #     self.donated_residues=[]
    #     for residue in source_seg.residues:
    #         if residue>f'{self.source_resseqnum1}{self.source_insertion1}' and residue<=f'{self.source_resseqnum2}{self.source_insertion2}':
    #             residue.set_resseqnum(next_available_resid)
    #             residue.set_chainID(graft_chainID)
    #             next_available_resid+=1
    #             self.mover_residues.append(residue)
    #         # only injest links that are internal to this set of residues
    #     self.donated_links=[]
    #     for l in g_links:
    #         if l.residue1 in self.donated_residues and l.residue2 in self.donated_residues:
    #             links.append(l)
    #                 injested_links+=1

    def assign_residues(self,Residues):
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
    def assign_residues(self,Residues,Links):
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

