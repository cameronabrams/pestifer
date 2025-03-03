# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import logging
logger=logging.getLogger(__name__)
from functools import singledispatchmethod
from pidibble.pdbrecord import PDBRecord

from ..baseobj import AncestorAwareObj, AncestorAwareObjList
from ..util.cifutil import CIFdict
from ..util.coord import ic_reference_closest
from ..scriptwriters import Psfgen, Filewriter
from ..stringthings import split_ri

class Link(AncestorAwareObj):
    """A class for handling covalent bonds between residues where at least one residue is non-protein
    
    Attributes
    ----------

    req_attr: list
    * name1: name of partner-1 residue
    * chainID1: chain ID of partner-1 residue
    * resseqnum1: residue sequence number of partner-1 residue
    * insertion1: insertion code of partner-1 residue
    * name2: name of partner-2 residue
    * chainID2: chain ID of partner-2 residue
    * resseqnum2: residue sequence number of partner-2 residue
    * insertion2: insertion code of partner-2 residue
    
    """
    req_attr=AncestorAwareObj.req_attr+['chainID1','resseqnum1','insertion1','chainID2','resseqnum2','insertion2']
    opt_attr=AncestorAwareObj.opt_attr+['name1','name2','altloc1','altloc2','resname1','resname2','sym1','sym2','link_distance','segname1','segname2','residue1','residue2','atom1','atom2','empty','segtype1','segtype2','ptnr1_label_asym_id','ptnr2_label_asym_id','ptnr1_label_seq_id','ptnr2_label_seq_id','ptnr1_label_comp_id','ptnr2_label_comp_id','ptnr1_auth_asym_id','ptnr2_auth_asym_id','ptnr1_auth_seq_id','ptnr2_auth_seq_id','ptnr1_auth_comp_id','ptnr2_auth_comp_id']    
    yaml_header='links'
    PDB_keyword='LINK'
    objcat='topol'
    
    patch_atomnames={
        'NGLA':['ND2','C1'],
        'NGLB':['ND2','C1'],
        'SGPA':['OG','C1'],
        'SGPB':['OG','C1'],
        '11aa':['O1','C1'],
        '11ab':['O1','C1'],
        '11bb':['O1','C1'],
        '12aa':['O2','C1'],
        '12ab':['O2','C1'],
        '12ba':['O2','C1'],
        '12bb':['O2','C1'],
        '13aa':['O3','C1'],
        '13ab':['O3','C1'],
        '13ba':['O3','C1'],
        '13bb':['O3','C1'],
        '14aa':['O4','C1'],
        '14ab':['O4','C1'],
        '14ba':['O4','C1'],
        '14bb':['O4','C1'],
        '16AT':['O6','C1'],
        '16BT':['O6','C1'],
        'SA26AT':['O6','C2'],
        'SA28AA':['O8','C2'],
        'SA29AT':['O9','C2'],
        'ZNHD':['NE2','ZN']
    }
    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)
    
    @__init__.register(PDBRecord)
    def _from_pdbrecord(self,pdbrecord):
        input_dict={
            'name1': pdbrecord.name1,
            'resname1':pdbrecord.residue1.resName,
            'chainID1':pdbrecord.residue1.chainID,
            'resseqnum1':pdbrecord.residue1.seqNum,
            'insertion1':pdbrecord.residue1.iCode,
            'name2': pdbrecord.name2,
            'resname2':pdbrecord.residue2.resName,
            'chainID2':pdbrecord.residue2.chainID,
            'resseqnum2':pdbrecord.residue2.seqNum,
            'insertion2':pdbrecord.residue2.iCode,
            'altloc2':pdbrecord.altLoc2,
            'altloc1':pdbrecord.altLoc1,
            'sym1':pdbrecord.sym1,
            'sym2':pdbrecord.sym2,
            'link_distance':pdbrecord.length,
            'segname1':pdbrecord.residue1.chainID,
            'segname2':pdbrecord.residue2.chainID,
            }
        input_dict.update({
                'residue1':None,
                'residue2':None,
                'atom1':None,
                'atom2':None,
                'empty':False,
                'segtype1':'UNSET',
                'segtype2':'UNSET',
            })
        super().__init__(input_dict)
        logger.debug(f'parsed {self}')
    
    @__init__.register(CIFdict)
    def _from_cifdict(self,cd):
        input_dict={}
        input_dict['name1']=cd['ptnr1_label_atom_id']
        input_dict['altloc1']=cd['pdbx_ptnr1_label_alt_id']
        input_dict['resname1']=cd['ptnr1_label_comp_id']
        # if the seq id is a dot, then we will revert to using the author's labels for this atom
        asym=cd['ptnr1_label_asym_id']
        seq=cd['ptnr1_label_seq_id']
        aseq=cd['ptnr1_auth_seq_id']
        if seq=='.':
            # assert len(asym)==1
            seq=int(aseq)
        else:
            seq=int(seq)
        input_dict['chainID1']=asym
        input_dict['resseqnum1']=seq
        input_dict['insertion1']=cd['pdbx_ptnr1_pdb_ins_code']

        input_dict['name2']=cd['ptnr2_label_atom_id']
        input_dict['altloc2']=cd['pdbx_ptnr2_label_alt_id']
        input_dict['resname2']=cd['ptnr2_label_comp_id']

        asym=cd['ptnr2_label_asym_id']
        seq=cd['ptnr2_label_seq_id']
        aseq=cd['ptnr2_auth_seq_id']
        if seq=='.':
            # assert len(asym)==1
            seq=int(aseq)
        else:
            seq=int(seq)
        input_dict['chainID2']=asym
        input_dict['resseqnum2']=seq

        input_dict['insertion2']=cd['pdbx_ptnr2_pdb_ins_code']
        input_dict['sym1']=cd['ptnr1_symmetry']
        input_dict['sym2']=cd['ptnr2_symmetry']
        input_dict['link_distance']=float(cd['pdbx_dist_value'])
        input_dict['segname1']=input_dict['chainID1']
        input_dict['segname2']=input_dict['chainID2']
        input_dict.update({
            'ptnr1_label_asym_id':cd['ptnr1_label_asym_id'],
            'ptnr2_label_asym_id':cd['ptnr2_label_asym_id'],
            'ptnr1_label_comp_id':cd['ptnr1_label_comp_id'],
            'ptnr2_label_comp_id':cd['ptnr2_label_comp_id'],
            'ptnr1_label_seq_id':cd['ptnr1_label_seq_id'],
            'ptnr2_label_seq_id':cd['ptnr2_label_seq_id'],
            'ptnr1_auth_asym_id':cd['ptnr1_auth_asym_id'],
            'ptnr2_auth_asym_id':cd['ptnr2_auth_asym_id'],
            'ptnr1_auth_comp_id':cd['ptnr1_auth_comp_id'],
            'ptnr2_auth_comp_id':cd['ptnr2_auth_comp_id'],
            'ptnr1_auth_seq_id':int(cd['ptnr1_auth_seq_id']),
            'ptnr2_auth_seq_id':int(cd['ptnr2_auth_seq_id']),
            })
        input_dict.update({
                'residue1':None,
                'residue2':None,
                'atom1':None,
                'atom2':None,
                'empty':False,
                'segtype1':'UNSET',
                'segtype2':'UNSET',
            })
        super().__init__(input_dict)

    @__init__.register(list)
    def _from_patchlist(self,L):
        s1,ri1=L[1].split(':')
        s2,ri2=L[2].split(':')
        r1,i1=split_ri(ri1)
        r2,i2=split_ri(ri2)
        input_dict={
            'chainID1':s1,
            'resseqnum1':r1,
            'insertion1':i1,
            'chainID2':s2,
            'resseqnum2':r2,
            'insertion2':i2,
            'residue1':None,
            'residue2':None,
            'atom1':None,
            'atom2':None,
            'empty':False,
            'segtype1':'UNSET',
            'segtype2':'UNSET',
        }
        super().__init__(input_dict)
        self.patchname=L[0]
        self.patchorder=[1,2] # default
        self.name1=Link.patch_atomnames[self.patchname][0]
        self.name2=Link.patch_atomnames[self.patchname][1]
        self.altloc1=''
        self.altloc2=''

    @__init__.register(str)
    def _from_shortcode(self,shortcode):
        I,J=shortcode.split('-')
        s1,ri1,a1=I.split('_')
        r1,i1=split_ri(ri1)
        s2,ri2,a2=J.split('_')
        r2,i2=split_ri(ri2)
        input_dict={
            'chainID1':s1,
            'resseqnum1':r1,
            'insertion1':i1,
            'name1':a1,
            'chainID2':s2,
            'resseqnum2':r2,
            'insertion2':i2,
            'name2':a2,
            'residue1':None,
            'residue2':None,
            'atom1':None,
            'atom2':None,
            'empty':False,
            'segtype1':'UNSET',
            'segtype2':'UNSET',
        }
        super().__init__(input_dict)
        self.patchname=''
        self.patchorder=[1,2] # default
        self.altloc1=''
        self.altloc2=''

    def set_patchname(self,force=False):
        """Determine the charmff patch residue name for this link based on
        residue names, atom names, and when necessary, 3D geometry of a 
        particular dihedral angle around the link's bond
        
        New Attributes
        --------------
        patchname: str
           charmff pres name
        
        patchorder: list
            [1,2] if residue 1 appears first in the pres listing
            [2,1] if residue 2 appears first in the pres listing
        """
        if hasattr(self,'patchname') and len(self.patchname)>0 and not force:
            logger.debug(f'Patchname for {str(self)} already set to {self.patchname}')
            return
        self.patchname=''
        self.patchorder=[1,2]
        logger.debug(f'patch assignment for link {str(self)}')
        if not self.residue1 and not self.residue2:
            logger.debug(f'missing residue')
            logger.debug(f'1 {self.residue1}')
            logger.debug(f'2 {self.residue2}')
            return
        logger.debug(f'resname1 {self.resname1} segtype2 {self.segtype2}')
        my_res12=[self.residue1,self.residue2]
        if self.resname1=='ASN' and self.segtype2=='glycan':
            # N-linked glycosylation site (toppar_all36_carb_glycopeptide)
            ICmap=[
                {'ICatomnames':['1CG','1ND2','2C1','2O5'],
                 'mapping':{'NGLA':168.99,'NGLB':-70.91}}
            ]
            self.patchname=ic_reference_closest(my_res12,ICmap)
        elif self.resname1=='SER' and self.segtype2=='glycan':
            # O-linked to serine (toppar_all36_carb_glycopeptide)
            ICmap=[
                {'ICatomnames':['1CB','1OG','2C1','2O5'],
                 'mapping':{'SGPA':45.37,'SGPB':19.87}}
            ]
            self.patchname=ic_reference_closest(my_res12,ICmap)
        elif self.resname1=='THR' and self.segtype2=='glycan':
            # O-linked to serine (toppar_all36_carb_glycopeptide)
            ICmap=[
                {'ICatomnames':['1CB','1OG1','2C1','2O5'],
                 'mapping':{'SGPA':69.9,'SGPB':33.16}}
            ]
            self.patchname=ic_reference_closest(my_res12,ICmap)
        elif self.name2=='C1' and self.segtype2=='glycan' and self.segtype1=='glycan':
            # all taken from 1xyz pres in top_all36_carb.rtf
            # including PHI angles of ICs with atoms in both residues
            if self.name1=='O1': # 1->1 link
                ICmap=[
                    {'ICatomnames':'1O5  1C1  1O1  2C1'.split(),
                     'mapping':{'11aa':103.46,'11ab':121.75,'11bb':-56.58}
                    },
                    {'ICatomnames':'1C1  1O1  2C1  2O5'.split(),
                     'mapping':{'11aa':103.54,'11ab':51.80,'11bb':-79.64}
                    },
                    {'atomnames':'1O1  2C1  2O5  2C5'.split(),
                     'mapping':{'11aa':64.56,'11ab':167.51,'11bb':172.18}}
                ]
                self.patchname=ic_reference_closest(my_res12,ICmap)
            elif self.name1=='O2': # 1->2 link
                ICmap=[
                    {'ICatomnames':'1C1  1C2  1O2  2C1'.split(),
                     'mapping':{'12aa':-132.81,'12ab':115.32,'12ba':-133.78,'12bb':117.14}
                    },
                    {'ICatomnames':'1C2  1O2  2C1  2O5'.split(),
                     'mapping':{'12aa':47.16,'12ab':86.93,'12ba':168.07,'12bb':-168.07}
                    }
                ]
                self.patchname=ic_reference_closest(my_res12,ICmap)
            elif self.name1=='O3': # 1->3 link
                ICmap=[
                    {'ICatomnames':'1C2  1C3  1O3  2C1'.split(),
                     'mapping':{'13aa':113.19,'13ab':-141.32,'13ba':-131.68,'13bb':-141.32}
                    },
                    {'ICatomnames':'1C3  1O3  2C1  2O5'.split(),
                     'mapping':{'13aa':65.46,'13ab':65.46,'13ba':-100.16,'13bb':-130.16}
                    }
                ]
                self.patchname=ic_reference_closest(my_res12,ICmap)
            elif self.name1=='O4': # 1->4 link
                ICmap=[
                    {'ICatomnames':'1C3  1C4  1O4  2C1'.split(),
                     'mapping':{'14aa':-86.29,'14ab':72.71,'14ba':-86.3,'14bb':81.86}},
                    {'ICatomnames':'1C4  1O4  2C1  2O5'.split(),
                     'mapping':{'14aa':133.57,'14ab':48.64,'14ba':-130.97,'14bb':-130.97}}
                ]
                self.patchname=ic_reference_closest(my_res12,ICmap)
            elif self.name1=='O6': # 1->6 link
                ICmap=[
                    {'ICatomnames':'1C6  1O6  2C1  2O5'.split(),
                     'mapping':{'16AT':71.24,'16BT':-63.49}}
                ]
                self.patchname=ic_reference_closest(my_res12,ICmap)
        elif self.name2=='C2' and self.segtype2=='glycan' and self.segtype1=='glycan':
            if self.name1=='O6':
                self.patchname='SA26AT'
            elif self.name1=='O8':
                self.patchname='SA28AA'
            elif self.name1=='O9':
                self.patchname='SA29AT'

        elif self.name1=='O6' and self.name2=='C2':
                self.patchname='SA26AT' 
        elif 'ZN' in self.resname1 and 'HIS' in self.resname2:
                # links in PDB files involving Zn usually list the Zn atom first, but the patch has it second
                self.patchname='ZNHD'
                self.patchorder=[2,1]
        else:
            logger.warning(f'Could not identify patch for link: {str(self)}')
            self.patchname='UNFOUND'

    def write_TcL(self,W:Psfgen,transform):
        """Insert the appropriate TcL commands to add this link in a psfgen script
        
        Assumes that if one of the two residues is an asparagine and the other is a 
        from a glycan, this requires an CHARMM NGLB patch

        A 2->6 intraglycan linkage requies the SA26T patch

        Others' patches are named 1xij where x is the carbon number of the upstream monomer,
        and i and j indicate whether the bond is axial or equatorial.
        
        Parameters
        ----------
        W: Psfgen
            the psfgen scriptwriter object
        transform: BiomT
            the designated transform under which this link is operational; used for its chainIDmap
        """
        chainIDmap=transform.chainIDmap
        seg1=self.residue1.chainID
        seg1=chainIDmap.get(seg1,seg1)
        seg2=self.residue2.chainID
        seg2=chainIDmap.get(seg2,seg2)
        rsn1=self.residue1.resseqnum
        rsn2=self.residue2.resseqnum
        ins1=self.residue1.insertion
        ins2=self.residue2.insertion
        logger.debug(f'Link: {self.residue1.chainID}->{seg1}:{rsn1}{ins1} {self.residue2.chainID}->{seg2}:{rsn2}{ins2}')
        if not self.patchname=='UNFOUND':
            if self.patchorder==[1,2]:
                W.addline(f'patch {self.patchname} {seg1}:{rsn1}{ins1} {seg2}:{rsn2}{ins2}')
            elif self.patchorder==[2,1]:
                W.addline(f'patch {self.patchname} {seg2}:{rsn2}{ins2} {seg1}:{rsn1}{ins1}')
        else:
            logger.warning(f'Could not identify patch for link: {str(self)}')
            W.comment(f'No patch found for {str(self)}')
    
    def update_residue(self,idx,**fields):
        if idx==1:
            if 'chainID' in fields:
                new_chainID=fields['chainID']
                if self.chainID1!=new_chainID:
                    logger.debug(f'updating link residue1 chainID from {self.chainID1} to {new_chainID}')
                    self.chainID1=new_chainID
                    self.residue1.set(chainID=new_chainID)
        elif idx==2:
            if 'chainID' in fields:
                new_chainID=fields['chainID']
                if self.chainID2!=new_chainID:
                    logger.debug(f'updating link residue2 chainID from {self.chainID2} to {new_chainID}')
                    self.chainID2=new_chainID
                    self.residue2.set(chainID=new_chainID)

    def __str__(self):
        if not hasattr(self,'residue1') or not self.residue1 or not hasattr(self,'residue2') or not self.residue2:
            return f'{self.chainID1}_{self.resname1}{self.resseqnum1}{self.insertion1}-{self.chainID2}_{self.resname2}{self.resseqnum2}{self.insertion2}'
        return f'{str(self.residue1)}-{str(self.residue2)}'

class LinkList(AncestorAwareObjList):
    """A class for handling lists of Links

    Methods
    -------
    assign_resiudes(Residues)
        scans the provided list of residues to assign actual residues to the 'residue1' and 'residue2'
        attributes of each link; sets the 'atom1' and 'atom2' attributes to point to the actual
        atoms; sets up the 'up' and 'down' pointers for every link; sets the segtype1 and segtype2
        attributes of every link; returns list of residues not assigned to links and links to which
        no residues were assigned.

    write_TcL(W,transform)
        calls write_TcL for each link

    prune_mutations(Mutations,Segments)
        Given a list of mutations, removes any links that were declared such that any mutation
        would make the link chemically impossible.  E.g., mutating an ASN of an N-linked
        glycosylation site results in loss of the ASN-BGLNA link.
    
    apply_segtypes(map)
        using the resname-to-segtype map provided, set the 'segtype1' and 'seqtype2' attributes
        of every element.
                
    """
    def assign_residues(self,Residues):
        """Assigns residue and atom pointers to each link; sets up the up and down links of both
        residues so that linked residue objects can reference one another; flags residues from
        list of residues passed in that are not assigned to any links

        Arguments
        ---------
        Residues: ResidueList
            List of residues to search in order to make residue assignments
        
        Returns
        -------
        ResidueList: list of residues from Residues that are not used for any assignments

        """
        logger.debug(f'Links: Assigning residues from list of {len(Residues)} residues')
        ignored_by_ptnr1=self.assign_objs_to_attr('residue1',Residues,resseqnum='resseqnum1',chainID='chainID1',insertion='insertion1')
        ignored_by_ptnr2=self.assign_objs_to_attr('residue2',Residues,resseqnum='resseqnum2',chainID='chainID2',insertion='insertion2')
        if hasattr(self,'patchname') and len(self.patchname)>0:
            # this link was most likely created when reading in patch records from a set of REMARKS in a pre-built psf file.
            # we need to get the precise atom names for this patch
            pass
        for link in self:
            try:
                link.residue1.linkTo(link.residue2,link)
            except:
                raise ValueError(f'Bad residue in link')
            link.atom1=link.residue1.atoms.get(name=link.name1,altloc=link.altloc1)
            link.atom2=link.residue2.atoms.get(name=link.name2,altloc=link.altloc2)
            link.segtype1=link.residue1.segtype
            link.segtype2=link.residue2.segtype
            # shortcodes don't provide resnames, so set them here
            if not hasattr(link,'resname1'):
                link.resname1=link.residue1.resname
            if not hasattr(link,'resname2'):
                link.resname2=link.residue2.resname
            link.set_patchname()
        # do cross-assignment to find true orphan links and dangling links
        orphan_1=ignored_by_ptnr1.assign_objs_to_attr('residue2',Residues,resseqnum='resseqnum2',chainID='chainID2',insertion='insertion2')
        orphan_2=ignored_by_ptnr2.assign_objs_to_attr('residue1',Residues,resseqnum='resseqnum1',chainID='chainID1',insertion='insertion1')
        orphans=orphan_1+orphan_2
        rlist=[]
        for link in ignored_by_ptnr1:
            rlist,list=link.residue2.get_down_group()
            rlist.insert(0,link.residue2)
            for r in rlist:
                Residues.remove(r)
        return Residues.__class__(rlist),self.__class__(ignored_by_ptnr1+ignored_by_ptnr2)
    
    def write_TcL(self,W:Psfgen,transform):
        for l in self:
            l.write_TcL(W,transform)

    # def remove_orphan_residues(self,Links,Residues):
    #     for dl in self:
    #         reslist,lnlist=dl.residue2.get_down_group()
    #         reslist.insert(0,dl.residue2)
    #         for r in reslist:
    #             Residues.remove(r)
    #         for l in lnlist:
    #             Links.remove(l)
    #     return reslist,lnlist

    def prune_mutations(self,Mutations,Segments):
        """Prune off any links and associated objects as a result of mutations
        
        Arguments
        ---------
        Mutations: MutationList
            list of mutations
            
        Segments: SegmentList
            list of assembled segments; might need to be modified if pruning gets rid
            of a whole segment's worth of residues
        
        """
        pruned={'residues':[],'links':self.__class__([]),'segments':[]}
        for m in Mutations:
            rlist,llist=[],[]
            left=self.get(chainID1=m.chainID,resseqnum1=m.resseqnum,insertion1=m.insertion)
            right=self.get(chainID2=m.chainID,resseqnum2=m.resseqnum,insertion2=m.insertion)
            if left:  # this is a link in which this mutation is the left member
                self.remove(left) # get rid of this link
                # we need to remove residue2 and everything downstream
                # remove downstream residues!
                rlist,llist=left.residue2.get_down_group()
                rlist.insert(0,left.residue2)
            elif right: # this is a link in which this mutation is the right member (should be very rare)
                self.remove(right)
                rlist,llist=right.residue2.get_down_group()
                rlist.insert(0,right.residue2)
            if rlist and llist:
                # logger.debug(f'Deleting residues down from and including {str(rlist[0])} due to a mutation')
                S=Segments.get_segment_of_residue(rlist[0])
                for r in rlist:
                    # logger.debug(f'...{str(r)}')
                    pruned['residues'].append(S.residues.remove(r))
                if len(S.residues)==0:
                    # logger.debug(f'All residues of {S.psfgen_segname} are deleted; {S.psfgen_segname} is deleted')
                    Segments.remove(S)
                    pruned['segments'].append(S)
                for l in llist:
                    pruned['links'].append(self.remove(l))
        return pruned

    def apply_segtypes(self,map):
        """Apply segtype values to each of the two residues using the map
        
        Parameters
        ----------
        map: dict
            map of segtypes for given resnames
        """
        self.map_attr('segtype1','resname1',map)
        self.map_attr('segtype2','resname2',map)

    def report(self,W:Filewriter):
        for l in self:
            W.addline(str(l))
