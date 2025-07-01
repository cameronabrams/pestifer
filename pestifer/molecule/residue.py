#Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
Defines the EmptyResidue and Residue classes for handling residues in a molecular structure.
"""
import logging
logger=logging.getLogger(__name__)

from argparse import Namespace
from functools import singledispatchmethod

from pidibble.baserecord import BaseRecord
from pidibble.pdbrecord import PDBRecord

from .atom import Atom, AtomList, Hetatm
from ..core.baseobj import AncestorAwareObj, AncestorAwareObjList
from ..core.labels import Labels
from ..objs.seqadv import Seqadv, SeqadvList
from ..objs.deletion import DeletionList
from ..objs.substitution import SubstitutionList
from ..core.stringthings import join_ri, split_ri
from ..util.cifutil import CIFdict
from ..util.coord import positionN

class EmptyResidue(AncestorAwareObj):
    """
    A class for handling missing residues in a molecular structure.
    This class represents residues that are not present in the structure, such as those that are missing
    due to low resolution or other reasons. It is used to track residues that are expected to be present
    but are not resolved in the coordinate file.
    
    Parameters
    ----------
    input_obj : str, BaseRecord, CIFdict, or EmptyResidue
        The input object can be a string representing a residue shortcode, a BaseRecord or PDBRecord object,
        or a CIFdict object containing residue information. If it is an EmptyResidue, it will be initialized
        with the same attributes.
    """

    req_attr=AncestorAwareObj.req_attr+['resname','resseqnum','insertion','chainID','resolved','segtype']
    """
    Required attributes for EmptyResidue.
    
    Attributes
    ----------
    resname : str
        The residue name.
    resseqnum : int
        The residue sequence number.
    insertion : str
        The insertion code, if applicable.
    chainID : str
        The chain ID.
    resolved : bool
        Indicates whether the residue is resolved (True) or not (False).
    segtype : str
        The segment type.
    """
    
    opt_attr=AncestorAwareObj.opt_attr+['model','id','auth_asym_id','auth_comp_id','auth_seq_id']
    """
    Optional attributes for EmptyResidue.
    
    Attributes
    ----------
    model : int
        The model number, if applicable.
    id : str
        The residue ID.
    auth_asym_id : str
        The author asymmetry ID, if applicable.
    auth_comp_id : str
        The author component ID, if applicable.
    auth_seq_id : int
        The author sequence ID, if applicable.
    """

    yaml_header='missings'
    """
    YAML header for EmptyResidue.
    """
    
    PDB_keyword='REMARK.465'
    """
    PDB keyword for EmptyResidue.
    """

    mmCIF_name='pdbx_unobs_or_zero_occ_residues'
    """
    mmCIF name for EmptyResidue.
    """

    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)
    
    @__init__.register(BaseRecord)
    @__init__.register(PDBRecord)
    def _from_pdbrecord(self,pdbrecord):
        input_dict={
            'model':pdbrecord.modelNum,
            'resname':pdbrecord.resName,
            'chainID':pdbrecord.chainID,
            'resseqnum':pdbrecord.seqNum,
            'insertion':pdbrecord.iCode,
            'resolved':False,
            'segtype':'UNSET'
        }
        super().__init__(input_dict)
    
    @__init__.register(CIFdict)
    def _from_cifdict(self,cd):
        mn=cd['pdb_model_num']
        if type(mn)==str and mn.isdigit:
            nmn=int(mn)
        else:
            nmn=1
        input_dict={
            'model':nmn,
            'resname':cd['label_comp_id'],
            'chainID':cd['label_asym_id'],
            'resseqnum':int(cd['label_seq_id']),
            'insertion':cd['pdb_ins_code'],
            'resolved':False,
            'segtype':'UNSET',
            'auth_asym_id':cd['auth_asym_id'],
            'auth_comp_id':cd['auth_comp_id'],
            'auth_seq_id':int(cd['auth_seq_id']),
        }
        super().__init__(input_dict)

    @__init__.register(str)
    def _from_shortcode(self,shortcode):
        # i:C:RRR###A
        # i model number, optional, default 1
        # C chain ID, optional, default 'A'
        # RRR one or three-letter residue code
        # ### integer residue sequence number
        # A insertion code, optional, default ''
        tokens=shortcode.split(':')
        has_modelnum=len(tokens)>2
        model_idx=-1
        chain_idx=-1
        model=1
        chainID='A'
        if has_modelnum:
            model_idx=0
            model=int(tokens[model_idx])
        has_chainID=len(tokens)>1
        if has_chainID:
            chain_idx=model_idx+1
            chainID=tokens[chain_idx]
        res_idx=chain_idx+1
        resncode=tokens[res_idx]
        rescode=''
        ri=''
        okrescode=True
        for byte in resncode:
            if byte.isalpha() and okrescode:
                rescode=rescode+byte
            else:
                okrescode=False
                ri=ri+byte
        rescode=rescode.upper()
        if len(rescode)==1:
            rescode=Labels.res_123[rescode]
        r,i=split_ri(ri)
        input_dict={
            'model':model,
            'resname':rescode,
            'chainID':chainID,
            'resseqnum':r,
            'insertion':i,
            'resolved':False,
            'segtype':'UNSET',
        }
        super().__init__(input_dict)

    def __str__(self):
        return f'{self.chainID}_{self.resname}{self.resseqnum}{self.insertion}*'

    def pdb_line(self):
        """
        Returns a PDB line representation of the EmptyResidue.
        The line is formatted according to the PDB standard for missing residues.
        
        Returns
        -------
        str
            A string representing the EmptyResidue in PDB format.
        """
        record_name,code=EmptyResidue.PDB_keyword.split('.')
        return '{:6s}{:>4d}   {:1s} {:3s} {:1s} {:>5d}{:1s}'.format(record_name,
        code,self.model,self.resname,self.chainID,self.resseqnum,self.insertion)

class EmptyResidueList(AncestorAwareObjList):
    """
    A class for handling lists of EmptyResidue objects.
    This class is used to manage collections of residues that are not present in the molecular structure,
    such as missing residues in a PDB file.

    This class does not add anything beyond the ``AncestorAwareObjList`` class.
    """
    pass

class Residue(EmptyResidue):
    """
    A class for handling residues in a molecular structure.
    This class extends the EmptyResidue class to include additional functionality for managing residues.
    
    Parameters
    ----------
    input_obj : Atom, Hetatm, EmptyResidue, str, or Residue
        The input object can be an Atom or Hetatm object, an EmptyResidue object, a string representing a residue shortcode,
        or another ``Residue`` object. If it is an EmptyResidue, it will be initialized with the same attributes.    
    """

    req_attr=EmptyResidue.req_attr+['atoms']
    """
    Required attributes for ``Residue`` in addition to those defined for ``EmptyResidue``.

    Attributes
    ----------
    atoms : AtomList
        A list of Atom objects that belong to this residue.
    """

    opt_attr=EmptyResidue.opt_attr+['up','down','uplink','downlink']
    """
    Optional attributes for ``Residue`` in addition to those defined for ``EmptyResidue``.

    - ``up``: list
        A list of residues that are linked to this residue in an upstream direction.
    - ``down``: list
        A list of residues that are linked to this residue in a downstream direction.
    - ``uplink``: list
        A list of links to residues that are connected to this residue in an upstream direction.
    - ``downlink``: list
        A list of links to residues that are connected to this residue in a downstream direction.
    """

    ignore_attr=EmptyResidue.ignore_attr+['atoms','up','down','uplink','downlink']
    """
    Attributes to ignore when comparing Residue objects.
    This includes the attributes defined in ``EmptyResidue`` as well as the additional attributes defined for ``Residue``.
    
    - ``atoms``: AtomList
        A list of Atom objects that belong to this residue.
    - ``up``: list
        A list of residues that are linked to this residue in an upstream direction.
    - ``down``: list
        A list of residues that are linked to this residue in a downstream direction.
    - ``uplink``: list
        A list of links to residues that are connected to this residue in an upstream direction.
    - ``downlink``: list
        A list of links to residues that are connected to this residue in a downstream direction.
    """

    _counter=0

    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)
    
    @__init__.register(Atom)
    @__init__.register(Hetatm)
    def _from_atom(self,a):
        input_dict={
            'resseqnum':a.resseqnum,
            'insertion':a.insertion,
            'resname':a.resname,
            'chainID':a.chainID,
            'atoms':AtomList([a]),
            'segtype':'UNSET',
            'resolved':True
        }
        for cif_xtra in ['auth_seq_id','auth_comp_id','auth_asym_id']:
            if hasattr(a,cif_xtra):
                input_dict[cif_xtra]=a.__dict__[cif_xtra]
        super().__init__(input_dict)
        self.index=Residue._counter
        Residue._counter+=1
        self.down=[]
        self.downlink=[]
        self.up=[]
        self.uplink=[]

    @__init__.register(EmptyResidue)
    def _from_emptyresidue(self,m):
        input_dict={
            'model':m.model,
            'resseqnum':m.resseqnum,
            'insertion':m.insertion,
            'resname':m.resname,
            'chainID':m.chainID,
            'atoms':AtomList([]),
            'segtype':'UNSET',
            'resolved':False
        }
        for cif_xtra in ['auth_asym_id','auth_comp_id','auth_seq_id']:
            if hasattr(m,cif_xtra):
                input_dict[cif_xtra]=m.__dict__[cif_xtra]
        super().__init__(input_dict)
        self.index=Residue._counter
        Residue._counter+=1
        self.down=[]
        self.downlink=[]
        self.up=[]
        self.uplink=[]

    @__init__.register(str)
    def _from_shortcode(self,shortcode):
        self.__init__(EmptyResidue(shortcode))

    def __str__(self):
        return super().__str__()[0:-1] # strip off the "*"

    def __lt__(self,other):
        """
        Compare this residue with another residue or a string representation of a residue to determine if this residue is "less than" the other; i.e., N-terminal to.  This is used for sorting residues in a sequence.
        The comparison is based on the residue sequence number and insertion code.
        If the other object is a string, it is split into residue sequence number and insertion code
        using the ``split_ri`` function.
        
        Parameters
        ----------
        other : Residue, str
            The other residue or string to compare with.
        
        Returns
        -------
        bool
            True if this residue is less than the other residue, False otherwise.
        """
        if hasattr(other,'resseqnum') and hasattr(other,'insertion'):
            o_resseqnum=other.resseqnum
            o_insertion=other.insertion
        elif type(other)==str:
            o_resseqnum,o_insertion=split_ri(other)
        else:
            raise Exception(f'I do not know how to make something of type {type(other)} a residue')
        if self.resseqnum<o_resseqnum:
            return True
        elif self.resseqnum==o_resseqnum:
            return self.insertion<o_insertion
        return False
    
    def __gt__(self,other):
        """
        Compare this residue with another residue or a string representation of a residue to determine if this residue is "greater than" the other; i.e., C-terminal to.  This is used for sorting residues in a sequence.

        Parameters
        ----------
        other : Residue, str
            The other residue or string to compare with.

        Returns
        -------
        bool
            True if this residue is greater than the other residue, False otherwise.
        """
        if hasattr(other,'resseqnum') and hasattr(other,'insertion'):
            o_resseqnum=other.resseqnum
            o_insertion=other.insertion
        elif type(other)==str:
            o_resseqnum,o_insertion=split_ri(other)
        else:
            raise Exception(f'I do not know how to make something of type {type(other)} a residue')
        if self.resseqnum>o_resseqnum:
            return True
        elif self.resseqnum==o_resseqnum:
            return self.insertion>o_insertion
        return False
    
    def __le__(self,other):
        """
        Compare this residue with another residue or a string representation of a residue to determine if this residue is "less than or equal to" the other; i.e., N-terminal to or same as.  This is used for sorting residues in a sequence.

        Parameters
        ----------
        other : Residue, str
            The other residue or string to compare with.

        Returns
        -------
        bool
            True if this residue is less than or equal to the other residue, False otherwise.
        """
        if self<other:
            return True
        return self.same_resid(other)
    
    def __ge__(self,other):
        """
        Compare this residue with another residue or a string representation of a residue to determine if this residue is "greater than or equal to" the other; i.e., C-terminal to or same as.  This is used for sorting residues in a sequence.

        Parameters
        ----------
        other : Residue, str
            The other residue or string to compare with.

        Returns
        -------
        bool
            True if this residue is greater than or equal to the other residue, False otherwise.
        """
        if self>other:
            return True
        return self.same_resid(other)
    
    def same_resid(self,other):
        """
        Check if this residue has the same residue sequence number and insertion code as another residue or a
        string representation of a residue. This is used to determine if two residues are the same in terms of their sequence position.
        
        Parameters
        ----------
        other : Residue, str
            The other residue or string to compare with.
            
        Returns
        -------
        bool
            True if this residue has the same residue sequence number and insertion code as the other residue,
            False otherwise."""
        if type(other)==type(self):
            o_resseqnum=other.resseqnum
            o_insertion=other.insertion
        elif type(other)==str:
            o_resseqnum,o_insertion=split_ri(other)
        return self.resseqnum==o_resseqnum and self.insertion==o_insertion
        
    def add_atom(self,a:Atom):
        """
        Add an atom to this residue if it matches the residue's sequence number, residue name,
        chain ID, and insertion code. This method is used to build a residue from its constituent
        atoms, ensuring that all atoms in the residue share the same sequence number, residue name,
        chain ID, and insertion code.
        
        Parameters
        ----------
        a : Atom
            The atom to be added to the residue.
        """
        if self.resseqnum==a.resseqnum and self.resname==a.resname and self.chainID==a.chainID and self.insertion==a.insertion:
            self.atoms.append(a)
            return True
        return False

    def set_chainID(self,chainID):
        """
        Set the chain ID for this residue and all its constituent atoms.
        This method updates the chain ID of the residue and all atoms within it to ensure consistency
        across the residue's structure. 
        
        Parameters
        ----------
        chainID : str
            The new chain ID to be set for the residue and its atoms.
        """
        self.chainID=chainID
        for a in self.atoms:
            a.chainID=chainID

    def linkTo(self,other,link):
        """
        Link this residue to another residue in a downstream direction.
        This method establishes a connection between this residue and another residue, allowing for
        traversal of the molecular structure in a downstream direction. It also updates the upstream
        and downstream links for both residues to maintain the connectivity information.
        
        Parameters
        ----------
        other : Residue
            The other residue to link to.
        link : str
            A string representing the link between the two residues. This could be a description of the type
            of link or a unique identifier for the link.
        """
        assert type(other)==type(self),f'type of other is {type(other)}; expected {type(self)}'
        self.down.append(other)
        self.downlink.append(link)
        other.up.append(self)
        other.uplink.append(link)

    def unlink(self,other,link):
        """
        Unlink this residue from another residue in a downstream direction.
        This method removes the connection between this residue and another residue, effectively
        severing the link between them. It updates the upstream and downstream links for both residues
        to reflect the removal of the connection.
        
        Parameters
        ----------
        other : Residue
            The other residue to unlink from.
        link : str
            A string representing the link between the two residues. This could be a description of the type
            of link or a unique identifier for the link.
        """
        self.down.remove(other)
        self.downlink.remove(link)
        other.up.remove(self)
        other.uplink.remove(link)

    def ri(self):
        """
        Returns the residue sequence number and insertion code as a string.
        """
        ins0='' if self.insertion==' ' else self.insertion
        return f'{self.resseqnum}{ins0}'
    
    def get_down_group(self):
        """
        Get the downstream group of residues.
        This method traverses the downstream links of this residue and collects all residues
        in the downstream direction.
        """
        res=[]
        lin=[]
        for d,dl in zip(self.down,self.downlink):
            res.append(d)
            lin.append(dl)
            # logger.debug(f'{str(self)}->{str(d)}')
            tres,tlin=d.get_down_group()
            res.extend(tres)
            lin.extend(tlin)
        return res,lin
    
class ResidueList(AncestorAwareObjList):
    """
    A class for handling lists of Residue objects.
    This class extends the AncestorAwareObjList to manage collections of residues in a molecular structure.
    It provides methods for initializing the list from various input types, indexing residues, and performing operations
    such as mapping chain IDs, retrieving residues and atoms, and handling residue ranges.
    """
    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)

    @__init__.register(AtomList)
    def _from_atomlist(self,atoms):
        R=[]
        for a in atoms:
            for r in R[::-1]:
                if r.add_atom(a):
                    break
            else:
                R.append(Residue(a))
        super().__init__(R)
    
    @__init__.register(EmptyResidueList)
    def _from_emptyresiduelist(self,input_list):
        R=[Residue(m) for m in input_list]
        super().__init__(R)

    def index(self,R):
        """
        Get the index of a residue in the list.
        
        Parameters
        ----------
        R : Residue
            The residue to find in the list.
        """
        for i,r in enumerate(self):
            if r is R:
                return i
        else:
            raise ValueError(f'Residue not found')

    def map_chainIDs_label_to_auth(self):
        """
        Create a mapping from chain IDs in the label (e.g., PDB format) to the author chain IDs.
        """
        self.chainIDmap_label_to_auth={}
        for r in self:
            if hasattr(r,'auth_asym_id'):
                label_Cid=r.chainID
                auth_Cid=r.auth_asym_id
                if not label_Cid in self.chainIDmap_label_to_auth:
                    self.chainIDmap_label_to_auth[label_Cid]=auth_Cid

    def get_residue(self,**fields):
        """
        Get a residue from the list based on specified fields.

        Parameters
        ----------
        **fields : keyword arguments
            The fields to match against residues in the list.

        Returns
        -------
        Residue
            The matching residue, or None if not found.
        """
        return self.get(**fields)
    
    def get_atom(self,atname,**fields):
        """
        Get an atom from the list based on its name and specified fields.
        
        Parameters
        ----------
        atname : str
            The name of the atom to retrieve.
        **fields : keyword arguments
            Additional fields to match against the atom's attributes.
        
        Returns
        -------
        Atom
            The matching atom, or None if not found."""
        S=('atoms',{'name':atname})
        return self.get_attr(S,**fields)
    
    def atom_serials(self,as_type=str):
        """
        Get a list of atom serial numbers for all atoms in the residues.
        
        Parameters
        ----------
        as_type : type, optional
            The type to which the serial numbers should be converted. Default is str.

        Returns
        -------
        list
            A list of atom serial numbers."""
        serlist=[]
        for res in self:
            for a in res.atoms:
                serlist.append(as_type(a.serial))
        return serlist
    
    def atom_resseqnums(self,as_type=str):
        """
        Get a list of residue sequence numbers for all atoms in the residues.
        
        Parameters
        ----------
        as_type : type, optional
            The type to which the residue sequence numbers should be converted. Default is str.
        
        Returns
        -------
        list
            A list of residue sequence numbers.
        """
        rlist=[]
        for res in self:
            for a in res.atoms:
                rlist.append(as_type(a.resseqnum))
        return rlist
    
    def caco_str(self,upstream_reslist,seglabel,molid_varname,tmat):
        """
        Generate a string representation of the ``caco`` command to position
        the N atom of a model-built residue based on the coordinates of the previous residue.
        
        Parameters
        ----------
        upstream_reslist : list
            A list of residues that are upstream of the current residue.
        seglabel : str
            The segment label for the residue.
        molid_varname : str
            The variable name for the molecule ID.
        tmat : numpy.ndarray
            A transformation matrix to apply to the coordinates of the residue.

        Returns
        -------
        str
            A string representing the ``caco`` command for positioning the residue.
        """
        r0=self[0]
        ur=upstream_reslist[-1]
        rN=positionN(ur,tmat)
        logger.debug(f'caco {rN}')
        return f'coord {seglabel} {r0.resseqnum}{r0.insertion} N {{{rN[0]:.5f} {rN[1]:.5f} {rN[2]:.5f}}}'
    
    def resrange(self,rngrec):
        """
        Get a range of residues based on a ResidueRange record.
        
        Parameters
        ----------
        rngrec : ResidueRange
            A record defining the range of residues to retrieve. It should contain the chain ID, residue sequence numbers, and insertion codes for the start and end of the range.

        Yields
        ------
        Residue
            Yields residues within the specified range.
        """
        subR=self.get(chainID=rngrec.chainID)
        subR.sort()
        r1=rngrec.resseqnum1
        i1=rngrec.insertion1
        r2=rngrec.resseqnum2
        i2=rngrec.insertion2
        R1=subR.get(resseqnum=r1,insertion=i1)
        if R1:
            R2=subR.get(resseqnum=r2,insertion=i2)
            if R2:
                idx1=subR.index(R1)
                idx2=subR.index(R2)
                assert idx2>=idx1
                for j in range(idx1,idx2+1):
                    yield self[j]
        return []
    
    def do_deletions(self,Deletions):
        """
        Apply a list of deletions to the residue list.
        
        Parameters
        ----------
        Deletions : list of ResidueRange
            A list of deletion ranges to apply. Each deletion range should contain the chain ID, residue sequence numbers, and insertion codes for the start and end of the deletion range.
        """
        delete_us=[]
        for d in Deletions:
            for dr in self.resrange(d):
                delete_us.append(dr)
        for d in delete_us:
            self.remove(d)

    def apply_segtypes(self):
        """
        Apply segment types to residues based on their residue names.
        This method uses the ``Labels.segtype_of_resname`` mapping to assign segment types to
        residues based on their residue names. It updates the `segtype` attribute of each residue
        in the residue list.
        """
        # logger.debug(f'residuelist:apply_segtypes {segtype_of_resname}')
        self.map_attr('segtype','resname',Labels.segtype_of_resname)
    
    def deletion(self,DL:DeletionList):
        """
        Remove residues from the residue list based on a DeletionList.
        This method iterates through the DeletionList, retrieves the residues to be deleted based on their
        chain ID, residue sequence numbers, and insertion codes, and removes them from the residue list
        
        Parameters
        ----------
        DL : DeletionList
            A list of deletion ranges to apply. Each deletion range should contain the chain ID, residue sequence numbers, and insertion codes for the start and end of the deletion range.
        """
        excised=[]
        for d in DL:
            chain=self.get(chainID=d.chainID)
            r1=chain.get(resseqnum=d.resseqnum1,insertion=d.insertion1)
            r2=chain.get(resseqnum=d.resseqnum2,insertion=d.insertion2)
            for r in chain:
                if r1<=r<=r2:
                    excised.append(r)
        for x in excised:
            self.remove(x)
            assert not x in self
        return excised

    def substitutions(self,SL:SubstitutionList):
        """
        Apply a list of substitutions to the residue list.
        This method iterates through the SubstitutionList, retrieves the residues to be substituted based on their
        chain ID, residue sequence numbers, and insertion codes, and replaces their residue names with the
        corresponding residue names from the substitution list. It also creates a new SeqadvList for any
        resolved residues that are substituted.
        
        Parameters
        ----------
        SL : SubstitutionList
            A list of substitutions to apply. Each substitution should contain the chain ID, residue sequence numbers,
            insertion codes, and the new residue sequence to substitute.
        
        Returns
        -------
        tuple
            A tuple containing a SeqadvList of new sequence advancements for resolved residues and a list of residues that were deleted.
        """
        delete_us=[]
        newseqadv=SeqadvList([]) # for holding single-residue changes for resolved residues
        for s in SL:
            subseq=s.subseq
            currsubidx=0
            chain=self.get(chainID=s.chainID)
            assert chain!=None,f'Error: no chain {s.chainID}'
            r1=chain.get(resseqnum=s.resseqnum1,insertion=s.insertion1)
            assert r1!=None,f'Error: no resseqnum {s.resseqnum1} insertion [{s.insertion1}]'
            r2=chain.get(resseqnum=s.resseqnum2,insertion=s.insertion2)
            assert r2!=None,f'Error: no resseqnum {s.resseqnum2} insertion [{s.insertion2}]'
            for r in chain:
                if r1<=r<=r2:
                    if currsubidx<len(subseq):
                        resname=Labels.res_123[subseq[currsubidx].upper()]
                        if r.resolved: # make a new seqadv for this mutation
                            input_dict={
                                'idCode':'I doubt I ever use this',
                                'typekey':'user',
                                'resname':r.resname,
                                'chainID':r.chainID,
                                'resseqnum':r.resseqnum,
                                'insertion':r.insertion,
                                'dbRes':resname
                            }
                            newseqadv.append(Seqadv(input_dict))
                        else:  # just change the residue name
                            r.name=resname
                        currsubidx+=1
                    else:
                        delete_us.append(r)
            if currsubidx<len(subseq):
                # we have unsubsituted residue(s) left that must be inserted
                pass
        for r in delete_us:
            self.remove(r)
        return newseqadv,delete_us

    def cif_residue_map(self):
        """
        Create a mapping of residues by chain ID and residue sequence number.
        This method iterates through the residue list and creates a dictionary where the keys are chain IDs
        and the values are dictionaries mapping residue sequence numbers to Namespace objects containing
        the residue information.
        """
        result={}
        for r in self:
            if hasattr(r,'label_asym_id'):
                if not r.chainID in result:
                    result[r.chainID]={}
                if not r.resseqnum in result[r.chainID]:
                    result[r.chainID][r.resseqnum]=Namespace(resseqnum=r.label_seq_id,chainID=r.label_asym_id,insertion=r.insertion)
        return result
    
    def apply_insertions(self,insertions):
        """
        Apply a list of insertions to the residue list.
        
        Parameters
        ----------
        insertions : list of Insertion
            A list of insertions to apply. Each insertion should contain the chain ID, residue sequence numbers,
            insertion codes, and the sequence of residues to insert.
        """
        for ins in insertions:
            c,r,i=ins.chainID,ins.resseqnum,ins.insertion
            inc_code=ins.integer_increment
            idx=self.iget(chainID=c,resseqnum=r,insertion=i)
            chainID=self[idx].chainID
            logger.debug(f'insertion begins after {r}{i} which is index {idx} in reslist, chain {chainID}')
            # add residues to residue list
            idx+=1
            i='A' if i in [' ',''] else chr(ord(i)+1)
            for olc in ins.sequence:
                if inc_code:
                    r+=1
                    i=''
                shortcode=f'{chainID}:{olc}{r}{i}'
                new_residue=EmptyResidue(shortcode)
                new_residue.atoms=AtomList([])
                new_residue.segtype='protein'
                self.insert(idx,new_residue)
                logger.debug(f'insertion of new empty residue {shortcode}')
                idx+=1
                i='A' if i in [' ',''] else chr(ord(i)+1)
    
    def renumber(self,links):
        """
        The possibility exists that empty residues added have resseqnums that conflict with existing resseqnums on the same chain if those resseqnums are in a different segtype (e.g., glycan).  This method will privilege protein residues in such conflicts, and it will renumber non-protein residues, updating any resseqnum records in links
        to match the new resseqnums.
        
        Parameters
        ----------
        links : LinkList
            A list of links that may contain residue sequence numbers that need to be updated.
        """
        protein_residues=self.get(segtype='protein')
        if len(protein_residues)==0: return
        min_protein_resseqnum=min([x.resseqnum for x in protein_residues])
        max_protein_resseqnum=max([x.resseqnum for x in protein_residues])
        non_protein_residues=ResidueList([])
        for p in self:
            if p.segtype!='protein':
                non_protein_residues.append(p)
        logger.debug(f'There are {len(protein_residues)} (resseqnums {min_protein_resseqnum} to {max_protein_resseqnum}) protein residues and {len(non_protein_residues)} non-protein residues')
        assert len(self)==(len(protein_residues)+len(non_protein_residues))
        non_protein_residues_in_conflict=ResidueList([])
        for np in non_protein_residues:
            # logger.debug(f'looking for conflict among protein for chain [{np.chainID}] resseqnum [{np.resseqnum}] insertion [{np.insertion}]')
            tst=protein_residues.get(chainID=np.chainID,resseqnum=np.resseqnum,insertion=np.insertion)
            if tst:
                # logger.debug(f'found it!')
                non_protein_residues_in_conflict.append(np)
        for npc in non_protein_residues_in_conflict:
            non_protein_residues.remove(npc)
        logger.debug(f'There are {len(non_protein_residues_in_conflict)} non-protein residues with resseqnums that conflict with protein residues')

        max_unused_resseqnum=max([max([x.resseqnum for x in protein_residues]),0 if len(non_protein_residues)==0 else max([x.resseqnum for x in non_protein_residues])])+1
        newtst=self.get(resseqnum=max_unused_resseqnum)
        assert newtst==[] #None
        mapper_by_chain={}
        for npc in non_protein_residues_in_conflict:
            c=npc.chainID
            if not c in mapper_by_chain:
                mapper_by_chain[c]={}
            old=join_ri(npc.resseqnum,npc.insertion)
            new=max_unused_resseqnum
            mapper_by_chain[c][old]=new
            max_unused_resseqnum+=1
            npc.resseqnum=new
            npc.insertion=''
            for a in npc.atoms:
                a.resseqnum=new
                a.insertion=''
            logger.debug(f'New resid: {c} {old} -> {new}')
        if mapper_by_chain:
            logger.debug(f'Remapping resids in links')
            for l in links:
                if l.chainID1 in mapper_by_chain:
                    old=join_ri(l.resseqnum1,l.insertion1)
                    if old in mapper_by_chain[l.chainID1]:
                        l.resseqnum1=mapper_by_chain[l.chainID1][old]
                        l.insertion1=''
                        logger.debug(f' remapped 1eft: {l.chainID1} {old} -> {l.resseqnum1}')
                if l.chainID2 in mapper_by_chain:
                    old=join_ri(l.resseqnum2,l.insertion2)
                    if old in mapper_by_chain[l.chainID2]:
                        l.resseqnum2=mapper_by_chain[l.chainID2][old]
                        l.insertion2=''
                        logger.debug(f' remapped right: {l.chainID2} {old} -> {l.resseqnum2}')

    def set_chainIDs(self,chainID):
        """
        Set the chain ID for all residues in the list to a specified value.
        
        Parameters
        ----------
        chainID : str
            The chain ID to set for all residues.
        """
        for r in self:
            r.set_chainID(chainID)

    def remap_chainIDs(self,the_map):
        """
        Remap the chain IDs of residues in the list according to a provided mapping.

        Parameters
        ----------
        the_map : dict
            A dictionary mapping old chain IDs to new chain IDs.
        """
        for r in self:
            if r.chainID in the_map:
                r.set_chainID(the_map[r.chainID])
