#Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
Defines the EmptyResidue and Residue classes for handling residues in a molecular structure.
"""

from __future__ import annotations

import logging
logger=logging.getLogger(__name__)

from argparse import Namespace
from functools import singledispatchmethod

from pydantic import Field, ConfigDict
from typing import Union, Any, ClassVar, Optional

from pidibble.baserecord import BaseRecord
from pidibble.pdbrecord import PDBRecord, PDBRecordDict

from .atom import Atom, AtomList, Hetatm
from ..objs.link import LinkList
from ..core.baseobj_new import BaseObj, BaseObjList
from ..core.labels import Labels
from ..objs.seqadv import Seqadv, SeqadvList
from ..objs.deletion import DeletionList
from ..objs.substitution import SubstitutionList
from ..core.stringthings import join_ri, split_ri
from ..util.cifutil import CIFdict
from ..util.coord import positionN
from .stateinterval import StateInterval, StateIntervalList
from mmcif.api.PdbxContainers import DataContainer

class EmptyResidue(BaseObj):
    """
    A class for handling missing residues in a molecular structure.
    This class represents residues that are not present in the structure, such as those that are missing
    due to low resolution or other reasons. It is used to track residues that are expected to be present
    but are not resolved in the coordinate file.
    """

    _required_fields = {'resname','resseqnum','insertion','chainID','resolved','segtype'}
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

    _optional_fields = {'model','obj_id','auth_asym_id','auth_comp_id','auth_seq_id'}
    """
    Optional attributes for EmptyResidue.
    
    Attributes
    ----------
    model : int
        The model number, if applicable.
    obj_id : str
        The residue ID.
    auth_asym_id : str
        The author asymmetry ID, if applicable.
    auth_comp_id : str
        The author component ID, if applicable.
    auth_seq_id : int
        The author sequence ID, if applicable.
    """

    _ignore_fields = {'empty','link','recordname'}
    """
    Fields that are ignored when comparing EmptyResidue objects.
    """
    
    resname: str = Field(..., description="The residue name.")
    resseqnum: int = Field(..., description="The residue sequence number.")
    insertion: str = Field(..., description="The insertion code, if applicable.")
    chainID: str = Field(..., description="The chain ID.")
    resolved: bool = Field(..., description="Indicates whether the residue is resolved (True) or not (False).")
    segtype: str = Field(..., description="The segment type.")
    model: int = Field(1, description="The model number, if applicable.")
    obj_id: int = Field(None, description="The residue ID.")
    auth_asym_id: str = Field(None, description="The author asymmetry ID, if applicable.")
    auth_comp_id: str = Field(None, description="The author component ID, if applicable.")
    auth_seq_id: int = Field(None, description="The author sequence ID, if applicable.")
    empty: bool = Field(False, description="Indicates whether the residue is empty (True) or not (False).")
    link: Any = Field(None, description="The link to the residue, if applicable.")
    recordname: str = Field('REMARK.465', description="The PDB record name for the residue.")

    _yaml_header: ClassVar[str] = 'missings'
    """
    YAML header for EmptyResidue.
    """

    _PDB_keyword: ClassVar[str] = 'REMARK.465'
    """
    PDB keyword for EmptyResidue.
    """
    _PDB_table_of_keyword: ClassVar[str] = 'MISSING'
    """
    PDB table keyword for EmptyResidue.
    """

    _CIF_CategoryName: ClassVar[str] = 'pdbx_unobs_or_zero_occ_residues'
    """
    mmCIF category name for EmptyResidue.
    """

    def describe(self):
        return f"<EmptyResidue {self.chainID}_{self.resname}{self.resseqnum}{self.insertion} ({'resolved' if self.resolved else 'unresolved'})>"
    
    class Adapter:
        """
        Adapter class for converting between different representations of EmptyResidue.
        This class provides methods to create an EmptyResidue from a PDBRecord, CIFdict, or a shortcode.
        """
        def __init__(self, resname=None, resseqnum=None, insertion=None, chainID=None, model=1, resolved=False, segtype='UNSET', auth_asym_id=None, auth_comp_id=None, auth_seq_id=None, obj_id=None):
            self.resname = resname
            self.resseqnum = resseqnum
            self.insertion = insertion
            self.chainID = chainID
            self.model = model
            self.resolved = resolved
            self.segtype = segtype
            self.auth_asym_id = auth_asym_id
            self.auth_comp_id = auth_comp_id
            self.auth_seq_id = auth_seq_id
            self.obj_id = obj_id

        def to_dict(self):
            return {k:v for k,v in self.__dict__.items() if v is not None}

        def to_input_str(self):
            return f"{self.model}:{self.chainID}:{self.resname}{self.resseqnum}{self.insertion}"
        
        @classmethod
        def from_pdbrecord(cls, pdbrecord: Union[PDBRecord, BaseRecord]) -> "EmptyResidue.Adapter":
            if pdbrecord.modelNum is None:
                pdbrecord.modelNum = 1
            elif type(pdbrecord.modelNum) == str and pdbrecord.modelNum.isdigit():
                pdbrecord.modelNum = int(pdbrecord.modelNum)
            elif type(pdbrecord.modelNum) == str and len(pdbrecord.modelNum) == 0:
                pdbrecord.modelNum = 1
            input_dict={
                'model':pdbrecord.modelNum,
                'resname':pdbrecord.resName,
                'chainID':pdbrecord.chainID,
                'resseqnum':pdbrecord.seqNum,
                'insertion':pdbrecord.iCode,
                'resolved':False,
                'segtype':'UNSET'
            }
            return cls(**input_dict)

        @classmethod
        def from_cifdict(cls, cd: CIFdict) -> "EmptyResidue.Adapter":
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
            return cls(**input_dict)

        @classmethod
        def from_shortcode(cls, shortcode: str) -> "EmptyResidue.Adapter":
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
            return cls(**input_dict)
    
    @BaseObj.from_input.register(Adapter)
    @classmethod
    def _from_adapter(cls, adapter: Adapter):
        """
        Create an EmptyResidue from an Adapter instance.
        
        Parameters
        ----------
        adapter : Adapter
            An instance of the Adapter class containing residue information.
        
        Returns
        -------
        EmptyResidue
            An instance of EmptyResidue initialized with the data from the adapter.
        """
        input_dict = adapter.to_dict()
        return cls(**input_dict)

    @singledispatchmethod
    @classmethod
    def new(cls, input_obj: Any):
        raise TypeError(f"Cannot create EmptyResidue from {type(input_obj)}. Use EmptyResidue.Adapter.from_pdbrecord or EmptyResidue.Adapter.from_cifdict instead.")
    
    @new.register(BaseRecord|PDBRecord)
    @classmethod
    def _from_pdbrecord(cls, pdbrecord: Union[BaseRecord, PDBRecord]) -> "EmptyResidue":
        logger.debug(f'Creating EmptyResidue from PDBRecord: {pdbrecord}')
        adapter = cls.Adapter.from_pdbrecord(pdbrecord)
        return cls._from_adapter(adapter)

    @new.register(CIFdict)
    @classmethod
    def _from_cifdict(cls, cd: CIFdict) -> "EmptyResidue":
        adapter = cls.Adapter.from_cifdict(cd)
        return cls._from_adapter(adapter)
    
    @new.register(str)
    @classmethod
    def _from_shortcode(cls, shortcode: str) -> "EmptyResidue":
        adapter = cls.Adapter.from_shortcode(shortcode)
        return cls._from_adapter(adapter)

    def __str__(self):
        return f'{self.chainID}_{self.resname}{self.resseqnum}{self.insertion}*'

    def pdb_line(self):
        """
        Returns a PDB line representation of the :class:`EmptyResidue`.
        The line is formatted according to the PDB standard for missing residues.
        
        Returns
        -------
        str
            A string representing the :class:`EmptyResidue` in PDB format.
        """
        record_name,code=EmptyResidue.PDB_keyword.split('.')
        return '{:6s}{:>4d}   {:1s} {:3s} {:1s} {:>5d}{:1s}'.format(record_name,
        code,self.model,self.resname,self.chainID,self.resseqnum,self.insertion)

class EmptyResidueList(BaseObjList[EmptyResidue]):
    """
    A class for handling lists of :class:`EmptyResidue` objects.
    This class is used to manage collections of residues that are not present in the molecular structure,
    such as missing residues in a PDB file.

    This class does not add anything beyond the :class:`~pestifer.core.baseobj.BaseObjList` class.
    """
    def describe(self):
        return f'<EmptyResidueList with {len(self)} residues>'
    
    @classmethod
    def from_pdb(cls, parsed: PDBRecordDict) -> "EmptyResidueList":
        """
        Create an EmptyResidueList from parsed PDB data.
        
        Parameters
        ----------
        parsed : PDBRecordDict
            A dictionary containing parsed PDB records.
        
        Returns
        -------
        EmptyResidueList
            A new instance of EmptyResidueList containing the missing residues.
        """
        if EmptyResidue._PDB_keyword not in parsed:
            return cls([])
        return cls([EmptyResidue.new(p) for p in parsed[EmptyResidue._PDB_keyword].tables[EmptyResidue._PDB_table_of_keyword]])
    
    @classmethod
    def from_cif(cls, cif_data: DataContainer) -> "EmptyResidueList":
        """
        Create an EmptyResidueList from CIF data.
        
        Parameters
        ----------
        cif_data : DataContainer
            A DataContainer instance containing CIF data.
        
        Returns
        -------
        EmptyResidueList
            A new instance of EmptyResidueList containing the missing residues.
        """
        obj = cif_data.getObj(EmptyResidue._CIF_CategoryName)
        if obj is None:
            return cls([])
        return cls([EmptyResidue.new(CIFdict(obj, i)) for i in range(len(obj))])

class Residue(EmptyResidue):
    """
    A class for handling residues in a molecular structure.
    This class extends the :class:`EmptyResidue` class to include additional functionality for managing residues.
    """
    model_config = ConfigDict(arbitrary_types_allowed=True)
    _required_fields = EmptyResidue._required_fields | {'atoms'}
    """
    Required attributes for :class:`~pestifer.molecule.residue.Residue` in addition to those defined for :class:`~pestifer.molecule.residue.EmptyResidue`.

    Attributes
    ----------
    atoms : :class:`~pestifer.molecule.atom.AtomList`
        A list of :class:`~pestifer.molecule.atom.Atom` objects that belong to this residue.
    """
    atoms: AtomList = Field(default_factory=AtomList, description="A list of atoms that belong to this residue.")

    # _optional_fields = EmptyResidue._optional_fields | {'uplink', 'downlink'}
    _optional_fields = EmptyResidue._optional_fields | {'up', 'down', 'uplink', 'downlink'}
    """
    Optional attributes for :class:`~pestifer.molecule.residue.Residue` in addition to those defined for :class:`~pestifer.molecule.residue.EmptyResidue`.

    - ``up``: list
        A list of residues that are linked to this residue in an upstream direction.
    - ``down``: list
        A list of residues that are linked to this residue in a downstream direction.
    - ``uplink``: list
        A list of links to residues that are connected to this residue in an upstream direction.
    - ``downlink``: list
        A list of links to residues that are connected to this residue in a downstream direction.
    """
    up: Optional['ResidueList'] = Field(default_factory=lambda: ResidueList())  # = Field(default=[], description="A list of residues linked to this residue in an upstream direction.")
    down: Optional['ResidueList'] = Field(default_factory=lambda: ResidueList())  # = Field(default=[], description="A list of residues linked to this residue in a downstream direction.")
    uplink: LinkList = Field(default_factory=LinkList, description="A list of links to residues connected to this residue in an upstream direction.")
    downlink: LinkList = Field(default_factory=LinkList, description="A list of links to residues connected to this residue in a downstream direction.")

    _ignore_fields = EmptyResidue._ignore_fields | {'atoms', 'up', 'down', 'uplink', 'downlink'}
    """
    Attributes to ignore when comparing Residue objects.
    This includes the attributes defined in :class:`~pestifer.molecule.residue.EmptyResidue` as well as the additional attributes defined for :class:`~pestifer.molecule.residue.Residue`.

    - ``atoms``: :class:`~pestifer.molecule.atom.AtomList`
        A list of :class:`~pestifer.molecule.atom.Atom` objects that belong to this residue.
    - ``up``: list
        A list of residues that are linked to this residue in an upstream direction.
    - ``down``: list
        A list of residues that are linked to this residue in a downstream direction.
    - ``uplink``: list
        A list of links to residues that are connected to this residue in an upstream direction.
    - ``downlink``: list
        A list of links to residues that are connected to this residue in a downstream direction.
    """

    _counter: ClassVar[int] = 0

    _chainIDmap_label_to_auth = {}

    def describe(self):
        """
        Returns a string description of the Residue object, including its chain ID, residue name, sequence number, and insertion code.
        
        Returns
        -------
        str
            A string representation of the Residue object.
        """
        return f'<Residue {self.chainID}_{self.resname}{self.resseqnum}{self.insertion} ({len(self.atoms)} atoms)>'

    class Adapter(EmptyResidue.Adapter):
        """
        Adapter class for converting between different representations of Residue.
        This class provides methods to create a Residue from a PDBRecord, CIFdict, or a shortcode.
        It extends the :class:`EmptyResidue.Adapter` class to include additional attributes specific to Residue.
        """
        def __init__(self, resname=None, resseqnum=None, insertion=None, chainID=None, model=1, resolved=False, segtype='UNSET', auth_asym_id=None, auth_comp_id=None, auth_seq_id=None, obj_id=None, atoms: AtomList=None):
            super().__init__(resname=resname, resseqnum=resseqnum, insertion=insertion, chainID=chainID, model=model, resolved=resolved, segtype=segtype, auth_asym_id=auth_asym_id, auth_comp_id=auth_comp_id, auth_seq_id=auth_seq_id, obj_id=obj_id)
            self.atoms = atoms

        @classmethod
        def from_atom(cls, a: Union[Atom, Hetatm]) -> "Residue.Adapter":
            """
            Create a Residue.Adapter from an Atom or Hetatm object.
            
            Parameters
            ----------
            a : :class:`~pestifer.molecule.atom.Atom` or :class:`~pestifer.molecule.atom.Hetatm`
                The atom or hetatm object to create the residue adapter from.
            
            Returns
            -------
            Residue.Adapter
                An instance of Residue.Adapter initialized with the atom's attributes.
            """
            return cls(resname=a.resname, resseqnum=a.resseqnum, insertion=a.insertion, chainID=a.chainID,
                       resolved=True, segtype='UNSET', atoms = AtomList([a]),
                       auth_asym_id=getattr(a, 'auth_asym_id', None),
                       auth_comp_id=getattr(a, 'auth_comp_id', None),
                       auth_seq_id=getattr(a, 'auth_seq_id', None))
        

        @classmethod
        def from_emptyresidue(cls, m: "EmptyResidue") -> "Residue.Adapter":
            """
            Create a Residue.Adapter from an EmptyResidue object.
            
            Parameters
            ----------
            m : :class:`~pestifer.molecule.residue.EmptyResidue`
                The empty residue object to create the residue adapter from.
            
            Returns
            -------
            Residue.Adapter
                An instance of Residue.Adapter initialized with the empty residue's attributes.
            """
            return cls(resname=m.resname, resseqnum=m.resseqnum, insertion=m.insertion, chainID=m.chainID,
                       model=m.model, segtype=m.segtype, resolved=False, atoms=AtomList([]),
                       auth_asym_id=getattr(m, 'auth_asym_id', None),
                       auth_comp_id=getattr(m, 'auth_comp_id', None),
                       auth_seq_id=getattr(m, 'auth_seq_id', None))

    @EmptyResidue.new.register(Atom|Hetatm)
    @classmethod
    def _from_atom(cls, a: Union[Atom, Hetatm]) -> "Residue":
        """
        Create a Residue from an Atom or Hetatm object.
        
        Parameters
        ----------
        a : :class:`~pestifer.molecule.atom.Atom` or :class:`~pestifer.molecule.atom.Hetatm`
            The atom or hetatm object to create the residue from.
        
        Returns
        -------
        Residue
            An instance of Residue initialized with the atom's attributes.
        """
        adapter = cls.Adapter.from_atom(a)
        return cls._from_adapter(adapter)

    @EmptyResidue.new.register(EmptyResidue)
    @classmethod
    def _from_emptyresidue(cls, m: "EmptyResidue") -> "Residue":
        """
        Create a Residue from an EmptyResidue object.
        
        Parameters
        ----------
        m : :class:`~pestifer.molecule.residue.EmptyResidue`
            The empty residue object to create the residue from.
        
        Returns
        -------
        Residue
            An instance of Residue initialized with the empty residue's attributes.
        """
        adapter = cls.Adapter.from_emptyresidue(m)
        return cls._from_adapter(adapter)

    def __str__(self):
        return super().__str__()[0:-1] # strip off the "*"

    def __lt__(self,other):
        """
        Compare this residue with another residue or a string representation of a residue to determine if this residue is "less than" the other; i.e., N-terminal to.  This is used for sorting residues in a sequence.
        The comparison is based on the residue sequence number and insertion code.
        If the other object is a string, it is split into residue sequence number and insertion code
        using the :func:`~pestifer.core.stringthings.split_ri` function.

        Parameters
        ----------
        other : :class:`~pestifer.molecule.residue.Residue`, str
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
        Compare this residue with another residue or a string representation of a residue to determine if this residue is "greater than" the other; i.e., C-terminal to.  This is used for sorting residues in a sequence. If the other object is a string, it is split into residue sequence number and insertion code using the :func:`~pestifer.core.stringthings.split_ri` function.

        Parameters
        ----------
        other : :class:`~pestifer.molecule.residue.Residue`, str
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
        other : :class:`~pestifer.molecule.residue.Residue`, str
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
        other : :class:`~pestifer.molecule.residue.Residue`, str
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
        other : :class:`~pestifer.molecule.residue.Residue`, str
            The other residue or string to compare with.
            
        Returns
        -------
        bool
            True if this residue has the same residue sequence number and insertion code as the other residue,
            False otherwise.
        """
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
        a : :class:`~pestifer.molecule.atom.Atom`
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
        self.set(chainID=chainID)
        for a in self.atoms:
            a.chainID=chainID

    def link_to(self,other,link):
        """
        Link this residue to another residue in a downstream direction.
        This method establishes a connection between this residue and another residue, allowing for
        traversal of the molecular structure in a downstream direction. It also updates the upstream
        and downstream links for both residues to maintain the connectivity information.
        
        Parameters
        ----------
        other : :class:`~pestifer.molecule.residue.Residue`
            The other residue to link to.
        link : :class:`~pestifer.objs.link.Link`
            A link object representing the connection between the two residues.
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
        other : :class:`~pestifer.molecule.residue.Residue`
            The other residue to unlink from.
        link : :class:`~pestifer.objs.link.Link`
            A link object representing the connection between the two residues.
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
        res=ResidueList()
        lin=LinkList()
        for d,dl in zip(self.down,self.downlink):
            res.append(d)
            lin.append(dl)
            # logger.debug(f'{str(self)}->{str(d)}')
            tres,tlin=d.get_down_group()
            res.extend(tres)
            lin.extend(tlin)
        return res,lin


class ResidueList(BaseObjList[Residue]):
    """
    A class for handling lists of :class:`~pestifer.molecule.residue.Residue` objects.
    This class extends the :class:`~pestifer.core.baseobj.AncestorAwareObjList` to manage collections of residues in a molecular structure.
    It provides methods for initializing the list from various input types, indexing residues, and performing operations
    such as mapping chain IDs, retrieving residues and atoms, and handling residue ranges.
    """
    def describe(self):
        """
        Returns a string description of the ResidueList object, including the number of residues it contains.
        
        Returns
        -------
        str
            A string representation of the ResidueList object.
        """
        return f'<ResidueList with {len(self)} residues>'
    
    @classmethod
    def from_atomlist(cls, atoms: AtomList):
        R = []
        for a in atoms:
            for r in R[::-1]:
                if r.add_atom(a):
                    break
            else:
                R.append(Residue.new(a))
        return cls(R)
    
    @classmethod
    def from_emptyresiduelist(cls, input_list: EmptyResidueList):
        """
        Create a ResidueList from an EmptyResidueList.
        
        Parameters
        ----------
        input_list : EmptyResidueList
            A list of EmptyResidue objects to convert into a ResidueList.
        
        Returns
        -------
        ResidueList
            An instance of ResidueList initialized with the residues created from the empty residues.
        """
        R = [Residue.new(m) for m in input_list]
        return cls(R)

    # @BaseObjList.__init__.register(EmptyResidueList)
    # def _from_emptyresiduelist(self, input_list: EmptyResidueList):
    #     R = [Residue.new(m) for m in input_list]
    #     super().__init__(R)

    # def index(self, R: Residue):
    #     """
    #     Get the index of a residue in the list.
        
    #     Parameters
    #     ----------
    #     R : :class:`~pestifer.molecule.residue.Residue`
    #         The residue to find in the list.
    #     """
    #     for i,r in enumerate(self):
    #         if r is R:
    #             return i
    #     else:
    #         raise ValueError(f'Residue not found')

    def map_chainIDs_label_to_auth(self):
        """
        Create a mapping from chain IDs in the label (e.g., PDB format) to the author chain IDs.
        """
        self._chainIDmap_label_to_auth = {}
        for r in self:
            if hasattr(r,'auth_asym_id'):
                label_Cid=r.chainID
                auth_Cid=r.auth_asym_id
                if not label_Cid in self._chainIDmap_label_to_auth:
                    self._chainIDmap_label_to_auth[label_Cid]=auth_Cid

    # def get_residue(self,**fields):
    #     """
    #     Get a residue from the list based on specified fields.

    #     Parameters
    #     ----------
    #     fields : keyword arguments
    #         The fields to match against residues in the list.

    #     Returns
    #     -------
    #     :class:`~pestifer.molecule.residue.Residue`
    #         The matching residue, or None if not found.
    #     """
    #     return self.get(**fields)
    
    # def get_atom(self,atname,**fields):
    #     """
    #     Get an atom from the residue's list of atoms based on its name and specified fields.
        
    #     Parameters
    #     ----------
    #     atname : str
    #         The name of the atom to retrieve.
    #     fields : keyword arguments
    #         Additional fields to match against the atom's attributes.
        
    #     Returns
    #     -------
    #     :class:`~pestifer.molecule.atom.Atom`
    #         The matching atom, or None if not found."""
    #     S=('atoms',{'name':atname})
    #     return self.get_attr(S,**fields)
    
    def atom_serials(self, as_type=str):
        """
        Get a list of atom serial numbers for all atoms in the residues.
        
        Parameters
        ----------
        as_type : type, optional
            The type to which the serial numbers should be converted. Default is str.

        Returns
        -------
        list
            A list of atom serial numbers.
        """
        serlist=[]
        for res in self:
            for a in res.atoms:
                serlist.append(as_type(a.serial))
        return serlist

    def atom_resseqnums(self, as_type=str):
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
        the N atom of a model-built residue based on the coordinates of the previous residue.  Uses the :func:`~pestifer.util.coord.positionN` function to calculate the position of the N atom based on the upstream residue's coordinates and the transformation matrix.
        
        Parameters
        ----------
        upstream_reslist : list
            A list of residues that are upstream of the current residue.
        seglabel : str
            The segment label for the residue.
        molid_varname : str
            The variable name for the molecule ID. (unused for now)
        tmat : :class:`numpy.ndarray`
            A transformation matrix to apply to the coordinates of the residue.

        Returns
        -------
        str
            A string representing the VMD/psfgen ``caco`` command for positioning the residue.
        """
        r0=self[0]
        ur=upstream_reslist[-1]
        rN=positionN(ur,tmat)
        logger.debug(f'caco {rN}')
        return f'coord {seglabel} {r0.resseqnum}{r0.insertion} N {{{rN[0]:.5f} {rN[1]:.5f} {rN[2]:.5f}}}'
    
    def resrange(self,rngrec):
        """
        Yield residues from specified residue range.
        
        Parameters
        ----------
        rngrec : :class:`argparse.Namespace`
            A record defining the range of residues to retrieve. It must have the following attributes:

            - `chainID`: The chain ID of the residues.
            - `resseqnum1`: The starting residue sequence number.
            - `insertion1`: The starting insertion code.
            - `resseqnum2`: The ending residue sequence number.
            - `insertion2`: The ending insertion code.

            For example, a :class:`~pestifer.objs.deletion.Deletion` object can be used to define the range

        Yields
        ------
        :class:`~pestifer.molecule.residue.Residue`
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
        Deletions : list of :class:`~pestifer.objs.deletion.Deletion`
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
        This method uses the :attr:`~pestifer.core.labels.Labels.segtype_of_resname` mapping to assign segment types to
        residues based on their residue names. It updates the :attr:`~pestifer.molecule.residue.Residue.segtype` attribute of each residue
        in the residue list.
        """
        # logger.debug(f'residuelist:apply_segtypes {segtype_of_resname}')
        self.map_attr('segtype','resname',Labels.segtype_of_resname)
    
    def deletion(self,DL:DeletionList):
        """
        Remove residues from the residue list based on a :class:`~pestifer.objs.deletion.DeletionList`.
        This method iterates through the DeletionList, retrieves the residues to be deleted based on their
        chain ID, residue sequence numbers, and insertion codes, and removes them from the residue list
        
        Parameters
        ----------
        DL : :class:`~pestifer.objs.deletion.DeletionList`
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
        Apply a list of substitutions to the residue list. This method iterates through the :class:`~pestifer.objs.substitution.SubstitutionList`, retrieves the residues to be substituted based on their chain ID, residue sequence numbers, and insertion codes, and replaces their residue names with the corresponding residue names from the substitution list. It also creates a new :class:`~pestifer.objs.seqadv.SeqadvList` for any resolved residues that are substituted.  Any residues that are deleted as a result of the substitutions are also returned in a list.

        Parameters
        ----------
        SL : :class:`~pestifer.objs.substitution.SubstitutionList`
            A list of substitutions to apply. Each substitution should contain the chain ID, residue sequence numbers,
            insertion codes, and the new residue sequence to substitute.
        
        Returns
        -------
        tuple : (:class:`~pestifer.objs.seqadv.SeqadvList`, list of :class:`~pestifer.molecule.residue.Residue`)
            A tuple containing a :class:`~pestifer.objs.seqadv.SeqadvList` of new sequence advancements for resolved residues and a list of residues that were deleted.
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
                            newseqadv.append(Seqadv(
                                idCode='I doubt I ever use this',
                                typekey='user',
                                resname=r.resname,
                                chainID=r.chainID,
                                resseqnum=r.resseqnum,
                                insertion=r.insertion,
                                dbRes=resname
                            ))
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
        insertions : list of :class:`~pestifer.objs.insertion.Insertion`
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
                new_empty_residue=EmptyResidue.new(shortcode)
                new_residue=Residue.new(new_empty_residue)
                new_residue.atoms=AtomList([])
                new_residue.segtype='protein'
                self.insert(idx,new_residue)
                logger.debug(f'insertion of new empty residue {shortcode}')
                idx+=1
                i='A' if i in [' ',''] else chr(ord(i)+1)

    def renumber(self, links: LinkList):
        """
        The possibility exists that empty residues added have resseqnums that conflict with existing resseqnums on the same chain if those resseqnums are in a different segtype (e.g., glycan).  This method will privilege protein residues in such conflicts, and it will renumber non-protein residues, updating any resseqnum records in links
        to match the new resseqnums.
        
        Parameters
        ----------
        links : :class:`~pestifer.objs.link.LinkList`
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
        assert newtst in [None, []] #None
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

    def state_bounds(self, state_func):
        """
        Get the state bounds for each residue in the list based on a state function.
        
        Parameters
        ----------
        state_func : function
            A function that takes a residue and returns its state bounds.
        
        Returns
        -------
        StateIntervalList
            A list of state intervals for each residue in the residue list.
        """

        return StateIntervalList.process_itemlist(self, state_func=state_func)

Residue.model_rebuild()
