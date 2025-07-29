# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Disufulfide bonds are covalent linkages between two CYS residues in a protein.
These are represented in PDB files as SSBOND records, and in mmCIF files
as struct_conn records.
"""
import logging
logger=logging.getLogger(__name__)

from functools import singledispatchmethod
from pidibble.pdbrecord import PDBRecord, PDBRecordDict

from pydantic import Field
from typing import ClassVar, Optional, Any, Dict, Set

from ..core.baseobj_new import BaseObj, BaseObjList
from ..util.cifutil import CIFdict
from ..core.scripters import PsfgenScripter
from ..core.stringthings import split_ri, join_ri

from mmcif.api.PdbxContainers import DataContainer

class SSBond(BaseObj):
    """
    A class for handling disulfide bonds between two CYS residues in a protein.
    """

    _required_fields = {'chainID1', 'resseqnum1', 'insertion1', 'chainID2', 'resseqnum2', 'insertion2'}
    """
    Required attributes for an SSBond object.
    These attributes must be provided when creating an SSBond object.
    
    - ``chainID1``: The chain ID of the first CYS residue.
    - ``resseqnum1``: The residue sequence number of the first CYS residue.
    - ``insertion1``: The insertion code of the first CYS residue.
    - ``chainID2``: The chain ID of the second CYS residue.
    - ``resseqnum2``: The residue sequence number of the second CYS residue.
    - ``insertion2``: The insertion code of the second CYS residue.
    """

    _optional_fields = {'serial_number', 'residue1', 'residue2', 'resname1', 'resname2', 'sym1', 'sym2', 'length', 'ptnr1_auth_asym_id', 'ptnr2_auth_asym_id', 'ptnr1_auth_seq_id', 'ptnr2_auth_seq_id'}
    """
    Optional attributes for an SSBond object.
    These attributes are not required but can be provided to enhance the SSBond object.
    
    - ``serial_number``: The serial number of the SSBond.
    - ``residue1``: The first CYS residue object associated with the SSBond.
    - ``residue2``: The second CYS residue object associated with the SSBond.
    - ``resname1``: The residue name of the first CYS residue.
    - ``resname2``: The residue name of the second CYS residue.
    - ``sym1``: The symmetry operator for the first CYS residue.
    - ``sym2``: The symmetry operator for the second CYS residue.
    - ``length``: The length of the SSBond.
    - ``ptnr1_auth_asym_id``: The author asym ID of the first CYS residue in mmCIF format.
    - ``ptnr2_auth_asym_id``: The author asym ID of the second CYS residue in mmCIF format.
    - ``ptnr1_auth_seq_id``: The author sequence ID of the first CYS residue in mmCIF format.
    - ``ptnr2_auth_seq_id``: The author sequence ID of the second CYS residue in mmCIF format.
    """

    chainID1: str = Field(..., description="Chain ID of the first CYS residue")
    resseqnum1: int = Field(..., description="Residue sequence number of the first CYS residue")
    insertion1: str = Field(..., description="Insertion code of the first CYS residue")
    chainID2: str = Field(..., description="Chain ID of the second CYS residue")
    resseqnum2: int = Field(..., description="Residue sequence number of the second CYS residue")
    insertion2: str = Field(..., description="Insertion code of the second CYS residue")

    serial_number: int = Field(0, description="Serial number of the SSBond")
    residue1: Optional[Any] = Field(None, description="First CYS residue object associated with the SSBond")
    residue2: Optional[Any] = Field(None, description="Second CYS residue object associated with the SSBond")
    resname1: str = Field('CYS', description="Residue name of the first CYS residue")
    resname2: str = Field('CYS', description="Residue name of the second CYS residue")
    sym1: str = Field('', description="Symmetry operator for the first CYS residue")
    sym2: str = Field('', description="Symmetry operator for the second CYS residue")
    length: float = Field(0.0, description="Length of the SSBond")
    ptnr1_auth_asym_id: str = Field('', description="Author asym ID of the first CYS residue in mmCIF format")
    ptnr2_auth_asym_id: str = Field('', description="Author asym ID of the second CYS residue in mmCIF format")
    ptnr1_auth_seq_id: str = Field('', description="Author sequence ID of the first CYS residue in mmCIF format")
    ptnr2_auth_seq_id: str = Field('', description="Author sequence ID of the second CYS residue in mmCIF format")

    _yaml_header: ClassVar[str] = 'ssbonds'
    """
    YAML header for SSBond objects.
    This header is used to identify SSBond objects in YAML files.
    """

    _objcat: ClassVar[str] = 'topol'
    """
    Category of the SSBond object.
    This categorization is used to group SSBond objects in the object manager.
    """

    _PDB_keyword: ClassVar[str] = 'SSBOND'
    """
    Keyword used to identify SSBond records in PDB files.
    """

    _CIF_CategoryName: ClassVar[str] = 'struct_conn'
    """
    Name of the CIF category that contains information for SSBond objects.
    """
    _CIF_CategoryElementTypes: ClassVar[Dict[str, Set]] = {'conn_type_id': {'disulf'}}
    """
    CIF category element types for SSBond objects; a CIF Category is effectively a list of dictionaries, and _CIF_CategoryElementTypes[keyname] is a set of values that are valid for that key, indicating the element is to be interpreted as a SSBond.
    """

    def describe(self):
        """
        Describe the SSBond object.
        
        Returns
        -------
        str
            A string description of the SSBond object, including chain IDs, residue sequence numbers, and insertion codes.
        """
        return f"SSBond(chainID1={self.chainID1}, resseqnum1={self.resseqnum1}, insertion1={self.insertion1}, chainID2={self.chainID2}, resseqnum2={self.resseqnum2}, insertion2={self.insertion2})"
    
    class Adapter:
        def __init__(self, chainID1: str, resseqnum1: int, insertion1: str, chainID2: str, resseqnum2: int, insertion2: str, serial_number: int = 0, residue1: Any = None, residue2: Any = None, resname1: str = 'CYS', resname2: str = 'CYS', sym1: str = '', sym2: str = '', length: float = 0.0, ptnr1_auth_asym_id: str = '', ptnr2_auth_asym_id: str = '', ptnr1_auth_seq_id: str = '', ptnr2_auth_seq_id: str = ''):
            self.chainID1 = chainID1
            self.resseqnum1 = resseqnum1
            self.insertion1 = insertion1
            self.chainID2 = chainID2
            self.resseqnum2 = resseqnum2
            self.insertion2 = insertion2
            self.serial_number = serial_number
            self.residue1 = residue1
            self.residue2 = residue2
            self.resname1 = resname1
            self.resname2 = resname2
            self.sym1 = sym1
            self.sym2 = sym2
            self.length = length
            self.ptnr1_auth_asym_id = ptnr1_auth_asym_id
            self.ptnr2_auth_asym_id = ptnr2_auth_asym_id
            self.ptnr1_auth_seq_id = ptnr1_auth_seq_id
            self.ptnr2_auth_seq_id = ptnr2_auth_seq_id

        @classmethod
        def from_string(cls, raw: str):
            # shortcode format: C_RRR-D_SSS
            # C, D chainIDs
            # RRR, SSS resids
            s1, s2 = raw.split('-')
            r1, i1 = split_ri(s1.split('_')[1])
            r2, i2 = split_ri(s2.split('_')[1])
            return cls(chainID1=s1.split('_')[0], resseqnum1=r1, insertion1=i1, chainID2=s2.split('_')[0], resseqnum2=r2, insertion2=i2)
        
        @classmethod
        def from_pdbrecord(cls, pdbrecord: PDBRecord):
            """
            Create an Adapter instance from a PDBRecord.
            
            Parameters
            ----------
            pdbrecord : PDBRecord
                The PDBRecord object containing the SSBond information.
            
            Returns
            -------
            Adapter
                A new Adapter instance with the SSBond information from the PDBRecord.
            """
            return cls(
                chainID1=pdbrecord.residue1.chainID,
                resseqnum1=pdbrecord.residue1.seqNum,
                insertion1=pdbrecord.residue1.iCode,
                chainID2=pdbrecord.residue2.chainID,
                resseqnum2=pdbrecord.residue2.seqNum,
                insertion2=pdbrecord.residue2.iCode,
                serial_number=pdbrecord.serNum,
                resname1='CYS',
                resname2='CYS',
                sym1=pdbrecord.sym1,
                sym2=pdbrecord.sym2,
                length=pdbrecord.length
            )

        @classmethod
        def from_cifdict(cls, cd: CIFdict):
            """
            Create an Adapter instance from a CIFdict.
            
            Parameters
            ----------
            cd : CIFdict
                The CIFdict object containing the SSBond information.
            
            Returns
            -------
            Adapter
                A new Adapter instance with the SSBond information from the CIFdict.
            """
            return cls(
                chainID1=cd['ptnr1_label_asym_id'],
                resseqnum1=int(cd['ptnr1_label_seq_id']),
                insertion1=cd['pdbx_ptnr1_pdb_ins_code'],
                chainID2=cd['ptnr2_label_asym_id'],
                resseqnum2=int(cd['ptnr2_label_seq_id']),
                insertion2=cd['pdbx_ptnr2_pdb_ins_code'],
                serial_number=int(cd['id'].strip('disulf')),
                resname1='CYS',
                resname2='CYS',
                sym1=cd['ptnr1_symmetry'],
                sym2=cd['ptnr2_symmetry'],
                length=float(cd['pdbx_dist_value'])
            )
    
        @classmethod
        def from_patchlist(cls, L: list):
            s1,ri1=L[0].split(':')
            s2,ri2=L[1].split(':')
            r1,i1=split_ri(ri1)
            r2,i2=split_ri(ri2)
            input_dict={
                'serial_number':0,
                'resname1':'CYS',
                'resname2':'CYS',
                'chainID1':s1,
                'chainID2':s2,
                'resseqnum1':r1,
                'resseqnum2':r2,
                'insertion1':i1,
                'insertion2':i2,
                'length':0.0,
                'sym1':'',
                'sym2':'',
                'residue1':None,
                'residue2':None
            }
            return cls(**input_dict)

        def to_string(self) -> str:
            """
            Convert the SSBond object to a shortcode string representation.
            
            Returns
            -------
            str
                A shortcode string representation of the SSBond object in the format C_RRR-D_SSS.
                Where C and D are chain IDs, and RRR and SSS are residue sequence numbers with insertion codes.
            """
            return f"{self.chainID1}_{join_ri(self.resseqnum1, self.insertion1)}-{self.chainID2}_{join_ri(self.resseqnum2, self.insertion2)}"
    
        def to_dict(self) -> dict:
            """
            Convert the SSBond object to a dictionary representation.
            
            Returns
            -------
            dict
                A dictionary representation of the SSBond object, including all attributes.
            """
            return {k: v for k, v in self.__dict__.items() if not v is None}
        
    @BaseObj.from_input.register(Adapter)
    @classmethod
    def _from_adapter(cls, adapter: Adapter):
        return cls(**adapter.to_dict())

    @singledispatchmethod
    @classmethod
    def new(cls, raw: Any) -> "SSBond":
        raise TypeError(f"Cannot create SSBond from {type(raw)}: {raw}")
    
    @new.register(str)
    @classmethod
    def _from_string(cls, raw: str) -> "SSBond":
        """
        Create a new SSBond instance from a shortcode string.
        
        Parameters
        ----------
        raw : str
            The shortcode string in the format C_RRR-D_SSS.
        
        Returns
        -------
        SSBond
            A new instance of SSBond.
        """
        adapter = cls.Adapter.from_string(raw)
        instance = cls._from_adapter(adapter)
        return instance
    
    @new.register(PDBRecord)
    @classmethod
    def _from_pdbrecord(cls, pdbrecord: PDBRecord) -> "SSBond":
        """
        Create a new SSBond instance from a PDBRecord.
        
        Parameters
        ----------
        pdbrecord : PDBRecord
            The PDBRecord object containing the SSBond information.
        
        Returns
        -------
        SSBond
            A new instance of SSBond.
        """
        adapter = cls.Adapter.from_pdbrecord(pdbrecord)
        instance = cls._from_adapter(adapter)
        return instance
    
    @new.register(CIFdict)
    @classmethod
    def _from_cifdict(cls, cd: CIFdict) -> "SSBond":
        """
        Create a new SSBond instance from a CIFdict.
        
        Parameters
        ----------
        cd : CIFdict
            The CIFdict object containing the SSBond information.
        
        Returns
        -------
        SSBond
            A new instance of SSBond.
        """
        adapter = cls.Adapter.from_cifdict(cd)
        instance = cls._from_adapter(adapter)
        return instance
    
    @new.register(list)
    @classmethod
    def _from_patchlist(cls, L: list) -> "SSBond":
        """
        Create a new SSBond instance from a patch list.
        
        Parameters
        ----------
        L : list
            A list containing two elements, each representing a residue in the format 'C:nnn' or 'D:nnn'.
        
        Returns
        -------
        SSBond
            A new instance of SSBond.
        """
        adapter = cls.Adapter.from_patchlist(L)
        instance = cls._from_adapter(adapter)
        return instance

    def to_input_string(self) -> str:
        """
        Convert the SSBond object to a string representation for input.
        
        Returns
        -------
        str
            A string representation of the SSBond object in the format C_RRR-D_SSS.
            Where C and D are chain IDs, and RRR and SSS are residue sequence numbers with insertion codes.
        """
        return self.Adapter(**(self.model_dump())).to_string()

    def __str__(self):
        """
        Return a shortcode representation of the SSBond.
        The shortcode representation includes the chain IDs, residue sequence numbers, and insertion codes
        of both residues involved in the SSBond.
        If the residue attributes are set, it uses their chain IDs, sequence numbers, and insertion codes instead.
        """
        return self.to_input_string()

    def pdb_line(self):
        """
        Write the SSBond as it would appear in a PDB file.

        Returns
        -------
        str
            The PDB line for the SSBOND record.
        """
        pdbline='SSBOND'+\
                '{:4d}'.format(self.serial_number)+\
                '{:>4s}'.format(self.resname1)+\
                '{:>2s}'.format(self.chainID1)+\
                '{:5d}'.format(self.resseqnum1)+\
                '{:1s}'.format(self.insertion1)+\
                '{:>6s}'.format(self.resname2)+\
                '{:>2s}'.format(self.chainID2)+\
                '{:5d}'.format(self.resseqnum2)+\
                '{:1s}'.format(self.insertion2)+\
                ' '*23+\
                '{:>6s}'.format(self.sym1)+\
                '{:>7s}'.format(self.sym2)+\
                '{:6.2f}'.format(self.length)
        return pdbline
    
    def write_TcL(self,W:PsfgenScripter,transform):
        """
        Writes the Tcl command to create a disulfide bond in a Psfgen script using the DISU patch. This method generates the Tcl command to create a disulfide bond between two CYS residues
        in a protein. It uses the chainID mapping from the Transform object to ensure that the
        correct chain IDs are used in the command.
        
        Parameters
        ----------
        W : PsfgenScripter
            The Psfgen script writer object to which the Tcl command will be written.
        transform : Transform
            A Transform object that contains the chainID mapping for the SSBond.
        """
        chainIDmap=transform.chainIDmap 
        # ok since these are only going to reference protein segments; protein segment names are the chain IDs
        logger.debug(f'writing patch for {self}')
        c1=chainIDmap.get(self.chainID1,self.chainID1)
        c2=chainIDmap.get(self.chainID2,self.chainID2)
        r1=self.resseqnum1
        r2=self.resseqnum2
        W.addline(f'patch DISU {c1}:{r1} {c2}:{r2}')

class SSBondList(BaseObjList[SSBond]):
    """
    A class for handling a list of SSBond objects.
    This class inherits from BaseObjList and provides methods to manage
    a list of SSBond objects.
    """

    def describe(self):
        return f'SSBondList with {len(self)} SSBonds'

    @classmethod
    def from_pdb(cls, pdb: PDBRecordDict) -> 'SSBondList':
        """
        Create a SSBondList from a PDBRecordDict.
        """
        if SSBond._PDB_keyword not in pdb:
            return cls([])
        return cls([SSBond.new(x) for x in pdb[SSBond._PDB_keyword]])

    @classmethod
    def from_cif(cls, dc: DataContainer) -> 'SSBondList':
        """
        Create a SSBondList from a CIF DataContainer.

        Parameters
        ----------
        dc : DataContainer
            A CIF DataContainer containing the necessary fields to create SSBond objects.

        Returns
        -------
        SSBondList
            An instance of SSBondList created from the CIF DataContainer.
        """
        L = []
        cif_category = dc.getObj(SSBond._CIF_CategoryName)
        if cif_category is None:
            return cls([])
        for i in range(len(cif_category)):
            for key, valset in SSBond._CIF_CategoryElementTypes.items():
                objTypeid = cif_category.getValue(key, i)
                if objTypeid in valset:
                    this_link = SSBond.new(CIFdict(cif_category, i))
                    L.append(this_link)
        return cls(L)

    def assign_residues(self,Residues):
        """
        Assigns a list of Residue objects to the SSBond residues.
        
        Parameters
        ----------
        Residues : list of Residue
            A list of Residue objects to be assigned to the SSBond residues.
        
        Returns
        -------
        SSBondList
            A new SSBondList object containing the residues that were not assigned because no applicable residues were found.
        """
        logger.debug(f'SSBonds: Assigning residues from list of {len(Residues)} residues')
        ignored_by_ptnr1=self.assign_objs_to_attr('residue1',Residues,resseqnum='resseqnum1',chainID='chainID1',insertion='insertion1')
        ignored_by_ptnr2=self.assign_objs_to_attr('residue2',Residues,resseqnum='resseqnum2',chainID='chainID2',insertion='insertion2')
        return self.__class__(ignored_by_ptnr1+ignored_by_ptnr2)
    
    def write_TcL(self,W:PsfgenScripter,transform):
        """
        Writes the Tcl commands to create disulfide bonds in a Psfgen script.
        This method iterates over each SSBond in the list and writes the Tcl command to create
        the disulfide bond using the DISU patch. It uses the chainID mapping from the Transform
        object to ensure that the correct chain IDs are used in the command.
        
        Parameters
        ----------
        W : PsfgenScripter
            The Psfgen script writer object to which the Tcl commands will be written.
        transform : Transform
            A Transform object that contains the chainID mapping for the SSBonds.
        """
        for s in self:
            s.write_TcL(W,transform)
    def prune_mutations(self,Mutations):
        """
        Prunes the SSBondList by removing SSBonds that are associated with mutations. This method iterates over each Mutation in the provided list and checks if there are
        any SSBonds that match the chain ID, residue sequence number, and insertion code of
        the Mutation. If a match is found, the SSBond is removed from the list 
        and added to a new SSBondList object. The method returns this new SSBondList object.
        
        Parameters
        ----------
        Mutations : list of Mutation
            A list of Mutation objects that contain information about the mutations to be pruned.
        
        Returns
        -------
        SSBondList
            A new SSBondList object containing the SSBonds that were pruned.

        """
        pruned=self.__class__([])
        for m in Mutations:
            left=self.get(chainID1=m.chainID,resseqnum1=m.resseqnum,insertion1=m.insertion)
            if left:
                pruned.append(self.remove(left))
            right=self.get(chainID2=m.chainID,resseqnum2=m.resseqnum,insertion2=m.insertion)
            if right:
                pruned.append(self.remove(right))
        return pruned
    
    def map_attr(self,mapped_attr,key_attr,map):
        """
        Maps an attribute of the SSBondList to a new value using a mapping dictionary.  A dictionary that maps the values of the key_attr to new values for the mapped_attr.
        This method iterates over each SSBond in the list and applies the mapping to the
        specified attribute. It logs the mapping operation and the number of SSBonds being processed.
        
        Parameters
        ----------
        mapped_attr : str
            The attribute of the SSBondList to be mapped.
        key_attr : str
            The attribute of the SSBond that will be used as the key for the mapping.
        map : dict
        """
        logger.debug(f'Mapping {mapped_attr} {key_attr} in {len(self)} SSBonds using map {map}')
        super().map_attr(mapped_attr,key_attr,map)
