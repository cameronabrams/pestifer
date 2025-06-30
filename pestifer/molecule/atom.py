#Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Class for handling atoms in molecular structures.
"""
import logging
logger=logging.getLogger(__name__)

from functools import singledispatchmethod

from pidibble.baserecord import BaseRecord
from pidibble.pdbrecord import PDBRecord

from ..core.baseobj import AncestorAwareObj, AncestorAwareObjList
from ..util.cifutil import CIFdict
from ..util.util import reduce_intlist

class Atom(AncestorAwareObj):
    """
    A class for handling atoms in molecular structures.
    This class represents an atom with various attributes such as serial number, name, residue name,
    chain ID, residue sequence number, insertion code, coordinates (x, y, z), occupancy, beta factor,
    element symbol, charge, and optional attributes like segment name, empty status, link status,
    record name, and author sequence ID, component ID, asym ID, and atom ID.
    """

    req_attr=AncestorAwareObj.req_attr+['serial','name','altloc','resname','chainID','resseqnum','insertion','x','y','z','occ','beta','elem','charge']
    """
    Required attributes for the Atom class.
    These attributes must be provided when creating an Atom instance.
    
    Attributes
    ----------
    serial : int
        Serial number of the atom.
    name : str 
        Name of the atom.
    altloc : str
        Alternate location identifier for the atom.
    resname : str
        Residue name to which the atom belongs.
    chainID : str
        Chain identifier for the atom.
    resseqnum : int
        Residue sequence number for the atom.
    insertion : str
        Insertion code for the residue to which the atom belongs.
    x : float
        X coordinate of the atom.
    y : float
        Y coordinate of the atom.
    z : float
        Z coordinate of the atom.
    occ : float
        Occupancy of the atom.
    beta : float
        Beta factor of the atom.
    elem : str
        Element symbol of the atom.
    charge : str
        Charge of the atom.
    """
    
    opt_attr=AncestorAwareObj.opt_attr+['segname','empty','link','recordname','auth_seq_id','auth_comp_id','auth_asym_id','auth_atom_id']
    """
    Optional attributes for the Atom class.
    These attributes can be provided when creating an Atom instance, but are not required.
    
    Attributes
    ----------
    segname : str
        Segment name to which the atom belongs. Defaults to the chain ID.
    empty : bool
        Indicates whether the atom is empty. Defaults to False.
    link : str
        Link status of the atom. Defaults to 'None'.
    recordname : str
        Record name for the atom. Defaults to 'ATOM'.
    auth_seq_id : int
        Author sequence ID for the atom.
    auth_comp_id : str
        Author component ID for the atom.
    auth_asym_id : str
        Author asym ID for the atom.
    auth_atom_id : str
        Author atom ID for the atom.
    """
    
    yaml_header='atoms'
    """
    Header for YAML serialization of Atom objects.
    """

    PDB_keyword='ATOM'
    """
    Keyword used in PDB files to identify atom records.
    """
    
    mmCIF_name='atom_site'
    """
    Name used in mmCIF files to identify atom site records.
    """

    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)

    @__init__.register(BaseRecord)
    @__init__.register(PDBRecord)
    def _from_pdbrecord(self,pdbrecord):
        input_dict={
            'serial':pdbrecord.serial,
            'name':pdbrecord.name,
            'altloc':pdbrecord.altLoc,
            'resname':pdbrecord.residue.resName,
            'chainID':pdbrecord.residue.chainID,
            'resseqnum':pdbrecord.residue.seqNum,
            'insertion':pdbrecord.residue.iCode,
            'x':pdbrecord.x,
            'y':pdbrecord.y,
            'z':pdbrecord.z,
            'occ':pdbrecord.occupancy,
            'beta':pdbrecord.tempFactor,
            'elem':pdbrecord.element,
            'charge':pdbrecord.charge
        }
        input_dict['segname']=input_dict['chainID']
        input_dict['link']='None'
        input_dict['empty']=False
        super().__init__(input_dict)

    @__init__.register(CIFdict)
    def _from_cifdict(self,cifdict):
        input_dict={
            'recordname':'ATOM',
            'serial':int(cifdict['id']),
            'name':cifdict['label_atom_id'],
            'altloc':cifdict['label_alt_id'],
            'resname':cifdict['label_comp_id'],
            'chainID':cifdict['label_asym_id'],
            'resseqnum':cifdict['label_seq_id'],
            'insertion':cifdict['pdbx_pdb_ins_code'],
            'x':float(cifdict['cartn_x']),
            'y':float(cifdict['cartn_y']),
            'z':float(cifdict['cartn_z']),
            'occ':float(cifdict['occupancy']),
            'beta':float(cifdict['b_iso_or_equiv']),
            'elem':cifdict['type_symbol'],
            'charge':cifdict['pdbx_formal_charge'],
            'auth_seq_id':cifdict['auth_seq_id'],
            'auth_comp_id':cifdict['auth_comp_id'],
            'auth_asym_id':cifdict['auth_asym_id'],
            'auth_atom_id':cifdict['auth_atom_id']
        }
        # if the seq id is a dot, we revert to the author designations for seq id and asym id
        if input_dict['resseqnum']=='.':
            # logger.debug(f'dot-resseqnum detected in {cifdict}')
            input_dict['resseqnum']=input_dict['auth_seq_id']
            # input_dict['chainID']=input_dict['auth_asym_id']
        input_dict['resseqnum']=int(input_dict['resseqnum'])
        if input_dict['auth_seq_id'].isdigit():
            input_dict['auth_seq_id']=int(input_dict['auth_seq_id'])
        input_dict['segname']=input_dict['chainID']
        input_dict['link']='None'
        input_dict['empty']=False
        super().__init__(input_dict)
    
    def pdb_line(self):
        """
        Returns a string representation of the atom in PDB format.
        This method formats the atom's attributes into a PDB line string.
        The line includes the record name, serial number, atom name, alternate location,
        residue name, chain ID, residue sequence number, insertion code, coordinates (x, y, z),
        occupancy, beta factor, element symbol, and charge.

        Returns
        -------
        str
            A formatted string representing the atom in PDB format.
        """
        pdbline='{:<6s}'.format(self.recordname)+\
                '{:5d}'.format(self.serial)+' '+\
                '{:<4s}'.format(' '+self.resname if len(self.resname)<4 else self.resname)+\
                '{:1s}'.format(self.altloc)+\
                '{:<4s}'.format(self.resname)+\
                '{:1s}'.format(self.chainID)+\
                '{:4d}'.format(self.resseqnum)+\
                '{:1s}'.format(self.insertion)+'   '+\
                '{:8.3f}'.format(self.x)+\
                '{:8.3f}'.format(self.y)+\
                '{:8.3f}'.format(self.z)+\
                '{:6.2f}'.format(self.occ)+\
                '{:6.2f}'.format(self.beta)+\
                10*' '+'{:>2s}'.format(self.elem)+'{:2s}'.format(self.charge)
        return pdbline
    
    def overwritePosition(self,other):
        """
        Overwrites the position of this atom with the position of another atom.
        This method updates the x, y, and z coordinates of this atom to match those of another atom.
        
        Parameters
        ----------
        other : Atom
            The atom whose position will be used to overwrite this atom's position.
        """
        self.x=other.x
        self.y=other.y
        self.z=other.z

class AtomList(AncestorAwareObjList):
    """
    A class for handling lists of Atom objects.
    This class inherits from AncestorAwareObjList and provides methods to manage
    a list of Atom objects, including serialization, reserialization, and position overwriting.
    """

    def reserialize(self):
        """
        Reserializes the AtomList by updating the serial numbers of each atom.
        This method assigns a new serial number to each atom in the list, starting from 1
        and incrementing for each atom. It also stores the original serial number in the
        `_ORIGINAL_` dictionary of each atom for reference.
        """
        serial=1
        for a in self:
            if not '_ORIGINAL_' in a.__dict__:
                a._ORIGINAL_={}
            a._ORIGINAL_['serial']=a.serial
            a.serial=serial
            serial+=1
            
    def adjustSerials(self,Ters):
        """
        Adjusts the serial numbers of atoms in the AtomList based on the provided TerList.
        This method reduces the serial numbers of atoms in the AtomList by the number of
        ignored serials in the TerList. It updates the `_ORIGINAL_` dictionary of each atom
        to store the original serial number before adjustment.

        Parameters
        ----------
        Ters : TerList
            A list of Ter objects containing serial numbers to be ignored (TER records in old-timey PDB files)
        """
        ignored_serials=[x.serial for x in Ters]
        if not ignored_serials:
            return
        logger.debug(f'These serials must be deleted: {ignored_serials}')
        ril=reduce_intlist([x.serial for x in self])
        logger.debug(f'Prior to ignore, serials populate {ril}')
        for a in self:
            try:
                n=next(x[0] for x in enumerate(ignored_serials) if x[1] > a.serial)
            except StopIteration:
                pass
            if n>0:
                if not '_ORIGINAL_' in a.__dict__:
                    a._ORIGINAL_={}
                a._ORIGINAL_['serial']=a.serial
                a.serial-=n
                logger.debug(f'Atom orig serial {a._ORIGINAL_["serial"]} to {a.serial}')

    def overwritePositions(self,other):
        """
        Overwrites the positions of atoms in this AtomList with the positions of atoms in another AtomList.
        This method iterates through both AtomLists, ensuring they are of equal length,
        and updates the position (x, y, z) of each atom in this AtomList to match the corresponding atom in the other AtomList.
        
        Parameters
        ----------
        other : AtomList
            The AtomList whose atom positions will be used to overwrite the positions in this AtomList.
        
        Raises
        ------
        AssertionError
            If the lengths of the two AtomLists are not equal, an assertion error is raised.
        """
        assert len(self)==len(other),'Error: atom lists not equal length'
        for sa,oa in zip(self,other):
            sa.overwritePosition(oa)
    
    def apply_psf_resnames(self,psfatoms):
        """
        Applies residue names from a PSF atom list to the corresponding atoms in this AtomList.
        This method iterates through both AtomLists, ensuring they are of equal length,
        and updates the `resname` attribute of each atom in this AtomList to match the corresponding atom in the PSF atom list.
        
        Parameters
        ----------
        psfatoms : AtomList
            The PSF atom list containing residue names to be applied.
        """
        for myatom,psfatom in zip(self,psfatoms):
            myatom.resname=psfatom.resname

class Hetatm(Atom):
    """
    A class for handling heteroatoms in molecular structures.
    This class inherits from the Atom class and represents heteroatoms with additional attributes.
    It includes the same attributes as Atom, but is specifically used for heteroatoms in PDB files.
    """

    PDB_keyword='HETATM'
    """
    Keyword used in PDB files to identify heteroatom records.
    """
    
    yaml_header='hetatoms'
    """
    Header for YAML serialization of Hetatm objects.
    """