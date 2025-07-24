# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A C-rotation is a transformation in which atoms are rotated around a given bond by a given amount.  The "C" 
designation means that only the "downstream" atoms of the bond are moved; upstream atoms, along with the atoms
of the bond itself, are not moved.  The set of upstream atoms and the set of downstream atoms must naturally
have no topological connection *other* than the bond itself.  Typically, this can be used to execute rotation
of a backbone loop in a C-terminal loop, a side-chain angle, usually in the service
of reducing steric clashes.  The primary job of this class is to translate the C-rotation shortcodes
specified by the user into TcL commands to be incorporated in a psfgen script.
"""
import logging
logger=logging.getLogger(__name__)
from typing import ClassVar
from pydantic import Field
from ..core.baseobj_new import BaseObj, BaseObjList
from ..core.scripters import VMDScripter, PsfgenScripter
from ..core.stringthings import split_ri, join_ri

class Crot(BaseObj):
    """
    A class for managing so-called "C-rotations" in a molecular structure.
    """
    _required_fields = ['angle']
    _optional_fields = ['chainID', 'resseqnum1', 'insertion1', 'resseqnum2', 'insertion2', 
                        'resseqnum3', 'insertion3', 'segname',
                        'atom1', 'atom2',
                        'segname1', 'segname2', 'segnamei', 'resseqnumi', 'insertioni','atomi',
                        'segnamejk', 'resseqnumj', 'insertionj', 'atomj', 'resseqnumk', 
                        'insertionk', 'atomk', 'degrees']
    """
    Optional attributes for a Crot object.
    These attributes are used to specify additional parameters for the C-rotation:
    
    - ``chainID``: The chain ID of the segment to which the rotation applies.
    - ``resseqnum1``: The residue sequence number of the first residue involved in the rotation.
    - ``insertion1``: The insertion code of the first residue involved in the rotation.
    - ``resseqnum2``: The residue sequence number of the second residue involved in the rotation.
    - ``insertion2``: The insertion code of the second residue involved in the rotation.
    - ``resseqnum3``: The residue sequence number of the third residue involved in the rotation (optional, used for ALPHA rotations).
    - ``segname``: The segment name of the segment to which the rotation applies.
    - ``atom1``: The name of the first atom involved in the rotation.
    - ``atom2``: The name of the second atom involved in the rotation.
    - ``segname1``: The segment name of the first atom involved in the rotation.
    - ``segname2``: The segment name of the second atom involved in the rotation.
    - ``segnamei``: The segment name of the first atom in the ANGLEIJK rotation.
    - ``resseqnumi``: The residue sequence number of the first atom in the ANGLEIJK rotation.
    - ``insertioni``: The insertion code of the first atom in the ANGLEIJK rotation.
    - ``atomi``: The name of the first atom in the ANGLEIJK rotation.
    - ``segnamejk``: The segment name of the second atom in the ANGLEIJK rotation.
    - ``resseqnumj``: The residue sequence number of the second atom in the ANGLEIJK rotation.
    - ``insertionj``: The insertion code of the second atom in the ANGLEIJK rotation.
    - ``atomj``: The name of the second atom in the ANGLEIJK rotation.
    - ``resseqnumk``: The residue sequence number of the third atom in the ANGLEIJK rotation.
    - ``insertionk``: The insertion code of the third atom in the ANGLEIJK rotation.
    - ``atomk``: The name of the third atom in the ANGLEIJK rotation.
    - ``degrees``: The angle in degrees by which the atoms will be rotated.
    """
    _attr_choices = {'angle': ['PHI', 'phi', 'PSI', 'psi', 'OMEGA', 'omega', 'CHI1', 'chi1', 'CHI2', 'chi2', 'ANGLEIJK', 'angleijk', 'ALPHA', 'alpha']}
    """
    Optional attribute dependencies for the Crot object.
    This dictionary defines the dependencies for optional attributes based on the value of the ``angle`` attribute.
    The keys are the possible values of the ``angle`` attribute, and the values are lists
    of required optional attributes for that angle type.
    The dependencies are as follows:
    
    - ``PHI``: Requires ``chainID``, ``resseqnum1``, and ``resseqnum2``.
    - ``PSI``: Requires ``chainID``, ``resseqnum1``, and ``resseqnum2``.
    - ``OMEGA``: Requires ``chainID``, ``resseqnum1``, and ``resseqnum2``.
    - ``CHI1``: Requires ``chainID`` and ``resseqnum1``.
    - ``CHI2``: Requires ``chainID`` and ``resseqnum1``.
    - ``ANGLEIJK``: Requires ``segnamei``, ``resseqnumi``, ``atomi``, ``segnamejk``, ``resseqnumj``, ``atomj``, ``resseqnumk``, and ``atomk``.
    - ``ALPHA``: Requires ``chainID``, ``resseqnum1``, ``resseqnum2``, and ``resseqnum3``.
    """
    _attr_dependencies = {'angle': {
                            'PHI':      ['chainID', 'resseqnum1', 'insertion1', 'resseqnum2', 'insertion2'],
                            'PSI':      ['chainID', 'resseqnum1', 'insertion1', 'resseqnum2', 'insertion2'],
                            'OMEGA':    ['chainID', 'resseqnum1', 'insertion1', 'resseqnum2', 'insertion2'],
                            'CHI1':     ['chainID', 'resseqnum1', 'insertion1'],
                            'CHI2':     ['chainID', 'resseqnum1', 'insertion1'],
                            'ANGLEIJK': ['segnamei', 'resseqnumi', 'atomi',
                                         'segnamejk', 'resseqnumj', 'atomj', 'resseqnumk', 'atomk'],
                            'ALPHA':    ['chainID', 'resseqnum1', 'resseqnum2', 'resseqnum3']}}
    angle: str = Field(..., description="Type of angle to be rotated (e.g., PHI, PSI, OMEGA, CHI1, CHI2, ANGLEIJK, ALPHA)")
    chainID: str = Field(None, description="Chain ID of the segment to which the rotation applies")
    resseqnum1: int = Field(None, description="Residue sequence number of the first residue involved in the rotation")
    resseqnum2: int = Field(None, description="Residue sequence number of the second residue involved in the rotation")
    resseqnum3: int = Field(None, description="Residue sequence number of the third residue involved in the rotation")
    insertion1: str = Field(None, description="Insertion code of the first residue involved in the rotation")
    insertion2: str = Field(None, description="Insertion code of the second residue involved in the rotation")
    insertion3: str = Field(None, description="Insertion code of the third residue involved in the rotation")
    segname: str = Field(None, description="Segment name of the segment to which the rotation applies")
    atom1: str = Field(None, description="Name of the first atom involved in the rotation")
    atom2: str = Field(None, description="Name of the second atom involved in the rotation")
    segname1: str = Field(None, description="Segment name of the first atom involved in the rotation")
    segname2: str = Field(None, description="Segment name of the second atom involved in the rotation")
    segnamei: str = Field(None, description="Segment name of the first atom involved in the rotation")
    resseqnumi: int = Field(None, description="Residue sequence number of the first atom in the ANGLEIJK rotation")
    insertioni: str = Field(None, description="Insertion code of the first atom in the ANGLEIJK rotation")
    atomi: str = Field(None, description="Name of the first atom in the ANGLEIJK rotation")
    segnamejk: str = Field(None, description="Segment name of the second atom in the ANGLEIJK rotation")
    resseqnumj: int = Field(None, description="Residue sequence number of the second atom in the ANGLEIJK rotation")
    insertionj: str = Field(None, description="Insertion code of the second atom in the ANGLEIJK rotation")
    atomj: str = Field(None, description="Name of the second atom in the ANGLEIJK rotation")
    resseqnumk: int = Field(None, description="Residue sequence number of the third atom in the ANGLEIJK rotation")
    insertionk: str = Field(None, description="Insertion code of the third atom in the ANGLEIJK rotation")
    atomk: str = Field(None, description="Name of the third atom in the ANGLEIJK rotation")
    degrees: float = Field(0.0, description="Angle in degrees by which the atoms will be rotated")

    _yaml_header: ClassVar[str] = 'crotations'
    """
    YAML header for Crot objects.
    This header is used to identify Crot objects in YAML files.
    """

    _objcat: ClassVar[str] = 'coord'
    """
    Category of the Crot object.
    This categorization is used to group Crot objects in the object manager. 
    """
    
    class Adapter:
        """
        A class to represent the shortcode format for Crot, so that we can register to BaseObj.from_input rather than defining a local from_input.

        The shortcode format is:
        - For PHI, PSI, OMEGA: angle,chainID,resseqnum1,insertion1,resseqnum2,insertion2,degrees
        - For CHI1, CHI2: angle,chainID,resseqnum1,insertion1,degrees
        - For ANGLEIJK: angle,segnamei,resseqnumi,insertioni,atomi,segnamejk,resseqnumj,insertionj,atomj,resseqnumk,insertionk,atomk,degrees
        - For ALPHA: angle,chainID,resseqnum1,insertion1,resseqnum2,insertion2,[resseqnum3,insertion3]
        """
        def __init__(self, angle: str, chainID: str = None, resseqnum1: int = None, insertion1: str = None,
                     resseqnum2: int = None, insertion2: str = None, resseqnum3: int = None, insertion3: str = None,
                     segname: str = None, atom1: str = None, atom2: str = None,
                     segname1: str = None, segname2: str = None,
                     segnamei: str = None, resseqnumi: int = None, insertioni: str = None, atomi: str = None,
                     segnamejk: str = None, resseqnumj: int = None, insertionj: str = None, atomj: str = None,
                     resseqnumk: int = None, insertionk: str = None, atomk: str = None,
                     degrees: float = 0.0):
            self.angle = angle.upper()
            self.chainID = chainID
            self.resseqnum1 = resseqnum1
            self.insertion1 = insertion1
            self.resseqnum2 = resseqnum2
            self.insertion2 = insertion2
            self.resseqnum3 = resseqnum3
            self.insertion3 = insertion3
            self.segname = segname
            self.atom1 = atom1
            self.atom2 = atom2
            self.segname1 = segname1
            self.segname2 = segname2
            self.segnamei = segnamei
            self.resseqnumi = resseqnumi
            self.insertioni = insertioni
            self.atomi = atomi
            self.segnamejk = segnamejk
            self.resseqnumj = resseqnumj
            self.insertionj = insertionj
            self.atomj = atomj
            self.resseqnumk = resseqnumk
            self.insertionk = insertionk
            self.atomk = atomk
            self.degrees=degrees
        
        @classmethod
        def from_string(cls, raw: str):
            dat=raw.split(',')
            angle=dat[0].upper()
            match angle:
                case 'PHI' | 'phi' | 'PSI' | 'psi' | 'OMEGA' | 'omega':
                    chainID=dat[1]
                    resseqnum1,insertion1=split_ri(dat[2])
                    resseqnum2,insertion2=split_ri(dat[3])
                    degrees=float(dat[4])
                    return cls(angle, chainID=chainID, resseqnum1=resseqnum1, insertion1=insertion1, resseqnum2=resseqnum2, insertion2=insertion2, degrees=degrees)
                case 'CHI1' | 'chi1' | 'CHI2' | 'chi2':
                    chainID=dat[1]
                    resseqnum1,insertion1=split_ri(dat[2])
                    degrees=float(dat[3])
                    return cls(angle, chainID=chainID, resseqnum1=resseqnum1, insertion1=insertion1, degrees=degrees)
                case 'ANGLEIJK' | 'angleijk':
                    segnamei=dat[1]
                    resseqnumi,insertioni=split_ri(dat[2])
                    atomi=dat[3]
                    segnamejk=dat[4]
                    resseqnumj,insertionj=split_ri(dat[5])
                    atomj=dat[6]
                    resseqnumk,insertionk=split_ri(dat[7])
                    atomk=dat[8]
                    degrees=float(dat[9])
                    return cls(angle, segnamei=segnamei, resseqnumi=resseqnumi, insertioni=insertioni, atomi=atomi,
                               segnamejk=segnamejk, resseqnumj=resseqnumj, insertionj=insertionj, atomj=atomj,
                               resseqnumk=resseqnumk, insertionk=insertionk, atomk=atomk, degrees=degrees)
                case 'ALPHA' | 'alpha':
                    chainID=dat[1]
                    resseqnum1,insertion1=split_ri(dat[2])
                    resseqnum2,insertion2=split_ri(dat[3])
                    if len(dat) < 5:
                        resseqnum3 = resseqnum2
                        insertion3 = insertion2
                    else:
                        resseqnum3,insertion3 = split_ri(dat[4])
                    return cls(angle, chainID=chainID, resseqnum1=resseqnum1, insertion1=insertion1,
                               resseqnum2=resseqnum2, insertion2=insertion2, resseqnum3=resseqnum3,
                               insertion3=insertion3)
                case _:
                    raise ValueError(f"Unknown angle type: {angle}")

        def to_string(self) -> str:
            match self.angle:
                case 'PHI' | 'phi' | 'PSI' | 'psi' | 'OMEGA' | 'omega':
                    return f'{self.angle.upper()},{self.chainID},{join_ri(self.resseqnum1,self.insertion1)},{join_ri(self.resseqnum2,self.insertion2)},{self.degrees:.4f}'
                case 'CHI1' | 'chi1' | 'CHI2' | 'chi2':
                    return f'{self.angle.upper()},{self.chainID},{join_ri(self.resseqnum1,self.insertion1)},{self.degrees:.4f}'
                case 'ANGLEIJK' | 'angleijk':
                    return f'{self.angle},{self.segnamei},{join_ri(self.resseqnumi,self.insertioni)},{self.atomi},{self.segnamejk},{join_ri(self.resseqnumj,self.insertionj)},{self.atomj},{join_ri(self.resseqnumk,self.insertionk)},{self.atomk},{self.degrees:.4f}'
                case 'ALPHA' | 'alpha':
                    return f'{self.angle},{self.chainID},{join_ri(self.resseqnum1,self.insertion1)},{join_ri(self.resseqnum2,self.insertion2)},{join_ri(self.resseqnum3,self.insertion3)}'

        def to_dict(self) -> dict:
            """
            Converts the Crot object to a dictionary representation.  Only attributes that are not None are included.
            
            Returns
            -------
            dict
                A dictionary representation of the Crot object.
            """
            return {k: v for k, v in self.__dict__.items() if v is not None}

    @BaseObj.from_input.register(Adapter)
    @classmethod
    def _from_shortcode(cls, shortcode: Adapter):
        """
        Converts a Crot.Adapter object to a Crot object.
        
        Parameters
        ----------
        shortcode : Adapter
            The Adapter object containing the C-rotation parameters.
        
        Returns
        -------
        Crot
            A Crot object initialized with the parameters from the Adapter.
        """
        return cls(**(shortcode.to_dict()))
    
    def to_input_string(self) -> str:
        """
        Converts the CleavageSite object to a string representation for input.

        Returns
        -------
        str
            A string representation of the CleavageSite object in the format:
            C:R1-R2
        """
        return self.Adapter(**(self.model_dump())).to_string()
    
    def to_input_string(self) -> str:
        """
        Converts the Crot object to a string representation for input.
        
        Returns
        -------
        str
            A string representation of the Crot object in the format:
            angle,chainID,resseqnum1,insertion1,resseqnum2,insertion2,degrees
        """
        return self.Adapter(**(self.model_dump())).to_string()

    def describe(self):
        return f'Crot({self.to_input_string()})'

    def write_TcL(self,W:VMDScripter,chainIDmap={},**kwargs):
        """
        Write the Tcl commands to perform the C-rotation in a VMD script.
        
        Parameters
        ----------
        W : VMDScripter
            The VMD script writer object to which the Tcl commands will be written.
        chainIDmap : dict, optional
            A dictionary mapping original chain IDs to new chain IDs. If not provided, the original chain ID will be used.
        **kwargs : dict, optional
            Additional keyword arguments that may be used in the future.
        """
        the_chainID=chainIDmap.get(self.chainID,self.chainID)
        molid_varname=W.molid_varname
        molid=f'${molid_varname}'
        # endIsCterm=kwargs.get('endIsCterm',True)
        if self.angle in ['PHI','PSI','OMEGA']:
            W.addline('set r1 [[atomselect {} "chain {} and resid {} and name CA"] get residue]'.format(molid,the_chainID,self.resseqnum1))
            W.addline('set r2 [[atomselect {} "chain {} and resid {} and name CA"] get residue]'.format(molid,the_chainID,self.resseqnum2))
            if self.resseqnum1<=self.resseqnum2:
                direction='C'
            else:
                direction='N'
            W.addline(f'brot {molid} $r1 $r2 {self.angle.lower()} {direction} {self.degrees}')
            # if endIsCterm:
            #     W.addline('Crot_{}_toCterm $r1 $r2 {} {} {}'.format(self.angle.lower(),the_chainID,molid,self.degrees))
            # else:
            #     W.addline('Crot_{} $r1 $r2 {} {} {}'.format(self.angle.lower(),the_chainID,molid,self.degrees))
        elif self.angle in ['CHI1','CHI2']:  # this is a side-chain bond
            W.addline('set r1 [[atomselect {} "chain {} and resid {} and name CA"] get residue]'.format(molid,the_chainID,self.resseqnum1))
            W.addline(f'brot {molid} $r1 -1 {self.angle[:-1].lower()} {self.angle[-1]} {self.degrees}')
        elif self.angle=='ANGLEIJK':
            W.addline('set rotsel [atomselect {} "segname {}"]'.format(molid,self.segnamejk))
            W.addline('set ri [lindex [[atomselect {} "segname {} and resid {} and name {}"] get {{x y z}}] 0]'.format(molid,self.segnamei,self.resseqnumi,self.atomi))
            W.addline('set rj [lindex [[atomselect {} "segname {} and resid {} and name {}"] get {{x y z}}] 0]'.format(molid,self.segnamejk,self.resseqnumj,self.atomj))
            W.addline('set rk [lindex [[atomselect {} "segname {} and resid {} and name {}"] get {{x y z}}] 0]'.format(molid,self.segnamejk,self.resseqnumk,self.atomk))
            W.addline('set rij [vecsub $ri $rj]')
            W.addline('set rjk [vecsub $rj $rk]')
            W.addline('set cijk [veccross $rij $rjk]')
            W.addline('$rotsel move [trans center $rj origin $rj axis $cijk {} degrees]'.format(self.degrees))
        elif self.angle=='ALPHA':
            W.addline('set r1 [[atomselect {} "chain {} and resid {} and name CA"] get residue]'.format(molid,the_chainID,self.resseqnum1))
            W.addline('set r2 [[atomselect {} "chain {} and resid {} and name CA"] get residue]'.format(molid,the_chainID,self.resseqnum2))
            W.addline('set rterm [[atomselect {} "chain {} and resid {} and name CA"] get residue]'.format(molid,the_chainID,self.resseqnum3))
            W.addline('fold_alpha $r1 $r2 $rterm {}'.format(molid))

class CrotList(BaseObjList[Crot]):
    """
    A class for managing lists of Crot objects.
    This class inherits from BaseObjList and provides methods to handle multiple Crot objects.
    It allows for writing Tcl commands for all Crot objects in the list.
    """

    def describe(self):
        """
        Returns a string description of the CrotList, including the number of Crot objects it contains.
        
        Returns
        -------
        str
            A string describing the CrotList.
        """
        return f"CrotList with {len(self)} Crot objects"
    
    def _validate_item(self, item: Crot) -> None:
        """
        Validate that the item is an instance of Crot.
        
        Parameters
        ----------
        item : Crot
            The item to validate.
        
        Raises
        ------
        TypeError
            If the item is not an instance of Crot.
        """
        if not isinstance(item, Crot):
            raise TypeError(f"Item must be an instance of Crot, got {type(item)}")

    def write_TcL(self,W:PsfgenScripter,chainIDmap={},**kwargs):
        """
        Write the Tcl commands for all Crot objects in the list.
        
        Parameters
        ----------
        W : PsfgenScripter
            The Psfgen script writer object to which the Tcl commands will be written.
        chainIDmap : dict, optional
            A dictionary mapping original chain IDs to new chain IDs. If not provided, the original chain ID will be used.
        **kwargs : dict, optional
            Additional keyword arguments that may be used in the future.
        """
        for c in self:
            c.write_TcL(W,chainIDmap=chainIDmap,**kwargs)
