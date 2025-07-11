# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Disufulfide bonds are covalent linkages between two CYS residues in a protein.
These are represented in PDB files as SSBOND records, and in mmCIF files
as struct_conn records."""
import logging
logger=logging.getLogger(__name__)

from functools import singledispatchmethod
from pidibble.pdbrecord import PDBRecord

from ..core.baseobj import AncestorAwareObj, AncestorAwareObjList
from ..util.cifutil import CIFdict
from ..core.scripters import PsfgenScripter
from ..core.stringthings import split_ri

class SSBond(AncestorAwareObj):
    """
    A class for handling disulfide bonds between two CYS residues in a protein.
    """

    req_attr=AncestorAwareObj.req_attr+['chainID1','resseqnum1','insertion1','chainID2','resseqnum2','insertion2']
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
    
    opt_attr=AncestorAwareObj.opt_attr+['serial_number','residue1','residue2','resname1','resname2','sym1','sym2','length','ptnr1_auth_asym_id','ptnr2_auth_asym_id','ptnr1_auth_seq_id','ptnr2_auth_seq_id']
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
    
    yaml_header='ssbonds'
    """
    YAML header for SSBond objects.
    This header is used to identify SSBond objects in YAML files.
    """
    
    objcat='topol'
    """
    Category of the SSBond object.
    This categorization is used to group SSBond objects in the object manager.
    """
    
    PDB_keyword='SSBOND'
    """
    Keyword used to identify SSBond records in PDB files.
    """

    mmCIF_name='struct_conn'
    """
    Name of the SSBond record in mmCIF files.
    This is used to identify SSBond records in mmCIF files.
    """
    
    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)

    @__init__.register(PDBRecord)
    def _from_pdbrecord(self,pdbrecord):
        input_dict={
            'serial_number':pdbrecord.serNum,
            'resname1':pdbrecord.residue1.resName,
            'chainID1':pdbrecord.residue1.chainID,
            'resseqnum1':pdbrecord.residue1.seqNum,
            'insertion1':pdbrecord.residue1.iCode,
            'resname2':pdbrecord.residue2.resName,
            'chainID2':pdbrecord.residue2.chainID,
            'resseqnum2':pdbrecord.residue2.seqNum,
            'insertion2':pdbrecord.residue2.iCode,
            'sym1':pdbrecord.sym1,
            'sym2':pdbrecord.sym2,
            'length':pdbrecord.length,
            'residue1':None,
            'residue2':None
        }
        super().__init__(input_dict)
        logger.debug(f'parsed {self}')
    
    @__init__.register(CIFdict)
    def _from_cifdict(self,cd):
        input_dict={
            'serial_number':int(cd['id'].strip('disulf')),
            'resname1':'CYS',
            'resname2':'CYS',
            'chainID1':cd['ptnr1_label_asym_id'],
            'chainID2':cd['ptnr2_label_asym_id'],
            'resseqnum1':int(cd['ptnr1_label_seq_id']),
            'resseqnum2':int(cd['ptnr2_label_seq_id']),
            'insertion1':cd['pdbx_ptnr1_pdb_ins_code'],
            'insertion2':cd['pdbx_ptnr2_pdb_ins_code'],
            'sym1':cd['ptnr1_symmetry'],
            'sym2':cd['ptnr2_symmetry'],
            'length':float(cd['pdbx_dist_value']),
            'residue1':None,
            'residue2':None,
            'ptnr1_auth_asym_id':cd['ptnr1_auth_asym_id'],
            'ptnr2_auth_asym_id':cd['ptnr2_auth_asym_id'],
            'ptnr1_auth_seq_id':int(cd['ptnr1_auth_seq_id']),
            'ptnr2_auth_seq_id':int(cd['ptnr2_auth_seq_id'])
        }
        super().__init__(input_dict)
    
    @__init__.register(str)
    def _from_shortcode(self,shortcode):
        # shortcode format: C_RRR-D_SSS
        # C, D chainIDs
        # RRR, SSS resids
        s1=shortcode.split('-')
        r1=s1[0].split('_')
        r2=s1[1].split('_')
        input_dict={
            'serial_number':0,
            'resname1':'CYS',
            'resname2':'CYS',
            'chainID1':r1[0],
            'chainID2':r2[0],
            'resseqnum1':int(r1[1]),
            'resseqnum2':int(r2[1]),
            'insertion1':'',
            'insertion2':'',
            'length':0.0,
            'sym1':'',
            'sym2':'',
            'residue1':None,
            'residue2':None
        }
        super().__init__(input_dict)

    @__init__.register(list)
    def _from_patchlist(self,L):
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
        super().__init__(input_dict)

    def __str__(self):
        """
        Return a shortcode representation of the SSBond.
        The shortcode representation includes the chain IDs, residue sequence numbers, and insertion codes
        of both residues involved in the SSBond.
        If the residue attributes are set, it uses their chain IDs, sequence numbers, and insertion codes instead.
        """
        c1=self.chainID1
        r1=self.resseqnum1
        i1=self.insertion1
        if self.residue1:
            c1=self.residue1.chainID
            r1=self.residue1.resseqnum
            i1=self.residue1.insertion
        c2=self.chainID2
        r2=self.resseqnum2
        i2=self.insertion2
        if self.residue2:
            c2=self.residue2.chainID
            r2=self.residue2.resseqnum
            i2=self.residue2.insertion
        # return f'{self.chainID1}_{self.resseqnum1}{self.insertion1}-{self.chainID2}_{self.resseqnum2}{self.insertion2}'
        return f'{c1}_{r1}{i1}-{c2}_{r2}{i2}'

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

class SSBondList(AncestorAwareObjList):
    """
    A class for handling a list of SSBond objects.
    This class inherits from AncestorAwareObjList and provides methods to manage
    a list of SSBond objects.
    """

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
