# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import logging
logger=logging.getLogger(__name__)

from functools import singledispatchmethod
from pidibble.pdbrecord import PDBRecord

from ..core.baseobj import AncestorAwareObj, AncestorAwareObjList
from ..util.cifutil import CIFdict
from ..core.scriptwriters import Psfgen
from ..core.stringthings import split_ri

class SSBond(AncestorAwareObj):
    """A class for handling disulfide bonds between two CYS residues in a protein.
    Attributes
    ----------
    req_attr : list
        * chainID1 : str
            chain ID of the first CYS residue in the disulfide bond
        * resseqnum1 : int
            resid of the first CYS residue in the disulfide bond
        * insertion1 : str
            insertion code of the first CYS residue in the disulfide bond
        * chainID2 : str
            chain ID of the second CYS residue in the disulfide bond
        * resseqnum2 : int
            resid of the second CYS residue in the disulfide bond
        * insertion2 : str
            insertion code of the second CYS residue in the disulfide bond
    opt_attr : list
        * serial_number : int
            serial number of the disulfide bond in the PDB file
        * residue1 : Residue
            the first CYS residue in the disulfide bond (optional, can be None)
        * residue2 : Residue
            the second CYS residue in the disulfide bond (optional, can be None)
        * resname1 : str
            residue name of the first CYS residue in the disulfide bond (default 'CYS')
        * resname2 : str
            residue name of the second CYS residue in the disulfide bond (default 'CYS')
        * sym1 : str
            symmetry operator for the first CYS residue (optional, can be empty)
        * sym2 : str
            symmetry operator for the second CYS residue (optional, can be empty)
        * length : float
            length of the disulfide bond (optional, can be 0.0)
        * ptnr1_auth_asym_id : str
            author chain ID of the first CYS residue in the disulfide bond (optional, can be empty)
        * ptnr2_auth_asym_id : str
            author chain ID of the second CYS residue in the disulfide bond (optional, can be empty)
        * ptnr1_auth_seq_id : int
            author sequence ID of the first CYS residue in the disulfide bond (optional, can be 0)
        * ptnr2_auth_seq_id : int
            author sequence ID of the second CYS residue in the disulfide bond (optional, can be 0)
    yaml_header : str
        'ssbonds'
    objcat : str
        'topol'
    PDB_keyword : str
        'SSBOND'
    mmCIF_name : str
        'struct_conn'
    """
    req_attr=AncestorAwareObj.req_attr+['chainID1','resseqnum1','insertion1','chainID2','resseqnum2','insertion2']
    opt_attr=AncestorAwareObj.opt_attr+['serial_number','residue1','residue2','resname1','resname2','sym1','sym2','length','ptnr1_auth_asym_id','ptnr2_auth_asym_id','ptnr1_auth_seq_id','ptnr2_auth_seq_id']
    yaml_header='ssbonds'
    objcat='topol'
    PDB_keyword='SSBOND'
    mmCIF_name='struct_conn'

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
        """Write the SSBond as it would appear in a PDB file
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
    
    def write_TcL(self,W:Psfgen,transform):
        """Writes the Tcl command to create a disulfide bond in a Psfgen script using the DISU patch.
        Parameters
        ----------
        W : Psfgen
            The Psfgen script writer object to which the Tcl command will be written.
        transform : Transform
            A Transform object that contains the chainID mapping for the SSBond.
        This method generates the Tcl command to create a disulfide bond between two CYS residues
        in a protein. It uses the chainID mapping from the Transform object to ensure that the
        correct chain IDs are used in the command.
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
    """A class for handling a list of SSBond objects.
    This class inherits from AncestorAwareObjList and provides methods to manage
    a list of SSBond objects.
    """
    def assign_residues(self,Residues):
        """Assigns a list of Residue objects to the SSBond residues.
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
    def write_TcL(self,W:Psfgen,transform):
        """Writes the Tcl commands to create disulfide bonds in a Psfgen script.
        This method iterates over each SSBond in the list and writes the Tcl command to create
        the disulfide bond using the DISU patch. It uses the chainID mapping from the Transform
        object to ensure that the correct chain IDs are used in the command.
        Parameters
        ----------
        W : Psfgen
            The Psfgen script writer object to which the Tcl commands will be written.
        transform : Transform
            A Transform object that contains the chainID mapping for the SSBonds.
        """
        for s in self:
            s.write_TcL(W,transform)
    def prune_mutations(self,Mutations):
        """Prunes the SSBondList by removing SSBonds that are associated with mutations.
        Parameters
        ----------
        Mutations : list of Mutation
            A list of Mutation objects that contain information about the mutations to be pruned.
        Returns
        -------
        SSBondList
            A new SSBondList object containing the SSBonds that were pruned.
        This method iterates over each Mutation in the provided list and checks if there are
        any SSBonds that match the chain ID, residue sequence number, and insertion code of
        the Mutation. If a match is found, the SSBond is removed from the list 
        and added to a new SSBondList object. The method returns this new SSBondList object.
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
        """Maps an attribute of the SSBondList to a new value using a mapping dictionary.
        Parameters
        ----------
        mapped_attr : str
            The attribute of the SSBondList to be mapped.
        key_attr : str
            The attribute of the SSBond that will be used as the key for the mapping.
        map : dict
            A dictionary that maps the values of the key_attr to new values for the mapped_attr.
        This method iterates over each SSBond in the list and applies the mapping to the
        specified attribute. It logs the mapping operation and the number of SSBonds being processed.
        """
        logger.debug(f'Mapping {mapped_attr} {key_attr} in {len(self)} SSBonds using map {map}')
        super().map_attr(mapped_attr,key_attr,map)
