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
from functools import singledispatchmethod

from ..core.baseobj import AncestorAwareObj, AncestorAwareObjList
from ..core.scripters import VMDScripter, PsfgenScripter

class Crot(AncestorAwareObj):
    """
    A class for managing so-called "C-rotations" in a molecular structure.
    """
    req_attr=AncestorAwareObj.req_attr+['angle']
    """
    Required attributes for a Crot object.
    
    - ``angle``: The type of angle to be rotated. This can be one of the following:

        - ``PHI``: Backbone phi torsion angle.
        - ``PSI``: Backbone psi torsion angle.
        - ``OMEGA``: Backbone omega torsion angle.
        - ``CHI1``: Side-chain chi1 torsion angle.
        - ``CHI2``: Side-chain chi2 torsion angle.
        - ``ANGLEIJK``: Angle defined by atoms i, j, and k.
        - ``ALPHA``: Set of rotations that define the alpha helix.
    """

    opt_attr=AncestorAwareObj.opt_attr+['chainID','resseqnum1','resseqnum2','resseqnum3','segname','atom1','atom2','segname1','segname2','segnamei','resseqnumi','atomi','segnamejk','resseqnumj','atomj','resseqnumk','atomk','degrees']
    """
    Optional attributes for a Crot object.
    These attributes are used to specify additional parameters for the C-rotation:
    
    - ``chainID``: The chain ID of the segment to which the rotation applies.
    - ``resseqnum1``: The residue sequence number of the first residue involved in the rotation.
    - ``resseqnum2``: The residue sequence number of the second residue involved in the rotation.
    - ``resseqnum3``: The residue sequence number of the third residue involved in the rotation (optional, used for ALPHA rotations).
    - ``segname``: The segment name of the segment to which the rotation applies.
    - ``atom1``: The name of the first atom involved in the rotation.
    - ``atom2``: The name of the second atom involved in the rotation.
    - ``segname1``: The segment name of the first atom involved in the rotation.
    - ``segname2``: The segment name of the second atom involved in the rotation.
    - ``segnamei``: The segment name of the first atom in the ANGLEIJK rotation.
    - ``resseqnumi``: The residue sequence number of the first atom in the ANGLEIJK rotation.
    - ``atomi``: The name of the first atom in the ANGLEIJK rotation.
    - ``segnamejk``: The segment name of the second atom in the ANGLEIJK rotation.
    - ``resseqnumj``: The residue sequence number of the second atom in the ANGLEIJK rotation.
    - ``atomj``: The name of the second atom in the ANGLEIJK rotation.
    - ``resseqnumk``: The residue sequence number of the third atom in the ANGLEIJK rotation.
    - ``atomk``: The name of the third atom in the ANGLEIJK rotation.
    - ``degrees``: The angle in degrees by which the atoms will be rotated.
    """

    attr_choices=AncestorAwareObj.attr_choices.copy()
    """
    Attribute choices for the Crot object.
    This dictionary defines the valid choices for the ``angle`` attribute.
    The choices are:
    
    - ``PHI``: Backbone phi torsion angle.
    - ``PSI``: Backbone psi torsion angle.
    - ``OMEGA``: Backbone omega torsion angle.
    - ``CHI1``: Side-chain chi1 torsion angle.
    - ``CHI2``: Side-chain chi2 torsion angle.
    - ``ANGLEIJK``: Angle defined by atoms i, j, and k.
    - ``ALPHA``: Set of rotations that define the alpha helix.
    """
    attr_choices.update({'angle':['PHI','PSI','OMEGA','CHI1','CHI2','ANGLEIJK','ALPHA']})


    opt_attr_deps=AncestorAwareObj.opt_attr_deps.copy()
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
    opt_attr_deps.update({
        'PHI':['chainID','resseqnum1','resseqnum2'],
        'PSI':['chainID','resseqnum1','resseqnum2'],
        'OMEGA':['chainID','resseqnum1','resseqnum2'],
        'CHI1':['chainID','resseqnum1'],
        'CHI2':['chainID','resseqnum1'],
        'ANGLEIJK':['segnamei','resseqnumi','atomi','segnamejk','resseqnumj','atomj','resseqnumk','atomk'],
        'ALPHA':['chainID','resseqnum1','resseqnum2','resseqnum3']
        })

    yaml_header='crotations'
    """
    YAML header for Crot objects.
    This header is used to identify Crot objects in YAML files.
    """
    
    objcat='coord'
    """
    Category of the Crot object.
    This categorization is used to group Crot objects in the object manager. 
    """
    
    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)
    
    @__init__.register(str)
    def _from_shortcode(self,shortcode):
        dat=shortcode.split(',')
        input_dict={}
        input_dict['angle']=dat[0].upper()
        if input_dict['angle']=='PHI' or input_dict['angle']=='PSI' or input_dict['angle']=='OMEGA':
            # this is a backbone torsion, so we need both an owner
            # residue and a residue indicating the end of the 
            # set of residues that will be reoriented by the
            # rotation
            input_dict['chainID']=dat[1]
            input_dict['resseqnum1']=int(dat[2])
            input_dict['resseqnum2']=int(dat[3])
            input_dict['degrees']=float(dat[4])
        elif input_dict['angle']=='CHI1' or input_dict['angle']=='CHI2':
            input_dict['chainID']=dat[1]
            input_dict['resseqnum1']=int(dat[2])
            input_dict['resseqnum2']=-1
            input_dict['degrees']=float(dat[3])
        elif input_dict['angle']=='ANGLEIJK':
            input_dict['segnamei']=dat[1]
            input_dict['resseqnumi']=int(dat[2])
            input_dict['atomi']=dat[3]
            input_dict['segnamejk']=dat[4]
            input_dict['resseqnumj']=int(dat[5])
            input_dict['atomj']=dat[6]
            input_dict['resseqnumk']=int(dat[7])
            input_dict['atomk']=dat[8]
            input_dict['degrees']=float(dat[9])
        elif input_dict['angle']=='ALPHA':
            input_dict['chainID']=dat[1]
            input_dict['resseqnum1']=int(dat[2])
            input_dict['resseqnum2']=int(dat[3])
            if len(dat)<5:
                input_dict['resseqnum3']=input_dict['resseqnum2']
            else:
                input_dict['resseqnum3']=int(dat[4])

        super().__init__(input_dict)
    
    def to_shortcode(self):
        """
        Convert the Crot object to a shortcode representation.
        This method generates a shortcode representation of the Crot object.
        The shortcode is a string that encodes the attributes of the Crot object in a specific
        format. The format depends on the type of angle (e.g., ``PHI``, ``PSI``, ``OMEGA``, ``CHI1``, ``CHI2``, ``ANGLEIJK``, ``ALPHA``).
        The shortcode is stored in the ``shortcode`` attribute of the Crot object.
        If the ``shortcode`` attribute already exists, the method returns without modifying it.
        """
        if 'shortcode' in self.__dict__:
            return
        ret=[f'{self.angle}']
        if self.angle=='PHI' or self.angle=='PSI' or self.angle=='OMEGA':
            # this is a backbone torsion, so we need both an owner
            # residue and a residue indicating the end of the 
            # set of residues that will be reoriented by the
            # rotation
            ret.append(f'{self.chainID}')
            ret.append(f'{self.resseqnum1}')
            ret.append(f'{self.resseqnum2}')
            ret.append(f'{self.degrees:.4f}')
        elif self.angle=='CHI1' or self.angle=='CHI2':
            ret.append(f'{self.chainID}')
            ret.append(f'{self.resseqnum1}')
            ret.append(f'-1')
            ret.append(f'{self.degrees:.4f}')
        elif self.angle=='ANGLEIJK':
            ret.append(f'{self.segnamei}')
            ret.append(f'{self.resseqnumi}')
            ret.append(f'{self.atomi}')
            ret.append(f'{self.segnamejk}')
            ret.append(f'{self.resseqnumj}')
            ret.append(f'{self.atomj}')
            ret.append(f'{self.resseqnumk}')
            ret.append(f'{self.atomk}')
            ret.append(f'{self.degrees:.4f}')
        elif self.angle=='ALPHA':
            ret.append(f'{self.chainID}')
            ret.append(f'{self.resseqnum1}')
            ret.append(f'{self.resseqnum2}')
            ret.append(f'{self.resseqnum3}')
        self.shortcode=','.join(ret)
    
    def __str__(self):
        self.to_shortcode()
        return self.shortcode
    
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

class CrotList(AncestorAwareObjList):
    """
    A class for managing lists of Crot objects.
    This class inherits from AncestorAwareObjList and provides methods to handle multiple Crot objects.
    It allows for writing Tcl commands for all Crot objects in the list.
    """

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
