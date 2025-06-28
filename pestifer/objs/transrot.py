# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import logging
logger=logging.getLogger(__name__)
from functools import singledispatchmethod

from ..core.baseobj import AncestorAwareObj, AncestorAwareObjList
from ..core.scriptwriters import VMD

class RotTrans(AncestorAwareObj):
    """A class for handling translation and rotation of segments in a molecular structure.
    This class represents a move operation that can either be a translation or a rotation.
    Attributes
    ----------
    req_attr : list
        * movetype : str
            The type of move operation, either 'TRANS' for translation or 'ROT' for rotation.
    opt_attr : list
        * x : float
            The x-coordinate for translation (required if movetype is TRANS).
        * y : float
            The y-coordinate for translation (required if movetype is TRANS).
        * z : float
            The z-coordinate for translation (required if movetype is TRANS).
        * axis : str
            The axis of rotation (required if movetype is ROT).
        * angle : float
            The angle of rotation in degrees (required if movetype is ROT).
    yaml_header : str
        'transrot'
    objcat : str
        'coord'
    """
    
    yaml_header='transrot'
    objcat='coord'
    req_attr=AncestorAwareObj.req_attr+['movetype']
    opt_attr=AncestorAwareObj.opt_attr+['x','y','z','axis','angle']
    attr_choices=AncestorAwareObj.attr_choices.copy()
    attr_choices.update({'movetype':['trans','rot','TRANS','ROT']})
    opt_attr_deps=AncestorAwareObj.opt_attr_deps.copy()
    opt_attr_deps.update({'TRANS':['x','y','z'],'ROT':['axis','angle']})

    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)
    
    @__init__.register(str)
    def _from_shortcode(self,shortcode):
        dat=shortcode.split(',')
        input_dict=dict(movetype=dat[0].upper())
        if input_dict['movetype']=='TRANS':
            input_dict['x']=float(dat[1])
            input_dict['y']=float(dat[2])
            input_dict['z']=float(dat[3])
        elif input_dict['movetype']=='ROT':
            input_dict['axis']=str(dat[1])
            input_dict['angle']=float(dat[2])
        else:
            logger.warning(f'move type {input_dict["movetype"]} not recognized')
        super().__init__(input_dict)

    def to_shortcode(self):
        """Convert the RotTrans object to a shortcode representation.
        This method generates a shortcode that represents the move type and its parameters.
        """
        if 'shortcode' in self.__dict__:
            return
        ret=[f'{self.movetype}']
        if self.movetype=='TRANS':
            ret.append(f'{self.x:.2f}')
            ret.append(f'{self.y:.2f}')
            ret.append(f'{self.z:.2f}')
        elif self.movetype=='ROT':
            ret.append(f'{self.axis}')
            ret.append(f'{self.angle:.2f}')
        self.shortcode=','.join(ret)
    
    def __str__(self):
        self.to_shortcode()
        return self.shortcode
    
    def write_TcL(self,W:VMD,**kwargs):
        """Write the Tcl commands to perform the move operation in VMD.
        This method generates the Tcl commands to either translate or rotate a segment in VMD.
        Parameters
        ----------
        W : VMD
            The VMD script writer object to which the Tcl commands will be written.
        **kwargs : dict
            Additional keyword arguments (not used in this method).
        """
        molid_varname=W.molid_varname
        molid=f'${molid_varname}'
        W.addline(f'set mover [atomselect {molid} all]')
        if self.movetype=='TRANS':
            W.addline(f'$mover moveby [list {self.x} {self.y} {self.z}]')
        elif self.movetype=='ROT':
            W.addline(f'set COM [measure center $mover weight mass]')
            W.addline(f'$mover move [trans origin $COM axis {self.axis} {self.angle}]')

class RotTransList(AncestorAwareObjList):
    """A class for handling lists of RotTrans objects.
    This class inherits from AncestorAwareObjList and provides methods to manage
    a list of RotTrans objects.
    """
    def write_TcL(self,W:VMD,**kwargs):
        """Write the Tcl commands for each RotTrans object in the list.
        This method iterates over each RotTrans object in the list and calls its `write_TcL` method
        to generate the Tcl commands. The commands are written to the provided VMD script writer object.
        Parameters
        ----------
        W : VMD
            The VMD script writer object to which the Tcl commands will be written.
        **kwargs : dict
            Additional keyword arguments (not used in this method).
        """
        for c in self:
            c.write_TcL(W,**kwargs)