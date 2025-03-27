# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import logging
logger=logging.getLogger(__name__)
from functools import singledispatchmethod

from ..baseobj import AncestorAwareObj, AncestorAwareObjList
from ..scriptwriters import VMD

class RotTrans(AncestorAwareObj):
    """A class for managing global rotations and translations """
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
        molid_varname=W.molid_varname
        molid=f'${molid_varname}'
        W.addline(f'set mover [atomselect {molid} all]')
        if self.movetype=='TRANS':
            W.addline(f'$mover moveby [list {self.x} {self.y} {self.z}]')
        elif self.movetype=='ROT':
            W.addline(f'set COM [measure center $mover weight mass]')
            W.addline(f'$mover move [trans origin $COM axis {self.axis} {self.angle}]')

class RotTransList(AncestorAwareObjList):
    def write_TcL(self,W:VMD,**kwargs):
        for c in self:
            c.write_TcL(W,**kwargs)