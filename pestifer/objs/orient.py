# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import logging
logger=logging.getLogger(__name__)
from functools import singledispatchmethod

from ..baseobj import AncestorAwareObj, AncestorAwareObjList
from ..scriptwriters import VMD

class Orient(AncestorAwareObj):
    req_attr=AncestorAwareObj.req_attr+['axis']
    opt_attr=AncestorAwareObj.opt_attr+['refatom']
    yaml_header='orient'
    objcat='coord'

    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)
    
    @__init__.register(str)
    def _from_shortcode(self,shortcode):
        logger.debug(f'orient shortcode {shortcode}')
        dat=shortcode.split(',')
        input_dict=dict(axis=dat[0])
        if len(dat)>1:
            input_dict['refatom']=dat[1]
        super().__init__(input_dict)

    def write_TcL(self,W:VMD):
        W.addline('set a [atomselect top all]')
        W.addline('set I [draw principalaxes $a]')
        adict=dict(x=r'{1 0 0}',y=r'{0 1 0}',z=r'{0 0 1}')
        W.addline(f'set A [orient $a [lindex $I 2] {adict[self.axis]}]')
        W.addline(r'$a move $A')
        if hasattr(self,'refatom'):
            W.addline(r'set com [measure center $a]')
            W.addline(r'$a moveby [vecscale $com -1]')
            W.addline(f'set z [[atomselect top "name {self.refatom}"] get z]')
            W.addline(r'if { $z < 0.0 } {')
            W.addline(r'   $a move [transaxis x 180 degrees]')
            W.addline(r'}')
            
class OrientList(AncestorAwareObjList):
    def write_TcL(self,W:VMD):
        for c in self:
            c.write_TcL(W)
