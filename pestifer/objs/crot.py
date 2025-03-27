# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import logging
logger=logging.getLogger(__name__)
from functools import singledispatchmethod

from ..baseobj import AncestorAwareObj, AncestorAwareObjList
from ..scriptwriters import VMD, Psfgen

class Crot(AncestorAwareObj):
    """A class for managing so-called 'C-rotations'
    
    A C-rotation is a transformation in which atoms are rotated around a given bond by a given amount.  The "C" 
    designation means that only the "downstream" atoms of the bond are moved; upstream atoms, along with the atoms
    of the bond itself, are not moved.  The set of upstream atoms and the set of downstream atoms must naturally
    have no topological connection *other* than the bond itself.  Typically, this can be used to execute rotation
    of a backbone look in a C-terminal loop, a side-chain angle, usually in the service
    of reducing steric clashses.  The primary job of this class is to translate the C-rotation shortcodes
    specified by the user into TcL commands to be incorporated in a psfgen script.

    NOTE: This is currently implemented in the cfapdbparse (2020) format, and has not been thoroughly tested.

    """
    req_attr=AncestorAwareObj.req_attr+['angle']
    opt_attr=AncestorAwareObj.opt_attr+['chainID','resseqnum1','resseqnum2','resseqnum3','segname','atom1','atom2','segname1','segname2','segnamei','resseqnumi','atomi','segnamejk','resseqnumj','atomj','resseqnumk','atomk','degrees']
    attr_choices=AncestorAwareObj.attr_choices.copy()
    attr_choices.update({'angle':['PHI','PSI','OMEGA','CHI1','CHI2','ANGLEIJK','ALPHA']})
    opt_attr_deps=AncestorAwareObj.opt_attr_deps.copy()
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
    objcat='coord'

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
    
    def write_TcL(self,W:VMD,chainIDmap={},**kwargs):
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
    def write_TcL(self,W:Psfgen,chainIDmap={},**kwargs):
        for c in self:
            c.write_TcL(W,chainIDmap=chainIDmap,**kwargs)
