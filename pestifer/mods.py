"""

.. module:: mods
   :synopsis: Manages all modifications that can be encoded in PDB/mmCIF files or user-specified
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import os
import logging
logger=logging.getLogger(__name__)
from pidibble.pdbrecord import PDBRecord
from pidibble.baserecord import BaseRecord
from .cifutil import CIFdict
from .basemod import AncestorAwareMod,AncestorAwareModList
from .stringthings import ByteCollector, split_ri

class Seqadv(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['idCode','resname','chainID','resseqnum','insertion']
    opt_attr=AncestorAwareMod.opt_attr+['database','dbAccession','dbRes','dbSeq','conflict']
    attr_choices=AncestorAwareMod.attr_choices.copy()
    attr_choices.update({'conflict':['ENGINEERED','ENGINEERED MUTATION','CONFLICT']})
    yaml_header='Seqadv'
    PDB_keyword='SEQADV'
    ''' These are possible conflict keywords
        - Cloning artifact
        - Expression tag
        - Conflict
        - Engineered
        - Variant 
        - Insertion
        - Deletion
        - Microheterogeneity
        - Chromophore
    '''    
    def __init__(self,input_obj):
        if type(input_obj)==dict:
            super().__init__(input_obj)
        elif type(input_obj)==PDBRecord:
            pdbrecord=input_obj
            input_dict={
                'idCode':pdbrecord.idCode,
                'resname':pdbrecord.residue.resName,
                'chainID':pdbrecord.residue.chainID,
                'resseqnum':pdbrecord.residue.seqNum,
                'insertion':pdbrecord.residue.iCode,
                'database':pdbrecord.database,
                'dbAccession':pdbrecord.dbAccession,
                'dbRes':pdbrecord.dbRes,
                'dbSeq':pdbrecord.dbSeq,
                'conflict':pdbrecord.conflict
            }
            super().__init__(input_dict)
        else:
            logger.error(f'Cannot initialize {self.__class__} from object type {type(input_obj)}')

    def pdb_line(self):
        return f'SEQADV {self.idCode:3s} {self.resname:>3s} {self.chainID:1s} {self.resseqnum:>4d}{self.insertion:1s} {self.database:>4s} {self.dbAccession:9s} {self.dbRes:3s} {self.dbSeq:>5d} {self.conflict:21s}          '

class SeqAdvList(AncestorAwareModList):
    pass

class Mutation(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['chainID','origresname','resseqnum','insertion','newresname','source']
    opt_attr=AncestorAwareMod.opt_attr+['resseqnumi']
    yaml_header='Mutations'
    def __init__(self,input_obj):
        if type(input_obj)==dict:
            if not 'insertion' in input_obj:
                input_obj['insertion']=''
            if not 'source' in input_obj:
                input_obj['source']='SEQADV'
            super().__init__(input_obj)
        elif type(input_obj)==Seqadv:
            sq=input_obj
            input_dict={
                'chainID':sq.chainID,
                'origresname':sq.resname,
                'newresname':sq.dbRes,
                'resseqnum':sq.resseqnum,
                'insertion':sq.insertion,
                'resseqnumi':f'{sq.resseqnum}{sq.insertion}',
                'source':'SEQADV'
            }
            super().__init__(input_dict)
        elif type(input_obj)==str:
            ### shortcode format: c:nnn,rrr,mmm
            ### c: chainID
            ### nnn: old resname
            ### rrr: resseqnum
            ### mmm: new resname
            s1=input_obj.split(':')
            s2=s1[1].split(',')
            input_dict={
                'chainID':s1[0],
                'origresname':s2[0],
                'resseqnumi':s2[1],
                'newresname':s2[2],
                'source':'USER'
            }
            if not input_dict['resseqnumi'][-1].isdigit():
                input_dict['insertion']=input_dict['resseqnumi'][-1]
                input_dict['resseqnum']=int(input_dict['resseqnumi'][:-1])
            else:
                input_dict['insertion']=''
                input_dict['resseqnum']=int(input_dict['resseqnumi'])
            super().__init__(input_dict)
        else:
            logger.error(f'Cannot initialize {self.__class__} from object type {type(input_obj)}')
    def __str__(self):
        return f'{self.chainID}:{self.origresname}{self.resseqnum}{self.insertion}{self.newresname}'

    def write_TcL(self):
         return f'    mutate {self.resseqnumi} {self.newresname}'

class MutationList(AncestorAwareModList):
    pass

class Deletion(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['chainID','resseqnum1','insertion1','resseqnum2','insertion2']
    opt_attr=AncestorAwareMod.opt_attr+['model']
    yaml_header='Deletions'

    def __init__(self,input_obj):
        if type(input_obj)==dict:
            super().__init__(input_obj)
        elif type(input_obj)==str:
            p1=input_obj.split('_')
            input_dict={}
            input_dict['chainID']=p1[0]
            p2=p1[1].split('-')
            ri1=p2[0]
            if ri1[-1].isalpha():
                input_dict['resseqnum1']=int(ri1[:-1])
                input_dict['insertion1']=ri1[-1]
            else:
                input_dict['resseqnum1']=int(ri1)
                input_dict['insertion1']=''
            if len(p2)>1:
                ri2=p2[1]
                if ri2[-1].isalpha():
                    input_dict['resseqnum2']=int(ri2[:-1])
                    input_dict['insertion2']=ri2[-1]
                else:
                    input_dict['resseqnum2']=int(ri2)
                    input_dict['insertion2']=''
            else:
                input_dict['resseqnum2']=input_dict['resseqnum1']
                input_dict['insertion2']=input_dict['insertion1']
            super().__init__(input_dict)
        else:
            logger.error(f'Cannot initialize {self.__class__} from object type {type(input_obj)}')

class DeletionList(AncestorAwareModList):
    pass

class Missing(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['resname','resseqnum','insertion','chainID']
    opt_attr=AncestorAwareMod.opt_attr+['model']
    yaml_header='Missing'
    PDB_keyword='REMARK.465'

    def __init__(self,input_obj):
        if type(input_obj)==dict:
            super().__init__(input_obj)
            self.model=''
        elif type(input_obj) in [PDBRecord,BaseRecord]:
            pdbrecord=input_obj
            input_dict={
                'model':pdbrecord.modelNum,
                'resname':pdbrecord.resName,
                'chainID':pdbrecord.chainID,
                'resseqnum':pdbrecord.seqNum,
                'insertion':pdbrecord.iCode
            }
            super().__init__(input_dict)
        elif type(input_obj)==CIFdict:
            cifdict=input_obj
            ic=cifdict['pdb_ins_code']
            mn=cifdict['pdb_model_num']
            if type(mn)==str and mn.isdigit:
                nmn=int(mn)
            else:
                nmn=0
            input_dict={
                'model':nmn,
                'resname':cifdict['auth_comp_id'],
                'chainID':cifdict['auth_asym_id'],
                'resseqnum':int(cifdict['auth_seq_id']),
                'insertion':' ' if ic=='?' else ic
            }
            super().__init__(input_dict)
        else:
            logger.error(f'Cannot initialize {self.__class__} from object type {type(input_obj)}')
            
    def pdb_line(self):
        record_name,code=Missing.PDB_keyword.split(',')
        return '{:6s}{:>4d}   {:1s} {:3s} {:1s} {:>5d}{:1s}'.format(record_name,
        code,self.model,self.resname,self.chainID,self.resseqnum,self.insertion)

class MissingList(AncestorAwareModList):
    pass

class Crot(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['angle','degrees']
    opt_attr=AncestorAwareMod.opt_attr+['chainID','resseqnum1','resseqnum2','segname','atom1','atom2','segname1','segname2','segnamei','resseqnumi','atomi','segnamejk','resseqnumj','atomj','resseqnumk','atomk']
    attr_choices=AncestorAwareMod.attr_choices.copy()
    attr_choices.update({'angle':['PHI','PSI','OMEGA','CHI1','CHI2','GLYCAN','LINK','ANGLEIJK']})
    opt_attr_deps=AncestorAwareMod.opt_attr_deps.copy()
    opt_attr_deps.update({
        'PHI':['chainID','resseqnum1','resseqnum2'],
        'PSI':['chainID','resseqnum1','resseqnum2'],
        'OMEGA':['chainID','resseqnum1','resseqnum2'],
        'CHI1':['chainID','resseqnum1'],
        'CHI2':['chainID','resseqnum1'],
        'GLYCAN':['segname','resseqnum1','atom1','resseqnum2','atom2'],
        'LINK':['segname1','segname2','resseqnum1','atom1','resseqnum2','atom2'],
        'ANGLEIJK':['segnamei','resseqnumi','atomi','segnamejk','resseqnumj','atomj','resseqnumk','atomk']
        })
    yaml_header='Crotations'
    @classmethod
    def from_shortcode(cls,shortcode):
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
        elif input_dict['angle']=='GLYCAN':
            input_dict['segname']=dat[1]
            input_dict['resseqnum1']=int(dat[2])
            input_dict['atom1']=dat[3]
            input_dict['resseqnum2']=int(dat[4])
            input_dict['atom2']=dat[5]
            input_dict['degrees']=float(dat[6])
        elif input_dict['angle']=='LINK':
            input_dict['segname1']=dat[1]
            input_dict['resseqnum1']=int(dat[2])
            input_dict['atom1']=dat[3]
            input_dict['segname2']=dat[4]
            input_dict['resseqnum2']=int(dat[5])
            input_dict['atom2']=dat[6]
            input_dict['degrees']=float(dat[7])
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
        inst=cls(input_dict)
        return inst
    
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
        elif self.angle=='GLYCAN':
            ret.append(f'{self.segname}')
            ret.append(f'{self.resseqnum1}')
            ret.append(f'{self.atom1}')
            ret.append(f'{self.resseqnum2}')
            ret.append(f'{self.atom2}')
            ret.append(f'{self.degrees:.4f}')
        elif self.angle=='LINK':
            ret.append(f'{self.segname1}')
            ret.append(f'{self.resseqnum1}')
            ret.append(f'{self.atom1}')
            ret.append(f'{self.segname2}')
            ret.append(f'{self.resseqnum2}')
            ret.append(f'{self.atom2}')
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
        self.shortcode=','.join(ret)
    
    def __str__(self):
        self.to_shortcode()
        return self.shortcode
        
    def psfgen_lines(self,**kwargs):
        molid=kwargs.get('molid','top')
        endIsCterm=kwargs.get('endIsCterm',False)
        retlines=[]
        if self.angle in ['PHI','PSI','OMEGA']:  # this is a backbone bond
            retlines.append('set r1 [[atomselect {} "chain {} and resid {} and name CA"] get residue]'.format(molid,self.chainID,self.resseqnum1))
            retlines.append('set r2 [[atomselect {} "chain {} and resid {} and name CA"] get residue]'.format(molid,self.chainID,self.resseqnum2))
            if endIsCterm:
                retlines.append('Crot_{}_toCterm $r1 $r2 {} {} {}'.format(self.angle.lower(),self.chainID,molid,self.degrees))
            else:
                retlines.append('Crot_{} $r1 $r2 {} {} {}'.format(self.angle.lower(),self.chainID,molid,self.degrees))
        elif self.angle in ['CHI1','CHI2']:  # this is a side-chain bond
            retlines.append('set r1 [[atomselect {} "chain {} and resid {} and name CA"] get residue]'.format(molid,self.chainID,self.resseqnum1))
            retlines.append('SCrot_{} $r1 {} {} {}'.format(self.angle.lower(),self.chainID,molid,self.degrees))
        elif self.angle=='GLYCAN':  # intra-glycan rotation
            retlines.append('set sel [atomselect {} "segname {}"]'.format(molid,self.segname))
            retlines.append('set i [[atomselect {} "segname {} and resid {} and name {}"] get index]'.format(molid,self.segname,self.resseqnum1,self.atom1))
            retlines.append('set j [[atomselect {} "segname {} and resid {} and name {}"] get index]'.format(molid,self.segname,self.resseqnum2,self.atom2))
            retlines.append('genbondrot {} $sel $i $j {}'.format(molid,self.degrees))
        elif self.angle=='LINK': # ASN-GLYcan rotation
            retlines.append('set sel [atomselect {} "segname {} {}"]'.format(molid,self.segname1,self.segname2))
            retlines.append('set i [[atomselect {} "segname {} and resid {} and name {}"] get index]'.format(molid,self.segname1,self.resseqnum1,self.atom1))
            retlines.append('set j [[atomselect {} "segname {} and resid {} and name {}"] get index]'.format(molid,self.segname2,self.resseqnum2,self.atom2))
            retlines.append('genbondrot {} $sel $i $j {}'.format(molid,self.degrees))
        elif self.angle=='ANGLEIJK':
            retlines.extend([
                'set rotsel [atomselect {} "segname {}"]'.format(molid,self.segnamejk),
                'set ri [lindex [[atomselect {} "segname {} and resid {} and name {}"] get {{x y z}}] 0]'.format(molid,self.segnamei,self.resseqnumi,self.atomi),
                'set rj [lindex [[atomselect {} "segname {} and resid {} and name {}"] get {{x y z}}] 0]'.format(molid,self.segnamejk,self.resseqnumj,self.atomj),
                'set rk [lindex [[atomselect {} "segname {} and resid {} and name {}"] get {{x y z}}] 0]'.format(molid,self.segnamejk,self.resseqnumk,self.atomk),
                'set rij [vecsub $ri $rj]',
                'set rjk [vecsub $rj $rk]',
                'set cijk [veccross $rij $rjk]',
                '$rotsel move [trans center $rj origin $rj axis $cijk {} degrees]'.format(self.degrees)
            ])
        return retlines

class CrotList(AncestorAwareModList):
    pass

class SSBond(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['chainID1','resseqnum1','insertion1','chainID2','resseqnum2','insertion2']
    opt_attr=AncestorAwareMod.opt_attr+['serial_number','resname1','resname2','sym1','sym2','length']
    yaml_header='SSBonds'
    PDB_keyword='SSBOND'

    def __init__(self,input_obj):
        if type(input_obj)==dict:
            super().__init__(input_obj)
        elif type(input_obj)==PDBRecord:
            pdbrecord=input_obj
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
                'length':pdbrecord.length
            }
            super().__init__(input_dict)
        elif type(input_obj)==CIFdict:
            d=input_obj
            input_dict={
                'serial_number':int(d['id'].strip('disulf')),
                'resname1':'CYS',
                'resname2':'CYS',
                'chainID1':d['ptnr1_auth_asym_id'],
                'chainID2':d['ptnr2_auth_asym_id'],
                'resseqnum1':int(d['ptnr1_auth_seq_id']),
                'resseqnum2':int(d['ptnr2_auth_seq_id']),
                'insertion1':d['pdbx_ptnr1_pdb_ins_code'],
                'insertion2':d['pdbx_ptnr2_pdb_ins_code'],
                'sym1':d['ptnr1_symmetry'],
                'sym2':d['ptnr2_symmetry'],
                'length':float(d['pdbx_dist_value'])
            }
            super().__init__(input_dict)
        elif type(input_obj)==str:
            #C_RRR-C_RRR
            s1=input_obj.split('-')
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
                'sym2':''
            }
            super().__init__(input_dict)
        else:
            logger.error(f'Cannot initialize {self.__class__} from object type {type(input_obj)}')

    def pdb_line(self):
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
    
    def write_TcL(self,B:ByteCollector,transform):
        chainIDmap=transform.chainIDmap 
        # ok since these are only going to reference protein segments; protein segment names are the chain IDs
        c1=chainIDmap.get(self.chainID1,self.chainID1)
        c2=chainIDmap.get(self.chainID2,self.chainID2)
        r1=self.resseqnum1
        r2=self.resseqnum2
        B.addline(f'patch DISU {c1}:{r1} {c2}:{r2}')

class SSBondList(AncestorAwareModList):
    def write_TcL(self,B:ByteCollector,transform,mods):
        undo_SSBonds=mods.get('SSBondDelete',SSBondDeleteList([]))
        for s in self:
            if undo_SSBonds.is_deleted(s,transform):
                continue
            s.write_TcL(B,transform)
        for s in mods.get('SSBonds',[]):
            s.writeTcL(B,transform)
    def prune_mutations(self,Mutations):
        for m in Mutations:
            left=self.get(chainID1=m.chainID,resseqnum1=m.resseqnum,insertion1=m.insertion)
            if left:
                self.remove(left)
            right=self.get(chainID2=m.chainID,resseqnum2=m.resseqnum,insertion2=m.insertion)
            if right:
                self.remove(right)

class SSBondDelete(SSBond):
    yaml_header='SSBondsDelete'

class SSBondDeleteList(SSBondList):
    def is_deleted(self,a_SSBond,transform):
        if self.get(
            chainID1=transform.chainIDmap[a_SSBond.chainID1],
            chainID2=transform.chainIDmap[a_SSBond.chainID2],
            resseqnum1=a_SSBond.resseqnum1,
            resseqnum2=a_SSBond.resseqnum2):
            return True
        return False

class Link(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['name1','chainID1','resseqnum1','insertion1','name2','chainID2','resseqnum2','insertion2']
    opt_attr=AncestorAwareMod.opt_attr+['altloc1','altloc2','resname1','resname2','sym1','sym2','link_distance','segname1','segname2','residue1','residue2','atom1','atom2','empty','segtype1','segtype2']    
    yaml_header='Links'
    PDB_keyword='LINK'

    def __init__(self,input_obj):
        if type(input_obj)==PDBRecord:
            pdbrecord=input_obj
            input_dict={
                'name1': pdbrecord.name1,
                'resname1':pdbrecord.residue1.resName,
                'chainID1':pdbrecord.residue1.chainID,
                'resseqnum1':pdbrecord.residue1.seqNum,
                'insertion1':pdbrecord.residue1.iCode,
                'name2': pdbrecord.name2,
                'resname2':pdbrecord.residue2.resName,
                'chainID2':pdbrecord.residue2.chainID,
                'resseqnum2':pdbrecord.residue2.seqNum,
                'insertion2':pdbrecord.residue2.iCode,
                'altloc2':pdbrecord.altLoc2,
                'altloc1':pdbrecord.altLoc1,
                'sym1':pdbrecord.sym1,
                'sym2':pdbrecord.sym2,
                'link_distance':pdbrecord.length,
                'segname1':pdbrecord.residue1.chainID,
                'segname2':pdbrecord.residue2.chainID,
                }
        elif type(input_obj)==CIFdict:
            d=input_obj
            input_dict={}
            input_dict['name1']=d['ptnr1_label_atom_id']
            al=d['pdbx_ptnr1_label_alt_id']
            input_dict['altloc1']=' ' if al=='?' else al
            input_dict['resname1']=d['ptnr1_auth_comp_id']
            input_dict['chainID1']=d['ptnr1_auth_asym_id']
            input_dict['resseqnum1']=int(d['ptnr1_auth_seq_id'])
            ic=d['pdbx_ptnr1_pdb_ins_code']
            input_dict['insertion1']=' ' if ic=='?' else ic
            input_dict['name2']=d['ptnr2_label_atom_id']
            al=d['pdbx_ptnr2_label_alt_id']
            input_dict['altloc2']=' ' if al=='?' else al
            input_dict['resname2']=d['ptnr2_auth_comp_id']
            input_dict['chainID2']=d['ptnr2_auth_asym_id']
            input_dict['resseqnum2']=int(d['ptnr2_auth_seq_id'])
            ic=d['pdbx_ptnr2_pdb_ins_code']
            input_dict['insertion2']=' ' if ic=='?' else ic
            input_dict['sym1']=d['ptnr1_symmetry']
            input_dict['sym2']=d['ptnr2_symmetry']
            input_dict['link_distance']=float(d['pdbx_dist_value'])
            input_dict['segname1']=input_dict['chainID1']
            input_dict['segname2']=input_dict['chainID2']
        elif type(input_obj)==dict:
            input_dict=input_obj
        else:
            logger.error(f'Cannot initialize {self.__class__} from object type {type(input_obj)}')
        input_dict.update({
                'residue1':None,
                'residue2':None,
                'atom1':None,
                'atom2':None,
                'empty':False,
                'segtype1':'UNSET',
                'segtype2':'UNSET',
            })
        super().__init__(input_dict)

    def write_TcL(self,B:ByteCollector,transform,mods):
        if self.resname1=='ASN' and self.segtype2=='GLYCAN':
            chainIDmap=transform.chainIDmap
            seg1=chainIDmap[self.segname1]
            seg2=transform.segname_by_type_map['GLYCAN'][self.segname2]
            B.addline(f'patch NGLB {seg1}:{self.resseqnum1}{self.insertion1} {seg2}:{self.resseqnum2}{self.insertion2}')
        else:
            # this is likely an intra-glycan linkage
            if self.name2=='C1' and self.segtype1=='GLYCAN':
                seg1=transform.segname_by_type_map['GLYCAN'][self.segname1]
                seg2=transform.segname_by_type_map['GLYCAN'][self.segname2]
                B.addline(f'set cn {self.name1[1]}')
                B.addline(f'set abi [axeq {self.resseqnum2} 0 {self.chainID2} {self.name2} {self.resseqnum1}]')
                B.addline(f'set abj [axeq {self.resseqnum1} 0 {self.chainID1} {self.name1} -1]')
                if self.name1=='O6':
                    B.addline('if { $abi == "a" } { set abi A }')
                    B.addline('if { $abi == "b" } { set abi B }')
                    B.addline('if { $abj == "b" } { set abj T }')
                B.addline('set pres "1$cn$abi$abj"')
                B.addline(f'patch $pres {seg1}:{self.resseqnum1}{self.insertion1} {seg2}:{self.resseqnum1}{self.insertion2}')
            elif self.name1=='C1' and self.segtype2=='GLYCAN':
                seg1=transform.segname_by_type_map['GLYCAN'][self.segname1]
                seg2=transform.segname_by_type_map['GLYCAN'][self.segname2]
                cmdj=f'[axeq {self.resseqnum2} 0 {self.chainID2} {self.name2} {self.resseqnum1}]'
                cmdi=f'[axeq {self.resseqnum1} 0 {self.chainID1} {self.name1} -1]'
                B.addline(f'patch 1{self.name2[1]:1s}{cmdi}{cmdj} {seg2}:{self.resseqnum2}{self.insertion2} {seg1}:{self.resseqnum1}{self.insertion1}')
            elif self.name1=='O6' and self.name2=='C2':
                seg1=transform.segname_by_type_map['GLYCAN'][self.segname1]
                seg2=transform.segname_by_type_map['GLYCAN'][self.segname2]
                B.addline(f'patch SA26AT {seg1}:{self.resseqnum1}{self.insertion1} {seg2}:{self.resseqnum2}{self.insertion2}')
            else:
                logger.warning(f'Could not identify patch for link: {str(self)}')
                B.comment(f'No patch found for {str(self)}')
        
    def __str__(self):
        return f'{self.chainID1}{self.resname1}{self.resseqnum1}{self.insertion1}-{self.chainID2}{self.resname2}{self.resseqnum2}{self.insertion2}'
    


class LinkList(AncestorAwareModList):
    def write_TcL(self,B:ByteCollector,transform,mods):
        for l in self:
            l.write_TcL(B,transform,mods)

    def prune_mutations(self,Mutations,Segments):
        for m in Mutations:
            rlist,llist=[],[]
            left=self.get(chainID1=m.chainID,resseqnum1=m.resseqnum,insertion1=m.insertion)
            right=self.get(chainID2=m.chainID,resseqnum2=m.resseqnum,insertion2=m.insertion)
            if left:  # this is a link in which this mutation is the left member
                self.remove(left) # get rid of this link
                # we need to remove residue2 and everything downstream
                # remove downstream residues!
                rlist,llist=left.residue2.get_down_group()
                rlist.insert(0,left.residue2)
            elif right: # this is a link in which this mutation is the right member (should be very rare)
                self.remove(right)
                rlist,llist=right.residue2.get_down_group()
                rlist.insert(0,right.residue2)
            if rlist and llist:
                # logger.debug(f'Deleting residues down from and including {str(rlist[0])} due to a mutation')
                S=Segments.get_segment_of_residue(rlist[0])
                for r in rlist:
                    # logger.debug(f'...{str(r)}')
                    S.residues.remove(r)
                if len(S.residues)==0:
                    # logger.debug(f'All residues of {S.psfgen_segname} are deleted; {S.psfgen_segname} is deleted')
                    Segments.remove(S)
                for l in llist:
                    self.remove(l)

    def apply_segtypes(self,map):
        self.map_attr('segtype1','resname1',map)
        self.map_attr('segtype2','resname2',map)

class Graft(AncestorAwareMod):
    yaml_header='Grafts'
    pass

class GraftList(AncestorAwareModList):
    pass

class Insertion(AncestorAwareMod):
    yaml_header='Insertions'
    pass

class InsertionList(AncestorAwareModList):
    pass

class Cleavage(AncestorAwareMod):
    yaml_header='Cleavages'
    pass

class CleavageList(AncestorAwareModList):
    pass

def apply_psf_info(p_struct,psf):
    if os.path.exists(psf):
        with open(psf,'r') as f:
            psf_lines=f.read().split('\n')
    metadata_lines=[x for x in psf_lines if x.startswith('REMARKS')]
    for ml in metadata_lines:
        words=ml.split()
        if ' '.join(words[1:3])=='patch DISU':
            c1,r1i=words[3].split(':')
            c2,r2i=words[4].split(':')
            r1,i1=split_ri(r1i)
            r2,i2=split_ri(r2i)
            p_struct['SSBONDS'].append(SSBond({'chainID1':c1,'chainID2':c2,'resseqnum1':r1,'resseqnum2':r2,'insertion1':i1,'insertion2':i2}))
        # elif ' '.join(words[1:3])=='patch '
            # add a new SSBOND record to p_struct
        # elif ...