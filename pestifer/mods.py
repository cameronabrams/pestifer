"""

.. module:: mods
   :synopsis: Manages all modifications
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
# from pestifer.mutation import Mutation
# from pestifer.ssbond import SSBond
# from pestifer.graft import Graft
# from pestifer.crot import Crot
# from pestifer.attach import Attach
# from pestifer.link import Link
# from pestifer.deletion import Deletion
# from pestifer.cleavage import Cleavage
# from pestifer.missing import Missing
from pidibble.pdbrecord import PDBRecord
from .basemod import AncestorAwareMod,AncestorAwareModList

class Seqadv(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['idCode','resname','chainID','resseqnum','insertion']
    opt_attr=AncestorAwareMod.opt_attr+['database','dbAccession','dbRes','dbSeq','conflict']
    yaml_header='Seqadv'
    PDB_keyword='SEQADV'

    @classmethod
    def from_pdbrecord(cls,pdbrecord:PDBRecord):
        """from_pdbrecord generates a AncestorAwareMod Seqadv from a SEQADV PDBRecord

        :param pdbrecord: a SEQADV PDBRecord
        :type pdbrecord: PDBRecord
        :return: a Seqadv AncestorAwareMod
        :rtype: Seqadv

        SEQADV
              idCode: 1GC1
             residue: resName: GLY; chainID: G; seqNum: 79; iCode: 
            database: UNP
         dbAccession: P04578
               dbRes: 
               dbSeq: 
            conflict: EXPRESSION TAG

        """
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
        inst=cls(input_dict)
        return inst
    def pdb_line(self):
        return f'SEQADV {self.idCode:3s} {self.resname:>3s} {self.chainID:1s} {self.resseqnum:>4d}{self.insertion:1s} {self.database:>4s} {self.dbAccession:9s} {self.dbRes:3s} {self.dbSeq:>5d} {self.conflict:21s}          '

class SeqAdvList(AncestorAwareModList):
    pass

class Mutation(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['chainID','origresname','newresname']
    opt_attr=AncestorAwareMod.opt_attr+['resseqnum','resseqnumi','insertion']
    yaml_header='Mutations'
    @classmethod
    def from_seqadv(cls,sq:Seqadv):
        if not 'MUTATION' in sq.conflict:
            return None
        input_dict={}
        input_dict['chainID']=sq.chainID
        input_dict['origresname']=sq.resname
        input_dict['newresname']=sq.dbRes
        input_dict['resseqnum']=sq.resseqnum
        input_dict['insertion']=sq.insertion
        input_dict['resseqnumi']=f'{sq.resseqnum}{sq.insertion}'
        inst=cls(input_dict)
        return inst
    
    @classmethod
    def from_shortcode(cls,shortcode):
        ### c:nnn,rrr,nnn
        s1=shortcode.split(':')
        s2=s1[1].split(',')
        input_dict={
            'chainID':s1[0],
            'origresname':s2[0],
            'resseqnumi':s2[1],
            'newresname':s2[2]
        }
        if not input_dict['resseqnumi'][-1].isdigit():
            input_dict['insertion']=input_dict['resseqnumi'][-1]
            input_dict['resseqnum']=int(input_dict['resseqnumi'][:-1])
        else:
            input_dict['insertion']=' '
            input_dict['resseqnum']=int(input_dict['resseqnumi'])
        inst=cls(input_dict)
        return inst

    def psfgen_lines(self):
         lines=['    mutate {self.resseqnumi} {self.newresname}']
         return lines

class MutationList(AncestorAwareModList):
    pass

class Missing(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['resname','resseqnum','insertion','chainID']
    opt_attr=AncestorAwareMod.opt_attr+['model']
    yaml_header='Missing'
    PDB_keyword='REMARK.465'

    @classmethod
    def from_pdbrecord(cls,pdbrecord:PDBRecord):
        input_dict={
            'model':pdbrecord.modelNum,
            'resname':pdbrecord.resName,
            'chainID':pdbrecord.chainID,
            'resseqnum':pdbrecord.seqNum,
            'insertion':pdbrecord.iCode
        }
        inst=cls(input_dict)
        return inst

    @classmethod
    def from_cifdict(cls,cifdict):
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
        inst=cls(input_dict)
        return inst

    def pdb_line(self):
        record_name,code=Missing.PDB_keyword.split(',')
        return '{:6s}{:>4d}   {:1s} {:3s} {:1s} {:>5d}{:1s}'.format(record_name,
        code,self.model,self.resname,self.chainID,self.resseqnum,self.insertion)

    def psfgen_lines(self):
        return ['     residue {} {}{} {}'.format(self.resname,self.resseqnum,self.insertion,self.chainID)]

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

class Atom(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['serial','name','altloc','resname','chainID','resseqnum','insertion','x','y','z','occ','beta','elem','charge']
    opt_attr=AncestorAwareMod.opt_attr+['segname','empty','link']
    yaml_header='Atoms'
    PDB_keyword='ATOM'

    @classmethod
    def from_pdbrecord(cls,pdbrecord:PDBRecord):
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
        inst=cls(input_dict)
        return inst
    @classmethod
    def from_cifdict(cls,cifdict):
        al=cifdict['label_alt_id']
        ic=cifdict['pdbx_pdb_ins_code']
        c=cifdict['pdbx_formal_charge']
        input_dict={
            'recordname':'ATOM',
            'serial':int(cifdict['id']),
            'name':cifdict['auth_atom_id'],
            'altloc':' ' if al=='.' else al,
            'resname':cifdict['auth_comp_id'],
            'chainID':cifdict['auth_asym_id'],
            'resseqnum':int(cifdict['auth_seq_id']),
            'insertion':' ' if ic=='?' else ic,
            'x':float(cifdict['cartn_x']),
            'y':float(cifdict['cartn_y']),
            'z':float(cifdict['cartn_z']),
            'occ':float(cifdict['occupancy']),
            'beta':float(cifdict['b_iso_or_equiv']),
            'elem':cifdict['type_symbol'],
            'charge':0.0 if c=='?' else float(c)
        }
        input_dict['segname']=input_dict['chainID']
        input_dict['link']='None'
        input_dict['empty']=False
        inst=cls(input_dict)
        return inst
    
    def pdb_line(self):
        pdbline='{:<6s}'.format(self.recordname)+\
                '{:5d}'.format(self.serial)+' '+\
                '{:<4s}'.format(' '+self.name if len(self.name)<4 else self.name)+\
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

class AtomList(AncestorAwareModList):
    pass

class Hetatm(Atom):
    PDB_keyword='HETATM'

class SSBond(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['chainID1','resseqnum1','insertion1','chainID2','resseqnum2','insertion2']
    opt_attr=AncestorAwareMod.opt_attr+['serial_number','resname1','resname2','sym1','sym2','length']
    yaml_header='SSBonds'
    PDB_keyword='SSBOND'

    @classmethod
    def from_pdbrecord(cls,pdbrecord:PDBRecord):
        """from_pdbrecord generates a AncestorAwareMod SSBond from a PDBRecord

        :param pdbrecord: an SSBOND PDBRecord
        :type pdbrecord: PDBRecord
        :return: a new SSBond instance
        :rtype: SSBond

        SSBOND
              serNum: 1
            residue1: resName: CYS; chainID: G; seqNum: 119; iCode: 
            residue2: resName: CYS; chainID: G; seqNum: 205; iCode: 
                sym1: 1555
                sym2: 1555
              length: 2.04

        """
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
        inst=cls(input_dict)
        return inst

    @classmethod
    def from_cifdict(cls,cifdict):
        d=cifdict
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
        inst=cls(input_dict)
        return inst

    @classmethod
    def from_shortcode(cls,shortcode):
        #C_RRR-C_RRR
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
            'sym2':''
        }
        inst=cls(input_dict)
        return inst
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
    def psfgen_lines(self):
        return ['patch DISU {}:{} {}:{}'.format(self.chainID1,self.resseqnum1,self.chainID2,self.resseqnum2)]

class SSBondList(AncestorAwareModList):
    pass

class Link(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['name1','chainID1','resseqnum1','iCode1','name2','chainID2','resseqnum2','iCode2']
    opt_attr=AncestorAwareMod.opt_attr+['altloc1','altloc2','resname1','resname2','sym1','sym2','link_distance']    
    yaml_header='Links'
    PDB_keyword='LINK'
    @classmethod
    def from_pdbrecord(cls,pdbrecord:PDBRecord):
        input_dict={
            'name1': pdbrecord.name1,
            'altloc1':pdbrecord.altLoc1,
            'resname1':pdbrecord.residue1.resName,
            'chainID1':pdbrecord.residue1.chainID,
            'resseqnum1':pdbrecord.residue1.seqNum,
            'iCode1':pdbrecord.residue1.iCode,
            'name2': pdbrecord.name2,
            'altloc2':pdbrecord.altLoc2,
            'resname2':pdbrecord.residue2.resName,
            'chainID2':pdbrecord.residue2.chainID,
            'resseqnum2':pdbrecord.residue2.seqNum,
            'iCode2':pdbrecord.residue2.iCode,
            'sym1':pdbrecord.sym1,
            'sym2':pdbrecord.sym2,
            'link_distance':pdbrecord.length
        }
        inst=cls(input_dict)
        inst.segname1=inst.chainID1
        inst.segname2=inst.chainID2
        inst.residue1=None
        inst.residue2=None
        inst.atom1=None
        inst.atom2=None
        inst.empty=False
        return inst
    @classmethod
    def from_cifdict(cls,cifdict):
        d=cifdict
        input_dict={}
        input_dict['name1']=d['ptnr1_label_atom_id']
        al=d['pdbx_ptnr1_label_alt_id']
        input_dict['altloc1']=' ' if al=='?' else al
        input_dict['resname1']=d['ptnr1_auth_comp_id']
        input_dict['chainID1']=d['ptnr1_auth_asym_id']
        input_dict['resseqnum1']=int(d['ptnr1_auth_seq_id'])
        ic=d['pdbx_ptnr1_pdb_ins_code']
        input_dict['icode1']=' ' if ic=='?' else ic
        input_dict['name2']=d['ptnr2_label_atom_id']
        al=d['pdbx_ptnr2_label_alt_id']
        input_dict['altloc2']=' ' if al=='?' else al
        input_dict['resname2']=d['ptnr2_auth_comp_id']
        input_dict['chainID2']=d['ptnr2_auth_asym_id']
        input_dict['resseqnum2']=int(d['ptnr2_auth_seq_id'])
        ic=d['pdbx_ptnr2_pdb_ins_code']
        input_dict['icode2']=' ' if ic=='?' else ic
        input_dict['sym1']=d['ptnr1_symmetry']
        input_dict['sym2']=d['ptnr2_symmetry']
        input_dict['link_distance']=float(d['pdbx_dist_value'])
        input_dict['segname1']=input_dict['chainID1']
        input_dict['segname2']=input_dict['chainID2']
        input_dict['residue1']=None
        input_dict['residue2']=None
        input_dict['atom1']=None
        input_dict['atom2']=None
        input_dict['empty']=False
        return cls(input_dict)

    def __str__(self):
        return f'{self.chainID1}{self.resname1}{self.resseqnum1}{self.iCode1}-{self.chainID2}{self.resname2}{self.resseqnum2}{self.iCode2}'
