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
import yaml
# from .seqadv import Seqadv
class BaseMod:
    # list of required attribute labels
    req_attr=[]
    # list of optional attribute labels
    opt_attr=[]
    # list of 2-tuples of attribute labels for mutually exclusive attributes
    alt_attr=[]
    # dictionary of lists of value choices for any attribute, keyed on attribute label
    attr_choices={}
    # dictionary of lists of labels of attributes that are required if the optional attributes whose labels key this dictionary are present 
    opt_attr_deps={}
    def __init__(self,input_dict={}):
        # mutually exclusive attribute labels should not appear in the list of required attributes
        assert all([(not x in self.req_attr and not y in self.req_attr) for x,y in self.alt_attr]),"Mutually exclusive attributes should not appear in the required list"
        for k,v in input_dict.items():
            if k in self.req_attr or k in self.opt_attr:
                self.__dict__[k]=v
        # all required attributes have values
        assert all([att in self.__dict__.keys() for att in self.req_attr]),f"Not all required attributes have values {self.__dict__}"
        # all mutual-exclusivity requirements are met
        assert all([
            ((x[0] in self.__dict__.keys() and not x[1] in self.__dict__.keys()) or 
             (x[1] in self.__dict__.keys() and not x[0] in self.__dict__.keys()) or not (x[0] in self.__dict__.keys() or x[1] in self.__dict__.keys()))
            for x in self.alt_attr]),"Mutual exclusivity requirements unmet"
        # all attributes with limited set of possible values have valid values
        assert all([self.__dict__[x] in choices for x,choices in self.attr_choices.items()]),"Invalid choices in one or more attributes"
        # for each optional attribute present, all required dependent attributes are also present
        for oa,dp in self.opt_attr_deps.items():
            if oa in self.__dict__.keys():
                assert all([x in self.__dict__.keys() for x in dp]),"Dependent required attributes of optional attributes not present"
    def __eq__(self,other):
        test_list=[self.__dict__[k]==other.__dict__[k] for k in self.req_attr]
        for k in self.opt_attr:
            if k in self.__dict__ and k in other.__dict__:
                test_list.append(self.__dict__[k]==other.__dict__[k])
        return all(test_list)
    def dump(self):
        retdict={}
        retdict['instanceOf']=type(self).__name__
        retdict.update(self.__dict__)
        return yaml.dump(retdict)
    def inlist(self,a_list):
        for s in a_list:
            if s==self:
                return True
        return False

class CloneableMod(BaseMod):
    @classmethod
    def clone(cls,inst,**options):
        input_dict={k:v for k,v in inst.__dict__.items() if k in inst.req_attr}
        input_dict.update({k:v for k,v in inst.__dict__.items() if k in inst.opt_attr})
        input_dict.update(options)
        return cls(input_dict)
    
class Seqadv(CloneableMod):
    req_attr=['idCode','resname','chainID','resseqnum','insertion']
    opt_attr=['database','dbAccession','dbRes','dbSeq','conflict']
    yaml_header='Seqadv'
    PDB_keyword='SEQADV'

    @classmethod
    def from_pdbrecord(cls,pdbrecord):
        dbseqraw=pdbrecord[43:48].strip()
        input_dict={
            'idCode':pdbrecord[7:11],
            'resname':pdbrecord[12:15].strip(),
            'chainID':pdbrecord[16:17],
            'resseqnum':int(pdbrecord[18:22]),
            'insertion':pdbrecord[22:23],
            'database':pdbrecord[24:28],
            'dbAccession':pdbrecord[29:38],
            'dbRes':pdbrecord[39:42].strip(),
            'dbSeq':int(dbseqraw) if dbseqraw.isdigit() else -1,
            'conflict':pdbrecord[49:70].strip()
        }
        inst=cls(input_dict)
        return inst
    def pdb_line(self):
        return f'SEQADV {self.idCode:3s} {self.resname:>3s} {self.chainID:1s} {self.resseqnum:>4d}{self.insertion:1s} {self.database:>4s} {self.dbAccession:9s} {self.dbRes:3s} {self.dbSeq:>5d} {self.conflict:21s}          '
    
class Mutation(CloneableMod):
    req_attr=['chainID','origresname','newresname']
    opt_attr=['resseqnum','resseqnumi','insertion']
    yaml_header='Mutations'
    @classmethod
    def from_seqdav(cls,sq:Seqadv):
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

class Missing(CloneableMod):
    req_attr=['resname','resseqnum','insertion','chainID']
    opt_attr=['model']
    yaml_header='Missing'
    PDB_keyword='REMARK,465'

    @classmethod
    def from_pdbrecord(cls,pdbrecord):
        input_dict={
            'model':pdbrecord[13:14].strip(),
            'resname':pdbrecord[15:18].strip(),
            'chainID':pdbrecord[19:20],
            'resseqnum':int(pdbrecord[21:26]),
            'insertion':pdbrecord[26:27]
        }
        inst=cls(input_dict)
        return inst

    @classmethod
    def from_cifdict(cls,cifdict):
        ic=cifdict['pdb_ins_code']
        input_dict={
            'model':int(cifdict['pdb_model_num']),
            'resname':cifdict['auth_comp_id'],
            'chainID':cifdict['auth_asym_id'],
            'resseqnum':int(cifdict['auth_seq_id']),
            'insertion':' ' if ic=='?' else ic
        }
        inst=cls(input_dict)
        return inst

    def pdb_line(self):
        return '{:6s}{:>4d}   {:1s} {:3s} {:1s} {:>5d}{:1s}'.format(self.record_name,
        self.code,self.rawmodel,self.rawresname,self.chainID,self.resseqnum,self.insertion)

    def psfgen_lines(self):
        return ['     residue {} {}{} {}'.format(self.resname,self.resseqnum,self.insertion,self.chainID)]

class Crot(CloneableMod):
    req_attr=['angle','degrees']
    opt_attr=['chainID','resseqnum1','resseqnum2','segname','atom1','atom2','segname1','segname2','segnamei','resseqnumi','atomi','segnamejk','resseqnumj','atomj','resseqnumk','atomk']
    attr_choices={'angle':['PHI','PSI','OMEGA','CHI1','CHI2','GLYCAN','LINK','ANGLEIJK']}
    opt_attr_deps={
        'PHI':['chainID','resseqnum1','resseqnum2'],
        'PSI':['chainID','resseqnum1','resseqnum2'],
        'OMEGA':['chainID','resseqnum1','resseqnum2'],
        'CHI1':['chainID','resseqnum1'],
        'CHI2':['chainID','resseqnum1'],
        'GLYCAN':['segname','resseqnum1','atom1','resseqnum2','atom2'],
        'LINK':['segname1','segname2','resseqnum1','atom1','resseqnum2','atom2'],
        'ANGLEIJK':['segnamei','resseqnumi','atomi','segnamejk','resseqnumj','atomj','resseqnumk','atomk']
        }
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
    
class Atom(CloneableMod):
    req_attr=['recordname','serial','name','altloc','resname','chainID','resseqnum','insertion','x','y','z','occ','beta','elem','charge']
    opt_attr=['segname','empty','link']
    yaml_header='Atoms'
    PDB_keyword='ATOM'

    @classmethod
    def from_pdbrecord(cls,pdbrecord):
        input_dict={
            'recordname':pdbrecord[:6],
            'serial':int(pdbrecord[6:11]),
            'name':pdbrecord[12:16].strip(),
            'altloc':pdbrecord[16:17],
            'resname':pdbrecord[17:21].strip(),
            'chainID':pdbrecord[21:22],
            'resseqnum':int(pdbrecord[22:26]),
            'insertion':pdbrecord[26:27],
            'x':float(pdbrecord[30:38]),
            'y':float(pdbrecord[38:46]),
            'z':float(pdbrecord[46:54]),
            'occ':float(pdbrecord[54:60]),
            'beta':float(pdbrecord[60:66]),
            'elem':pdbrecord[76:78].strip(),
            'charge':pdbrecord[78:80].strip()
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

class Hetatm(Atom):
    PDB_keyword='HETATM'

class SSBond(CloneableMod):
    req_attr=['chainID1','resseqnum1','insertion1','chainID2','resseqnum2','insertion2']
    opt_attr=['serial_number','resname1','resname2','sym1','sym2','length']
    yaml_header='SSBonds'
    PDB_keyword='SSBOND'

    @classmethod
    def from_pdbrecord(cls,pdbrecord):
        input_dict={
            'serial_number':int(pdbrecord[7:10]),
            'resname1':pdbrecord[11:14].strip(),
            'chainID1':pdbrecord[15:16].strip(),
            'resseqnum1':int(pdbrecord[17:21]),
            'insertion1':pdbrecord[21:22],
            'resname2':pdbrecord[25:28].strip(),
            'chainID2':pdbrecord[29:30].strip(),
            'resseqnum2':int(pdbrecord[31:35]),
            'insertion2':pdbrecord[35:36],
            'sym1':pdbrecord[59:65].strip(),
            'sym2':pdbrecord[66:72].strip(),
            'length':float(pdbrecord[73:78])
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

class Link(CloneableMod):
    yaml_header='Links'
    PDB_keyword='LINK'
    pass

ModOtherKeywords=['Title','Author','Notes']
ModTypesContainer=[
    Atom, Hetatm, SSBond, Link, Seqadv, Missing, Mutation, Crot
]
class ModTypesRegistry:
    by_yaml={}
    by_pdb_key={}
    for t in ModTypesContainer:
        if 'yaml_header' in t.__dict__:
            by_yaml[t.yaml_header]=t
        if 'PDB_keyword' in t.__dict__:
            by_pdb_key[t.PDB_keyword]=t
    @classmethod
    def modtype(cls,hdr):
        if hdr in cls.by_yaml:
            return cls.by_yaml[hdr]
        if hdr in cls.by_pdb_key:
            return cls.by_pdb_key[hdr]
        return None 
# class ModsContainer:
#     def __init__(self,input_dict):
#         self.M=[]
#         for k,v in input_dict.items():
#             if k in ModTypes:
#                 self.M.append(ModTypes[k](v))
#     def update(self,input_dict):
#         for k,v in input_dict.items():
#             if k in ModTypes:
#                 self.M.append(ModTypes[k](v))
#     def dump(self):
#         retstr=''
#         for m in self.M:
#             retstr+=m.dump()
#         return retstr
