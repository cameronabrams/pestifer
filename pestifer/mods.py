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
    req_attr=[]
    opt_attr=[]
    alt_attr=[]
    def __init__(self,input_dict={}):
        for k,v in input_dict.items():
            if k in self.req_attr or k in self.opt_attr:
                self.__dict__[k]=v
        assert all([att in self.__dict__.keys() for att in self.req_attr])
        assert all([(x[0] in self.__dict__.keys() or x[1] in self.__dict__.keys()) for x in self.alt_attr])
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
    alt_attr=[('resseqnum','insertion')]

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

class Atom(CloneableMod):
    req_attr=['recordname','serial','name','altloc','resname','chainID','resseqnum','insertion','x','y','z','occ','beta','elem','charge']
    opt_attr=['segname','empty','link']
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

class SSBond(CloneableMod):
    req_attr=['chainID1','resseqnum1','insertion1','chainID2','resseqnum2','insertion2']
    opt_attr=['serial_number','resname1','resname2','sym1','sym2','length']

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
    pass

ModTypes={
    'Mutations':Mutation,
    'Seqadv':Seqadv,
    # 'Grafts':Graft,
    # 'Deletions':Deletion,
    # 'Crotations':Crot,
    # 'Attachments':Attach,
    # 'Links':Link,
    'SSbonds':SSBond,
    # 'Cleavages':Cleavage,
    'Missing':Missing,
    'Atoms':Atom
}

class ModsContainer:
    def __init__(self,input_dict):
        self.M=[]
        for k,v in input_dict.items():
            if k in ModTypes:
                self.M.append(ModTypes[k](v))
    def update(self,input_dict):
        for k,v in input_dict.items():
            if k in ModTypes:
                self.M.append(ModTypes[k](v))
    def dump(self):
        retstr=''
        for m in self.M:
            retstr+=m.dump()
        return retstr
