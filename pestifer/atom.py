#Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Atoms
"""
from .mods import *
from pidibble.baserecord import BaseRecord
from functools import singledispatchmethod
from .util import reduce_intlist

class Atom(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['serial','name','altloc','resname','chainID','resseqnum','insertion','x','y','z','occ','beta','elem','charge']
    opt_attr=AncestorAwareMod.opt_attr+['segname','empty','link','recordname','label_seq_id','label_comp_id','label_asym_id','label_atom_id']
    yaml_header='atoms'
    PDB_keyword='ATOM'
    mmCIF_name='atom_site'

    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)

    @__init__.register(BaseRecord)
    @__init__.register(PDBRecord)
    def _from_pdbrecord(self,pdbrecord):
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
        super().__init__(input_dict)

    @__init__.register(CIFdict)
    def _from_cifdict(self,cifdict):
        input_dict={
            'recordname':'ATOM',
            'serial':int(cifdict['id']),
            'name':cifdict['auth_atom_id'],
            'altloc':cifdict['label_alt_id'],
            'resname':cifdict['auth_comp_id'],
            'chainID':cifdict['auth_asym_id'],
            'resseqnum':int(cifdict['auth_seq_id']),
            'insertion':cifdict['pdbx_pdb_ins_code'],
            'x':float(cifdict['cartn_x']),
            'y':float(cifdict['cartn_y']),
            'z':float(cifdict['cartn_z']),
            'occ':float(cifdict['occupancy']),
            'beta':float(cifdict['b_iso_or_equiv']),
            'elem':cifdict['type_symbol'],
            'charge':cifdict['pdbx_formal_charge'],
            'label_seq_id':int(cifdict['label_seq_id']) if cifdict['label_seq_id'].isdigit() else cifdict['label_seq_id'],
            'label_comp_id':cifdict['label_comp_id'],
            'label_asym_id':cifdict['label_asym_id'],
            'label_atom_id':cifdict['label_atom_id']
        }
        input_dict['segname']=input_dict['chainID']
        input_dict['link']='None'
        input_dict['empty']=False
        super().__init__(input_dict)
    
    def pdb_line(self):
        pdbline='{:<6s}'.format(self.recordname)+\
                '{:5d}'.format(self.serial)+' '+\
                '{:<4s}'.format(' '+self.resname if len(self.resname)<4 else self.resname)+\
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
    
    def overwritePosition(self,other):
        self.x=other.x
        self.y=other.y
        self.z=other.z

class AtomList(AncestorAwareModList):
    def reserialize(self):
        serial=1
        for a in self:
            if not '_ORIGINAL_' in a.__dict__:
                a._ORIGINAL_={}
            a._ORIGINAL_['serial']=a.serial
            a.serial=serial
            serial+=1
            
    def adjustSerials(self,Ters):
        ignored_serials=[x.serial for x in Ters]
        if not ignored_serials:
            return
        logger.debug(f'These serials must be deleted: {ignored_serials}')
        ril=reduce_intlist([x.serial for x in self])
        logger.debug(f'Prior to ignore, serials populate {ril}')
        for a in self:
            try:
                n=next(x[0] for x in enumerate(ignored_serials) if x[1] > a.serial)
            except StopIteration:
                pass
            if n>0:
                if not '_ORIGINAL_' in a.__dict__:
                    a._ORIGINAL_={}
                a._ORIGINAL_['serial']=a.serial
                a.serial-=n
                logger.debug(f'Atom orig serial {a._ORIGINAL_["serial"]} to {a.serial}')

    def overwritePositions(self,other):
        assert len(self)==len(other),'Error: atom lists not equal length'
        for sa,oa in zip(self,other):
            sa.overwritePosition(oa)
    
    def apply_psf_resnames(self,psfatoms):
        for myatom,psfatom in zip(self,psfatoms):
            myatom.resname=psfatom.resname

class Hetatm(Atom):
    PDB_keyword='HETATM'
    yaml_header='hetatoms'
