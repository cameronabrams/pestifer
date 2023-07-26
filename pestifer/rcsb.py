"""

.. module:: rcsb
   :synopsis: Manages all downloading and parsing of data from RCSB
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import urllib.request
import os
import logging
import yaml

from .mods import ModTypesRegistry as mtr
from pestifer import Resources

logger=logging.getLogger(__name__)

BASE_URL='https://files.rcsb.org/download'


def PDB_List(arg,d=','):
    return [x.strip() for x in arg.split(d)]

def func(obj, param, param_values):
    setattr(obj, param, param_values)
    return obj

class PDBRecord:

    def __init__(self,input_dict):
        self.__dict__.update(input_dict)

    @classmethod
    def base_parse(cls,pdbrecord,record_format,typemap):
        input_dict={}
        fields=record_format.get('fields',{})
        for k,v in fields.items():
            typestring,byte_range=v
            typ=typemap[typestring]
            # print(typestring,typ)
            input_dict[k]=typ(pdbrecord[byte_range[0]-1:byte_range[1]])
            if typ==str:
                input_dict[k]=input_dict[k].strip()
        return input_dict

    @classmethod
    def new(cls,pdbrecord,record_format,typemap):
        input_dict=cls.base_parse(pdbrecord,record_format,typemap)
        if input_dict:
            inst=cls(input_dict)
            return inst
        return None

    def update_2(self,pdbrecord,record_format,typemap):
        """ type-2 update of existing parsed record 
            if there is only one field in addition to
            'continuation', then it is assumed this field is a strict continuation; if not, each line is a distinct record """
        if not 'continuation' in record_format['fields']:
            logging.fatal(f'A type-2 pdb record must have a continuation field')
        input_dict=PDBRecord.base_parse(pdbrecord,record_format,typemap)
        for k in input_dict.keys():
            if not type(self.__dict__[k])==list:
                self.__dict__[k]=[self.__dict__[k],input_dict[k]]
            else:
                self.__dict__[k].append(input_dict[k])



class PDBParser:
    previous_key=None
    # keys_encountered=[]
    parsed={}
    mappers={'Integer':int,'String':str,'Float':float,'List':PDB_List}
    def __init__(self,**options):
        self.pdb_code=options.get('PDBcode','')
        self.overwrite=options.get('overwrite',False)
        pdb_format_file=os.path.join(
            os.path.dirname(Resources.__file__),
            'config','pdb_format.yaml')
        if os.path.exists(pdb_format_file):
            with open(pdb_format_file,'r') as f:
                self.pdb_format_dict=yaml.safe_load(f)
        for map,d in self.pdb_format_dict['delimiters'].items():
            if not map in self.mappers:
                self.mappers[map]=func(PDB_List,'d',d)
        # print(self.mappers)
    # def get_record_format(self,key):
    #     for rf in self.pdb_format_dict['record_formats']:
    #         if rf['key']==key:
    #             return rf
    #     return {}
            
    def fetch(self):
        self.filename=f'{self.pdb_code}.pdb'
        target_url=os.path.join(BASE_URL,self.filename)
        if not os.path.exists(self.filename) or self.overwrite:
            try:
                urllib.request.urlretrieve(target_url,self.filename)
            except:
                logger.warning(f'Could not fetch {self.filename}')

    def read(self):
        self.pdb_lines=[]
        with open(self.filename,'r') as f:
            self.pdb_lines=f.read().split('\n')
            if self.pdb_lines[-1]=='':
                self.pdb_lines=self.pdb_lines[:-1]
    def close_record(self,key):
        pass
    def parseline(self,line):
        key=line[:6].strip()
        record_format=self.pdb_format_dict['record_formats'][key]

    def parse(self,**options):
        for l in self.pdb_lines:
            key=l[:6].strip()
            record_format=self.pdb_format_dict['record_formats'][key]
            # print(key,record_format)
            record_type=record_format['type']
            if record_type==1: # one-time-single-line
                # print(key,record_format)
                if key in self.parsed:
                    logger.fatal(f'PDB record type {key} can appear only once in a PDB file.')
                self.parsed[key]=PDBRecord.new(l,record_format,self.mappers)
            elif record_type==2:
                if key!=self.previous_key:
                    self.parsed[key]=PDBRecord.new(l,record_format,self.mappers)
                    self.previous_key=key
                else:
                    self.parsed[key].update_2(l,record_format,self.mappers)
            

class MolData:
    captures=['ATOM','LINK','SSBOND','SEQADV','REMARK,465']
    _cloneables=['ATOM','LINK','SSBOND','SEQADV']
    def __init__(self,**options):
        self.Atoms=[]
        self.Seqadv=[]
        self.Links=[]
        self.SSBonds=[]
        self.BrokenSSBonds=[]
        self.Missing=[]
        self.Mutations=[]

        # initialize all mod lists with any user-input options
        if 'mods' in options:
            # mods that can be specified OUTSIDE of PDB files
            for mod in options['mods']:
                for k,v in mod.items():
                    # one or more keys specifies a list of mods
                    T=mtr.modtype(k)
                    if T:
                        for tm in v:
                            if type(tm)==str:
                                self.__dict__[k].append(T.from_shortcode(tm))
                            elif type(tm)==dict:
                                self.__dict__[k].append(T(tm))

    def parse_line(self,key,pdbrecord):
        Mod=mtr.modtype(key)
        if Mod:
            self.__dict__[Mod.yaml_header]=Mod.from_pdbrecord(pdbrecord)

class MolMeta:
    captures=['TITLE','KEYWDS','MASTER','HEADER','MDLTYP','DBREF','SEQRES','REVDAT','EXPDTA']
    def __init__(self):
        self.Title=''
        self.titlelines=[]
        self.keywords=[]
    
    def parse_line(self,key,pdbrecord):
        if key=='TITLE':
            self.parse_title(pdbrecord)
        elif key == 'KEYWDS':
            self.parse_keywords(pdbrecord)
        elif key == 'MASTER':
            self.ParseMasterRecord(pdbrecord)
        elif key == 'HEADER':
            self.ParseHeader(pdbrecord)
        elif key == 'MDLTYP':
            self.ParseModelType(pdbrecord)
        elif key == 'DBREF':
            self.ParseDBRef(pdbrecord)
        elif key == 'SEQRES':
            self.ParseSeqRes(pdbrecord)
        elif key == 'REVDAT':
            self.ParseRevisionDate(pdbrecord)
        elif key == 'EXPDTA':
            self.ParseExpDta(pdbrecord)
    def parse_title(self,pdbrecord):
        short=pdbrecord[10:80].strip()
        self.titlelines.append(short)
        self.Title=short if len(self.Title)==0 else (self.Title+short)
    def parse_keywords(self,pdbrecord):
        ctok=pdbrecord[9:10]
        kwdlist=pdbrecord[10:80].split(',')
        for i,k in enumerate(kwdlist):
            if i==0 and ctok.isdigit():
                pk=self.keywords[-1]
                nk=pk+' '+k.strip()
                self.keywords[-1]=nk
            else:
                self.keywords.append(k.strip())
    def generate_keywords_record(self):
        retstr=''
        if len(self.keywords)>0:
           cpst=', '.join(self.keywords)+' '
           i=0
           sp=[]
           lnlim=69
           while i<len(cpst):
               if cpst[i]==' ':
                   sp.append(i)
               i = i + 1
           nln=len(cpst)//lnlim+1
           beg=0
           for j in range(nln):
               end=-1
               # find rightmost ' ' such that length < lnlim
               for i in range(len(sp)-1,-1,-1):
                   if sp[i]-beg<69:
                      end=sp[i]
                      break
               retstr+='KEYWDS  {:>2s} {}\n'.format(' ' if j==0 else str(j+1),cpst[beg:end])
               beg=end+1
        return retstr

# class PDBRecord:
#     lines=[]
#     def __init__(self,pdbrecord:str):
#         key=pdbrecord[:6].strip()
#         subkey=pdbrecord[7:10].strip()
#         assert (subkey=='' or subkey.isdigit())




        