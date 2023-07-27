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

class Listparser:
    def __init__(self,d=','):
        self.d=d
    def parse(self,string):
        if self.d==None:
            return [x for x in string.split() if x.strip()!='']
        else:
            return [x.strip() for x in string.split(self.d) if x.strip()!='']
    
def list_parse(obj,d):
    return obj(d).parse

ListParsers={
    'CList':list_parse(Listparser,','),
    'SList':list_parse(Listparser,';'),
    'WList':list_parse(Listparser,None),
    'DList':list_parse(Listparser,':')
}

class PDBRecord:

    def __init__(self,input_dict):
        self.__dict__.update(input_dict)

    @classmethod
    def base_parse(cls,pdbrecord,record_format,typemap):
        while len(pdbrecord)<80:
            pdbrecord+=' '
        input_dict={}
        fields=record_format.get('fields',{})
        allowed_values=record_format.get('allowed',{})
        for k,v in fields.items():
            typestring,byte_range=v
            typ=typemap[typestring]
            assert byte_range[1]<=len(pdbrecord)
            # using columns beginning with "1" not "0"
            fieldstring=pdbrecord[byte_range[0]-1:byte_range[1]].strip()
            # print(typestring,typ)
            input_dict[k]='' if fieldstring=='' else typ(fieldstring)
            if typ==str:
                input_dict[k]=input_dict[k].strip()
            if typestring in allowed_values:
                assert input_dict[k] in allowed_values[typestring]
        return input_dict

    @classmethod
    def new(cls,pdbrecord,record_format,typemap,**kwargs):
        input_dict=cls.base_parse(pdbrecord,record_format,typemap)
        if input_dict:
            continue_on=kwargs.get('continue_on',None)
            if continue_on:
                continuation=input_dict.get('continuation','')
                if continuation!='':
                    # this is a continuation record
                    for attr in continue_on.__dict__.keys():
                        cont_attr=input_dict.get(attr,'')
                        if cont_attr!='':
                            if type(continue_on.__dict__[attr])!=list:
                                continue_on.__dict__[attr]=[continue_on.__dict__[attr]]
                            if type(cont_attr)==list:
                                continue_on.__dict__[attr].extend(cont_attr)
                            else:
                                continue_on.__dict__[attr].append(cont_attr)
                    return continue_on
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
            # just add to the end of the string
            if type(self.__dict__[k])==str:
                self.__dict__[k]+=' '+input_dict[k]
            else:
                if not type(self.__dict__[k])==list:
                    self.__dict__[k]=[self.__dict__[k]]
                if type(input_dict[k])==list:
                    self.__dict__[k].extend(input_dict[k])
                else:
                    self.__dict__[k].append(input_dict[k])

    def parse_tokens(self,record_format,typemap):
        if not 'tokens' in record_format:
            return
        attr_w_tokens=record_format['tokens']
        # print(attr_w_tokens)
        self.tokens={}
        current_parent_toknames={}
        for a in attr_w_tokens.keys():
            obj=self.__dict__[a] # expect to be a list
            assert type(obj)==list
            self.tokens[a]={}
            tdict=attr_w_tokens[a]
            for pt in self.__dict__[a]:
                toks=[x.strip() for x in pt.split(':')]
                assert len(toks)==2,f'Malformed token string {pt}'
                tokname,tokvalue=[x.strip() for x in pt.split(':')]
                assert tokname in tdict.keys(),f'Unrecognized token {tokname}'
                typ=typemap[tdict[tokname]['type']]
                if 'associated_to' in tdict[tokname]:
                    parent_tokname=tdict[tokname]['associated_to']
                    asstokname=current_parent_toknames[parent_tokname]
                    self.tokens[a][asstokname][tokname]=typ(tokvalue)

                else:
                    tokvalue=typ(tokvalue)
                    ntokname=f'{tokname}.{tokvalue}'
                    self.tokens[a][ntokname]={}
                    current_parent_toknames[tokname]=ntokname

class PDBParser:
    previous_key=None
    previous_record_format={}
    last_parsed_3=None
    # keys_encountered=[]
    parsed={}
    mappers={'Integer':int,'String':str,'Float':float}
    mappers.update(ListParsers)
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
                self.mappers[map]=Listparser(d).parse
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
                    if self.previous_key:
                        # print(self.previous_key)
                        self.parsed[self.previous_key].parse_tokens(self.previous_record_format,self.mappers)
                    self.previous_key=key
                    self.previous_record_format=record_format
                else:
                    self.parsed[key].update_2(l,record_format,self.mappers)
            elif record_type==3:
                if not key in self.parsed:
                    self.parsed[key]=[]
                    self.last_parsed_3=None
                parsed_record=PDBRecord.new(l,record_format,self.mappers,continue_on=self.last_parsed_3)
                if parsed_record==self.last_parsed_3:
                    pass
                else:
                    self.parsed[key].append(parsed_record)
                self.last_parsed_3=self.parsed[key][-1]
            

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




        