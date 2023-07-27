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

class Stringparser:
    def __init__(self,fmtdict,typemap):
        self.typemap=typemap
        self.fields={k:v for k,v in fmtdict.items()}
    def parse(self,record):
        input_dict={}
        total_bytes=0
        for k,v in self.fields.items():
            typestring,byte_range=v
            total_bytes+=byte_range[1]-byte_range[0]+1
        while len(record)<=total_bytes:
                record+=' '
        # print(f'({record})',len(record))
        for k,v in self.fields.items():
            typestring,byte_range=v
            typ=self.typemap[typestring]
            assert byte_range[1]<=len(record)
            # using columns beginning with "1" not "0"
            fieldstring=record[byte_range[0]-1:byte_range[1]]
            # print(k,f'({fieldstring})')
            fieldstring=fieldstring.strip()
            # print(typestring,typ)
            input_dict[k]='' if fieldstring=='' else typ(fieldstring)
            if typ==str:
                input_dict[k]=input_dict[k].strip()
        return BaseRecord(input_dict)

# def isempty(adict):
#     isempty=True
#     for v in adict.values():
#         isempty&=(v=='')
#     return isempty

def stringparse(fmtdict,typemap):
    return Stringparser(fmtdict,typemap).parse

class BaseRecord:
    def __init__(self,input_dict):
        self.__dict__.update(input_dict)
    def empty(self):
        isempty=True
        for v in self.__dict__.values():
            isempty&=(v=='')
        return isempty
    def __str__(self):
        return ';'.join([f'{k}:[{v}]' for k,v in self.__dict__.items()])

class PDBRecord(BaseRecord):
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
            fieldstring=pdbrecord[byte_range[0]-1:byte_range[1]]
            # print(typestring,typ)
            input_dict[k]='' if fieldstring.strip()=='' else typ(fieldstring)
            if type(input_dict[k])==str:
                input_dict[k]=input_dict[k].strip()
            # if typ==str:
            #     input_dict[k]=input_dict[k].strip()
            if typestring in allowed_values:
                assert input_dict[k] in allowed_values[typestring]

        return input_dict

    @classmethod
    def new(cls,pdbrecord,record_format,typemap,**kwargs):
        input_dict=cls.base_parse(pdbrecord,record_format,typemap)
        concats=record_format.get('concatenate',{})
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
            for cfield,subf in concats.items():
                if not cfield in input_dict:
                    input_dict[cfield]=[]
                    for f in subf:
                        # print(f,input_dict[f])
                        if input_dict[f]:
                            input_dict[cfield].append(input_dict[f])
                        # del input_dict[f]
            subrecords=record_format.get('subrecords',{})
            for srf,srecfmtwkey in subrecords.items():
                assert type(input_dict[srf])==str
                typestring,byte_range=srecfmtwkey['key']
                typ=typemap[typestring]
                fieldstring=pdbrecord[byte_range[0]-1:byte_range[1]].strip()
                key=typ(fieldstring)
                srf=srecfmtwkey['subrecord_formats'][key]
                print(key,srf.get('continues',[]))
                if not key in input_dict:
                    input_dict[key]=PDBRecord.new(pdbrecord,srf,typemap)
                else:
                    newrec=PDBRecord.new(pdbrecord,srf,typemap)
                    input_dict[key].merge(newrec,srf)
            inst=cls(input_dict)
            return inst
        return None

    def update_2(self,pdbrecord,record_format,typemap):
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
    def merge(self,other,record_format):
        continues=self.__dict__ if (not 'continues' in record_format or not record_format['continues']) else record_format['continues']
        # print(continues)
        # print(self.__dict__)
        for cfield in continues:
            if type(self.__dict__[cfield])==list:
                if type(other.__dict__[cfield])==list:
                    self.__dict__[cfield].extend(other.__dict__[cfield])
                else:
                    self.__dict__[cfield].append(other.__dict__[cfield])
            elif type(self.__dict__[cfield])==str: # assume a string
                # print(self.__dict__[cfield],other.__dict__[cfield])
                self.__dict__[cfield]+=' '+other.__dict__[cfield]
            elif type(self.__dict__[cfield])==int:
                self.__dict__[cfield]=[self.__dict__[cfield],other.__dict__[cfield]]
            else:
                self.__dict__[cfield]=[self.__dict__[cfield],other.__dict__[cfield]]
        # concats=record_format.get('concatenate',{})
        # for concatfield,subfields in concats.items():
        #     if not concatfield in self.__dict__.keys():
        #         self.__dict__[concatfield]=[]
            
class PDBParser:
    previous_key=None
    previous_record_format={}
    last_parsed_3=None
    last_parsed_5=None
    # keys_encountered=[]
    parsed={}
    mappers={'Integer':int,'String':str,'Float':float}
    mappers.update(ListParsers)
    comment_lines=[]
    comment_chars=['#']
    def __init__(self,**options):
        self.pdb_code=options.get('PDBcode','')
        # print(self.pdb_code)
        self.overwrite=options.get('overwrite',False)
        pdb_format_file=options.get('pdb_format_file',os.path.join(
            os.path.dirname(Resources.__file__),
            'config','pdb_format.yaml'))
        if os.path.exists(pdb_format_file):
            with open(pdb_format_file,'r') as f:
                self.pdb_format_dict=yaml.safe_load(f)
        delimiter_dict=self.pdb_format_dict.get('delimiters',{})
        for map,d in delimiter_dict.items():
            if not map in self.mappers:
                self.mappers[map]=Listparser(d).parse
        cformat_dict=self.pdb_format_dict.get('custom_formats',{})
        for cname,cformat in cformat_dict.items():
            if not cname in self.mappers:
                self.mappers[cname]=Stringparser(cformat,PDBParser.mappers).parse
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

    def parse(self,**options):
        for i,l in enumerate(self.pdb_lines):
            key=l[:6].strip()
            if key[0] in PDBParser.comment_chars:
                self.comment_lines.append([i,l])
                continue
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
                if parsed_record!=self.last_parsed_3:
                    self.parsed[key].append(parsed_record)
                self.last_parsed_3=self.parsed[key][-1]
            elif record_type==4:
                # - continuation requires which field continues
                # - serNum requires matches
                # - seqNum is added to keyname
                parsed_record=PDBRecord.new(l,record_format,self.mappers)
                if not key in self.parsed:
                    self.parsed[key]=[parsed_record]
                else:
                    if 'determinants' in record_format:
                        det=[parsed_record.__dict__[x] for x in     record_format['determinants']]
                        present_record=None
                        for r in self.parsed[key]:
                            testdet=[r.__dict__[x] for x in record_format['determinants']]
                            if det==testdet:
                                present_record=r
                                # print(present_record.__dict__)
                                break
                        if not present_record:
                            self.parsed[key].append(parsed_record)
                        else:
                            present_record.merge(parsed_record,record_format)
                    else:
                        self.parsed[key].append(parsed_record)
            elif record_type==5:
                if not key in self.parsed:
                    self.parsed[key]=[]
                    self.last_parsed_5=None
                parsed_record=PDBRecord.new(l,record_format,self.mappers,continue_on=self.last_parsed_5)
                if parsed_record!=self.last_parsed_5:
                    self.parsed[key].append(parsed_record)
                self.last_parsed_5=self.parsed[key][-1]
            elif record_type==6:
                parsed_record=PDBRecord.new(l,record_format,self.mappers)
                if not key in self.parsed:
                    self.parsed[key]=[]
                self.parsed[key].append(parsed_record)
            

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




        