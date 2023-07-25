"""

.. module:: rcsb
   :synopsis: Manages all downloading and parsing of data from RCSB
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import urllib.request
import os
import logging

from .mods import ModTypesRegistry as mtr

logger=logging.getLogger(__name__)

BASE_URL='https://files.rcsb.org/download'

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

class PDBParser:
    parseables=[]
    def __init__(self,**options):
        self.pdb_code=options.get('PDBcode','')
        self.pdb_lines=[]

    def fetch(self):
        filename=f'{self.pdb_code}.pdb'
        target_url=os.path.join(BASE_URL,filename)
        self.pdb_lines=[]
        try:
            urllib.request.urlretrieve(target_url,filename)
            with open(filename,'r') as f:
                self.pdb_lines=f.read().split('\n')
                if self.pdb_lines[-1]=='':
                    self.pdb_lines=self.pdb_lines[:-1]
        except:
            logger.warning(f'Could not fetch {filename}')

        return filename
    
    def parse(self,**options):
        d=MolData(mods=options.get('mods',{}))
        m=MolMeta()
        for l in self.pdb_lines:
            key=l[:6].strip()
            if key=="REMARK":
                token=l[7:10].strip()
                if token.isdigit():
                    key=f'{key},{token}'
            if key in MolData.captures:
                d.parse_line(key,l)
            elif key in MolMeta.captures:
                m.parse_line(key,l)
            else:
                logger.info(f'pdbrecord "{key}" is not parseable')

        return m,d

        