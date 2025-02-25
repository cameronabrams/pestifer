# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import logging
import networkx as nx

from .psftopoelement import PSFTopoElement,PSFTopoElementList

from ..config import Config, segtype_of_resname
from ..stringthings import split_ri

logger=logging.getLogger(__name__)

class PSFAtom(PSFTopoElement):
    def __init__(self,atomline):
        tokens=[x.strip() for x in atomline.split()]
        assert len(tokens)==9,f'cannot parse psf atomline: {atomline}'
        idx=int(tokens[0])
        super().__init__([idx])
        r,i=split_ri(tokens[2])
        self.serial=idx
        self.chainID=tokens[1]
        self.resseqnum=r
        self.insertion=i
        self.resname=tokens[3]
        self.name=tokens[4]
        self.type=tokens[5]
        self.charge=float(tokens[6])
        self.atomicwt=float(tokens[7])
        if len(segtype_of_resname)==0:
            c=Config()
        self.segtype=segtype_of_resname[self.resname]
        self.ligands=[]
    
    def __hash__(self):
        return self.serial
    
    def __eq__(self,other):
        return self.serial==other.serial

    def __str__(self):
        return f'{self.chainID}_{self.resname}_{self.resseqnum}{self.insertion}-{self.name}({self.serial})'

    def isH(self):
        return self.atomicwt<1.1
    
    def add_ligand(self,other):
        if not other in self.ligands:
            self.ligands.append(other)

    def is_pep(self,other):
        return (self.type=='C' and other.type=='NH1') or (self.type=='NH1' and other.type=='C')

class PSFAtomList(PSFTopoElementList):

    def __init__(self,atoms):
        super().__init__(atoms)

    def graph(self):
        g=nx.Graph()
        for a in self:
            for l in a.ligands:
                g.add_edge(a,l)
        logger.debug(f'atomlist graph has {len(g)} nodes')
        return g
