"""

.. module:: mods
   :synopsis: Manages all modifications
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import urllib.request
import os
import logging

logger=logging.getLogger(__name__)

BASE_URL='https://files.rcsb.org/download'

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
    

        