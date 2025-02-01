# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#
# Manages the collection of PDB files used as inputs for packmol
#
# Collections are subdivided by charmmff 'streams'
# 
import os
class PDBCollection:
    def __init__(self,basepath=''):
        self.basepath=basepath
        if os.path.isdir(self.basepath):
            pass

    def add_usercollection(self,user_path=''):
        pass

    
