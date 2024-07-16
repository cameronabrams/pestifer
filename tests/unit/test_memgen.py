# Author: Cameron F. Abrams, <cfa22@drexel.edu>

from pestifer.config import Config
import os
from pestifer.util import pdb_charmify

def test_memgen_charmify():
    c=Config(memgen_test=True)
    pdbfile='test.pdb'
    assert os.path.exists(pdbfile)
    df=c['user']['ambertools']['charmmlipid2amber_df']
    assert df.shape[0]==7162
    mdf={'PSM':df[df['search'].str.contains('PSM')]}
    assert len(mdf['PSM'])==127
    pdb_charmify(pdbfile,mdf,outfile='test-out.pdb')

