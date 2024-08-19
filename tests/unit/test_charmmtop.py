# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import unittest
from pestifer.config import ResourceManager
from pestifer.charmmtop import getResis
import glob
import os
import logging
logger=logging.getLogger(__name__)

class TestCHARMMtop(unittest.TestCase):

    def test_lipid_charmmtop(self):
        r=ResourceManager()
        toppardir=r['charmmff']['toppar']
        lipidstreamdir=os.path.join(toppardir,'stream/lipid')
        strlipids=glob.glob(os.path.join(lipidstreamdir,'*.str'))
        idx=['list' in x for x in strlipids].index(True)
        strlipids.remove(strlipids[idx])
        # logger.debug(f'strlipids in {", ".join([os.path.basename(x) for x in strlipids])}')
        tops=[os.path.join(toppardir,'top_all36_lipid.rtf')]+glob.glob(os.path.join(lipidstreamdir,'*.str'))

        resis=[]
        for t in tops:
            sublist=getResis(t)
            # logger.debug(f'{os.path.basename(t)} {len(sublist)} RESIs')
            resis.extend(sublist)
        logger.debug(f'{len(resis)} RESIs')
        for r in resis:
            logger.debug(f'{r.resname}, {r.synonym}, {r.num_atoms()}, {os.path.basename(r.topfile)}')
