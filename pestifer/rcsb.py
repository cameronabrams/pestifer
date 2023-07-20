import urllib.request
import os
import logging

logger=logging.getLogger(__name__)

BASE_URL='https://files.rcsb.org/download'

def get_pdb_file(pdb_code):
    filename=f'{pdb_code}.pdb'
    target_url=os.path.join(BASE_URL,filename)
    try:
        urllib.request.urlretrieve(target_url, filename)
    except:
        logger.warning()
    return filename
    