"""

.. module:: banner
   :synopsis: Defines the banner method for printing the banner to a logging channel
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
from .stringthings import my_logger
from .util import get_version

banner_message="""
    Pestifer {:s}
    https://github.com/cameronabrams/pestifer

    Cameron F. Abrams
    cfa22@drexel.edu

    Supported in part by Grants GM100472, AI15407, 
    and AI178833 from the NIH

    CHARMM force field files from the MacKerell Lab
    July 22 update
    """.format(get_version())
def banner(logf):
    my_logger(banner_message,logf,fill=' ',just='<')