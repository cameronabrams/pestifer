"""

.. module:: stringthings
   :synopsis: defines a byte-collector for convenient script writing and a file-collector for convenient clean-up
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import pandas as pd
import logging
logger=logging.getLogger(__name__)
import os
# _ANGSTROM_='Ångström'
class ByteCollector:
    def __init__(self,comment_char='#'):
        self.comment_char=comment_char
        self.byte_collector=''
    def reset(self):
        self.byte_collector=''
    def write(self,msg):
        self.byte_collector+=msg
    def addline(self,msg,end='\n'):
        self.byte_collector+=f'{msg}{end}'
    def comment(self,msg,end='\n'):
        comment_line=f'{self.comment_char} {msg}'
        self.addline(comment_line,end=end)
    def log(self,msg):
        my_logger(msg,self.addline)
    def banner(self,msg):
        my_logger(msg,self.addline,fill='#',width=80)
    def __str__(self):
        return self.byte_collector

def my_logger(msg,logf,width=67,fill='*',sep=', ',just='^'):
    fmt=r'{'+r':'+fill+just+f'{width}'+r'}'
    ll=' ' if just in ['^','>'] else ''
    rr=' ' if just in ['^','<'] else ''
    if type(msg)==list:
        rr=' ' if ' ' not in sep else ''
        lnlst=[]
        for tok in msg:
            ll_sofar=sum([len(x) for x in lnlst])
            test_ll=ll_sofar+len(tok)+len(sep)*(len(lnlst)+1)
            if test_ll>width-2*(1+len(sep)):
                outstr=ll+sep.join(lnlst)+sep+rr
                logf(fmt.format(outstr))
                lnlst=[tok]
            else:
                lnlst.append(tok)
        outstr=ll+sep.join(lnlst)+' '
        logf(fmt.format(outstr))
    elif type(msg)==pd.DataFrame:
        for ln in msg.to_string().split('\n'):
            outstr=ll+ln+rr
            logf(fmt.format(outstr))
    else:
        lns=msg.split('\n')
        for ln in lns:
            outstr=ll+ln+rr
            logf(fmt.format(outstr))
        
class FileCollector(list):
    def __init__(self):
        self.file_collection=[]
    def __iter__(self):
        for f in self.file_collection:
            yield f
    def append(self,item):
        self.file_collection.append(item)
    def extend(self,a_list):
        self.file_collection.extend(a_list)
    def flush(self):
        logger.debug(f'Flushing file collector: {len(self.file_collection)} files.')
        for f in self:
            if os.path.exists(f):
                os.remove(f)