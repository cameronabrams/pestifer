# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" Defines the ByteCollector and FileCollector classes
"""
import pandas as pd
from collections import UserList
import logging
logger=logging.getLogger(__name__)
import os
from .command import Command
from itertools import product
# _ANGSTROM_='Ångström'

import importlib.metadata

__pestifer_version__ = importlib.metadata.version("pestifer")

banner_message="""
    {} v. {}
    https://pestifer.readthedocs.io/en/latest/

    Cameron F. Abrams <cfa22@drexel.edu>

    Supported in part by Grants GM100472, AI154071, 
    and AI178833 from the NIH

    CHARMM force field files (July 22) from the 
    MacKerell Lab
    """.format(__package__.title(),__pestifer_version__)

def banner(logf):
    my_logger(banner_message,logf,fill=' ',just='<')

class ByteCollector:
    """A simple string manager
    
    The main object in a ByteCollector instance is a string of bytes (byte_collector).
    The string can be appended to by anoter string or the contents of a file.
    The string can have "comments" written to it.

    Attributes
    ----------
    byte_collector: str
       the string
    
    line_length: int
       number of bytes that represents the maximum length of one line of text
    
    comment_char: str(1)
       beginning-of-line character that signals the line is a comment

    Methods
    -------
    reset()
        blanks the string
    
    write(msg)
        appends msg to string

    addline(msg)
        appends msg to string with a line-end byte
        
    injest_file(filename)
        appends the contents of filename to the string
    
    comment(msg)
        appends the msg as a comment (or multiple comment
        lines) to the string

    log(msg)
        uses the auxiliary function my_logger to write
        a log line to the string
    
    banner(msg)
        uses the auxiliary function my_logger to write
        a banner line to the string
    
    """
    def __init__(self,comment_char='#',line_length=80):
        self.line_length=line_length
        self.comment_char=comment_char
        self.byte_collector=''
    
    def reset(self):
        """Resets the string"""
        self.byte_collector=''

    def write(self,msg):
        """Appends msg to the string
        
        Parameters
        ----------
        msg: str
           the message
        """
        self.byte_collector+=msg

    def addline(self,msg,end='\n'):
        """Appends msg to the string as a line
        
        Parameters
        ----------
        msg: str
           the message
        
        end: str, optional
            end-of-line byte
        """
        self.byte_collector+=f'{msg}{end}'

    def lastline(self,end='\n',exclude='#'):
        """Returns last line in the string
        
        Parameters
        ----------
        end: str, optional
            end-of-line byte
        exclude: str, optional
            comment byte
        """
        lines=[x for x in self.byte_collector.split(end) if (len(x)>0 and not x.startswith(exclude))]
        if len(lines)>0:
            return lines[-1]
        else:
            return None
    
    def has_statement(self,statement,end='\n',exclude='#'):
        """Determines if a particular statement is on at least one non-comment line
        
        Parameters
        ----------
        statement: str
            the statement; e.g., 'exit'
        end: str, optional
            end-of-line byte
        exclude: str, optional
            comment byte
        """
        lines=[x for x in self.byte_collector.split(end) if (len(x)>0 and not x.startswith(exclude))]
        if len(lines)>0:
            for l in lines:
                if statement in l:
                    return True
        return False
    
    def injest_file(self,filename):
        """Appends contents of file 'filename' to the string
        
        Parameters
        ----------
        filename: str
           the name of the file
        """
        with open(filename,'r') as f:
            self.byte_collector+=f.read()

    def comment(self,msg,end='\n'):
        """Appends msg as a comment to the string
        
        Parameters
        ----------
        msg: str
           the message
        
        end: str, optional
            end-of-line byte
        """
        comment_line=f'{self.comment_char} {msg}'
        comment_words=comment_line.split()
        comment_lines=['']
        current_line_idx=0
        for word in comment_words:
            test_line=' '.join(comment_lines[current_line_idx].split()+[word])
            if len(test_line)>self.line_length:
                comment_lines.append(f'{self.comment_char} {word}')
                current_line_idx+=1
            else:
                comment_lines[current_line_idx]=test_line
        for line in comment_lines:
            self.addline(line,end=end)

    def log(self,msg):
        my_logger(msg,self.addline)

    def banner(self,msg):
        my_logger(msg,self.addline,fill='#',width=80)

    def __str__(self):
        return self.byte_collector

def my_logger(msg,logf,width=67,fill='*',sep=', ',just='^',frame=''):
    """A fancy logger
    
    Parameters
    ----------
    msg: str, list
       the message to be logged, either as a single string or a list of strings
    
    logf: function
       writer; e.g., print, f.write, etc.

    width: int, optional
       linelength in bytes

    fill: str, optional
       single character used to fill blank spaces

    sep: str, optional
       single character used in join calls

    just: str, optional
       format character
    """
    fmt=r'{'+r':'+fill+just+f'{width}'+r'}'
    if frame:
        ffmt=r'{'+r':'+frame+just+f'{width}'+r'}'
    ll=' ' if just in ['^','>'] else ''
    rr=' ' if just in ['^','<'] else ''
    if frame:
        logf(ffmt.format(frame))
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
    if frame:
        logf(ffmt.format(frame))
            
class FileCollector(UserList):
    """A class for handling collections of files
    
    Methods
    -------
    flush()
       remove all files in the collection
       
    tarball()
       make a tarball of the collection
    """
    def flush(self):
        logger.debug(f'Flushing file collector: {len(self)} files.')
        for f in self:
            if os.path.exists(f):
                os.remove(f)
            else:
                logger.debug(f'{f}: not found.')
        self.clear()
    def tarball(self,basename):
        """Makes a tarball of the files in the collection
        
        Parameters
        ----------
        basename: str
            basename of the resulting tarball
        """
        filelist=' '.join([x for x in self])
        c=Command(f'tar zvcf {basename}.tgz {filelist}')
        c.run()
        logger.debug(f'generated tarball {basename}.tgz')

def split_ri(ri):
    """A simple utility function for splitting the integer resid and
    1-byte insertion code out of a string resid-insertion code
    concatenation
    
    Parameters
    ----------
    ri: the supposed resid-insertion concatenation
    
    Returns
    -------
    tuple(int, str): the integer resid and the 1-byte insertion code or '' if none
    """
    if ri[-1].isdigit(): # there is no insertion code
        r=int(ri)
        i=''
    else:
        r=int(ri[:-1])
        i=ri[-1]
    return r,i

def join_ri(resseqnum,insertion):
    if insertion=='':
        return resseqnum
    return f'{resseqnum}{insertion}'

def ri_range(val,split_chars=['-','#']):
    the_split=[val]
    for c in split_chars:
        the_splits=[x.split(c) for x in the_split]
        the_split=[]
        for s in the_splits:
            the_split.extend(s)
    return [split_ri(x) for x in the_split]

def get_boxsize_from_packmolmemgen(logname='packmol-memgen.log'):
    res=''
    extract_keys=[''.join(x) for x in product(['x_','y_','z_'],['min','max','len'])]
    boxinfo={}
    with open(logname,'r') as f:
        log=f.read().split('\n')
        for l in log:
            tok=l.split()
            if len(tok)>0:
                if tok[0] in extract_keys:
                    boxinfo[tok[0]]=float(tok[2])
    if boxinfo:
        ox=0.5*(boxinfo["x_max"]-boxinfo["x_min"])
        oy=0.5*(boxinfo["y_max"]-boxinfo["y_min"])
        oz=0.5*(boxinfo["z_max"]-boxinfo["z_min"])
        res =f'cellBasisVector1 {boxinfo["x_len"]} 0 0\n'
        res+=f'cellBasisVector2 0 {boxinfo["y_len"]} 0\n'
        res+=f'cellBasisVector3 0 0 {boxinfo["z_len"]}\n'
        res+=f'cellOrigin 0 0 0'
        # res+=f'cellOrigin {ox} {oy} {oz}'
    return res,boxinfo