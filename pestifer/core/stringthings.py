# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
Defines the :class:`ByteCollector` and :class:`FileCollector` classes, and defines the :func:`banner` and :func:`my_logger` functions.
"""
import logging
import os
import re
import subprocess
import shutil

import pandas as pd

from collections import UserList
from io import StringIO

logger=logging.getLogger(__name__)

import importlib.metadata

__pestifer_version__ = importlib.metadata.version("pestifer")

_banner_message="""
    pestifer v. {}
    https://pestifer.readthedocs.io/en/latest/

    Cameron F. Abrams <cfa22@drexel.edu>

    Supported in part by Grants GM100472, AI154071, 
    and AI178833 from the NIH

    CHARMM force field files (July 24) from the 
    MacKerell Lab
    """.format(__pestifer_version__)

_enhanced_banner_message="""
░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░░░░░░░░░░░░░▒███████████▓░░░░░░░░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░░░░░░░░░░░█████████████████▒░░░░░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░░░░░░░░░▓████████████████████░░░░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░░░░░░░░▓▓▓████████▒████████▓▓█░░░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░░░░░░░▓▒▒▒░▒████▒▒█▓░████▓░░▒▒█░░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░░░░░░▒█▒░▒██▓▓▓██▒░░███▓▓███░▒█▓░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░░░░░░▓▒░░█▓░░░░░░█▓▓▓░░░░░░▓▓░▒█░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░░░░░░▓▒░█▒░░░░░░░███▒░░░░░░░█░▒▓░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░░░░░░▒░░█▒░░░░░▒█▒░░██░░░░░░█░░▓░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░░░░░░░░░██████▒▓▒░░░░█▒██████░░░░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░░░░░░░█▓▒▓██░▒██░░▒░░███░▓█▓░▓█░░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░░░░░░░░▒░░░░░░▒███▓████░░░░░░▒░░░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░▒░▒░▒░▒░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░▒░▒░▒░▒░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░█████▒░▒█████▒░▓█████░▓██████░██▒▒█████░██████░██████░░░░░░░░
░░░░░░░░░██░░██▒░██░░█▒▒█▓░░▒█░█░░██░▓░██░░██░░▓▒▓█▒░░█░██▒░▒██░░░░░░░
░░░░░░░░░██░░▒██░██░░░▒▒██░░░▒▒░░░██░░░██░░██░░░▒▓█▒░░░░██▒░░██▒░░░░░░
░░░░░░░░░██░░▒██░██░▒▓░░▓███░░░░░░██░░░██░░██░▒█░▓█▒░▒▓░██▒░░██░░░░░░░
░░░░░░░░░██░░██▒░████▓░░░░▓███░░░░██░░░██░░██▓██░▓████▓░██▒▒██▒░░░░░░░
░░░░░░░░░██▓██░░▒██░░▓░▓░░░░▒██░░░██░░░██░░██░░▓░▓█▒░░█░█████░░░░░░░░░
░░░░░░░░░██▒░░░░▒██░░░░▒█░░░░▓█▒░░██░░░██░░██░░░░▓█▒░░░░██▒▓█▒░░░░░░░░
░░░░░░░░░██░░░░░▒█▓░░▓▓░██░░▓██░░░██░░░██░░██░░░░▓█▒░░▒▒██▒▒██░░░░░░░░
░░░░░░░░░██░░░░░▒██▒▓█▒░▒████▒░░░▓▓▓▒░▒▓▓▒▓██▒░░░██▓▒▓█░▓█▒░██▒░░░░░░░
░░░░░░░░░██░░░░░▓███▓▒░░░░░░░░░░░░░░░░░░░░░░░░░░░░▒▓███░██▒░▒█▓░░░░░░░
░░░░░░░░░▓█▒░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░▒▓█░░██▒░░░░░░
░░░░░░░░▓▓▒░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░▓█░░░░░░
░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░ https://pestifer.readthedocs.io/en/latest/ ░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░ Cameron F. Abrams, cfa22@drexel.edu ░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
"""

def banner(logf):
    """
    Writes a banner message to the log file

    Parameters
    ----------
    logf: file-like object
        The log file to which the banner message will be written.
    """
    my_logger(_banner_message,logf,fill=' ',just='<')

class ByteCollector:
    """
    A simple string manager.
    The main object in a ``ByteCollector`` instance is a string of bytes (``byte_collector``).
    The string can be appended to by another string or the contents of a file.
    The string can have comments written to it.
    The string can be written to a file.

    Parameters
    ----------
    comment_char: str, optional
        The character used to denote comments in the string.
        Defaults to '#'.
    line_length: int, optional
        The maximum length of a line in the string.
        Defaults to 80.
    """
    def __init__(self,comment_char='#',line_length=80):
        self.line_length=line_length
        self.comment_char=comment_char
        self.byte_collector=''
    
    # def __len__(self):
    #     """
    #     Returns the length of the byte collector string
    #     """
    #     return len(self.byte_collector)

    def reset(self):
        """
        Resets the string
        """
        self.byte_collector=''

    def write(self,msg):
        """
        Appends msg to the string
        
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
           the end-of-line byte
        """
        self.byte_collector+=f'{msg}{end}'

    def lastline(self,end='\n',exclude='#'):
        """
        Returns last line in the string

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
        """
        Determines if a particular statement is on at least one non-comment line

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

    def reassign(self,varname,varvalue,style='bash',end='\n'):
        """
        Reassigns variable varname to value varvalue if an assignment to varname already
        exists in the byte_collector 
        
        Parameters
        ----------
        varname: str
            the name of the variable to reassign
        varvalue: str
            the value to assign to the variable
        style: str, optional
            the style of the assignment (e.g., ``bash``, 'SLURM')
        end: str, optional
            end-of-line byte
        """
        lines=self.byte_collector.split(end)
        logger.debug(lines)
        if len(lines)>0:
            for i,l in enumerate(lines):
                if len(l) > 0:
                    if style=='bash':
                        if l.startswith('export'):
                            logger.debug(f'{varname} {varvalue} : export assignment {l}')
                            tokens=l.split()
                            assignment=tokens[1]
                            evarname,evarvalue=assignment.split('=')
                            if evarname==varname:
                                lines[i]=f'export {varname}={varvalue}'
                        elif l.startswith(varname) and '=' in l:
                            lines[i]=f'{varname}={varvalue}'
                    elif style=='SLURM': # will only replace SLURM variables!
                        if l.startswith('#SBATCH'):
                            tokens=l.split()
                            logger.debug(f'SLURM statement {l}')
                            if tokens[1].startswith('--'):
                                assignment=tokens[1][2:]
                                if not '=' in assignment:
                                    logger.debug(f'cannot parse SLURM directive {l}')
                                else:
                                    svarname,svarvalue=assignment.split('=')
                                    if svarname==varname:
                                        lines[i]=f'#SBATCH --{svarname}={varvalue}'
                            elif tokens[1].startswith('-'):
                                svarname=tokens[1][1:]
                                svarvalue=tokens[2]
                                if svarname==varname:
                                    lines[i]=f'#SBATCH -{svarname} {varvalue}'
        self.byte_collector=end.join(lines)
        logger.debug(self.byte_collector)

    # def addline(self,msg,end='\n'):
    #     """Appends msg to the string as a line
        
    #     Parameters
    #     ----------
    #     msg: str
    #        the message
        
    #     end: str, optional
    #         end-of-line byte
    #     """
    #     self.byte_collector+=f'{msg}{end}'

    # def lastline(self,end='\n',exclude='#'):
    #     """Returns last line in the string
        
    #     Parameters
    #     ----------
    #     end: str, optional
    #         end-of-line byte
    #     exclude: str, optional
    #         comment byte
    #     """
    #     lines=[x for x in self.byte_collector.split(end) if (len(x)>0 and not x.startswith(exclude))]
    #     if len(lines)>0:
    #         return lines[-1]
    #     else:
    #         return None
    
    # def has_statement(self,statement,end='\n',exclude='#'):
    #     """
    #     Determines if a particular statement is on at least one non-comment line
        
    #     Parameters
    #     ----------
    #     statement: str
    #         the statement; e.g., ``exit``
    #     end: str, optional
    #         end-of-line byte
    #     exclude: str, optional
    #         comment byte
    #     """
    #     lines=[x for x in self.byte_collector.split(end) if (len(x)>0 and not x.startswith(exclude))]
    #     if len(lines)>0:
    #         for l in lines:
    #             if statement in l:
    #                 return True
    #     return False

    def ingest_file(self,filename):
        """
        Appends contents of file ``filename`` to the string

        Parameters
        ----------
        filename: str
           the name of the file
        """
        with open(filename,'r') as f:
            self.byte_collector+=f.read()

    def write_file(self,filename):
        """
        Writes string to ``filename``

        Parameters
        ----------
        filename: str
           the name of the file
        """
        with open(filename,'w') as f:
            f.write(self.byte_collector)

    def comment(self,msg,end='\n'):
        """
        Appends ``msg`` as a comment to the string
        
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
        """
        Logs a message to the string
        
        Parameters
        ----------
        msg: str
           the message to log"""
        my_logger(msg,self.addline)

    def banner(self,msg):
        """
        Logs a banner message to the string
        
        Parameters
        ----------
        msg: str
           the message to log
        """
        my_logger(msg,self.addline,fill='#',width=80,just='^')

    def __str__(self):
        """
        Returns the string representation of the byte collector
        """
        return self.byte_collector

def my_logger(msg,logf,width=None,fill='',just='<',frame='',depth=0,**kwargs):
    """
    A fancy logger
    
    Parameters
    ----------
    msg: str, list
       the message to be logged, either as a single string, list, or dict
    
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
    if width is None:
        if logf is print:
            ts=shutil.get_terminal_size((80,20))
            width=ts.columns
        else:
            width=67
    fmt=r'{'+r':'+fill+just+f'{width}'+r'}'
    ll=' ' if just in '^>' else ''
    rr=' ' if just in '^<' else ''
    if frame:
        ffmt=r'{'+r':'+frame+just+f'{width}'+r'}'
        logf(ffmt.format(frame))
    if type(msg)==list:
        for tok in msg:
            my_logger(tok,logf,width=width,fill=fill,just=just,frame=False,depth=depth,kwargs=kwargs)
    elif type(msg)==dict:
        for key,value in msg.items():
            if type(value)==str or not hasattr(value,"__len__"):
                my_logger(f'{key}: {value}',logf,width=width,fill=fill,just=just,frame=False,depth=depth,kwargs=kwargs)
            else:
                my_logger(f'{key}:',logf,width=width,fill=fill,just=just,frame=False,depth=depth,kwargs=kwargs)
                my_logger(value,logf,width=width,fill=fill,just=just,frame=False,depth=depth+1,kwargs=kwargs)
    elif type(msg)==pd.DataFrame:
        dfoutmode=kwargs.get('dfoutmode','value')
        if dfoutmode=='value':
            my_logger([ll+x+rr for x in msg.to_string().split('\n')],logf,width=width,fill=fill,just=just,frame=False,depth=depth,kwargs=kwargs)
        elif dfoutmode=='info':
            buf=StringIO()
            msg.info(buf=buf)
            my_logger([ll+x+rr for x in buf.getvalue().split('\n')],logf,width=width,fill=fill,just=just,frame=False,depth=depth,kwargs=kwargs)
        else:
            return
    else:
        indent=f'{" "*depth*2}' if just=='<' and not kwargs.get('no_indent',False) else ''
        if type(msg)==str:
            lns=msg.split('\n')
            if len(lns)>1:
                my_logger(lns,logf,width=width,fill=fill,just=just,frame=False,depth=depth+1,kwargs=kwargs)
            else:
                outstr=indent+ll+f'{msg}'+rr
                logf(fmt.format(outstr))
        else:
            outstr=indent+ll+f'{msg}'+rr
            logf(fmt.format(outstr))
    if frame:
        logf(ffmt.format(frame))
            
class FileCollector(UserList):
    """
    A class for handling collections of files

    """
    def flush(self):
        """
        Deletes all files in the collection and clears the collection
        """
        logger.debug(f'Flushing file collector: {len(self)} files.')
        for f in self:
            if os.path.exists(f):
                os.remove(f)
            else:
                logger.debug(f'{f}: not found.')
        self.clear()
        
    def tarball(self,basename):
        """
        Makes a tarball of the files in the collection
        
        Parameters
        ----------
        basename: str
            basename of the resulting tarball
        """
        filelist=' '.join([x for x in self])
        subprocess.run(f'tar zvcf {basename}.tgz {filelist}',
                        shell=True, 
                        executable='/bin/bash',
                        check=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)
        logger.debug(f'generated tarball {basename}.tgz')

def split_ri(ri):
    """
    A simple utility function for splitting the integer resid and
    1-byte insertion code out of a string resid-insertion code
    concatenation

    Parameters
    ----------
    ri: str
        the string representation of a residue number, e.g., ``123A`` or ``123``

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
    """
    Joins a residue sequence number and an insertion code into a single string.

    Parameters
    ----------
    resseqnum: int
        The residue sequence number, e.g., 123.
    insertion: str
        The insertion code, e.g., ``A``. If there is no insertion code, this
        should be an empty string.
    """
    if insertion=='':
        return resseqnum
    return f'{resseqnum}{insertion}'

def ri_range(val,split_chars=['-','#']):
    """
    Splits a string representation of a range of residue numbers into a list of
    residue numbers. The string can contain multiple ranges separated by
    characters in ``split_chars``. The ranges can be specified as a single
    residue number, a range of residue numbers (e.g., ``123-456``),
    or a range of residue numbers with insertion codes (e.g., ``123A-456B``).

    Parameters
    ----------
    val: str
        The string representation of the residue number range.
    split_chars: list of str, optional
        A list of characters that can be used to split the string into multiple
        ranges. Defaults to ['-', '#']. 
    Returns
    -------
    list of str: A list of residue numbers in string format, e.g., ['123', '124', '125A', '126B'].
    """
    the_split=[val]
    for c in split_chars:
        the_splits=[x.split(c) for x in the_split]
        the_split=[]
        for s in the_splits:
            the_split.extend(s)
    return [split_ri(x) for x in the_split]

# def rescale_packmol_inp_box(fn,scale):
#     shutil.copy(fn,f'{fn}-bak')
#     with open(fn,'r') as f:
#         lines=f.read().split('\n')
#     for i in range(len(lines)):
#         l=lines[i]
#         tokens=[x.strip() for x in l.split()]
#         if tokens:
#             if tokens[0]=='inside':
#                 if tokens[1]=='box':
#                     llx,lly,llz,urx,ury,urz=list(map(float,tokens[2:]))
#                     llx*=scale[0]
#                     urx*=scale[0]
#                     lly*=scale[1]
#                     ury*=scale[1]
#                     llz*=scale[2]
#                     urz*=scale[2]
#                     newline=f'  inside box {llx:.2f} {lly:.2f} {llz:.2f} {urx:.2f} {ury:.2f} {urz:.2f}'
#                     lines[i]=newline
#     with open('tmp','w') as f:
#         for l in lines:
#             f.write(f'{l}\n')
#     shutil.move('tmp',fn)

def oxford(a_list,conjunction='or'):
    """ 
    Obey the Oxford comma when ending a list with a conjunction

    Parameters
    ----------
    a_list: list of str
        The list of items to join.
    conjunction: str, optional
        The conjunction to use when joining the list. Defaults to 'or'.

    Returns
    -------
    str: The Oxford-compliant string representation of the list.
    """
    if not a_list: return ''
    if len(a_list)==1:
        return a_list[0]
    elif len(a_list)==2:
        return f'{a_list[0]} {conjunction} {a_list[1]}'
    else:
        return ", ".join(a_list[:-1])+f', {conjunction} {a_list[-1]}'
    
def linesplit(line,cchar='!'):
    """
    Splits a line at the first occurrence of the character ``cchar``.
    """
    if not cchar in line:
        return line,''
    idx=line.index(cchar)
    if idx==0:
        return '',line[1:]
    return line[:idx],line[idx+1:]

def striplist(L):
    """
    Strips whitespace from each item in a list and removes empty strings.

    Parameters
    ----------
    L: list of str
        The list of strings to be stripped.
    Returns
    -------
    list of str: A new list with whitespace stripped from each item and empty strings removed.
    """
    l=[x.strip() for x in L]
    while '' in l:
        l.remove('')
    return l


def to_latex_math(s):
    """
    Converts a string into LaTeX math format.  Written by ChatGPT 4o.
    """
    def convert_token(token):
        # Match base, subscript, and superscript using regex
        match = re.fullmatch(r'\s*([a-zA-Z0-9]+)(?:_([a-zA-Z0-9^]+))?(?:\^([a-zA-Z0-9_]+))?\s*', token)
        if not match:
            return token  # return as-is if it doesn't match expected pattern
        base, sub, sup = match.groups()
        result = base
        if sub:
            if '^' in sub:
                # Handle cases like a_x^2 packed into sub
                sub, sup = sub.split('^', 1)
            result += f'$_{{{sub}}}$'
        if sup:
            result += f'$^{{{sup}}}$'
        return result

    parts = s.split(',')
    converted = [convert_token(part) for part in parts]
    return ', '.join(converted)
