# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" Defines the ByteCollector and FileCollector classes
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

banner_message="""
    pestifer v. {}
    https://pestifer.readthedocs.io/en/latest/

    Cameron F. Abrams <cfa22@drexel.edu>

    Supported in part by Grants GM100472, AI154071, 
    and AI178833 from the NIH

    CHARMM force field files (July 24) from the 
    MacKerell Lab
    """.format(__pestifer_version__)

enhanced_banner_message="""
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
        msg: str            for ln in:
                outstr=ll+ln+rr
                logf(fmt.format(outstr))

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
    
    def reassign(self,varname,varvalue,style='bash',end='\n'):
        """Reassigns variable varname to value varvalue if an assignment to varname already
           exists in the byte_collector """
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

    def injest_file(self,filename):
        """Appends contents of file 'filename' to the string
        
        Parameters
        ----------
        filename: str
           the name of the file
        """
        with open(filename,'r') as f:
            self.byte_collector+=f.read()

    def write_file(self,filename):
        """Writes string to 'filename' 
        
        Parameters
        ----------
        filename: str
           the name of the file
        """
        with open(filename,'w') as f:
            f.write(self.byte_collector)

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
        my_logger(msg,self.addline,fill='#',width=80,just='^')

    def __str__(self):
        return self.byte_collector

def my_logger(msg,logf,width=None,fill='',just='<',frame='',depth=0,**kwargs):
    """A fancy logger
    
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
        subprocess.run(f'tar zvcf {basename}.tgz {filelist}',
                        shell=True, 
                        executable='/bin/bash',
                        check=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)
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

def rescale_packmol_inp_box(fn,scale):
    shutil.copy(fn,f'{fn}-bak')
    with open(fn,'r') as f:
        lines=f.read().split('\n')
    for i in range(len(lines)):
        l=lines[i]
        tokens=[x.strip() for x in l.split()]
        if tokens:
            if tokens[0]=='inside':
                if tokens[1]=='box':
                    llx,lly,llz,urx,ury,urz=list(map(float,tokens[2:]))
                    llx*=scale[0]
                    urx*=scale[0]
                    lly*=scale[1]
                    ury*=scale[1]
                    llz*=scale[2]
                    urz*=scale[2]
                    newline=f'  inside box {llx:.2f} {lly:.2f} {llz:.2f} {urx:.2f} {ury:.2f} {urz:.2f}'
                    lines[i]=newline
    with open('tmp','w') as f:
        for l in lines:
            f.write(f'{l}\n')
    shutil.move('tmp',fn)

def oxford(a_list,conjunction='or'):
    """ Obey the Oxford comma when ending a list with a conjunction
    """
    if not a_list: return ''
    if len(a_list)==1:
        return a_list[0]
    elif len(a_list)==2:
        return f'{a_list[0]} {conjunction} {a_list[1]}'
    else:
        return ", ".join(a_list[:-1])+f', {conjunction} {a_list[-1]}'
    
def linesplit(line,cchar='!'):
    if not cchar in line:
        return line,''
    idx=line.index(cchar)
    if idx==0:
        return '',line[1:]
    return line[:idx],line[idx+1:]

def striplist(L):
    l=[x.strip() for x in L]
    while '' in l:
        l.remove('')
    return l


def to_latex_math(s):
    # written by ChatGPT
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
