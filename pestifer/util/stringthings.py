# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
Defines the :class:`ByteCollector` and :class:`FileCollector` classes, and defines the :func:`banner` and :func:`my_logger` functions.
"""
import ast
import logging
import operator
import os
import re
import subprocess
import shutil
import sys

import pandas as pd

from argparse import Namespace
from collections import UserList
from io import StringIO
from pathlib import Path
from typing import Callable

logger = logging.getLogger(__name__)

import importlib.metadata

__pestifer_version__ = importlib.metadata.version("pestifer")

_banner_message="""
    Pestifer v. {}
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

def banner(logf: Callable, args: Namespace):
    """
    Writes a banner message to the log file

    Parameters
    ----------
    logf: file-like object
        The log file to which the banner message will be written.
    """
    if args.banner:
        if args.kick_ass:
            my_logger(_enhanced_banner_message, logf, fill=' ', just='<')
        else:
            my_logger(_banner_message, logf, fill=' ', just='<')

class ByteCollector:
    """
    A simple string manager.

    Parameters
    ----------
    comment_char: str, optional
        The character used to denote comments in the string.
        Defaults to '#'.
    line_length: int, optional
        The maximum length of a line in the string.
        Defaults to 80.
    """
    def __init__(self, comment_char='#', line_length=80):
        self.line_length = line_length
        self.comment_char = comment_char
        self.byte_collector = ''

    def reset(self):
        """
        Resets the string
        """
        self.byte_collector = ''

    def write(self, msg):
        """
        Appends msg to the string
        
        Parameters
        ----------
        msg: str
           the message
        """
        self.byte_collector += msg

    def addline(self, msg, end='\n'):
        """Appends msg to the string as a line

        Parameters
        ----------
        msg: str
           the message

        end: str, optional
           the end-of-line byte
        """
        self.byte_collector += f'{msg}{end}'

    def lastline(self, end='\n', exclude='#'):
        """
        Returns last line in the string

        Parameters
        ----------
        end: str, optional
            end-of-line byte
        exclude: str, optional
            comment byte
        """
        lines = [x for x in self.byte_collector.split(end) if (len(x) > 0 and not x.startswith(exclude))]
        if len(lines) > 0:
            return lines[-1]
        else:
            return None

    def has_statement(self, statement: str, end='\n', exclude='#'):
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
        lines = [x for x in self.byte_collector.split(end) if (len(x) > 0 and not x.startswith(exclude))]
        if len(lines) > 0:
            for l in lines:
                if statement in l:
                    return True
        return False

    def reassign(self, varname: str, varvalue: str, style: str = 'bash', end: str = '\n'):
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
        lines = self.byte_collector.split(end)
        logger.debug(lines)
        if len(lines) > 0:
            for i, l in enumerate(lines):
                if len(l) > 0:
                    if style == 'bash':
                        if l.startswith('export'):
                            logger.debug(f'{varname} {varvalue} : export assignment {l}')
                            tokens = l.split()
                            assignment = tokens[1]
                            evarname, evarvalue = assignment.split('=')
                            if evarname == varname:
                                lines[i] = f'export {varname}={varvalue}'
                        elif l.startswith(varname) and '=' in l:
                            lines[i] = f'{varname}={varvalue}'
                    elif style == 'SLURM':  # will only replace SLURM variables!
                        if l.startswith('#SBATCH'):
                            tokens = l.split()
                            logger.debug(f'SLURM statement {l}')
                            if tokens[1].startswith('--'):
                                assignment = tokens[1][2:]
                                if not '=' in assignment:
                                    logger.debug(f'cannot parse SLURM directive {l}')
                                else:
                                    svarname, svarvalue = assignment.split('=')
                                    if svarname == varname:
                                        lines[i] = f'#SBATCH --{svarname}={varvalue}'
                            elif tokens[1].startswith('-'):
                                svarname = tokens[1][1:]
                                svarvalue = tokens[2]
                                if svarname == varname:
                                    lines[i] = f'#SBATCH -{svarname} {varvalue}'
        self.byte_collector = end.join(lines)
        logger.debug(self.byte_collector)

    def ingest_file(self, filename: Path | str):
        """
        Appends contents of file ``filename`` to the string

        Parameters
        ----------
        filename: str
           the name of the file
        """
        with open(filename, 'r') as f:
            self.byte_collector += f.read()

    def write_file(self, filename: Path | str):
        """
        Writes string to ``filename``

        Parameters
        ----------
        filename: str
           the name of the file
        """
        with open(filename, 'w') as f:
            f.write(self.byte_collector)

    def comment(self, msg: str, end: str = '\n'):
        """
        Appends ``msg`` as a comment to the string
        
        Parameters
        ----------
        msg: str
           the message
        
        end: str, optional
            end-of-line byte
        """
        comment_line = f'{self.comment_char} {msg}'
        comment_words = comment_line.split()
        comment_lines = ['']
        current_line_idx = 0
        for word in comment_words:
            test_line = ' '.join(comment_lines[current_line_idx].split() + [word])
            if len(test_line) > self.line_length:
                comment_lines.append(f'{self.comment_char} {word}')
                current_line_idx += 1
            else:
                comment_lines[current_line_idx] = test_line
        for line in comment_lines:
            self.addline(line, end=end)

    def log(self, msg: str):
        """
        Logs a message to the string
        
        Parameters
        ----------
        msg: str
           the message to log"""
        my_logger(msg, self.addline)

    def banner(self, msg: str):
        """
        Logs a banner message to the string
        
        Parameters
        ----------
        msg: str
           the message to log
        """
        my_logger(msg, self.addline, fill='#', width=80, just='^')

    def __str__(self):
        """
        Returns the string representation of the byte collector
        """
        return self.byte_collector

def my_logger(msg: str | list | dict | pd.DataFrame, logf: Callable, width=None, fill='', just='<', frame='', depth=0, **kwargs):
    """
    A fancy recursive logger
    
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
            ts = shutil.get_terminal_size((80, 20))
            width = ts.columns
        else:
            width = 67
    fmt = r'{'+r':'+fill+just+f'{width}'+r'}'
    ll = ' ' if just in '^>' else ''
    rr = ' ' if just in '^<' else ''
    if frame:
        ffmt = r'{'+r':'+frame+just+f'{width}'+r'}'
        logf(ffmt.format(frame))
    # logger.debug(f'Logging message of type {type(msg)} at depth {depth}; is list? {isinstance(msg, list)}; is dict? {isinstance(msg, dict)}; is DataFrame? {isinstance(msg, pd.DataFrame)}')
    if isinstance(msg, list):
        for tok in msg:
            my_logger(tok, logf, width=width, fill=fill, just=just, frame=False, depth=depth, kwargs=kwargs)
    elif isinstance(msg, dict):
        # logger.debug(f'logging dict of length {len(msg)}')
        for key, value in msg.items():
            if isinstance(value, str) or not hasattr(value, "__len__"):
                my_logger(f'{key}: {value}', logf, width=width, fill=fill, just=just, frame=False, depth=depth, kwargs=kwargs)
            elif isinstance(value, list):
                # logger.debug(f'key {key} is a list of length {len(value)}')
                for tok in value:
                    my_logger(tok, logf, width=width, fill=fill, just=just, frame=False, depth=depth+1, kwargs=kwargs)
            else:
                my_logger(f'{key}:', logf, width=width, fill=fill, just=just, frame=False, depth=depth, kwargs=kwargs)
                my_logger(value, logf, width=width, fill=fill, just=just, frame=False, depth=depth+1, kwargs=kwargs)
    elif isinstance(msg, pd.DataFrame):
        dfoutmode = kwargs.get('dfoutmode', 'value')
        if dfoutmode == 'value':
            my_logger([ll+x+rr for x in msg.to_string().split('\n')], logf, width=width, fill=fill, just=just, frame=False, depth=depth, kwargs=kwargs)
        elif dfoutmode == 'info':
            buf = StringIO()
            msg.info(buf=buf)
            my_logger([ll+x+rr for x in buf.getvalue().split('\n')], logf, width=width, fill=fill, just=just, frame=False, depth=depth, kwargs=kwargs)
        else:
            return
    else:
        indent = f'{" "*depth*2}' if just == '<' and not kwargs.get('no_indent', False) else ''
        if isinstance(msg, str):
            lns = msg.split('\n')
            if len(lns) > 1:
                my_logger(lns, logf, width=width, fill=fill, just=just, frame=False, depth=depth+1, kwargs=kwargs)
            else:
                outstr = indent + ll + f'{msg}' + rr
                logf(fmt.format(outstr))
        else:
            outstr = indent + ll + f'{msg}' + rr
            logf(fmt.format(outstr))
    if frame:
        logf(ffmt.format(frame))
            
class FileCollector(UserList):
    """
    A class for handling collections of files; inherits from :class:`~collections.UserList`.
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

    def tarball(self, basename, preferred_extension='tgz'):
        """
        Makes a tarball of the files in the collection
        
        Parameters
        ----------
        basename: str
            basename of the resulting tarball
        """
        filelist = ' '.join([x for x in self])
        subprocess.run(f'tar zvcf {basename}.{preferred_extension} {filelist}',
                        shell=True, 
                        executable='/bin/bash',
                        check=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)
        logger.debug(f'generated tarball {basename}.{preferred_extension}')

def oxford(a_list: list[str], conjunction: str = 'or') -> str | None:
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
    if not a_list: return None
    if len(a_list) == 1:
        return a_list[0]
    elif len(a_list) == 2:
        return f'{a_list[0]} {conjunction} {a_list[1]}'
    else:
        return ", ".join(a_list[:-1]) + f', {conjunction} {a_list[-1]}'

def linesplit(line: str, cchar: str = '!') -> tuple[str, str]:
    """
    Splits a line at the first occurrence of the character ``cchar``.
    """
    if not cchar in line:
        return line, ''
    idx = line.index(cchar)
    if idx == 0:
        return '', line[1:]
    return line[:idx], line[idx+1:]

def striplist(L: list[str]) -> list[str]:
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
    l = [x.strip() for x in L]
    while '' in l:
        l.remove('')
    return l

def to_latex_math(s: str) -> str:
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

def example_footer(author_name: str = '', author_email: str = '') -> str:
    """
    Returns a formatted footer string with author information.

    Parameters
    ----------
    author_name: str, optional
        The name of the author. Defaults to an empty string.
    author_email: str, optional
        The email of the author. Defaults to an empty string.

    Returns
    -------
    str: A formatted footer string.
    """
    footer = ''
    if author_name and author_email:
        footer = f""".. raw:: html

    <div class="autogen-footer">
        <p>Example author: {author_name} &nbsp;&nbsp;&nbsp; Contact: <a href="mailto:{author_email}">{author_email}</a></p>
    </div>"""
    return footer

def raise_clean(ErrorInstance: Exception):
    """
    Raises an error with a clean message showing no traceback.

    Parameters
    ----------
    ErrorInstance: Exception instance
        The exception instance to raise.
    """
    try:
        raise ErrorInstance
    except ErrorInstance.__class__ as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

# def special_update(dict, key, val, mode='strict'):
#     """
#     Special update function for dictionaries with different modes.

#     Parameters
#     ----------
#     dict: dict
#         The dictionary to update.
#     key: str
#         The key to update.
#     val: any
#         The value to set.
#     mode: str, optional
#         The mode of the update. Can be 'strict' or 'permissive'. Defaults to 'strict'.

#     Raises
#     ------
#     KeyError
#         If the mode is 'strict' and the key already exists.
#     """
#     if mode == 'strict' and key in dict:
#         raise KeyError(f'Key {key} already exists in dictionary.')
#     dict[key] = val

def plu(n: int, singform: str = '', plurform: str = 's') -> str:
    """
    Returns the singular or plural form of a word based on the count.

    Parameters
    ----------
    n: int
        The count to determine singular or plural.
    singform: str, optional
        The singular form of the word/word ending. Defaults to an empty string.
    plurform: str, optional
        The plural form of the word/word ending. Defaults to 's'.

    Returns
    -------
    str: The appropriate form of the word based on the count.

    Examples
    --------
    >>> for ncats in [1, 2]:
    ...     print(f'There {plu(ncats, "is", "are")} {ncats} cat{plu(ncats)}.')
    There is 1 cat.
    There are 2 cats.
    >>> for nberries in [1, 2]:
    ...     print(f'There {plu(nberries, "is", "are")} {nberries} berr{plu(nberries, "y", "ies")}.')
    There is 1 berry.
    There are 2 berries.
    >>> for noctopi in [1, 2]:
    ...     print(f'There {plu(noctopi, "is", "are")} {noctopi} octop{plu(noctopi, "us", "i")}.')
    There is 1 octopus.
    There are 2 octopi.

    """
    return singform if n == 1 else plurform

# Supported operators for comparisons
comp_ops = {
    ast.Eq: operator.eq,
    ast.NotEq: operator.ne,
    ast.Lt: operator.lt,
    ast.LtE: operator.le,
    ast.Gt: operator.gt,
    ast.GtE: operator.ge,
    ast.In: lambda a, b: a in b,
    ast.NotIn: lambda a, b: a not in b,
}

# Supported boolean ops
bool_ops = {
    ast.And: all,
    ast.Or: any,
}

def parse_filter_expression(expr: str):
    """
    Parses a logical expression string (e.g., "segtype == 'protein' and resid > 100")
    and returns a function f(obj) -> bool that evaluates the expression on obj's attributes.
    """
    expr_ast = ast.parse(expr, mode='eval')

    def _eval(node, obj):
        if isinstance(node, ast.Expression):
            return _eval(node.body, obj)

        elif isinstance(node, ast.BoolOp):
            op_func = bool_ops[type(node.op)]
            return op_func(_eval(value, obj) for value in node.values)

        elif isinstance(node, ast.Compare):
            # Left side can be Name or Attribute, only support simple Names for now
            if isinstance(node.left, ast.Name):
                left_val = getattr(obj, node.left.id)
            else:
                raise ValueError(f"Unsupported left operand type: {type(node.left)}")

            # Only support single comparator for simplicity
            op = node.ops[0]
            comparator = node.comparators[0]

            # Evaluate right side literals
            if isinstance(comparator, (ast.Constant)):
                right_val = comparator.value if hasattr(comparator, 'value') else comparator.s
            elif isinstance(comparator, ast.Name):
                right_val = getattr(obj, comparator.id)
            else:
                # For lists or tuples, recursively evaluate literals
                right_val = ast.literal_eval(comparator)

            op_func = comp_ops.get(type(op))
            if op_func is None:
                raise ValueError(f"Unsupported comparison operator: {type(op)}")

            return op_func(left_val, right_val)

        elif isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.Not):
            return not _eval(node.operand, obj)

        elif isinstance(node, ast.Constant):
            # True / False / None literals
            return node.value

        else:
            raise ValueError(f"Unsupported AST node: {type(node)}")

    def filter_func(obj):
        return _eval(expr_ast, obj)

    return filter_func
