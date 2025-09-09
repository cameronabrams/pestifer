# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import logging
from ..util.stringthings import ByteCollector

logger = logging.getLogger(__name__)

class GenericScripter:
    """
    A class for writing files with specific formats.
    
    Parameters
    ----------
    comment_char : str, optional
        The character used for comments in the script. Default is '#'.
    """
    def __init__(self, *args, **kwargs):
        comment_char = kwargs.pop('comment_char', '#')
        self.B = ByteCollector(comment_char=comment_char)
        self.is_written = False
        self.indent = ''

    def newfile(self, filename):
        """
        Initializes a new file for writing.
        
        Parameters
        ----------
        filename : str
            The name of the file to be created or overwritten.
        """
        self.filename = filename
        self.B.reset()
        self.is_written = False
        return self

    def ingest_file(self, filename):
        """
        Ingest a file into the script writer.

        Parameters
        ----------
        filename : str
            The name of the file to be ingested. The file's contents will be added to the ByteCollector.

        """
        self.B.ingest_file(filename)

    def addline(self, data, indents=0):
        """
        Add a line of code to the ByteCollector.
        
        Parameters
        ----------
        data : str
            The line of code to be added.
        indents : int, optional
            The number of indentation levels to apply. Default is 0.
        """
        self.B.addline(indents * self.indent + data)

    def banner(self, data):
        """
        Add a banner comment to the ByteCollector.

        Parameters
        ----------
        data : str
            The banner comment to be added.
        """
        self.B.banner(data)

    def comment(self, data):
        """
        Add a comment to the ByteCollector.

        Parameters
        ----------
        data : str
            The comment to be added.
        """
        self.B.comment(data)

    def has_statement(self, statement):
        """
        Check if a specific statement is present in the ByteCollector.

        Parameters
        ----------
        statement : str
            The statement to check for.
        """
        return self.B.has_statement(statement)

    def writefile(self, force=False):
        """
        Write the contents of the ByteCollector to a file.
        
        Parameters
        ----------
        force : bool, optional
            If True, overwrite the file even if it has already been written.
        """
        if not self.is_written or force:
            with open(self.filename, 'w') as f:
                f.write(str(self.B))
            self.is_written = True
        else:
            logger.debug(f'{self.filename} has already been written')

    def __repr__(self):
        return f"{self.__class__.__name__}"