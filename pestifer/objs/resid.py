
from __future__ import annotations
from pydantic import Field
from ..core.baseobj import BaseObj, BaseObjList

class ResID(BaseObj):
    """
    A class for handling residue numbers and insertion codes.
    This class is used to represent a residue number and an optional insertion code.
    It provides methods to split and join residue numbers and insertion codes.
    """

    _required_fields = {'resseqnum'}
    _optional_fields = {'insertion'}
    resseqnum: int = Field(..., description="Residue sequence number")
    insertion: str | None = Field(None, description="Residue insertion code")

    def __str__(self):
        """ 
        Returns the string representation of the residue number and insertion code.
        If there is no insertion code, it returns just the residue number.
        """
        if self.insertion is None or self.insertion == '':
            return str(self.resseqnum)
        return f'{self.resseqnum}{self.insertion}'

    @property
    def pdbresid(self) -> str:
        """
        Returns the residue number and insertion code in a format suitable for PDB files.
        If there is no insertion code, it returns just the residue number as a string.
        """
        if self.insertion is None or self.insertion == '':
            return str(self.resseqnum)+' '
        return f'{self.resseqnum}{self.insertion}'

    @property
    def resid(self) -> str | int:
        """
        Returns the residue number and insertion code as a string, or if there is no insertion code,
        returns just the residue number as an integer.
        This is equivalent to the string representation of the ResID object.
        """
        if self.insertion is None or self.insertion == '':
            return self.resseqnum
        return str(self)

    def __hash__(self):
        return hash((self.resseqnum, self.insertion))

    @classmethod
    def _adapt(cls, *args, **kwargs) -> dict:
        """
        Adapts the input to a dictionary format suitable for ResID instantiation.
        This method is used to convert various input types into a dictionary of parameters.
        """
        if args and len(args) == 2:
            # assume first arg is integer resnum and second is insertion code
            resseqnum, insertion = args
            if insertion is not None:
                return {'resseqnum': resseqnum, 'insertion': insertion}
            return {'resseqnum': resseqnum}
        elif args and isinstance(args[0], int):
            # assume first arg is integer resnum
            return {'resseqnum': args[0]}
        elif args and isinstance(args[0], str):
            resseqnum, insertion = ResID.split_ri(args[0])
            if insertion is not None:
                return {'resseqnum': resseqnum, 'insertion': insertion}
            return {'resseqnum': resseqnum}
        elif args and isinstance(args[0], ResID):
            if args[0].insertion is not None:
                return {'resseqnum': args[0].resseqnum, 'insertion': args[0].insertion}
            return {'resseqnum': args[0].resseqnum}
        return super()._adapt(*args, **kwargs)

    @staticmethod
    def split_ri(ri) -> tuple[int, str | None]:
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
    tuple[int, str | None]: the integer resid and the 1-byte insertion code or None if none
    """
        if ri[-1].isdigit(): # there is no insertion code
            r = int(ri)
            i = None
        else:
            r = int(ri[:-1])
            i = ri[-1]
        return r, i

    def __lt__(self, other: ResID) -> bool:
        """
        Compares two ResID objects based on their residue sequence number and insertion code.
        """
        if self.resseqnum < other.resseqnum:
            return True
        if self.resseqnum == other.resseqnum and (self.insertion or '') < (other.insertion or ''):
            return True
        return False

    def __le__(self, other: ResID) -> bool:
        """
        Checks if self is less than or equal to other
        """
        if self.resseqnum <= other.resseqnum:
            return True
        if self.resseqnum == other.resseqnum and (self.insertion or '') <= (other.insertion or ''):
            return True
        return False

    def __eq__(self, other: ResID | dict) -> bool:
        """
        Checks if two ResID objects are equal based on their residue sequence number and insertion code.
        """
        if isinstance(other, ResID):
            return self.resseqnum == other.resseqnum and (self.insertion or '') == (other.insertion or '')
        elif isinstance(other, dict):
            return self.resseqnum == other.get('resseqnum') and (self.insertion or '') == (other.get('insertion') or '')
        return False

    def __add__(self, other: ResID | int) -> ResID:
        """
        Adds an integer to the residue sequence number of this ResID object.
        If the other object is a ResID, it adds the residue sequence numbers.
        """
        if isinstance(other, ResID):
            return ResID(resseqnum=self.resseqnum + other.resseqnum, insertion=self.insertion)
        elif isinstance(other, int):
            return ResID(resseqnum=self.resseqnum + other, insertion=self.insertion)

    def increment(self, by_seqnum: bool = False) -> ResID:
        """
        Increments the insertion code of this ResID object by one.
        If there is no insertion code, it sets it to 'A'.
        """
        if by_seqnum:
            return ResID(resseqnum=self.resseqnum + 1, insertion=None)
        if self.insertion is None or self.insertion == '':
            return ResID(resseqnum=self.resseqnum, insertion='A')
        else:
            new_ord = ord(self.insertion) + 1
            if new_ord == 64:
                new_ord = 97
            elif new_ord > 122:
                raise ValueError(f"Insertion code cannot exceed '{chr(122)}'")
            new_insertion_code = chr(new_ord)
            return ResID(resseqnum=self.resseqnum, insertion=new_insertion_code)

    @staticmethod
    def join_ri(resseqnum: int, insertion: str | None = None) -> str:
        """
        Joins a residue sequence number and an insertion code into a single string.

        Parameters
        ----------
        resseqnum: int
            The residue sequence number, e.g., 123.
        insertion: str | None
            The insertion code, e.g., ``A``. If there is no insertion code, this
            should be an empty string or None.
        
        Returns
        -------
        str: The combined residue number and insertion code as a string.
        """
        if insertion is None or insertion == '':
            return str(resseqnum)
        return f'{resseqnum}{insertion}'

class ResIDList(BaseObjList[ResID]):
    """
    A list of ResID objects.
    This class is used to handle a list of residue numbers and insertion codes.
    It inherits from BaseObjList and provides methods to describe the list.
    """
    
    def describe(self) -> str:
        """
        Describe the ResIDList.

        Returns
        -------
        str
            A string description of the ResIDList, including the number of residues.
        """
        return f"<ResIDList with {len(self)} residues>"

    def __init__(self, *args, **kwargs):
        """
        Initializes the ResIDList with a list of ResID objects.
        If the input is a single ResID object, it is converted to a list.
        """
        if len(args) == 1 and isinstance(args[0], ResID):
            args = (args[0],)
        elif len(args) == 1 and isinstance(args[0], str):
            args = ResIDList.ri_range(args[0]),
        super().__init__(*args, **kwargs)

    @staticmethod
    def ri_range(val, split_chars: tuple[str] = ('-', '#')) -> tuple[ResID]:
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
        tuple of ResID: A tuple of residue instances, e.g., (ResID('123'), ResID('124'), ResID('125A'), ResID('126B')).
        """
        the_split = [val]
        for c in split_chars:
            the_splits = [x.split(c) for x in the_split]
            the_split = []
            for s in the_splits:
                the_split.extend(s)
        return tuple([ResID(x) for x in the_split])