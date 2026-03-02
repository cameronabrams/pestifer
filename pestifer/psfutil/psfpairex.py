# Author: Cameron F. Abrams, <cfa22@drexel.edu>

"""
Definition of the :class:`PSFPairEx` and :class:`PSFPairExList` classes for handling
non-bonded pair exclusions in PSF topology files.

In the dual-topology paradigm for alchemical free energy calculations in NAMD, atoms
belonging to one alchemical end-state must not interact non-bonded with atoms belonging
to the other end-state.  These cross-topology exclusions are recorded in the ``NNB``
section of the PSF file as explicit pairs of atom serial numbers, using the same
integer-pair line format as the ``BOND`` section.

The :class:`PSFPairEx` class represents a single such excluded pair, while the
:class:`PSFPairExList` class is a collection of such pairs.
These are descendants of the :class:`PSFTopoElement <.psftopoelement.PSFTopoElement>` and
:class:`PSFTopoElementList <.psftopoelement.PSFTopoElementList>` classes, respectively,
which provide a framework for representing and manipulating generic PSF topology elements.
"""

from functools import singledispatchmethod

from .psftopoelement import PSFTopoElement, PSFTopoElementList, LineList


class PSFPairEx(PSFTopoElement):
    """
    A class representing a cross-topology non-bonded pair exclusion in a PSF topology file,
    used in the dual-topology alchemical free energy paradigm to prevent atoms belonging to
    one end-state from interacting non-bonded with atoms belonging to the other end-state.
    """

    def __eq__(self, other):
        """
        Check if two PSFPairEx objects are equal by comparing their serial numbers.
        Two excluded pairs are considered equal if they have the same serial numbers
        in either forward or reverse order.

        Parameters
        ----------
        other : PSFPairEx
            The other excluded pair to compare against.
        """
        return (
            [self.serial1, self.serial2] == [other.serial1, other.serial2]
            or [self.serial1, self.serial2] == [other.serial2, other.serial1]
        )


class PSFPairExList(PSFTopoElementList):
    """
    A class representing a list of :class:`PSFPairEx` objects.
    This class inherits from :class:`PSFTopoElementList <.psftopoelement.PSFTopoElementList>`
    and provides methods for managing a collection of non-bonded pair exclusions in a PSF
    topology file.  It can be initialized from a list of lines or from a
    :class:`PSFTopoElementList <.psftopoelement.PSFTopoElementList>` object.
    """

    @singledispatchmethod
    def __init__(self, input_data, **kwargs):
        super().__init__(input_data)

    @__init__.register(PSFTopoElementList)
    def _from_tel(self, input_obj, **kwargs):
        super().__init__(input_obj)

    @__init__.register(LineList)
    def _from_list(self, lines, include_serials=[]):
        P = []
        for line in lines:
            tokens = [int(x) for x in line.split()]
            assert len(tokens) % 2 == 0, f'Poorly formatted NNB line in psffile?'
            if not include_serials:
                for l, r in zip(tokens[:-1:2], tokens[1::2]):
                    P.append(PSFPairEx([l, r]))
            else:
                for l, r in zip(tokens[:-1:2], tokens[1::2]):
                    if include_serials[l - 1] and include_serials[r - 1]:
                        P.append(PSFPairEx([l, r]))
        super().__init__(P)
