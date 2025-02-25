# Author: Cameron F. Abrams, <cfa22@drexel.edu>
from scipy.constants import physical_constants

_SYMBOLS_={
    'ANGSTROM':'Å',
    'CUBED':'³',
    'SQUARED':'²'
}
_UNITS_={
    'SQUARE-ANGSTROMS':f'{_SYMBOLS_["ANGSTROM"]}{_SYMBOLS_["SQUARED"]}',
    'CUBIC-ANGSTROMS':f'{_SYMBOLS_["ANGSTROM"]}{_SYMBOLS_["CUBED"]}',
    }
g_per_amu=physical_constants['atomic mass constant'][0]*1000
A_per_cm=1.e8
A3_per_cm3=A_per_cm**3