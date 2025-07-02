# Author: Cameron F. Abrams, <cfa22@drexel.edu>
from scipy.constants import physical_constants, Avogadro
import numpy as np

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
cm3_per_A3=1.0/A3_per_cm3
n_per_mol=Avogadro

def nmolec_in_cuA(MW_g,density_gcc,volume_A3):
    """
    Calculate the number of molecules in a given volume in cubic Angstroms.

    Parameters
    ----------
    MW_g : float
        Molecular weight in grams per mole.
    density_gcc : float
        Density in grams per cubic centimeter.
    volume_A3 : float
        Volume in cubic Angstroms.

    Returns
    -------
    int
        Number of molecules in the specified volume.
    """
    return int(np.floor(density_gcc/MW_g*cm3_per_A3*n_per_mol*volume_A3))

def cuA_of_nmolec(MW_g,density_gcc,nmolec):
    """
    Calculate the volume in cubic Angstroms of a given number of molecules.

    Parameters
    ----------
    MW_g : float
        Molecular weight in grams per mole.
    density_gcc : float
        Density in grams per cubic centimeter.
    nmolec : int
        Number of molecules.

    Returns
    -------
    float
        Volume in cubic Angstroms.
    """
    return nmolec/(density_gcc/MW_g*cm3_per_A3*n_per_mol)