# Author: Cameron F. Abrams <cfa22@drexel.edu>.
""" 
Some calculations based on atom coordinates
"""
import logging
import os

import numpy as np
import pandas as pd

from pidibble.pdbparse import PDBParser

from ..core.labels import Labels

logger = logging.getLogger(__name__)

def lawofcos(a: np.ndarray, b: np.ndarray) -> float:
    """
    return the cosine of the angle defined by vectors a and b if they share a vertex (the LAW OF COSINES)

    Parameters
    ----------
    a : np.ndarray
        First vector
    b : np.ndarray
        Second vector
    """
    return np.dot(a, b) / np.sqrt(np.dot(a, a) * np.dot(b, b))

def build_tmat(RotMat: np.ndarray, TransVec: np.ndarray) -> np.ndarray:
    """
    Builds a 4 x 4 homogeneous transformation matrix 
    
    Parameters
    ----------
    RotMat: numpy.ndarray
        3 x 3 rotation matrix
    TransVec: numpy.ndarray
        translation vector
    """
    tmat = np.identity(4, dtype=float)
    for i in range(3):
        for j in range(3):
            tmat[i][j] = RotMat[i][j]
        tmat[i][3] = TransVec[i]
    return tmat

def measure_dihedral(a1, a2, a3, a4):
    """
    Measure dihedral angle IN RADIANS of a1->a2--a3->a4
    
    Parameters
    ----------
    a1, a2, a3, a4: Atom
       Atom objects with x, y, z attributes

    Returns
    -------
    float
         Dihedral angle in radians between the four atoms a1, a2, a3, a4 in radians.
    """
    v1 = np.array([a1.x, a1.y, a1.z])
    v2 = np.array([a2.x, a2.y, a2.z])
    v3 = np.array([a3.x, a3.y, a3.z])
    v4 = np.array([a4.x, a4.y, a4.z])
    d12 = v1 - v2
    d23 = v2 - v3
    d34 = v3 - v4
    c12 = np.cross(d12, d23)
    c23 = np.cross(d23, d34)
    c123 = np.cross(d23, c12)
    c12 /= np.linalg.norm(c12)
    c23 /= np.linalg.norm(c23)
    c123 /= np.linalg.norm(c123)
    cosP = np.dot(c12, c23)
    sinP = np.dot(c23, c123)
    if np.isclose(cosP, 1.0):
        phi = np.pi
    elif np.isclose(cosP, 0.0):
        phi = 0.0
    else:
        phi = np.arccos(cosP)
    if sinP < 0:
        phi *= -1
    return phi


def positionN(res, tmat):
    """
    Given residue res, calculate the nominal position of the amide nitrogen
    in the next residue based on the positions of CA, C, and O.
       
    Parameters
    ----------
    res: Residue
       Residue object with atoms CA, C, and O (or OT1)
    tmat: numpy.ndarray
       4x4 homogeneous transformation matrix

    Returns
    -------
    rN: position of N atom (np.ndarray(3))
    """
    CA = res.atoms.get(lambda x: x.name == 'CA')
    C = res.atoms.get(lambda x: x.name == 'C')
    O = res.atoms.get(lambda x: x.name == 'O')
    if not O:
        logger.debug(f'Is this a C-terminus?')
        O = res.atoms.get(lambda x: x.name == 'OT1')
    rCA = np.dot(tmat, np.array([CA.x, CA.y, CA.z, 1.0]))[:3]
    rC = np.dot(tmat, np.array([C.x, C.y, C.z, 1.0]))[:3]
    rO = np.dot(tmat, np.array([O.x, O.y, O.z, 1.0]))[:3]
    R21 = rC - rCA
    r21 = R21 / np.linalg.norm(R21)
    R32 = rO - rC
    r32 = R32 / np.linalg.norm(R32)
    c = np.cross(r21, r32)
    mat = np.array([c, r21, r32])
    b = np.array([0, -np.cos(np.pi / 180.0 * 114.44), np.cos(np.pi / 180.0 * 123.04)])
    amat = np.linalg.inv(mat)
    rnd = np.dot(amat, b)
    rN = rC + rnd * 1.355
    R34 = rC - rN
    logger.debug(f'positionN: C-N bond length {np.linalg.norm(R34):.4f}')
    return rN

def coorddf_from_pdb(pdb, segtypes=False):
    """
    Read a PDB file and return a DataFrame with atom coordinates and other information.

    Parameters
    ----------
    pdb : str
        Path to the PDB file.
    segtypes : bool, optional
        If True, include segment types in the DataFrame. Default is False.

    Returns
    -------
    pd.DataFrame
        DataFrame with atom coordinates and other information.
    """
    p=PDBParser(filepath=pdb).parse()
    atlist=p.parsed['ATOM']+p.parsed.get('HETATM',[])
    serial = [x.serial for x in atlist]
    name = [x.name for x in atlist]
    x = [x.x for x in atlist]
    y = [x.y for x in atlist]
    z = [x.z for x in atlist]
    alt = [x.altLoc for x in atlist]
    resname = [x.residue.resName for x in atlist]
    resid = [x.residue.seqNum for x in atlist]
    chain = [x.residue.chainID for x in atlist]
    ins = [x.residue.iCode for x in atlist]
    basedict = {'name': name, 'x': x, 'y': y, 'z': z, 'resname': resname, 'resid': resid, 'chain': chain, 'altloc': alt, 'insertion': ins}
    if segtypes:
        basedict['segtype'] = [Labels.segtype_of_resname[x] for x in resname]
    return pd.DataFrame(basedict, index=serial)

def mic_shift(point, ref, box):
    """
    Given a point, a reference point, and a box, return the point shifted
    into the periodic box defined by the reference point.

    Parameters
    ----------
    point : np.ndarray
        The point to be shifted, as a 3-element array.
    ref : np.ndarray
        The reference point, as a 3-element array.
    box : np.ndarray
        The box dimensions, as a 3x3 array where each row is a basis vector
        defining the periodic box.
    """
    hbox = np.diagonal(box) / 2
    cpoint = point.copy()
    d = cpoint - ref
    boxlengths = np.zeros(3, dtype=int)
    for i in range(3):
        while d[i] < -hbox[i]:
            cpoint[i] += box[i][i]
            d = cpoint - ref
            boxlengths[i] += 1
        while d[i] >= hbox[i]:
            cpoint[i] -= box[i][i]
            d = cpoint - ref
            boxlengths[i] -= 1
    return cpoint, boxlengths