# Author: Cameron F. Abrams <cfa22@drexel.edu>.
""" Some calculations based on atom coordinates
"""
import logging
import numpy as np
logger=logging.getLogger(__name__)

def measure_dihedral(a1,a2,a3,a4):
    """Measure dihedral angle IN RADIANS of a1->a2--a3->a4"""
    v1=np.array([a1.x,a1.y,a1.z])
    v2=np.array([a2.x,a2.y,a2.z])
    v3=np.array([a3.x,a3.y,a3.z])
    v4=np.array([a4.x,a4.y,a4.z])
    d12=v1-v2
    d23=v2-v3
    d34=v3-v4
    c12=np.cross(d12,d23)
    c23=np.cross(d23,d34)
    c123=np.cross(d23,c12)
    c12/=np.linalg.norm(c12)
    c23/=np.linalg.norm(c23)
    c123/=np.linalg.norm(c123)
    cosP=np.dot(c12,c23)
    sinP=np.dot(c23,c123)
    if np.isclose(cosP,1.0):
        phi=np.pi
    elif np.isclose(cosP,0.0):
        phi=0.0
    else:
        phi=np.arccos(cosP)
    if sinP<0:
        phi*=-1
    return phi

def ic_reference_closest(res12,ICmaps):
    """Given the two Residues in res12 and the maps in ICmaps, 
    return the mapping key to which the given IC values are
    closest in a Euclidean sense
    
    Parameters
    ----------
    res12: list
       exactly two Residue objects which must have lists of atoms attributes
    
    ICMaps: list
       list of dicts of the following format
          key           desc
          ---           ----
          ICatomnames   list of atom names as they appear in the charmff IC
          mapping       dictionary of patch name to IC value

    This method will identify the four atoms of the IC and reference
    them directly when calling the measure_dihedral function. 
    
    The list of computed dihedral values from the set of atoms is a "point"
    in "IC-space", and each patch has its own "reference point" in this space.

    The reference point to which the point is closest is identified as the 
    desired result.

    """
    for ic in ICmaps:
        logger.debug(f'icmap {ic}')
        ic['atoms']=[]
        for n in ic['ICatomnames']:
            r=int(n[0])-1
            an=n[1:]
            at=res12[r].atoms.get(name=an)
            ic['atoms'].append(at)
            # logger.debug(f'Assigned atom {at.name} of {at.resname}{at.resseqnum}')
    map_points={}
    the_point=[]
    for ic in ICmaps:
        value=measure_dihedral(*(ic['atoms']))*180.0/np.pi
        # logger.debug(f'{ic["ICatomnames"]} value {value:.2f}')
        the_point.append(value)
        for m,v in ic['mapping'].items():
            if not m in map_points:
                map_points[m]=[]
            map_points[m].append(v)
    for k,v in map_points.items():
        map_points[k]=np.array(v)
    the_point=np.array(the_point)
    logger.debug(f'ic the point: {the_point}')
    # calculate Euclidean distance adhering to the periodicity
    # of dihedral-angle space
    displacements={k:(the_point-v) for k,v in map_points.items()}
    for n,d in displacements.items():
        for i in range(len(d)):
            if d[i]<180.0:
                d[i]+=180.0
            if d[i]>180.0:
                d[i]-=180.0
    norms={k:np.linalg.norm(d) for k,d in displacements.items()}
    logger.debug(f'norms {norms}')
    the_one=[k for (k,v) in sorted(norms.items(),key=lambda x: x[1])][0]
    logger.debug(f'returning {the_one}')
    return the_one

def positionN(res,tmat):
    """Given residue res, calculate the nominal position of the amide nitrogen
       in the next residue based on the positions of CA, C, and O.
       
    Parameters
    ----------
    res: Residue
    tmat: 4x4 homogeneous transformation matrix
    
    Returns
    -------
    rN: position of N atom (np.ndarray(3))
    """
    CA=res.atoms.get(name='CA')
    C=res.atoms.get(name='C')
    O=res.atoms.get(name='O')
    if not O:
        logger.debug(f'Is this a C-terminus?')
        O=res.atoms.get(name='OT1')
    rCA=np.dot(tmat,np.array([CA.x,CA.y,CA.z,1.0]))[:3]
    rC=np.dot(tmat,np.array([C.x,C.y,C.z,1.0]))[:3]
    rO=np.dot(tmat,np.array([O.x,O.y,O.z,1.0]))[:3]
    R21=rC-rCA
    r21=R21/np.linalg.norm(R21)
    R32=rO-rC
    r32=R32/np.linalg.norm(R32)
    c=np.cross(r21,r32)
    mat=np.array([c,r21,r32])
    b=np.array([0,-np.cos(np.pi/180.0*114.44),np.cos(np.pi/180.0*123.04)])
    amat=np.linalg.inv(mat)
    rnd=np.dot(amat,b)
    rN=rC+rnd*1.355
    R34=rC-rN
    logger.debug(f'positionN: C-N bond length {np.linalg.norm(R34):.4f}')
    return rN

def coorddf_from_pdb(pdb,segtypes=False):
    import pandas as pd
    import os
    from pidibble.pdbparse import PDBParser
    p=PDBParser(PDBcode=os.path.splitext(pdb)[0]).parse()
    atlist=p.parsed['ATOM']+p.parsed.get('HETATM',[])
    serial=[x.serial for x in atlist]
    name=[x.name for x in atlist]
    x=[x.x for x in atlist]
    y=[x.y for x in atlist]
    z=[x.z for x in atlist]
    alt=[x.altLoc for x in atlist]
    resname=[x.residue.resName for x in atlist]
    resid=[x.residue.seqNum for x in atlist]
    chain=[x.residue.chainID for x in atlist]
    ins=[x.residue.iCode for x in atlist]
    basedict={'name':name,'x':x,'y':y,'z':z,'resname':resname,'resid':resid,'chain':chain,'altloc':alt,'insertion':ins}
    if segtypes:
        from .config import segtype_of_resname, Config
        if len(segtype_of_resname)==0: c=Config()
        basedict['segtype']=[segtype_of_resname[x] for x in resname]
    return pd.DataFrame(basedict,index=serial)

def mic_shift(point,ref,box):
    hbox=np.diagonal(box)/2
    cpoint=point.copy()
    d=cpoint-ref
    boxlengths=np.zeros(3,dtype=int)
    for i in range(3):
        while d[i]<-hbox[i]:
            cpoint[i]+=box[i][i]
            d=cpoint-ref
            boxlengths[i]+=1
        while d[i]>=hbox[i]:
            cpoint[i]-=box[i][i]
            d=cpoint-ref
            boxlengths[i]-=1
    return cpoint,boxlengths