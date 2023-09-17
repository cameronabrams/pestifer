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