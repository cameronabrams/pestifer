#Author Cameron F. Abrams, <cfa22@drexel.edu>
"""
Functions for handling the colvars module of NAMD.
This module provides functions to write colvar specifications to a script writer.
It includes functions to declare distance collective variables, harmonic distance biases,
and to handle groups of atoms specified by their names or serial numbers in a PDB file.
It also includes functions to write colvar specifications for single harmonic distance biases.
The module uses the PDBParser to parse PDB files and extract atom names and serial numbers.
It ensures that atom names are unique and that the specified groups contain valid atoms.
It also provides functions to write colvar specifications for distances between groups of atoms.
"""
from pidibble.pdbparse import PDBParser
import os
import logging
from copy import deepcopy
from ..core.scripters import Filewriter
logger=logging.getLogger(__name__)

def colvar_writer(specs,scripter:Filewriter,pdb=''):
    """
    Writes the colvar specifications to a script writer.
    This function processes the specifications for groups of atoms, distances, and harmonic biases,
    and writes the corresponding colvar declarations to the provided script writer.
    
    Parameters
    ----------
    specs : dict
        A dictionary containing the colvar specifications. It should include ``groups``, ``distances``,
        and ``harmonics`` keys, each containing relevant specifications.
    scripter : Filewriter
        An instance of a script writer that will be used to write the colvar specifications.
    pdb : str, optional
        The path to a PDB file from which atom names and serial numbers will be extracted.
        If provided, the function will use this PDB file to resolve atom names and serial numbers
        for the groups specified in the colvar specifications.
        If not provided, it will assume that atom serial numbers are already specified in the groups.
        Default is an empty string.
    """
    logger.debug(f'{specs}')
    local_specs=deepcopy(specs)
    atom_names=[]
    atom_serials=[]
    if pdb:
        p=PDBParser(PDBcode=os.path.splitext(pdb)[0]).parse()
        atoms=p.parsed['ATOM']
        atom_names=[a.name for a in atoms]
        atom_serials=[a.serial for a in atoms]

    groups=local_specs.get('groups',{})
    for groupname,groupspecs in groups.items():
        atomnames=groupspecs.get('atomnames',[])
        serials=groupspecs.get('serials',[])
        if atomnames:
            for a in atomnames:
                assert atom_names.count(a)==1,f'Cannot find unique atom {a}; you should use serial numbers'
                serials.append(atom_serials[atom_names.index(a)])
        groupspecs['atomNumbers']=' '.join([f'{i}' for i in serials])
        assert len(groupspecs['atomNumbers'])>0,f'No atoms in group {groupname}'

    distances=local_specs.get('distances',{})
    logger.debug(f'distances {distances}')
    for distancename,distancespecs in distances.items():
        for i in range(len(distancespecs['groups'])):
            distancespecs['groups'][i]=groups[distancespecs['groups'][i]]
            distancespecs['name']=distancename
        declare_distance_cv(distancespecs,scripter)

    harmonics=local_specs.get('harmonics',{})
    for harmonicname,harmonicspecs in harmonics.items():
        harmonicspecs['name']=harmonicname
        declare_harmonic_distance_biases(harmonicspecs,scripter)

def declare_distance_cv(data,scripter:Filewriter):
    """
    Writes the colvar stanza describing a distance collective variable defined by groups to a colvar input script

    Parameters
    ----------
    data: dict
        dictionary containing colvar specifications
    scripter: Filewriter
        An instance of a script writer that will be used to write the colvar specifications.
    """
    name=data['name']
    scripter.addline( 'colvar {')
    scripter.addline(f'    name {name}')
    scripter.addline( '    distance {')
    scripter.addline( '        group1 {')
    scripter.addline( '            atomNumbers { '+f'{data["groups"][0]["atomNumbers"]}'+' }')
    scripter.addline( '        }')
    scripter.addline( '        group2 {')
    scripter.addline( '            atomNumbers { '+f'{data["groups"][1]["atomNumbers"]}'+' }')
    scripter.addline( '        }')
    scripter.addline( '    }')
    scripter.addline( '}')

def declare_distance_cv_atoms(data,scripter:Filewriter):
    """
    Writes the colvar stanza describing a distance collective variable defined by two atoms specified by their serial numbers to a colvar input script.

    Parameters
    ----------
    data: dict
        dictionary containing serial numbers of atoms i and j and a colvar name
    scripter: Filewriter
        An instance of a script writer that will be used to write the colvar specifications.
    """
    i=data['serial_i']
    j=data['serial_j']
    name=data['colvars']
    scripter.addline( 'colvar {')
    scripter.addline(f'    name {name}')
    scripter.addline( '    distance {')
    scripter.addline( '        group1 {')
    scripter.addline( '            atomNumbers { '+f'{i}'+' }')
    scripter.addline( '        }')
    scripter.addline( '        group2 {')
    scripter.addline( '            atomNumbers { '+f'{j}'+' }')
    scripter.addline( '        }')
    scripter.addline( '    }')
    scripter.addline( '}')

def declare_harmonic_distance_biases(data,scripter):
    """
    Writes the harmonic distance bias for a colvar defined by two groups of atoms to a colvar input script.

    Parameters
    ----------
    data: dict
        dictionary containing colvar specifications, including the name of the colvar, the groups of atoms
        involved, the force constant, the distance between the groups, and optional target distance and number of steps.
    scripter: Filewriter
        An instance of a script writer that will be used to write the colvar specifications.
    """
    name=data['name']
    colvars=' '.join(data['colvars'])
    k=data['forceConstant']
    distance=' '.join([f'{x:.5f}' for x in data['distance']])
    targ_distance=' '.join([f'{x:.5f}' for x in data.get('targ_distance',[])])
    targ_numsteps=data.get('targ_numsteps','')
    scripter.addline( 'harmonic {')
    scripter.addline(f'    name {name}')
    scripter.addline(f'    colvars {colvars}')
    scripter.addline(f'    forceConstant {k}')
    scripter.addline(f'    centers {distance}')
    if targ_distance: scripter.addline(f'    targetCenters {targ_distance}')
    if targ_numsteps: scripter.addline(f'    targetNumSteps {targ_numsteps}')
    scripter.addline( '}')

def declare_single_harmonic_distance_bias(data,scripter:Filewriter):
    """
    Writes the harmonic bias for a colvar defined by a single group of atoms to a colvar input script.
    
    Parameters
    ----------
    data: dict
        dictionary containing colvar specifications, including the name of the colvar, the force constant,
        the initial distance, and optional target distance and number of steps.
    scripter: Filewriter
        An instance of a script writer that will be used to write the colvar specifications.
    """
    name=data['colvars']
    k=data['forceConstant']
    init_distance=data['distance']
    targ_distance=data.get('targ_distance','')
    targ_numsteps=data.get('targ_numsteps','')
    scripter.addline( 'harmonic {')
    scripter.addline(f'    colvars {name}')
    scripter.addline(f'    forceConstant {k}')
    scripter.addline(f'    centers {init_distance}')
    if targ_distance: scripter.addline(f'    targetCenters {targ_distance}')
    if targ_numsteps: scripter.addline(f'    targetNumSteps {targ_numsteps}')
    scripter.addline( '}')