#Author Cameron F. Abrams, <cfa22@drexel.edu>
"""functions for handling the colvars module of namd
"""
from pidibble.pdbparse import PDBParser
import os
import logging
logger=logging.getLogger(__name__)

def colvar_writer(specs,writer,pdb=''):
    logger.debug(f'{specs}')
    atom_names=[]
    atom_serials=[]
    if pdb:
        p=PDBParser(PDBcode=os.path.splitext(pdb)[0]).parse()
        atoms=p.parsed['ATOM']
        atom_names=[a.name for a in atoms]
        atom_serials=[a.serial for a in atoms]

    groups=specs.get('groups',{})
    for groupname,groupspecs in groups.items():
        atomnames=groupspecs.get('atomnames',[])
        serials=groupspecs.get('serials',[])
        if atomnames:
            for a in atomnames:
                assert atom_names.count(a)==1,f'Cannot find unique atom {a}; you should use serial numbers'
                serials.append(atom_serials[atom_names.index(a)])
        groupspecs['atomNumbers']=' '.join([f'{i}' for i in serials])
        assert len(groupspecs['atomNumbers'])>0,f'No atoms in group {groupname}'

    distances=specs.get('distances',{})
    logger.debug(f'distances {distances}')
    for distancename,distancespecs in distances.items():
        for i in range(len(distancespecs['groups'])):
            distancespecs['groups'][i]=groups[distancespecs['groups'][i]]
            distancespecs['name']=distancename
        declare_distance_cv(distancespecs,writer)

    harmonics=specs.get('harmonics',{})
    for harmonicname,harmonicspecs in harmonics.items():
        harmonicspecs['name']=harmonicname
        declare_harmonic_distance_biases(harmonicspecs,writer)

def declare_distance_cv(data,writer):
    name=data['name']
    writer.addline( 'colvar {')
    writer.addline(f'    name {name}')
    writer.addline( '    distance {')
    writer.addline( '        group1 {')
    writer.addline( '            atomNumbers { '+f'{data["groups"][0]["atomNumbers"]}'+' }')
    writer.addline( '        }')
    writer.addline( '        group2 {')
    writer.addline( '            atomNumbers { '+f'{data["groups"][1]["atomNumbers"]}'+' }')
    writer.addline( '        }')
    writer.addline( '    }')
    writer.addline( '}')

def declare_distance_cv_atoms(data,writer):
    """writes the colvar stanza to a colvar input script
    
    Parameters
    ----------
    data: dict
        dictionary containing serial numbers of atoms i and j and a colvar name
    writer: scriptwriter
    """
    i=data['serial_i']
    j=data['serial_j']
    name=data['colvars']
    writer.addline( 'colvar {')
    writer.addline(f'    name {name}')
    writer.addline( '    distance {')
    writer.addline( '        group1 {')
    writer.addline( '            atomNumbers { '+f'{i}'+' }')
    writer.addline( '        }')
    writer.addline( '        group2 {')
    writer.addline( '            atomNumbers { '+f'{j}'+' }')
    writer.addline( '        }')
    writer.addline( '    }')
    writer.addline( '}')

def declare_harmonic_distance_biases(data,writer):
    name=data['name']
    colvars=' '.join(data['colvars'])
    k=data['forceConstant']
    distance=' '.join([f'{x:.5f}' for x in data['distance']])
    targ_distance=' '.join([f'{x:.5f}' for x in data.get('targ_distance',[])])
    targ_numsteps=data.get('targ_numsteps','')
    writer.addline( 'harmonic {')
    writer.addline(f'    name {name}')
    writer.addline(f'    colvars {colvars}')
    writer.addline(f'    forceConstant {k}')
    writer.addline(f'    centers {distance}')
    if targ_distance: writer.addline(f'    targetCenters {targ_distance}')
    if targ_numsteps: writer.addline(f'    targetNumSteps {targ_numsteps}')
    writer.addline( '}') 

def declare_single_harmonic_distance_bias(data,writer):
    """writes the harmonic bias for a colvar"""
    name=data['colvars']
    k=data['forceConstant']
    init_distance=data['distance']
    targ_distance=data.get('targ_distance','')
    targ_numsteps=data.get('targ_numsteps','')
    writer.addline( 'harmonic {')
    writer.addline(f'    colvars {name}')
    writer.addline(f'    forceConstant {k}')
    writer.addline(f'    centers {init_distance}')
    if targ_distance: writer.addline(f'    targetCenters {targ_distance}')
    if targ_numsteps: writer.addline(f'    targetNumSteps {targ_numsteps}')
    writer.addline( '}')