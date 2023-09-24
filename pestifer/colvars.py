#Author Cameron F. Abrams, <cfa22@drexel.edu>
"""functions for handling the colvars module of namd
"""

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
    name=data['name']
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

def declare_harmonic_distance_bias(data,writer):
    """writes the harmonic bias for a colvar"""
    name=data['name']
    k=data['k']
    init_distance=data['distance']
    targ_distance=data['targ_distance']
    targ_numsteps=data['targ_numsteps']
    writer.addline( 'harmonic {')
    writer.addline(f'    colvars {name}')
    writer.addline(f'    forceConstant {k}')
    writer.addline(f'    centers {init_distance}')
    writer.addline(f'    targetCenters {targ_distance}')
    writer.addline(f'    targetNumSteps {targ_numsteps}')
    writer.addline( '}')