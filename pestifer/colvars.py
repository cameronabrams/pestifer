def declare_distance_cv_atoms(data,writer):
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

def declare_harmonic_distance_bias(data,writer,**options):
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