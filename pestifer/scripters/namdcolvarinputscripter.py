# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import logging

from copy import deepcopy
from pathlib import Path
from pidibble.pdbparse import PDBParser

from .genericscripter import GenericScripter

from ..molecule.atom import AtomList

logger = logging.getLogger(__name__)

class NAMDColvarInputScripter(GenericScripter):  # TODO: refactor as a JSONScripter

    def __init__(self, *args, **kwargs):
        """
        Initialize a new NAMD colvar input script writer.
        
        Parameters
        ----------
        *args, **kwargs: Additional arguments and keyword arguments to pass to the GenericScripter base class.
        """
        super().__init__(*args, **kwargs)
        self.default_ext = '.in'
        self.specs = {}

    def construct_on_pdb(self, specs: dict = {}, pdb: Path | str = None, **kwargs):
        """
        Construct the NAMD colvar input script.

        Parameters
        ----------
        specs : dict
            A dictionary containing the specifications for the colvar input script.
        **kwargs : Additional keyword arguments to customize the script generation.
        """
        self.specs = deepcopy(specs)
        atom_names = []
        atom_serials = []
        if pdb:
            p = PDBParser(filepath=pdb).parse()
            atoms: AtomList = p.parsed['ATOM']
            atom_names = [a.name for a in atoms.data]
            atom_serials = [a.serial for a in atoms.data]

        groups = specs.get('groups', {})
        for groupname, groupspecs in groups.items():
            atomnames = groupspecs.get('atomnames', [])
            serials = groupspecs.get('serials', [])
            if atomnames:
                for a in atomnames:
                    assert atom_names.count(a) == 1, f'Cannot find unique atom {a}; you should use serial numbers'
                    serials.append(atom_serials[atom_names.index(a)])
            groupspecs['atomNumbers'] = ' '.join([f'{i}' for i in serials])
            assert len(groupspecs['atomNumbers']) > 0, f'No atoms in group {groupname}'

        distances = specs.get('distances', {})
        logger.debug(f'distances {distances}')
        for distancename, distancespecs in distances.items():
            for i in range(len(distancespecs['groups'])):
                distancespecs['groups'][i] = groups[distancespecs['groups'][i]]
                distancespecs['name'] = distancename
            self.declare_distance_cv(distancespecs)

        harmonics = specs.get('harmonics', {})
        for harmonicname, harmonicspecs in harmonics.items():
            harmonicspecs['name'] = harmonicname
            self.declare_harmonic_distance_biases(harmonicspecs)

    def declare_distance_cv(self, data):
        """
        Writes the colvar stanza describing a distance collective variable defined by groups to a colvar input script

        Parameters
        ----------
        data: dict
            dictionary containing colvar specifications
        """
        name = data['name']
        self.addline( 'colvar {')
        self.addline(f'    name {name}')
        self.addline( '    distance {')
        self.addline( '        group1 {')
        self.addline( '            atomNumbers { '+f'{data["groups"][0]["atomNumbers"]}'+' }')
        self.addline( '        }')
        self.addline( '        group2 {')
        self.addline( '            atomNumbers { '+f'{data["groups"][1]["atomNumbers"]}'+' }')
        self.addline( '        }')
        self.addline( '    }')
        self.addline( '}')

    def declare_distance_cv_atoms(self, data):
        """
        Writes the colvar stanza describing a distance collective variable defined by two atoms specified by their serial numbers to a colvar input script.

        Parameters
        ----------
        data: dict
            dictionary containing serial numbers of atoms i and j and a colvar name
        """
        i = data['serial_i']
        j = data['serial_j']
        name = data['colvars']
        self.addline( 'colvar {')
        self.addline(f'    name {name}')
        self.addline( '    distance {')
        self.addline( '        group1 {')
        self.addline( '            atomNumbers { '+f'{i}'+' }')
        self.addline( '        }')
        self.addline( '        group2 {')
        self.addline( '            atomNumbers { '+f'{j}'+' }')
        self.addline( '        }')
        self.addline( '    }')
        self.addline( '}')

    def declare_harmonic_distance_biases(self, data):
        """
        Writes the harmonic distance bias for a colvar defined by two groups of atoms to a colvar input script.

        Parameters
        ----------
        data: dict
            dictionary containing colvar specifications, including the name of the colvar, the groups of atoms
            involved, the force constant, the distance between the groups, and optional target distance and number of steps.
        """
        name = data['name']
        colvars = ' '.join(data['colvars'])
        k = data['forceConstant']
        distance = ' '.join([f'{x:.5f}' for x in data['distance']])
        targ_distance = ' '.join([f'{x:.5f}' for x in data.get('targ_distance', [])])
        targ_numsteps = data.get('targ_numsteps', '')
        self.addline( 'harmonic {')
        self.addline(f'    name {name}')
        self.addline(f'    colvars {colvars}')
        self.addline(f'    forceConstant {k}')
        self.addline(f'    centers {distance}')
        if targ_distance: self.addline(f'    targetCenters {targ_distance}')
        if targ_numsteps: self.addline(f'    targetNumSteps {targ_numsteps}')
        self.addline( '}')

    def declare_single_harmonic_distance_bias(self, data):
        """
        Writes the harmonic bias for a colvar defined by a single group of atoms to a colvar input script.
        
        Parameters
        ----------
        data: dict
            dictionary containing colvar specifications, including the name of the colvar, the force constant,
            the initial distance, and optional target distance and number of steps.
        """
        name = data['colvars']
        k = data['forceConstant']
        init_distance = data['distance']
        targ_distance = data.get('targ_distance', '')
        targ_numsteps = data.get('targ_numsteps', '')
        self.addline( 'harmonic {')
        self.addline(f'    colvars {name}')
        self.addline(f'    forceConstant {k}')
        self.addline(f'    centers {init_distance}')
        if targ_distance: self.addline(f'    targetCenters {targ_distance}')
        if targ_numsteps: self.addline(f'    targetNumSteps {targ_numsteps}')
        self.addline( '}')