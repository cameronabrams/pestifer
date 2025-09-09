# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A module for handling bilayer structures in molecular simulations.
"""

import logging

import numpy as np

from ..charmmff.charmmffcontent import CHARMMFFContent
from ..core.artifacts import ArtifactDict
from ..util.stringthings import my_logger
from ..util.units import _UNITS_, _SYMBOLS_, cuA_of_nmolec
from ..scripters import PackmolScripter

sA_  = _SYMBOLS_['ANGSTROM']
sA2_ = _UNITS_['SQUARE-ANGSTROMS']
sA3_ = _UNITS_['CUBIC-ANGSTROMS']

logger = logging.getLogger(__name__)

def orthohexagonal_cell(a: float, L_star: float, origin: tuple[float, float]=(0.0, 0.0)):
    m = round(L_star / a)  # enforce matching x-side
    cands_n = [max(1, round(m / np.sqrt(3)) + dn) for dn in (-1, 0, 1)]
    best = None
    for n in cands_n:
        Lx = m*a
        Ly = n*np.sqrt(3)*a
        diff = abs(Lx - Ly)
        if best is None or diff < best[0]:
            best = (diff, n, Lx, Ly)
    diff, n, Lx, Ly = best
    pts_R0 = [(r*a + origin[0], q*np.sqrt(3)*a + origin[1]) for q in range(n) for r in range(m)]
    pts_shift = [(r*a + 0.5*a + origin[0], q*np.sqrt(3)*a + 0.5*np.sqrt(3)*a + origin[1]) for q in range(n) for r in range(m)]
    return pts_R0 + pts_shift, Lx, Ly

def squareish_lattice(Lx: float, Ly: float, m: int):
    cands_n = [m + dn for dn in (-4, -3, -2, -1, 0, 1, 2, 3, 4)]
    best = None
    for n in cands_n:
        tLx = m*Lx
        tLy = n*Ly
        diff = abs(tLx - tLy)
        if best is None or diff < best[0]:
            best = (diff, n, tLx, tLy)
    return m, best[1]

def orthohexagonal_patch(SAPL: float, lipid_count_ceiling = 100):
    a = np.sqrt(2*SAPL/(3*np.sqrt(3)))  # area per point = SAPL
    best = None
    for nsites_1d in range(4,12):
        L_star = round(nsites_1d*a)
        pts_R0, Lx, Ly = orthohexagonal_cell(a, L_star)
        npts = len(pts_R0)
        if best is None or npts < lipid_count_ceiling:
            best = (npts, Lx, Ly, pts_R0)
    return best


class BilayerSpecString:
    """ 
    A class for handling bilayer specification strings in memgen format.
    The specification string is a string that describes the composition of the bilayer in terms of species and their fractions.
    
    Parameters
    ----------
    specstring : str, optional
        The packmol-memgen-format specification string for the bilayer.
    fracstring : str, optional
        The mole-fraction specification string for the bilayer.
    leaflet_delimiter : str, optional
        The delimiter used to separate the left and right leaflets in the specification string.
    species_delimiter : str, optional
        The delimiter used to separate species in the specification string.
    """
    def __init__(self, specstring='', fracstring='', leaflet_delimiter='//', species_delimiter=':'):
        self.specstring = specstring
        self.fracstring = fracstring
        self.leaflet_delimiter = leaflet_delimiter
        self.species_delimiter = species_delimiter
        self.extrastrings = {}
        if specstring:
            Sleft, Sright = (specstring.split(leaflet_delimiter) + [specstring])[:2]
            Lleft, Lright = Sleft.split(species_delimiter), Sright.split(species_delimiter)
            nSpleft, nSright = len(Lleft), len(Lright)
            if fracstring:
                FSleft, FSright = (fracstring.split(leaflet_delimiter) + [fracstring])[:2]
                Fleft, Fright = np.array([float(x) for x in FSleft.split(species_delimiter)]), np.array([float(x) for x in FSright.split(species_delimiter)])
                Fleft /= np.sum(Fleft)
                Fright /= np.sum(Fright)
                assert len(Fleft) == nSpleft, f'Left layer has {nSpleft} species but {len(Fleft)} fractions specified'
                assert len(Fright) == nSright, f'Right layer has {nSright} species but {len(Fright)} fractions specified'
            else:
                Fleft = np.array([1.0/nSpleft]*nSpleft)
                Fright = np.array([1.0/nSright]*nSright)

            self.left = [dict(name=n, frac=x) for n, x in zip(Lleft, Fleft)]
            self.right = [dict(name=n, frac=x) for n, x in zip(Lright, Fright)]

    def add_specstring(self, attr_name, specstring='', attr_type=str):
        """
        Adds a specification string to the BilayerSpecString object.

        This method updates the internal state of the object to include the new specification string.

        Parameters
        ----------
        attr_name : str
            The name of the attribute to which the specification string will be added.
        specstring : str, optional
            The specification string to be added.
        attr_type : type, optional
            The type to which the specification string will be converted.
        """
        if specstring:
            self.extrastrings[attr_name] = specstring
            Sleft, Sright = (specstring.split(self.leaflet_delimiter) + [specstring])[:2]
            Lleft, Lright = [attr_type(x) for x in Sleft.split(self.species_delimiter)], [attr_type(x) for x in Sright.split(self.species_delimiter)]
            nSpleft, nSright = len(Lleft), len(Lright)

            if nSpleft == len(self.left):
                for i in range(nSpleft):
                    self.left[i][attr_name] = Lleft[i]
            elif nSpleft == 1:
                for i in range(len(self.left)):
                    self.left[i][attr_name] = Lleft[0]
            if nSright == len(self.right):
                for i in range(nSright):
                    self.right[i][attr_name] = Lright[i]
            elif nSright == 1:
                for i in range(len(self.right)):
                    self.right[i][attr_name] = Lright[0]

def specstrings_builddict(lipid_specstring='', lipid_ratio_specstring='', lipid_conformers_specstring='0',
                          solvent_specstring='TIP3', solvent_ratio_specstring=''):
    """
    Builds a dictionary of bilayer specifications from the provided specification strings.
    
    Parameters
    ----------
    lipid_specstring : str, optional
        The specification string for the lipid bilayer.
    lipid_ratio_specstring : str, optional
        The mole-fraction specification string for the lipid bilayer.
    lipid_conformers_specstring : str, optional
        The conformer specification string for the lipid bilayer.
    solvent_specstring : str, optional
        The specification string for the solvent.
    solvent_ratio_specstring : str, optional
        The mole-fraction specification string for the solvent.

    Returns
    -------
    dict
        A dictionary containing the bilayer specifications.
    """
    L = BilayerSpecString(specstring=lipid_specstring, fracstring=lipid_ratio_specstring)
    L.add_specstring('conf', lipid_conformers_specstring, int)
    C = BilayerSpecString(specstring=solvent_specstring, fracstring=solvent_ratio_specstring)
    return {
        'upper_leaflet': L.left,
        'lower_leaflet': L.right,
        'upper_chamber': C.left,
        'lower_chamber': C.right
    }

class Bilayer:
    """
    A class for handling bilayer structures in molecular simulations.
    This class represents a bilayer composed of lipids and solvent, with specifications for each leaflet and chamber.

    Parameters
    ----------
    composition_dict : dict, optional
        A dictionary containing the composition of the bilayer, including leaflets and chambers.
    leaflet_nlipids : dict, optional
        A dictionary specifying the number of lipids per leaflet in a patch.
        Default is {'upper': 100, 'lower': 100}.
    solvent_to_key_lipid_ratio : float, optional
        The ratio of solvent molecules to key lipid molecules in the bilayer.
        Default is 32.0.
    neutralizing_salt : list, optional
        A list containing the names of the cation and anion used for neutralizing the bilayer. Default is ['POT', 'CLA'].
    salt_concentration : float, optional
        The concentration of salt in the bilayer solution, in molarity (M). Default is 0.0.
    solution_gcc : float, optional
        The density of the solution in grams per cubic centimeter (gcc). Default is 1.0.
    charmmffcontent: CHARMMFFContent object
        The CHARMM force field content object that allows access to PDB structures and RESI/PATCH information.
    solvent_specstring : str, optional
        The specification string for the solvent in the bilayer.
        Default is 'TIP3'.
    solvent_ratio_specstring : str, optional
        The mole-fraction specification string for the solvent.
        Default is '1.0'.
    """

    def __init__(self, composition_dict: dict = {}, 
                 leaflet_nlipids: dict[str, int] = dict(upper=100, lower=100), 
                 solvent_to_key_lipid_ratio: float = 32.0,
                 neutralizing_salt: list[str] = ['POT', 'CLA'], 
                 salt_concentration: float = 0.0, 
                 solution_gcc: float = 1.0, 
                 charmmffcontent: CHARMMFFContent = None, 
                 solvent_specstring: str = 'TIP3', 
                 solvent_ratio_specstring: str = '1.0'):

        # leaflet_nlipids is the number of lipids per leaflet in a patch
        self.leaflet_nlipids = leaflet_nlipids
        self.charmmffcontent = charmmffcontent
        # attributes that will be set later
        self.area = 0.0
        self.artifacts = ArtifactDict()

        if not composition_dict:
            logger.debug('Empty bilayer')
            return None
        
        # user need not have specified the solvent composition in the upper and lower chambers
        if 'upper_chamber' not in composition_dict or 'lower_chamber' not in composition_dict:
            logger.debug(f'Provided composition dictionary does not include solvent chamber specifications')
            logger.debug(f'Using specstrings \'{solvent_specstring}\' and \'{solvent_ratio_specstring}\' to build solvent composition')
            C = BilayerSpecString(specstring=solvent_specstring, fracstring=solvent_ratio_specstring)
            if 'upper_chamber' not in composition_dict:
                composition_dict['upper_chamber'] = C.left
            if 'lower_chamber' not in composition_dict:
                composition_dict['lower_chamber'] = C.right

        lipid_names = [x['name'] for x in composition_dict['upper_leaflet']]
        lipid_names += [x['name'] for x in composition_dict['lower_leaflet']]
        self.lipid_names = list(set(lipid_names))

        solvent_names = [x['name'] for x in composition_dict['upper_chamber']]
        solvent_names += [x['name'] for x in composition_dict['lower_chamber']]
        self.solvent_names = list(set(solvent_names + neutralizing_salt))
        self.species_names = self.lipid_names + self.solvent_names

        # Tell the CHARMMFFContent object to expose required PDBs and residue topologies
        self.charmmffcontent.provision(resnames=self.species_names)

        # complete leaflet entries in composition dictionary with species counts, charges, and MWs
        for l in ['upper_leaflet', 'lower_leaflet']:
            adjective = 'upper' if l == 'upper_leaflet' else 'lower'
            L = composition_dict[l]
            logger.debug(f'Leaflet {l} composition: {L}')
            for d in L:
                resi = self.charmmffcontent.get_resi(d['name'])
                if resi is None:
                    raise ValueError(f'Cannot find residue {d["name"]} in CHARMMFFContent')
                if not 'patn' in d:
                    d['patn'] = int(d['frac'] * leaflet_nlipids[adjective])
                if not 'charge' in d:
                    d['charge'] = resi.charge
                if not 'MW' in d:
                    d['MW'] = resi.mass


        # complete chamber entries in composition dictionary with species counts, charges, and MWs
        for c in ['upper_chamber', 'lower_chamber']:
            adjective = 'upper' if c == 'upper_chamber' else 'lower'
            L = composition_dict[c]
            Nsol = 0
            AMW = 0.0
            for d in L:
                resi = self.charmmffcontent.get_resi(d['name'])
                if not 'patn' in d:
                    d['patn'] = int(d['frac'] * leaflet_nlipids[adjective] * solvent_to_key_lipid_ratio)
                if not 'charge' in d:
                    d['charge'] = resi.charge
                if not 'MW' in d:
                    d['MW'] = resi.mass
                Nsol += d['patn']
                AMW += d['MW'] * d['patn']
            if Nsol > 0:
                AMW /= Nsol
            else:
                AMW = 0.0
            if salt_concentration > 0.0:
                Npm = int(AMW / 1000 * salt_concentration / solution_gcc * Nsol)
                if Npm > 0:
                    cation_name, anion_name = neutralizing_salt
                    logger.debug(f'Salting at {salt_concentration} M, soln density {solution_gcc} gcc')
                    logger.debug(f'-> adding {Npm} {cation_name} and {Npm} {anion_name} to {c}')
                    cation, anion = self.charmmffcontent.get_resi(cation_name), self.charmmffcontent.get_resi(anion_name)
                    n_cation = int(np.round(Npm / np.abs(cation.charge), 0))
                    n_anion = int(np.round(Npm / np.abs(anion.charge), 0))
                    composition_dict[c].append({'name': cation_name, 'patn': n_cation, 'charge': cation.charge, 'MW': cation.mass})
                    composition_dict[c].append({'name': anion_name, 'patn': n_anion, 'charge': anion.charge, 'MW': anion.mass})
        # set up some short-cut object labes
        self.slices = {'lower_chamber': {}, 'lower_leaflet': {}, 'upper_leaflet': {}, 'upper_chamber': {}}
        self.LC = self.slices['lower_chamber']
        self.LL = self.slices['lower_leaflet']
        self.UL = self.slices['upper_leaflet']
        self.UC = self.slices['upper_chamber']
        for layer, data in self.slices.items():
            data['composition'] = composition_dict[layer]
        # if the bilayer is asymmetric (each leaflet has a unique composition), we cannot assume a priori that
        # each leaflet in a patch has the same number of lipids.  We set a flag to indicate that the patch is 
        # asymmetric and that the number of lipids in each leaflet may need to be adjusted after equilibration
        # and measurment of the pressure profile.
        ul_lx, ll_lx = [(x['name'], x['frac']) for x in self.slices['upper_leaflet']['composition']], [(x['name'], x['frac']) for x in self.slices['lower_leaflet']['composition']]
        self.asymmetric = set(ul_lx) != set(ll_lx)

        self.species_data = {}
        self.addl_streamfiles = []
        pdbrepository = self.charmmffcontent.pdbrepository if self.charmmffcontent is not None else None
        if pdbrepository is not None:
            for l in self.species_names:
                logger.debug(f'Getting pdb for {l}')
                if not l in pdbrepository:
                    raise Exception(f'Cannot find {l} in PDB repository')
                pdbstruct = pdbrepository.checkout(l)
                self.species_data[l] = pdbstruct
                for p in self.species_data[l].get_parameters():
                    if p.endswith('.str') and not p in self.addl_streamfiles:
                        self.addl_streamfiles.append(p)
        logger.debug(f'Additional stream files:')
        my_logger(self.addl_streamfiles, logger.debug)

        self.total_charge = 0.0
        for layer, data in self.slices.items():
            data['charge'] = 0.0
            data['maxthickness'] = 0.0
            data['composition'] = composition_dict[layer]
            data['avgMW'] = 0.0
            data['patn'] = 0
            for species in data['composition']:
                logger.debug(f'Layer {layer} species {species["name"]} patn {species["patn"]} charge {species["charge"]} MW {species["MW"]}')
                data['patn'] += species['patn']
                data['charge'] += species['charge'] * species['patn']
                if 'chamber' in layer:
                    data['avgMW'] += species['MW'] * species['patn']
                elif 'leaflet' in layer:
                    for lipid in data['composition']:
                        lipid['reference_length'] = self.species_data[lipid['name']].get_head_tail_length(conformerID=lipid.get('conf', 0))
                        if lipid['reference_length'] > data['maxthickness']:
                            data['maxthickness'] = lipid['reference_length']
                    logger.debug(f'{layer} maxthickness {data["maxthickness"]:.3f}')
            logger.debug(f'Layer {layer} patn {data["patn"]} charge {data["charge"]:.3f} e avgMW*N {data["avgMW"]:.3f} g/mol')
            self.total_charge += data['charge']
            data['avgMW'] /= data['patn']

        if self.total_charge != 0.0:
            logger.debug(f'Total charge of bilayer is {self.total_charge:.3f} e')
            cation_name, anion_name = neutralizing_salt
            if self.total_charge > 0.0:  # need to include anion
                ion_name = anion_name
            else:
                ion_name = cation_name
            if ion_name not in self.species_names:
                self.species_names.append(ion_name)
                self.species_data[ion_name] = pdbrepository.checkout(ion_name)
            ion_resi = self.charmmffcontent.get_resi(ion_name)
            ion_q = ion_resi.charge
            logger.debug(f'Adding {ion_name} with charge {ion_q:.3f} e')
            ion_patn = int(np.round(np.abs(self.total_charge), 0) / np.abs(ion_q))
            if ion_patn % 2 == 0:
                lc_ion_patn = ion_patn // 2
                uc_ion_patn = ion_patn // 2
            else:
                lc_ion_patn = ion_patn // 2
                uc_ion_patn = ion_patn // 2 + 1
            for chamber, nions in zip(['upper_chamber', 'lower_chamber'], [uc_ion_patn, lc_ion_patn]):
                chamber_names = [x['name'] for x in self.slices[chamber]['composition']]
                if ion_name not in chamber_names:
                    self.slices[chamber]['composition'].append({'name': ion_name, 'patn': nions, 'charge': ion_q, 'MW': ion_resi.mass})
                else:
                    # if the ion is already in the chamber, just add to the number of ions
                    for species in self.slices[chamber]['composition']:
                        if species['name'] == ion_name:
                            species['patn'] += nions
                            break

        # finally, check out all required PDB input files for packmol
        self.register_species_pdbs = []
        for layer, data in self.slices.items():
            for species in data['composition']:
                species_name = species['name']
                conformerID = species.get('conf', 0)
                noh = species.get('noh', False)
                species['local_name'] = self.species_data[species_name].get_pdb(conformerID=conformerID, noh=noh)
                if species['local_name'] not in self.register_species_pdbs:
                    self.register_species_pdbs.append(species['local_name'])
                # logger.debug(f'Checked out {species_name} as {species["local_name"]}')

    def spec_out(self, SAPL=75.0, xy_aspect_ratio=1.0, half_mid_zgap=1.0, solution_gcc=1.0, rotation_pm=10.0):
        """
        Specs out a patch of the bilayer with specified parameters.

        Parameters
        ----------
        SAPL : float, optional
            The surface area per lipid in Å². Default is 75.0.
        xy_aspect_ratio : float, optional
            The aspect ratio of the patch in the x and y dimensions. Default is 1.0.
        half_mid_zgap : float, optional
            The half mid-plane gap in Å. Default is 1.0 Å.
        solution_gcc : float, optional
            The density of the solution in grams per cubic centimeter (gcc). Default is 1.0.
        rotation_pm : float, optional
            The rotation angle in degrees for the patch. Default is 10.0 degrees.
        """
        patch_area = SAPL * self.leaflet_nlipids['upper']  # assume symmetric
        Lx = np.sqrt(patch_area / xy_aspect_ratio)
        Ly = xy_aspect_ratio * Lx
        self.patch_area = patch_area
        lc_vol = cuA_of_nmolec(self.LC['avgMW'], solution_gcc, self.LC['patn'])
        uc_vol = cuA_of_nmolec(self.UC['avgMW'], solution_gcc, self.UC['patn'])
        lc_thickness = lc_vol / self.patch_area
        uc_thickness = uc_vol / self.patch_area
        ll_maxthickness = self.LL['maxthickness']
        ul_maxthickness = self.UL['maxthickness']
        ll_actthickness = np.cos(np.deg2rad(rotation_pm)) * ll_maxthickness
        ul_actthickness = np.cos(np.deg2rad(rotation_pm)) * ul_maxthickness
        # guarantees that the longest lipid is longer than the width of the leaflet
        zmin = 0.0
        zmax = lc_thickness + ll_actthickness + 2 * half_mid_zgap + ul_actthickness + uc_thickness
        self.LC['z-lo'] = zmin
        self.LC['z-hi'] = self.LC['z-lo'] + lc_thickness
        self.LL['z-lo'] = self.LC['z-hi']
        self.LL['z-hi'] = self.LL['z-lo'] + ll_actthickness
        self.UL['z-lo'] = self.LL['z-hi'] + 2 * half_mid_zgap
        self.midplane_z = self.LL['z-hi'] + half_mid_zgap
        self.UL['z-hi'] = self.UL['z-lo'] + ul_actthickness
        self.UC['z-lo'] = self.UL['z-hi']
        self.UC['z-hi'] = zmax
        logger.debug(f'++++++++++++++ {zmax:8.3f} ++++++++++++++')
        logger.debug(f'                                    ^  ')
        logger.debug(f'    U P P E R   C H A M B E R {zmax-self.UC["z-lo"]:8.3f}')
        logger.debug(f'                                    v  ')
        logger.debug(f'-------------- {self.UC["z-lo"]:8.3f} --------------')
        logger.debug(f'                                    ^  ')
        logger.debug(f'    U P P E R   L E A F L E T {self.UC["z-lo"]-self.UL["z-lo"]:8.3f}    ')
        logger.debug(f'                                    v  ')
        logger.debug(f'-------------- {self.UL["z-lo"]:8.3f} --------------')
        logger.debug(f'    M I D P L A N E  H A L F  G A P    ')
        logger.debug(f'-------------- {self.midplane_z:8.3f} --------------')
        logger.debug(f'    M I D P L A N E  H A L F  G A P    ')
        logger.debug(f'-------------- {self.LL["z-hi"]:8.3f} --------------')
        logger.debug(f'                                    ^  ')
        logger.debug(f'    L O W E R   L E A F L E T {self.LL["z-hi"]-self.LL["z-lo"]:8.3f}    ')
        logger.debug(f'                                    v  ')
        logger.debug(f'-------------- {self.LC["z-hi"]:8.3f} --------------')
        logger.debug(f'                                    ^  ')
        logger.debug(f'    L O W E R   C H A M B E R {self.LC["z-hi"]-zmin:8.3f}')
        logger.debug(f'                                    v  ')
        logger.debug(f'++++++++++++++ {zmin:8.3f} ++++++++++++++')
        self.patch_ll_corner = np.array([0, 0, zmin])
        self.patch_ur_corner = np.array([Lx, Ly, zmax])
        # box and origin
        self.box = np.array([[Lx, 0, 0], [0, Ly, 0], [0, 0, zmax - zmin]])
        self.origin = np.array([self.box[i][i] / 2 for i in range(3)])

    def write_packmol(self, pm: PackmolScripter, half_mid_zgap=2.0, rotation_pm=0.0, nloop=100):
        """
        Writes the packmol input for the bilayer patch to the provided Packmol object.

        Parameters
        ----------
        pm : Packmol
            The Packmol ScriptWriter object to which the bilayer patch specifications will be written.
        half_mid_zgap : float, optional
            The half mid-plane gap in Å. Default is 2.0 Å.
        rotation_pm : float, optional
            The rotation angle in degrees for the patch. Default is 0.0 degrees.
        nloop : int, optional
            The number of loops for packing the bilayer. Default is 100.
        """
        pm.addline(f'pbc {" ".join([f"{_:.3f}" for _ in self.patch_ll_corner])} {" ".join([f"{_:.3f}" for _ in self.patch_ur_corner])}')
        ll = self.patch_ll_corner
        ur = self.patch_ur_corner
        for leaflet in [self.LL, self.UL]:
            logger.debug(f'Leaflet species to pack:')
            my_logger(leaflet["composition"], logger.debug)
            for specs in leaflet['composition']:
                name = specs['name']
                logger.debug(f'Packing {name}')
                lipid_max_length = self.species_data[name].get_max_internal_length(conformerID=specs.get('conf', 0))
                lipid_headtail_length = self.species_data[name].get_head_tail_length(conformerID=specs.get('conf', 0))
                lipid_overhang = lipid_max_length - lipid_headtail_length
                leaflet_thickness = leaflet['z-hi'] - leaflet['z-lo']
                ref_atoms = self.species_data[name].get_ref_atoms()
                hs = ' '.join([f"{x['serial']}" for x in ref_atoms['heads']])
                ts = ' '.join([f"{x['serial']}" for x in ref_atoms['tails']])
                n = specs['patn']
                pm.addline(f'structure {specs["local_name"]}')
                pm.comment(f'  max int length {lipid_max_length:.3f}, head-tail length {lipid_headtail_length:.3f}, overhang {lipid_overhang:.3f}')
                pm.addline(f'number {n}', indents=1)
                # if the maximum length of the lipid is less than the desired leaflet thickness minus a margin,
                # we can pack directly into the leaflet slab using contstrained rotation to orient
                if lipid_max_length < leaflet_thickness:
                    logger.debug(f' -> lipid length {lipid_max_length:.3f} < leaflet thickness {leaflet_thickness:.3f}')
                    if leaflet is self.LL:
                        constrain_rotation = 180.0
                    elif leaflet is self.UL:
                        constrain_rotation = 0.0
                    inside_z_lo = leaflet['z-lo']
                    inside_z_hi = leaflet['z-hi']
                    pm.addline(f'inside box {ll[0]:.3f} {ll[1]:.3f} {inside_z_lo:.3f} {ur[0]:.3f} {ur[1]:.3f} {inside_z_hi:.3f}', indents=1)
                    pm.addline(f'constrain_rotation x {constrain_rotation} {rotation_pm}', indents=1)
                    pm.addline(f'constrain_rotation y {constrain_rotation} {rotation_pm}', indents=1)
                else:
                    # we need to pack by specifying some atoms above a plane and other atoms below a different
                    # plane, and the "inside box" should refer to the whole cell
                    # if the head-tail length is GREATER than the leaflet thickness, we can use the
                    # explicit leafleat boundaries as the planes
                    if lipid_headtail_length > leaflet_thickness:
                        below_plane_z = leaflet['z-lo']
                        above_plane_z = leaflet['z-hi'] - half_mid_zgap
                    else:
                        # if the head-tail length is less than the leaflet thickness, while the max internal
                        # length is still greater, then we can't use the leaflet boundaries as the planes
                        span = leaflet_thickness - lipid_headtail_length
                        below_plane_z = leaflet['z-lo'] + span / 2.0
                        above_plane_z = leaflet['z-hi'] - span / 2.0
                    pm.addline(f'inside box {ll[0]:.3f} {ll[1]:.3f} {ll[2]:.3f} {ur[0]:.3f} {ur[1]:.3f} {ur[2]:.3f}', indents=1)
                    pm.addline(f'atoms {hs}', indents=1)
                    if leaflet is self.LL: # heads are low
                        pm.addline(f'below plane 0. 0. 1. {below_plane_z:.3f}', indents=2)
                    elif leaflet is self.UL: # heads are high
                        pm.addline(f'above plane 0. 0. 1. {above_plane_z:.3f}', indents=2)
                    pm.addline('end atoms', indents=1)
                    pm.addline(f'atoms {ts}', indents=1)
                    if leaflet is self.LL: # tails are high
                        pm.addline(f'above plane 0. 0. 1. {above_plane_z:.3f}', indents=2)
                    elif leaflet is self.UL: # tails are low
                        pm.addline(f'below plane 0. 0. 1. {below_plane_z:.3f}', indents=2)
                    pm.addline('end atoms', indents=1)

                pm.addline(f'nloop {nloop}', indents=1)
                pm.addline(f'end structure')
        for chamber in [self.LC, self.UC]:
            logger.debug(f'Chamber species to pack:')
            my_logger(chamber["composition"], logger.debug)
            for specs in chamber['composition']:
                name = specs['name']
                n = specs['patn']
                pm.addline(f'structure {specs["local_name"]}')
                pm.addline(f'number {n}', indents=1)
                inside_z_lo = chamber['z-lo']
                inside_z_hi = chamber['z-hi']
                pm.addline(f'inside box {ll[0]:.3f} {ll[1]:.3f} {inside_z_lo:.3f} {ur[0]:.3f} {ur[1]:.3f} {inside_z_hi:.3f}', indents=1)
                pm.addline(f'nloop {nloop}', indents=1)
                pm.addline(f'end structure')
        pm.writefile()
