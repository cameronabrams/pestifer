# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""Species-resolved mass-density profiles :math:`\\rho(z)` along a bilayer normal.

Given a CHARMM/NAMD PSF (per-atom mass and residue name), a matching single-frame
coordinate set (a PDB or a NAMD binary ``.coor``), and an XSC cell file, this
module computes and plots the mass density of water, lipid, protein, and ions as
a function of :math:`z`.  The bilayer midplane is placed at :math:`z=0` and the
bulk solvent is made continuous through the periodic boundary (so it does not
spuriously read zero at the box edges).  For a multicomponent bilayer the total
lipid profile can be decomposed into one profile per lipid species.

The species of an atom is inferred from its residue name: explicit water, ion,
and amino-acid residue-name sets are recognized, and everything else is treated
as lipid.  This is appropriate for the membrane-protein systems pestifer builds;
the per-species sets can be overridden if needed.

This is a single-frame analysis.  For a smoother, trajectory-averaged profile,
average the per-frame results over a production DCD (e.g. via ``catdcd``/VMD).
"""
import logging
import struct

import numpy as np

logger = logging.getLogger(__name__)

# 1 amu/A^3 expressed in g/cm^3
AMU_PER_A3_TO_G_PER_CC = 1.66053906660

#: residue names treated as water
WATER_RESNAMES = frozenset(
    {'TIP3', 'TIP3P', 'TIP4', 'TIP5', 'TP3M', 'SPC', 'SPCE', 'SWM4',
     'WAT', 'HOH', 'OH2', 'OPC'})
#: residue names treated as monatomic ions
ION_RESNAMES = frozenset(
    {'POT', 'CLA', 'SOD', 'CES', 'MG', 'CAL', 'ZN2', 'ZN', 'CL', 'NA',
     'K', 'LIT', 'RUB', 'BAR', 'CD2', 'FE2'})
#: residue names treated as protein
PROTEIN_RESNAMES = frozenset(
    {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'HSD',
     'HSE', 'HSP', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR',
     'TRP', 'TYR', 'VAL', 'ACE', 'CT3', 'NME', 'NMA'})

# fixed styling for the canonical species curves
_SPECIES_STYLE = {
    'water':   dict(color='#1f77b4', lw=1.8),
    'lipid':   dict(color='#d62728', lw=1.8),
    'protein': dict(color='#2ca02c', lw=1.8),
    'ion':     dict(color='#7f7f7f', lw=1.4),
}


def _parse_psf(path):
    """Return ``(resnames, masses)`` in PSF atom order."""
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines) and '!NATOM' not in lines[i]:
        i += 1
    if i == len(lines):
        raise ValueError(f'{path}: no !NATOM record found; not a PSF?')
    natom = int(lines[i].split()[0])
    resnames = np.empty(natom, dtype=object)
    masses = np.empty(natom, dtype=float)
    for k, line in enumerate(lines[i + 1:i + 1 + natom]):
        p = line.split()
        # XPLOR/CHARMM PSF: id segname resid resname name type charge mass ...
        resnames[k] = p[3]
        masses[k] = float(p[7])
    return resnames, masses


def _read_pdb_z(path, natom):
    z = np.empty(natom, dtype=float)
    k = 0
    with open(path) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                if k < natom:
                    z[k] = float(line[46:54])
                k += 1
    if k != natom:
        raise ValueError(f'{path}: {k} atoms but PSF has {natom}')
    return z


def _read_namdbin_z(path, natom):
    """Read z from a NAMD binary coordinate file (int32 count + N*3 float64)."""
    with open(path, 'rb') as f:
        raw = f.read()
    for endian in ('<', '>'):
        n = struct.unpack(endian + 'i', raw[:4])[0]
        if n == natom and len(raw) >= 4 + n * 24:
            xyz = np.frombuffer(raw[4:4 + n * 24], dtype=endian + 'f8').reshape(n, 3)
            return np.ascontiguousarray(xyz[:, 2])
    raise ValueError(f'{path}: NAMD binary atom count does not match PSF ({natom})')


def _read_z(path, natom):
    if path.lower().endswith(('.pdb', '.ent')):
        return _read_pdb_z(path, natom)
    return _read_namdbin_z(path, natom)


def _parse_xsc_cell(path):
    """Return ``(lateral_area, c_z)`` from the last XSC data line (orthorhombic)."""
    with open(path) as f:
        data = [ln for ln in f if ln.strip() and not ln.startswith('#')]
    if not data:
        raise ValueError(f'{path}: no data lines')
    v = [float(x) for x in data[-1].split()]
    a_x, b_y, c_z = v[1], v[5], v[9]
    return a_x * b_y, c_z


def classify_species(resnames, water=WATER_RESNAMES, ions=ION_RESNAMES,
                     protein=PROTEIN_RESNAMES):
    """Map an array of residue names to ``water``/``ion``/``protein``/``lipid``."""
    cls = np.empty(len(resnames), dtype=object)
    for i, r in enumerate(resnames):
        if r in water:
            cls[i] = 'water'
        elif r in ions:
            cls[i] = 'ion'
        elif r in protein:
            cls[i] = 'protein'
        else:
            cls[i] = 'lipid'
    return cls


class DensityProfile:
    """Compute species-resolved mass-density profiles for one membrane frame."""

    def __init__(self, psf, coor, xsc, **species_sets):
        self.resnames, self.masses = _parse_psf(psf)
        self.z = _read_z(coor, len(self.resnames))
        self.area, self.c_z = _parse_xsc_cell(xsc)
        self.cls = classify_species(self.resnames, **species_sets)
        self.lipid_resnames = sorted(set(self.resnames[self.cls == 'lipid']))

    def compute(self, dz=1.0, lipid_components=False):
        """Return ``(z_centers, profiles)`` with the midplane at ``z=0``.

        ``profiles`` is an ordered dict ``label -> rho(z)`` in g/cm^3.  The
        canonical species present (water, lipid, protein, ion) are always
        included; when ``lipid_components`` is True the individual lipid
        residue names are added after the total lipid curve.
        """
        lip = self.cls == 'lipid'
        if not lip.any():
            raise ValueError('no lipid atoms found; cannot locate a midplane')
        z0 = np.average(self.z[lip], weights=self.masses[lip])
        # wrap into the periodic cell about the midplane so bulk solvent is
        # continuous through the z-PBC instead of reading zero at the box edges
        zc = self.z - z0
        zc -= self.c_z * np.round(zc / self.c_z)

        nbins = max(1, int(round(self.c_z / dz)))
        edges = np.linspace(-self.c_z / 2, self.c_z / 2, nbins + 1)
        centers = 0.5 * (edges[:-1] + edges[1:])
        slab_vol = self.area * (self.c_z / nbins)

        def density(mask):
            hist, _ = np.histogram(zc[mask], bins=edges, weights=self.masses[mask])
            return hist / slab_vol * AMU_PER_A3_TO_G_PER_CC

        profiles = {}
        for species in ('water', 'lipid', 'protein', 'ion'):
            mask = self.cls == species
            if mask.any():
                profiles[species] = density(mask)
        if lipid_components:
            for rn in self.lipid_resnames:
                profiles[f'lipid:{rn}'] = density(self.resnames == rn)
        return centers, profiles

    def plot(self, outfile, title='', dz=1.0, lipid_components=False,
             figsize=(6.4, 4.4), dpi=150):
        """Compute and render the profile to ``outfile``; returns ``outfile``."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        centers, profiles = self.compute(dz=dz, lipid_components=lipid_components)
        components = [k for k in profiles if k.startswith('lipid:')]
        comp_colors = {}
        if components:
            cmap = plt.get_cmap('tab10' if len(components) <= 10 else 'tab20')
            comp_colors = {k: cmap(i % cmap.N) for i, k in enumerate(components)}

        fig, ax = plt.subplots(figsize=figsize)
        for label, rho in profiles.items():
            if label in _SPECIES_STYLE:
                ax.plot(centers, rho, label=label, **_SPECIES_STYLE[label])
            else:  # a lipid component
                rn = label.split(':', 1)[1]
                ax.plot(centers, rho, label=rn, lw=1.1, color=comp_colors[label])

        ax.set_xlabel(r'$z$ relative to bilayer midplane (Å)')
        ax.set_ylabel(r'mass density (g cm$^{-3}$)')
        if title:
            ax.set_title(title)
        ax.set_xlim(centers.min(), centers.max())
        ax.set_ylim(bottom=0)
        ax.grid(alpha=0.25)
        # place the legend outside the axes so it never collides with the data,
        # however many lipid components are shown (bbox_inches='tight' below keeps
        # it from being clipped)
        ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5),
                  frameon=False, fontsize='small')
        fig.tight_layout()
        fig.savefig(outfile, dpi=dpi, bbox_inches='tight')
        plt.close(fig)
        logger.info(f'wrote {outfile}')
        return outfile
