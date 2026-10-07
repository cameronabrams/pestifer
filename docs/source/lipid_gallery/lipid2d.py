"""Draw a lipid's 2D structure in the conventional orientation: headgroup up, tails hanging down.

Open Babel lays out from the 3D conformer and offers no control over orientation, so a POPC came
out lying on its side.  Here obabel is used only for what it is good at -- perceiving bond orders
from 3D geometry (in POPC it finds the oleoyl C=C, both ester C=O and the P=O) -- and RDKit does
the depiction, which can then be rotated deliberately.

The orientation rule is derived from the molecule rather than tuned per species, so one routine
serves every non-sterol lipid in the collection.  Sterols are out of scope: they have no acyl
chains and no head/tail axis, and the convention does not apply to them.

Usage:  python lipid2d.py <conformer.pdb> <out.svg>     (also writes <out.png>)
"""
import os
import re
import subprocess
import sys
import tempfile

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

MIN_CHAIN_CARBONS = 8
_IDX = 'pdb_row'

# CHARMM names acyl carbons two different ways and both must be read:
#   C21/C318/CA16  -- glycerophospholipids and cardiolipins: chain tag, then position
#   C10S/C16F      -- sphingolipids: position, THEN chain tag (S = sphingoid base, F = N-acyl)
# A single pattern catches only the first, which is why sphingomyelin and the ceramides came back
# as "no acyl chains found".
_CHAIN_PATTERNS = (re.compile(r'C(?P<tag>[0-9A-Z])(?P<pos>[0-9]+)$'),
                   re.compile(r'C(?P<pos>[0-9]+)(?P<tag>[SF])$'))


def pdb_atom_names(pdb):
    return [ln[12:16].strip() for ln in open(pdb) if ln.startswith(('ATOM', 'HETATM'))]


def _element_of(name):
    """CHARMM atom name -> element symbol: the first alphabetic character, and only that.

    Deliberately not a longest-match against the periodic table.  A cardiolipin's arm carbons are
    named CA1..CA18, and two-letter matching reads `CA16` as CALCIUM -- which is exactly what
    happened when the element was left for Open Babel to guess, giving a structure drawn with
    `CaH3` groups hanging off it.  Lipid topologies use only C, H, O, N, P and S, so the first
    letter is unambiguous here.
    """
    for ch in name:
        if ch.isalpha():
            return ch.upper()
    return 'C'


def _with_element_column(pdb):
    """A copy of the conformer PDB with columns 77-78 filled in.

    pestifer's conformer PDBs stop at the occupancy/B columns, leaving the element field blank,
    so any reader has to guess the element from the atom NAME -- and guessing is what produced
    the calcium.  Writing the column makes the element explicit before bond perception runs,
    rather than correcting it afterwards when it has already influenced the perceived bonds.
    """
    fh = tempfile.NamedTemporaryFile('w', suffix='.pdb', delete=False)
    for ln in open(pdb):
        if ln.startswith(('ATOM', 'HETATM')):
            row = ln.rstrip('\n').ljust(78)
            ln = row[:76] + f'{_element_of(row[12:16].strip()):>2s}' + '\n'
        fh.write(ln)
    fh.close()
    return fh.name


def acyl_chains(names_by_index):
    """``[(pivot, tip, [moving])]`` per acyl chain, in the index space of ``names_by_index``.

    ``pivot`` is the chain's first carbon -- the ester or amide carbonyl -- and the moving set is
    every carbon beyond it.  Rotating that set about the pivot swings the tail rigidly: bond
    lengths and angles within the chain are untouched, and so is the carbonyl oxygen, which hangs
    off the pivot rather than moving with it.

    A short run is not a chain: a glycerophosphocholine's C11-C15 parse as a five-member "chain 1"
    whose tip sits up in the headgroup, so the length floor is what keeps the headgroup out.
    """
    groups = {}
    for idx, nm in names_by_index.items():
        for pat in _CHAIN_PATTERNS:
            m = pat.fullmatch(nm or '')
            if m:
                groups.setdefault(m.group('tag'), []).append((int(m.group('pos')), idx))
                break
    out = []
    for members in groups.values():
        if len(members) >= MIN_CHAIN_CARBONS:
            members.sort()
            out.append((members[0][1], members[-1][1], [i for _, i in members[1:]]))
    return out


def head_indices(mol):
    """Headgroup reference atoms, chosen by ELEMENT so no name list has to be maintained.

    Phosphorus first -- every glycerophospholipid and sphingomyelin has it, a cardiolipin has two
    -- then nitrogen, which carries the ceramides and other N-headed lipids, then oxygen for a
    neutral glycerolipid.  Taken from what the molecule has rather than from its resname, so a
    species nobody enumerated still stands up the right way.
    """
    phosphorus = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'P']
    if phosphorus:
        return phosphorus
    # No phosphate: the polar region IS the head, so take nitrogen and oxygen together.  Using
    # nitrogen alone put a ceramide's reference on its mid-backbone amide and left the molecule
    # leaning; the hydroxyls are as much of its headgroup as the amide is.
    polar = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in ('N', 'O')]
    return polar


def _rotate(xy, theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.atleast_2d(xy) @ np.array([[c, s], [-s, c]])


def _assign_formal_charges(mol):
    """Charge the atoms whose valence demands it, from valence alone rather than a species list.

    A mol block perceived from 3D carries no formal charges, so RDKit refuses a choline nitrogen
    (four bonds, no charge) and an ammonium one (one bond plus three explicit hydrogens).  Both
    are caught by counting TOTAL valence rather than heavy-atom degree -- an earlier version
    checked degree alone and so handled phosphatidylcholine while rejecting every
    phosphatidylethanolamine and phosphatidylserine.
    """
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N' and atom.GetDegree() + atom.GetNumExplicitHs() == 4:
            atom.SetFormalCharge(1)
        elif atom.GetSymbol() == 'O' and atom.GetDegree() == 1 and not atom.GetNumExplicitHs():
            nb = atom.GetNeighbors()[0]
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
            if nb.GetSymbol() == 'P' and bond.GetBondTypeAsDouble() == 1.0:
                atom.SetFormalCharge(-1)
                atom.SetNoImplicit(True)


def _read_molecule(pdb):
    """Parse the conformer, strip hydrogens, and keep each surviving atom's PDB row.

    obabel is run WITHOUT its hydrogen-delete flag so that mol-block atom order matches the PDB
    exactly; the hydrogens come off afterwards in RDKit.  Deleting them in obabel instead left
    polar hydrogens behind on some lipids -- a cardiolipin arrived with 102 atoms against 100
    heavy ones -- which silently shifted a positional name mapping.  Tagging every atom with its
    PDB row before the strip removes the need for the counts to agree at all.
    """
    tmp = _with_element_column(pdb)
    try:
        block = subprocess.run(['obabel', tmp, '-omol'], capture_output=True, text=True,
                               check=True).stdout
    finally:
        os.unlink(tmp)
    mol = Chem.MolFromMolBlock(block, sanitize=False, removeHs=False)
    if mol is None:
        raise RuntimeError(f'could not read a molecule from {pdb}')
    names = pdb_atom_names(pdb)
    if mol.GetNumAtoms() != len(names):
        raise RuntimeError(f'{pdb}: {len(names)} PDB atoms but {mol.GetNumAtoms()} in the mol '
                           f'block; the name mapping would be wrong')
    for atom in mol.GetAtoms():
        atom.SetIntProp(_IDX, atom.GetIdx())
    _assign_formal_charges(mol)
    Chem.SanitizeMol(mol)
    mol = Chem.RemoveAllHs(mol)
    return mol, {a.GetIdx(): names[a.GetIntProp(_IDX)] for a in mol.GetAtoms()}


def draw(pdb, out_svg, width=420, height=560, fan_degrees=12.0):
    mol, names_by_index = _read_molecule(pdb)
    chains, heads = acyl_chains(names_by_index), head_indices(mol)
    if not chains:
        raise RuntimeError(f'{pdb}: no acyl chains found -- a sterol, or an unknown naming scheme')
    if not heads:
        raise RuntimeError(f'{pdb}: no headgroup atom found')

    rdDepictor.Compute2DCoords(mol)
    conf = mol.GetConformer()
    xy = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y]
                   for i in range(mol.GetNumAtoms())])

    # 1. Stand it up: headgroup above the chain attachments.
    #
    # NOT on the head->tips vector, which was the first attempt and does not work: RDKit's
    # depiction sends the two tails in opposite directions, so their centroid lands almost on the
    # headgroup, the vector is near zero and its direction is noise.  The attachments are a stable
    # reference because they sit on the backbone.
    pivots = [p for p, _, _ in chains]
    up = xy[heads].mean(axis=0) - xy[pivots].mean(axis=0)
    # rotate `up` onto +y.  Written as "quarter turn minus the current angle" rather than packed
    # into one arctan2 call: the compact form had its sign inverted, which stood the headgroup
    # DOWN, and that is invisible in the expression but obvious the moment it is drawn.
    xy = _rotate(xy, np.pi / 2 - np.arctan2(up[1], up[0]))

    # 2. Swing each tail rigidly about its own carbonyl until it hangs downward, fanned slightly
    # so that two chains leaving nearby attachments do not draw on top of one another.
    fan = np.deg2rad(fan_degrees)
    order = np.argsort([xy[p][0] for p, _, _ in chains])
    for rank, ci in enumerate(order):
        pivot, tip, moving = chains[ci]
        lean = -fan + 2 * fan * (rank / max(len(chains) - 1, 1))
        want = -np.pi / 2 + lean
        have = np.arctan2(*(xy[tip] - xy[pivot])[::-1])
        xy[moving] = _rotate(xy[moving] - xy[pivot], want - have) + xy[pivot]

    for i, (x, y) in enumerate(xy):
        conf.SetAtomPosition(i, (float(x), float(y), 0.0))

    for drawer, path in ((rdMolDraw2D.MolDraw2DSVG(width, height), out_svg),
                         (rdMolDraw2D.MolDraw2DCairo(width, height),
                          out_svg[:-4] + '.png' if out_svg.endswith('.svg') else None)):
        if path is None:
            continue
        drawer.drawOptions().addStereoAnnotation = False
        drawer.drawOptions().explicitMethyl = False
        rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
        drawer.FinishDrawing()
        txt = drawer.GetDrawingText()
        with open(path, 'wb' if isinstance(txt, bytes) else 'w') as fh:
            fh.write(txt)
    return mol.GetNumAtoms(), len(chains), len(heads)


if __name__ == '__main__':
    n, c, h = draw(sys.argv[1], sys.argv[2])
    print(f'{sys.argv[1]}: {n} heavy atoms, {c} chain(s), {h} head atom(s) -> {sys.argv[2]}')
