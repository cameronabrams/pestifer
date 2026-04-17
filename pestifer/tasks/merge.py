# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`MergeTask` class for combining two or more pre-built
PSF/PDB systems into a single PSF/PDB system.

Usage is described in the :ref:`config_ref tasks merge` documentation.

The task uses a two-step VMD/psfgen approach executed in a single script:

1. **Segment renaming** (VMD ``mol load`` + ``atomselect``): any systems whose
   segment names collide with earlier systems get their conflicting segments
   renamed in a temporary copy.
2. **Topology merge** (psfgen ``readpsf``): the (possibly renamed) PSF/PDB
   files are merged by loading them sequentially into psfgen, which handles
   serial-number renumbering and topology-table reconstruction automatically.

Segment-name collisions are resolved according to ``collision_strategy``:

* ``'enumerate'`` *(default)*: strip the last character of the conflicting
  name and append a digit (``0``–``9``, then fall back to two digits) until
  a unique name is found.  E.g. ``PROA`` → ``PRO0``, ``PRO1``, …
* ``'error'``: raise :class:`ValueError` immediately on the first unresolved
  collision.  Use when every rename must be explicit.

Explicit renames take priority over the automatic strategy and are specified
per-system via ``segname_map``.
"""

import logging
import os

from .psfgen import PsfgenTask
from ..core.artifacts import *
from ..psfutil.psfcontents import PSFContents

logger = logging.getLogger(__name__)


class MergeTask(PsfgenTask):
    """
    Merge two or more pre-built PSF/PDB systems into a single PSF/PDB system.

    Each system is specified by a ``psf`` file plus either a ``pdb`` or
    ``coor`` coordinate file.  Segment name collisions are detected
    automatically and resolved according to the ``collision_strategy``.

    Parameters (from YAML ``specs``)
    ---------------------------------
    systems : list of dict
        Each entry must contain ``psf`` and either ``pdb`` or ``coor``.
        An optional ``segname_map`` dict provides explicit renames for that
        system (``{old_segname: new_segname}``).
    collision_strategy : str, optional
        ``'enumerate'`` *(default)* or ``'error'``.
    """

    _yaml_header = 'merge'

    def do(self) -> int:
        self.next_basename('merge')
        systems_specs: list[dict] = self.specs.get('systems', [])
        collision_strategy: str = self.specs.get('collision_strategy', 'enumerate')

        if len(systems_specs) < 2:
            raise ValueError('merge task requires at least 2 systems')

        # Normalise each spec: convert coor → pdb if needed, validate files.
        systems = []
        for i, spec in enumerate(systems_specs):
            psf = spec.get('psf', '')
            pdb = spec.get('pdb', '')
            coor = spec.get('coor', '')
            if not psf:
                raise ValueError(f'systems[{i}]: psf file must be specified')
            if not pdb and not coor:
                raise ValueError(f'systems[{i}]: either pdb or coor must be specified')
            if coor and not pdb:
                # Convert binary NAMD coordinate file to PDB so VMD can read it
                pdb = self.coor_to_pdb(coor, psf)
            if not os.path.exists(psf):
                raise FileNotFoundError(f'systems[{i}]: psf file not found: {psf}')
            if not os.path.exists(pdb):
                raise FileNotFoundError(f'systems[{i}]: pdb file not found: {pdb}')
            systems.append({
                'psf': psf,
                'pdb': pdb,
                'segname_map': dict(spec.get('segname_map', {})),
            })

        # Resolve segment name collisions across all systems.
        resolved = self._resolve_collisions(systems, collision_strategy)

        # Collect topology stream files from all input PSF REMARKS (union).
        CC = self.resource_manager.charmmff_content
        all_streamfiles: list[str] = []
        for sys in resolved:
            for sf in self._stream_files_from_psf(sys['psf']):
                if sf not in all_streamfiles:
                    all_streamfiles.append(sf)

        # Copy stream files to CWD so psfgen can find them.
        for sf in all_streamfiles:
            CC.copy_charmmfile_local(sf)

        # Build a single VMD/psfgen script that (a) renames segments where
        # needed, then (b) merges all systems via readpsf.
        pg = self.get_scripter('psfgen')
        pg.newscript(self.basename, additional_topologies=all_streamfiles)

        # Step A: rename segments in each system that requires it.
        for i, sys in enumerate(resolved):
            rename_map = sys['effective_segname_map']
            if rename_map:
                tmp_psf = f'_merge_tmp_{i}.psf'
                tmp_pdb = f'_merge_tmp_{i}.pdb'
                pg.comment(f'Rename segments in system {i}: {rename_map}')
                pg.addline(f'mol load psf {sys["psf"]} pdb {sys["pdb"]}')
                for old_name, new_name in rename_map.items():
                    pg.addline(f'set _ms [atomselect top "segname {old_name}"]')
                    pg.addline(f'$_ms set segname {new_name}')
                    pg.addline(f'$_ms delete')
                pg.addline(f'set _ma [atomselect top all]')
                pg.addline(f'$_ma writepsf {tmp_psf}')
                pg.addline(f'$_ma writepdb {tmp_pdb}')
                pg.addline(f'$_ma delete')
                pg.addline(f'mol delete top')
                sys['_use_psf'] = tmp_psf
                sys['_use_pdb'] = tmp_pdb
            else:
                sys['_use_psf'] = sys['psf']
                sys['_use_pdb'] = sys['pdb']

        # Step B: merge with psfgen readpsf (no coord guessing or regen needed
        # because all topology comes directly from the input PSF files).
        pg.comment('Merge all systems via psfgen readpsf')
        for sys in resolved:
            pg.load_project(sys['_use_psf'], sys['_use_pdb'])

        # Carry forward patch remarks from all input PSFs, applying any
        # segment renames so the provenance stays consistent.
        seen_patch_remarks: set[str] = set()
        for sys in resolved:
            rename_map = sys['effective_segname_map']
            for remark in self._patch_remarks_from_psf(sys['psf']):
                if rename_map:
                    remark = self._rename_patch_remark(remark, rename_map)
                if remark not in seen_patch_remarks:
                    pg.addline(f'remarks {remark}')
                    seen_patch_remarks.add(remark)

        pg.writescript(self.basename, guesscoord=False, regenerate=False)
        self.result = pg.runscript(keep_tempfiles=False)
        if self.result != 0:
            return self.result

        # Register pipeline artifacts.
        self.register(self.basename, key='tcl', artifact_type=PsfgenInputScriptArtifact)
        self.register(self.basename, key='log', artifact_type=PsfgenLogFileArtifact)
        self.register(
            dict(
                pdb=PDBFileArtifact(self.basename, pytestable=True),
                psf=PSFFileArtifact(self.basename, pytestable=True),
            ),
            key='state',
            artifact_type=StateArtifacts,
        )
        self.register(
            [CharmmffStreamFileArtifact(sf) for sf in all_streamfiles],
            key='charmmff_streamfiles',
            artifact_type=CharmmffStreamFileArtifacts,
        )

        # Log the final segment layout for the user.
        final_segnames = {}
        for sys in resolved:
            for orig in sys['segnames']:
                final = sys['effective_segname_map'].get(orig, orig)
                final_segnames[final] = sys['psf']
        logger.info(
            f'merge: produced {self.basename}.psf with '
            f'{len(final_segnames)} segments: {list(final_segnames.keys())}'
        )

        return self.result

    # ------------------------------------------------------------------
    # Collision resolution
    # ------------------------------------------------------------------

    def _resolve_collisions(
        self,
        systems: list[dict],
        strategy: str,
    ) -> list[dict]:
        """Compute effective segment rename maps for all systems.

        Returns an extended copy of *systems* where each entry gains two
        new keys:

        ``'segnames'``
            List of original segment names found in this system's PSF.
        ``'effective_segname_map'``
            Dict ``{old: new}`` containing only the renames that are actually
            needed (empty if no renames are required for this system).
        """
        used: set[str] = set()
        result = []

        for sys in systems:
            psf_contents = PSFContents(sys['psf'])
            seg_names = [s.segname for s in psf_contents.segments]
            user_map: dict[str, str] = sys.get('segname_map', {})
            effective_map: dict[str, str] = {}

            for name in seg_names:
                # Apply user-supplied rename first; default is identity.
                final = user_map.get(name, name)

                if final in used:
                    if strategy == 'error':
                        raise ValueError(
                            f"Segment name collision: '{final}' (from system "
                            f"'{sys['psf']}') conflicts with an earlier system. "
                            f"Provide an explicit segname_map or use "
                            f"collision_strategy: enumerate."
                        )
                    # Auto-enumerate: shorten to 3 chars + digit suffix.
                    final = self._unique_name(final, used)
                    effective_map[name] = final
                elif final != name:
                    # User-supplied rename — record it even if not a collision
                    # so the VMD script knows to apply it.
                    effective_map[name] = final

                used.add(final)

            result.append({
                **sys,
                'segnames': seg_names,
                'effective_segname_map': effective_map,
            })

        return result

    @staticmethod
    def _unique_name(name: str, used: set[str]) -> str:
        """Return a 4-char segment name that is not in *used*.

        Strips the last character and appends incrementing digits.
        Tries single digits (0–9), then two-digit suffixes (00–99).

        Raises
        ------
        ValueError
            If no unique name can be found within the search range.
        """
        base3 = name[:3]
        for i in range(10):
            candidate = f'{base3}{i}'
            if candidate not in used:
                return candidate
        base2 = name[:2]
        for i in range(100):
            candidate = f'{base2}{i:02d}'
            if candidate not in used:
                return candidate
        raise ValueError(
            f"Cannot find a unique segment name derived from '{name}'; "
            f"provide explicit segname_map entries."
        )

    # ------------------------------------------------------------------
    # PSF REMARKS parsing
    # ------------------------------------------------------------------

    @staticmethod
    def _patch_remarks_from_psf(psf_filename: str) -> list[str]:
        """Return patch remark texts (without the leading 'REMARKS ' prefix) from
        a PSF's REMARKS section.  Each entry looks like 'patch DISU A:5 A:55'."""
        remarks = []
        with open(psf_filename, 'r') as f:
            for line in f:
                line = line.strip()
                if '!NATOM' in line:
                    break
                if line.startswith('REMARKS patch'):
                    # Normalise whitespace while preserving seg:resid tokens.
                    remarks.append(' '.join(line[len('REMARKS '):].split()))
        return remarks

    @staticmethod
    def _rename_patch_remark(remark: str, rename_map: dict) -> str:
        """Apply a segment rename map to a patch remark text.

        Handles the psfgen colon format where segment and residue are joined
        as ``seg:resid`` (e.g. ``patch DISU A:5 A:55``).
        """
        tokens = remark.split()
        # tokens[0] == 'patch', tokens[1] == patchname, rest are seg:resid pairs
        result = tokens[:2]
        for token in tokens[2:]:
            if ':' in token:
                seg, resid = token.split(':', 1)
                result.append(f'{rename_map.get(seg, seg)}:{resid}')
            else:
                result.append(rename_map.get(token, token))
        return ' '.join(result)

    @staticmethod
    def _stream_files_from_psf(psf_filename: str) -> list[str]:
        """Return the list of ``.str`` stream-file basenames recorded in a
        PSF's REMARKS section.

        Mirrors the logic in :class:`ContinuationTask`.
        """
        stream_files = []
        with open(psf_filename, 'r') as f:
            for line in f:
                line = line.strip()
                if '!NATOM' in line:
                    break
                if line.startswith('REMARKS topology'):
                    tokens = line.split()
                    name = os.path.basename(tokens[-1])
                    if name.endswith('.str') and name not in stream_files:
                        stream_files.append(name)
        return stream_files
