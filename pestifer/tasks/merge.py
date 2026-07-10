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

* ``'enumerate'`` *(default)*: rename the conflicting segment to a fresh
  single-character segid whose leading character is not used by any other
  segment.  Because a PSF has no chain column, VMD re-derives each atom's
  chain ID from its segid's leading character every time a PDB is regenerated
  (solvate, ``coor``→``pdb``); a single-character segid therefore keeps the
  segment's chain ID unique *and* stable through the whole pipeline, so merged
  copies of the same structure do not collapse onto one chain in chain-based
  views.
* ``'error'``: raise :class:`ValueError` immediately on the first unresolved
  collision.  Use when every rename must be explicit.

Explicit renames take priority over the automatic strategy and are specified
per-system via ``segname_map``.
"""

import logging
import os

from .psfgen import PsfgenTask
from ..core.artifacts import *
from ..core.errors import PestiferBuildError
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
    pipeline_system : dict, optional
        Rename options applied to the pipeline state when it is automatically
        prepended.  Accepted keys: ``segname_map`` and ``chainID_map``.
        Ignored when there is no prior pipeline state.
    systems : list of dict
        Each entry must contain ``psf`` and either ``pdb`` or ``coor``.
        An optional ``segname_map`` dict provides explicit renames for that
        system's PSF segment names (``{old_segname: new_segname}``).
        An optional ``chainID_map`` dict renames PDB chain IDs for that
        system (``{old_chainID: new_chainID}``).
    collision_strategy : str, optional
        ``'enumerate'`` *(default)* or ``'error'``.  Under ``'enumerate'`` a colliding segment
        is renamed to a fresh single-character segid whose leading character is unused, so its
        VMD-derived chain ID stays unique through the PSF->PDB regenerations that follow (a PSF
        carries no chain column, so the chain is re-derived from the segid each time).
    """

    _yaml_header = 'merge'

    @classmethod
    def pipeline_contract(cls, specs):
        from .pipeline_contract import TaskContract, STATE
        # combines pre-built systems from its `systems` list, folding the current state in as
        # one of them if present -- so it can start a pipeline and does not discard prior work
        return TaskContract(requires=(), provides=(STATE,))

    def do(self) -> int:
        self.next_basename('merge')
        systems_specs: list[dict] = list(self.specs.get('systems', []))
        collision_strategy: str = self.specs.get('collision_strategy', 'enumerate')

        state: StateArtifacts = self.get_current_artifact('state')
        if state and state.psf and (state.pdb or state.coor):
            pipeline_overrides: dict = self.specs.get('pipeline_system', {})
            systems_specs.insert(0, {
                'psf': state.psf.name,
                'pdb': state.pdb.name if state.pdb else '',
                'coor': state.coor.name if state.coor else '',
                'segname_map': dict(pipeline_overrides.get('segname_map', {})),
                'chainID_map': dict(pipeline_overrides.get('chainID_map', {})),
            })

        if len(systems_specs) < 2:
            raise PestiferBuildError('merge task requires at least 2 systems')

        # Normalise each spec: convert coor → pdb if needed, validate files.
        systems = []
        for i, spec in enumerate(systems_specs):
            psf = spec.get('psf', '')
            pdb = spec.get('pdb', '')
            coor = spec.get('coor', '')
            if not psf:
                raise PestiferBuildError(f'systems[{i}]: psf file must be specified')
            if not pdb and not coor:
                raise PestiferBuildError(f'systems[{i}]: either pdb or coor must be specified')
            if coor and not pdb:
                # Convert binary NAMD coordinate file to PDB so VMD can read it
                pdb = self.coor_to_pdb(coor, psf)
            if not os.path.exists(psf):
                raise PestiferBuildError(f'systems[{i}]: psf file not found: {psf}')
            if not os.path.exists(pdb):
                raise PestiferBuildError(f'systems[{i}]: pdb file not found: {pdb}')
            systems.append({
                'psf': psf,
                'pdb': pdb,
                'segname_map': dict(spec.get('segname_map', {})),
                'chainID_map': dict(spec.get('chainID_map', {})),
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

        # Copy stream files to CWD so psfgen can find them, and drop any that cannot be
        # resolved.  A PSF built outside pestifer can record a topology file pestifer does
        # not ship -- e.g. VMD's NAMD-specific 'toppar_water_ions_namd.str' (the standard
        # 'toppar_water_ions.str', which is shipped and already in this list, is its
        # equivalent).  copy_charmmfile_local only warns on a miss and writes nothing, so
        # passing the name on anyway would emit a 'topology <missing>' line that aborts
        # psfgen with "Unable to open topology file".  The merge uses readpsf, so the full
        # structure comes from the input PSFs and a missing template is safe to skip.
        available_streamfiles: list[str] = []
        dropped_streamfiles: list[str] = []
        for sf in all_streamfiles:
            CC.copy_charmmfile_local(sf)
            if os.path.exists(sf):
                available_streamfiles.append(sf)
            else:
                dropped_streamfiles.append(sf)
                logger.warning(
                    f'merge: topology stream file {sf!r} recorded in an input PSF is not '
                    f'available in pestifer\'s force field and was not found locally; '
                    f'skipping it (the merged structure is read from the input PSFs).')
        all_streamfiles = available_streamfiles

        # Build a single VMD/psfgen script that (a) renames segments where
        # needed, then (b) merges all systems via readpsf.
        pg = self.get_scripter('psfgen')
        pg.newscript(self.basename, additional_topologies=all_streamfiles)

        # Step A: rename segments in each system that requires it.
        tmp_psfs = []
        tmp_pdbs = []
        for i, sys in enumerate(resolved):
            rename_map = sys['effective_segname_map']
            chain_map = sys.get('chainID_map', {})
            if rename_map or chain_map:
                tmp_psf = f'_merge_tmp_{i}.psf'
                tmp_pdb = f'_merge_tmp_{i}.pdb'
                pg.comment(f'Rename segments/chains in system {i}: segnames={rename_map} chains={chain_map}')
                pg.addline(f'mol load psf {sys["psf"]} pdb {sys["pdb"]}')
                for old_name, new_name in rename_map.items():
                    pg.addline(f'set _ms [atomselect top "segname {old_name}"]')
                    pg.addline(f'$_ms set segname {new_name}')
                    # new_name is a single character, so it is also a valid chainID; set it here
                    # so the merged PDB carries the distinct chain immediately (readpsf/coordpdb
                    # preserves input chains, so without this the renamed segment would keep its
                    # old chain until a later PSF round-trip re-derives it from the segid).
                    pg.addline(f'$_ms set chain {new_name}')
                    pg.addline(f'$_ms delete')
                for old_chain, new_chain in chain_map.items():
                    pg.addline(f'set _mc [atomselect top "chain {old_chain}"]')
                    pg.addline(f'$_mc set chain {new_chain}')
                    pg.addline(f'$_mc delete')
                pg.addline(f'set _ma [atomselect top all]')
                pg.addline(f'$_ma writepsf {tmp_psf}')
                pg.addline(f'$_ma writepdb {tmp_pdb}')
                pg.addline(f'$_ma delete')
                pg.addline(f'mol delete top')
                sys['_use_psf'] = tmp_psf
                sys['_use_pdb'] = tmp_pdb
                tmp_psfs.append(tmp_psf)
                tmp_pdbs.append(tmp_pdb)
            else:
                sys['_use_psf'] = sys['psf']
                sys['_use_pdb'] = sys['pdb']

        # Step B: merge with psfgen readpsf (no coord guessing or regen needed
        # because all topology comes directly from the input PSF files).
        pg.comment('Merge all systems via psfgen readpsf')
        for sys in resolved:
            pg.load_project(sys['_use_psf'], sys['_use_pdb'])

        pg.writescript(self.basename, guesscoord=False, regenerate=False)
        self.result = pg.runscript(keep_tempfiles=False)
        if self.result != 0:
            return self.result

        # psfgen's readpsf carries remarks from input PSFs, but only for
        # unmodified segment names.  Inject any missing renamed patch remarks
        # directly into the merged PSF now that it has been written.
        expected_remarks: list[str] = []
        seen_remarks: set[str] = set()
        for sys in resolved:
            rename_map = sys['effective_segname_map']
            for remark in self._patch_remarks_from_psf(sys['psf']):
                if rename_map:
                    remark = self._rename_patch_remark(remark, rename_map)
                if remark not in seen_remarks:
                    expected_remarks.append(remark)
                    seen_remarks.add(remark)
        if expected_remarks:
            self._inject_patch_remarks(f'{self.basename}.psf', expected_remarks)

        # readpsf carries the input PSFs' topology remarks verbatim, including any
        # unavailable stream file we dropped above; strip those so the merged PSF does not
        # advertise a topology file pestifer cannot open (which would re-trigger the same
        # failure in a downstream task that reads this PSF).
        if dropped_streamfiles:
            self._strip_topology_remarks(f'{self.basename}.psf', dropped_streamfiles)

        # Register temporary segment-rename intermediates (written by VMD in Step A).
        if tmp_psfs:
            self.register([PSFFileArtifact(p) for p in tmp_psfs], key='merge_tmp_psfs', artifact_type=PSFFileArtifactList)
        if tmp_pdbs:
            self.register([PDBFileArtifact(p) for p in tmp_pdbs], key='merge_tmp_pdbs', artifact_type=PDBFileArtifactList)

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
            [CharmmffTopFileArtifact(x) for x in pg.topologies if x.endswith('.rtf')],
            key='charmmff_topfiles',
            artifact_type=CharmmffTopFileArtifacts,
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
        used_segids: set[str] = set()
        # Leading character of every original segid across all systems.  A renamed segment must
        # get a leading character OUTSIDE this set (and outside those already assigned), because
        # a PSF carries no chainID -- every time the pipeline regenerates a PDB from the PSF
        # (solvate, coor->pdb) VMD re-derives the chainID from the segid's leading character
        # (segid 'H0' -> chain 'H').  If two merged copies keep segids that share a leading
        # character (H, H0), they collapse onto one chainID downstream and chain-based views
        # show one copy instead of several.  Giving each renamed segment a fresh single-character
        # segid makes its derived chainID unique and survives every regeneration.
        reserved_leads: set[str] = set()
        seg_name_lists: list[list[str]] = []
        for sys in systems:
            psf_contents = PSFContents(sys['psf'])
            seg_names = [s.segname for s in psf_contents.segments]
            seg_name_lists.append(seg_names)
            user_map: dict[str, str] = sys.get('segname_map', {})
            for name in seg_names:
                reserved_leads.add((user_map.get(name, name) or ' ')[0])

        assigned_leads: set[str] = set()
        result = []

        for sys, seg_names in zip(systems, seg_name_lists):
            user_map = sys.get('segname_map', {})
            effective_map: dict[str, str] = {}

            for name in seg_names:
                # Apply user-supplied rename first; default is identity.
                final = user_map.get(name, name)

                if final in used_segids:
                    if strategy == 'error':
                        raise PestiferBuildError(
                            f"Segment name collision: '{final}' (from system "
                            f"'{sys['psf']}') conflicts with an earlier system. "
                            f"Provide an explicit segname_map or use "
                            f"collision_strategy: enumerate."
                        )
                    # Auto-resolve: fresh single-character segid whose leading character (and
                    # thus its VMD-derived chainID) is not used by any other segment.
                    final = self._fresh_segid(used_segids, reserved_leads | assigned_leads)
                    assigned_leads.add(final[0])
                    effective_map[name] = final
                elif final != name:
                    # User-supplied rename — record it even if not a collision
                    # so the VMD script knows to apply it.
                    effective_map[name] = final

                used_segids.add(final)

            result.append({
                **sys,
                'segnames': seg_names,
                'effective_segname_map': effective_map,
            })

        return result

    @staticmethod
    def _fresh_segid(used_segids: set[str], forbidden_leads: set[str]) -> str:
        """Return a fresh single-character segid that is not already a segid and whose character
        is not among *forbidden_leads*.

        A single character guarantees the VMD-derived chainID (the segid's leading character)
        equals the segid itself, so the two stay consistent through every PSF->PDB regeneration.
        Since a PDB chainID is a single column, at most 62 (``A-Za-z0-9``) distinct chains are
        possible; exhausting the pool is a hard error.
        """
        import string
        for c in string.ascii_uppercase + string.ascii_lowercase + string.digits:
            if c not in forbidden_leads and c not in used_segids:
                return c
        raise PestiferBuildError(
            "merge: ran out of unique single-character chain IDs while resolving segment "
            "collisions (the merged system has more than ~62 distinct chains). Provide explicit "
            "segname_map entries to disambiguate."
        )

    # ------------------------------------------------------------------
    # PSF REMARKS parsing
    # ------------------------------------------------------------------

    @staticmethod
    def _inject_patch_remarks(psf_path: str, expected_remarks: list[str]) -> None:
        """Ensure every entry in *expected_remarks* appears in the REMARKS header
        of *psf_path*.

        PSFContents parses the TITLE token window as psflines[fl : l1-1], where
        l1 is the !NATOM line index, so the line at l1-1 is excluded.  PSF
        convention is that this excluded line is the single blank separator before
        each section.  We therefore insert new remarks after the last existing
        REMARKS line, then reconstruct exactly one blank line before !NATOM.
        """
        with open(psf_path) as f:
            lines = f.readlines()
        existing: set[str] = set()
        natom_idx: int | None = None
        for i, line in enumerate(lines):
            stripped = line.strip()
            if '!NATOM' in stripped:
                natom_idx = i
                break
            if stripped.startswith('REMARKS patch'):
                existing.add(' '.join(stripped[len('REMARKS '):].split()))
        if natom_idx is None:
            return
        to_add = [r for r in expected_remarks if r not in existing]
        if not to_add:
            return
        # Walk back from !NATOM to find the last non-blank line (the insertion point).
        insert_after = natom_idx - 1
        while insert_after > 0 and not lines[insert_after].strip():
            insert_after -= 1
        new_lines = lines[:insert_after + 1]
        for r in to_add:
            new_lines.append(f' REMARKS {r}\n')
        new_lines.append('\n')          # exactly one blank separator before !NATOM
        new_lines.extend(lines[natom_idx:])  # !NATOM and everything after
        MergeTask._fix_ntitle_count(new_lines)
        with open(psf_path, 'w') as f:
            f.writelines(new_lines)

    @staticmethod
    def _fix_ntitle_count(lines: list[str]) -> None:
        """Rewrite the ``!NTITLE`` count line in-place so it equals the actual number of title
        (REMARKS) lines between it and the blank separator before ``!NATOM``.

        Adding or removing header REMARKS (as :meth:`_inject_patch_remarks` and
        :meth:`_strip_topology_remarks` do) without correcting this count yields a PSF that NAMD
        rejects with ``DIDN'T FIND "NATOM"`` -- it reads the declared number of title lines and,
        if that is wrong, never lands on the ``!NATOM`` record.
        """
        ntitle_idx = next((i for i, l in enumerate(lines) if l.strip().endswith('!NTITLE')), None)
        natom_idx = next((i for i, l in enumerate(lines) if '!NATOM' in l), None)
        if ntitle_idx is None or natom_idx is None:
            return
        end = natom_idx
        while end > ntitle_idx + 1 and not lines[end - 1].strip():
            end -= 1
        n_title = end - (ntitle_idx + 1)
        lines[ntitle_idx] = f'{n_title:>8} !NTITLE\n'

    @staticmethod
    def _strip_topology_remarks(psf_path: str, dropped_basenames: list[str]) -> None:
        """Remove ``REMARKS topology <file>`` header lines whose file basename is in
        *dropped_basenames* from the PSF, so the merged structure does not advertise a
        topology stream file that pestifer cannot resolve."""
        dropped = set(dropped_basenames)
        with open(psf_path) as f:
            lines = f.readlines()
        out, natom_seen, removed = [], False, False
        for line in lines:
            if not natom_seen and '!NATOM' in line:
                natom_seen = True
            if (not natom_seen and 'REMARKS topology' in line
                    and line.split() and os.path.basename(line.split()[-1]) in dropped):
                removed = True
                continue
            out.append(line)
        if removed:
            MergeTask._fix_ntitle_count(out)
            with open(psf_path, 'w') as f:
                f.writelines(out)

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
