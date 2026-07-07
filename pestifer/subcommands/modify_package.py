# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
Subcommand to modify the pestifer package.  This subcommand is only exposed if the Pestifer installation includes the entire source tree (i.e., it was installed with pip -e from a git repository).

Modifications are organized under a required *category* positional and a *verb*:

- ``example`` -- manage the example library (``add``/``update``/``delete``/``rename``/``author``).
- ``pdb-repo`` -- contribute residue coordinates to the built-in PDB repository (``add-entry``).
- ``charmmff`` -- modify CHARMM force-field content: contribute a built-in custom residue
  (``add-residue``), or regenerate force-field-derived resources (``regenerate-segtypes``,
  ``update-atomselect-macros``).
- ``ledger`` -- inspect or reverse recorded modifications (``show``/``revert``).

Every mutating command records what it did in an append-only ledger
(``pestifer/resources/modifications.jsonl``); ``ledger show`` lists them and ``ledger revert``
reverses one by git-reverting its commit and curating the ledger.

Because every modification is meant to become a pull request, the git workflow is folded into
the command and **runs by default**: each invocation verifies your working tree is clean, makes
the change on a fresh branch (auto-named ``modpkg/<category>-<verb>-<detail>``), and commits
exactly the files it touched -- then prints the ``git push`` / ``gh pr create`` steps.  Use
``--branch NAME`` to name the branch yourself, or ``--no-branch`` to skip branching and commit
and just apply the change to the working tree.
"""
from dataclasses import dataclass

import argparse as ap
import logging
import os
import re

from . import Subcommand

from ..core import ledger
from ..core.resourcemanager import ResourceManager
from ..util import gitutil

logger = logging.getLogger(__name__)


def _add_residue(RM, args, out=print):
    """Install a custom residue; return (touched paths, summary)."""
    result = RM.add_custom_residue(args.file, segtype=args.segtype, force=args.force)
    out(f"Installed custom residue file: {result['dest']}")
    out(f"  RESI defined: {', '.join(result['resnames'])}")
    if result['patch_names']:
        out(f"  PRES defined: {', '.join(result['patch_names'])}")
    if result['segtype_added']:
        out(f"  classified {', '.join(result['segtype_added'])} as segtype '{result['segtype']}'")
    else:
        out(f"  (segtype already registered for {', '.join(result['resnames'])})")
    out('  resource cache cleared; it will rebuild on the next run.')
    return result['touched_paths'], result


def _add_pdb_entry(RM, args, out=print):
    """Install make-pdb-collection entries into the PDB repository; return (touched, summary).

    ``args.entry_dir`` may be a single entry directory or a directory of entry
    subdirectories (a whole stream/substream), so a batch install is one command.
    """
    result = RM.add_pdb_collection(args.entry_dir, collection=args.collection, force=args.force)
    entries = result['entries']
    for resname, coll, kind in entries:
        content = 'box' if kind == 'box' else 'conformers'
        out(f"  {resname} -> {coll} ({content})")
    ncoll = len(result['collections'])
    out(f"Installed {len(entries)} entr{'y' if len(entries) == 1 else 'ies'} into "
        f"{ncoll} collection{'s' if ncoll != 1 else ''}: {', '.join(result['collections'])}")
    if result['created_collections']:
        out(f"  created collection{'s' if len(result['created_collections']) != 1 else ''}: "
            f"{', '.join(result['created_collections'])}")
    out('  resource cache cleared; it will rebuild on the next run.')
    return result['touched_paths'], result


def _do_example(RM, args):
    """Dispatch an ``example`` verb onto the ResourceManager example-management methods."""
    verb = args.verb
    if verb == 'add':
        RM.add_example(args.scriptname, example_id=args.id, author_name=args.author_name,
                       author_email=args.author_email, title=args.title, db_id=args.db_id,
                       auxiliary_inputs=args.auxiliary_inputs, outputs=args.outputs)
    elif verb == 'update':
        RM.update_example(args.example_id, shortname=args.name, author_name=args.author_name,
                          author_email=args.author_email, title=args.title, db_id=args.db_id,
                          auxiliary_inputs=args.auxiliary_inputs, outputs=args.outputs)
    elif verb == 'delete':
        RM.delete_example(args.example_id)
    elif verb == 'rename':
        RM.rename_example(args.example_id, new_name=args.new_name)
    elif verb == 'author':
        RM.set_example_author(args.example_id, args.author_name, args.author_email)


def _slugify(text) -> str:
    """Lowercase, hyphenate, and trim a value for use in a branch name."""
    return re.sub(r'[^a-z0-9]+', '-', str(text).lower()).strip('-')[:40]


def _unique_branch(repo_root, base) -> str:
    """Return ``base``, or ``base-2``/``base-3``/... if a branch by that name already exists."""
    name, n = base, 2
    while gitutil.branch_exists(repo_root, name):
        name = f"{base}-{n}"
        n += 1
    return name


def _auto_branch_name(repo_root, category, verb, args) -> str:
    """Derive a unique, readable branch name for a modify-package contribution."""
    detail = ''
    if category == 'example':
        detail = _slugify(getattr(args, 'scriptname', '') or getattr(args, 'new_name', '')
                          or getattr(args, 'example_id', ''))
    elif category == 'pdb-repo':
        detail = _slugify(os.path.basename(os.path.normpath(args.entry_dir)))
    elif category == 'charmmff' and verb == 'add-residue':
        detail = _slugify(os.path.splitext(os.path.basename(args.file))[0])
    base = f"modpkg/{category}-{verb}" + (f"-{detail}" if detail else '')
    return _unique_branch(repo_root, base)


def _summary_for(category, verb, args, commit_summary) -> str:
    """One-line human description of a modification, used for the ledger and commit message."""
    if category == 'example':
        return {
            'add': f"add example from {os.path.basename(getattr(args, 'scriptname', ''))}",
            'rename': f"rename example {getattr(args, 'example_id', '')} to {getattr(args, 'new_name', '')}",
            'delete': f"delete example {getattr(args, 'example_id', '')}",
            'update': f"update example {getattr(args, 'example_id', '')}",
            'author': f"set author for example {getattr(args, 'example_id', '')}",
        }.get(verb, str(verb))
    if commit_summary is not None and 'resnames' in commit_summary:
        return (f"add {', '.join(commit_summary['resnames'])} to built-in custom "
                f"[segtype: {commit_summary['segtype']}]")
    if commit_summary is not None and 'entries' in commit_summary:
        entries, colls = commit_summary['entries'], commit_summary['collections']
        if len(entries) == 1:
            resname, coll, _kind = entries[0]
            return f"add {resname} coordinates to the {coll} collection"
        return (f"add {len(entries)} entries to {'the ' if len(colls) == 1 else ''}"
                f"{', '.join(colls)} collection{'s' if len(colls) != 1 else ''}")
    if verb == 'regenerate-segtypes':
        return 'regenerate derived segtype classification'
    if verb == 'update-atomselect-macros':
        return 'regenerate atomselect macros'
    return 'modify-package change'


def _ledger_show(RM, args):
    """Print the modification ledger (optionally filtered/limited)."""
    entries = ledger.read(RM.resources_path)
    print(ledger.format_entries(entries, limit=args.limit, category=args.category_filter))


def _ledger_revert(RM, args, repo_root):
    """
    Reverse a previously-recorded modification by git-reverting the commit that made it, then
    curate the ledger to keep the audit trail (the original entry is marked reverted and a new
    ``revert`` entry is appended).  Runs on a fresh branch by default (like the other verbs).
    """
    entry_id = args.id
    snapshot = ledger.read(RM.resources_path)          # read BEFORE reverting (has the entry)
    target = next((e for e in snapshot if e.get('id') == entry_id), None)
    if target is None:
        raise RuntimeError(f'no ledger entry #{entry_id}')
    if target.get('category') == 'ledger':
        raise RuntimeError(f'entry #{entry_id} is itself a ledger action and cannot be reverted')
    if target.get('reverted'):
        raise RuntimeError(f'ledger entry #{entry_id} is already reverted')

    ledger_path = ledger.ledger_file(RM.resources_path)
    sha = gitutil.commit_introducing(repo_root, ledger_path, f'"id": {entry_id},')
    if not sha:
        raise RuntimeError(
            f'could not find the commit that recorded entry #{entry_id}; it may have been made '
            'with --no-branch (uncommitted) or already rewritten -- revert it by hand')

    contribute = not getattr(args, 'no_branch', False)
    if contribute and not gitutil.worktree_is_clean(repo_root):
        raise RuntimeError('working tree is not clean; commit or stash your changes first, '
                           'or pass --no-branch')
    branch_name = getattr(args, 'branch', None) or _unique_branch(repo_root, f'modpkg/revert-{entry_id}')
    if contribute:
        if gitutil.branch_exists(repo_root, branch_name):
            raise gitutil.GitError(f"branch '{branch_name}' already exists")
        gitutil.create_and_checkout_branch(repo_root, branch_name)

    # reverse the recorded change into the working tree (this also removes the entry's own
    # ledger line, since that commit added it); then rewrite the ledger from the snapshot with
    # the original marked reverted and a new revert entry appended, so the audit trail survives
    gitutil.revert_into_worktree(repo_root, sha)
    new_id = max((e.get('id', 0) for e in snapshot), default=0) + 1
    for e in snapshot:
        if e.get('id') == entry_id:
            e['reverted'] = True
            e['reverted_by'] = new_id
    snapshot.append({
        'schema': ledger.SCHEMA_VERSION, 'id': new_id,
        'timestamp': ledger._now(), 'category': 'ledger', 'verb': 'revert',
        'summary': f"revert #{entry_id}: {target.get('summary', '')}",
        'files': target.get('files', []), 'author': gitutil.git_identity(repo_root),
        'branch': branch_name if contribute else None, 'reverted': False, 'reverts': entry_id,
    })
    ledger.write_all(RM.resources_path, snapshot)

    print(f"Reverted modification #{entry_id} ({target.get('summary', '')}); recorded as #{new_id}.")
    if contribute:
        gitutil.commit_all(repo_root, f"revert(modify-package): #{entry_id} {target.get('summary', '')}")
        print(f"\nCommitted the revert to new branch '{branch_name}':")
        print("\nReview it with `git show`, then push and open a pull request:")
        print(f"    git push -u origin {branch_name}")
        print("    gh pr create      # or open the PR on GitHub")


@dataclass
class ModifyPackageSubcommand(Subcommand):
    name: str = 'modify-package'
    short_help: str = "modify the pestifer package"
    long_help: str = "Modify the pestifer package: manage examples, contribute PDB-repository coordinates, or modify CHARMM force-field content."

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        """
        Modify the pestifer source package, dispatching on the ``category`` and ``verb`` positionals.

        .. note::

            The command ``pestifer modify-package`` can only be invoked if ``pestifer`` is installed as an editable source tree, i.e. the directory containing the ``pestifer`` package.  A simple pip installation of the pestifer package will not allow this command to run, as the pestifer package will not have a full source tree available.  If you want to modify the package, you must clone the repository from GitHub and then execute ``pip install -e .`` from the package root.
        """
        RM = ResourceManager()
        category = args.category
        verb = getattr(args, 'verb', None)
        repo_root = gitutil.package_repo_root()

        # ---- ledger category: read-only show, or the (self-contained) revert flow ---------
        if category == 'ledger':
            if verb == 'show':
                _ledger_show(RM, args)
            elif verb == 'revert':
                _ledger_revert(RM, args, repo_root)
            return True

        # a contribution branch is opened by default; --no-branch applies the change in place
        contribute = not getattr(args, 'no_branch', False)
        explicit_branch = getattr(args, 'branch', None)
        touched: list = []
        commit_summary = None
        branch_name = None

        if contribute:
            if not gitutil.worktree_is_clean(repo_root):
                raise RuntimeError(
                    'working tree is not clean; commit or stash your changes first, or pass '
                    '--no-branch to apply the change in place (a contribution branch must '
                    'contain only the contribution)')
            if explicit_branch is not None:
                if gitutil.branch_exists(repo_root, explicit_branch):
                    raise gitutil.GitError(
                        f"branch '{explicit_branch}' already exists; choose a different --branch name or delete it first")
                branch_name = explicit_branch
            else:
                branch_name = _auto_branch_name(repo_root, category, verb, args)

        # ---- apply the modification -------------------------------------------------------
        if category == 'example':
            _do_example(RM, args)
            # example management touches files it does not enumerate; with a clean tree
            # (contribute mode) every current change is attributable to this operation
            touched = gitutil.changed_paths(repo_root)

        elif category == 'pdb-repo':
            if verb == 'add-entry':
                paths, commit_summary = _add_pdb_entry(RM, args)
                touched += paths

        elif category == 'charmmff':
            if verb == 'add-residue':
                paths, commit_summary = _add_residue(RM, args)
                touched += paths
            elif verb == 'regenerate-segtypes':
                from ..core.labels import _DERIVED_SEGTYPES_PATH
                derived = RM.regenerate_derived_segtypes()
                n = sum(len(v) for v in derived.values())
                print(f"Regenerated force-field-derived segtype classification: {n} residues "
                      f"across {len(derived)} segtypes -> {os.path.relpath(_DERIVED_SEGTYPES_PATH, repo_root)}")
                touched.append(str(_DERIVED_SEGTYPES_PATH))
            elif verb == 'update-atomselect-macros':
                RM.update_atomselect_macros()
                touched.append(os.path.join(RM.resource_path['tcl'], 'macros.tcl'))

        summary = _summary_for(category, verb, args, commit_summary)

        # ---- record the modification in the ledger and fold it into the touched set -------
        if touched:
            rel_files = sorted(os.path.relpath(p, repo_root) for p in touched)
            entry = ledger.append(
                RM.resources_path, category=category, verb=verb, summary=summary,
                files=rel_files, author=gitutil.git_identity(repo_root), branch=branch_name)
            touched.append(str(ledger.ledger_file(RM.resources_path)))
            print(f"Recorded modification #{entry['id']} in the ledger.")

        # ---- commit on a fresh branch, or leave the change in the working tree -------------
        if contribute:
            if not touched:
                raise RuntimeError('no package files were modified; nothing to commit '
                                   '(pass --no-branch if you did not intend to open a branch)')
            gitutil.create_and_checkout_branch(repo_root, branch_name)
            gitutil.stage_and_commit(repo_root, touched, f"contrib({category}): {summary}")
            print(f"\nCommitted to new branch '{branch_name}':")
            for p in touched:
                print(f"    {os.path.relpath(p, repo_root)}")
            print("\nReview it with `git show`, then push and open a pull request:")
            print(f"    git push -u origin {branch_name}")
            print("    gh pr create      # or open the PR on GitHub")
        return True

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)

        # shared git-flow options, mixed into every verb.  By default modify-package opens a
        # fresh contribution branch and commits exactly the files it touched; --branch names
        # that branch explicitly, and --no-branch skips branching/committing entirely.
        contrib = ap.ArgumentParser(add_help=False)
        g = contrib.add_mutually_exclusive_group()
        g.add_argument('--branch', type=str, default=None, metavar='NAME',
                       help='name for the contribution branch (default: an auto-generated name); '
                            'the change is committed on a new branch off HEAD, requires a clean working tree')
        g.add_argument('--no-branch', dest='no_branch', action='store_true',
                       help='do not open a branch or commit; just apply the change to the working tree')

        cats = self.parser.add_subparsers(dest='category', required=True, metavar='CATEGORY',
                                          help='category of modification: example | pdb-repo | charmmff | ledger')

        # ---- example ----------------------------------------------------------------
        ex = cats.add_parser('example', help='manage the example library')
        exv = ex.add_subparsers(dest='verb', required=True, metavar='VERB',
                                help='add | update | delete | rename | author')
        p = exv.add_parser('add', parents=[contrib], help='add a new example from a YAML script')
        p.add_argument('scriptname', type=str, help='YAML file of the example')
        p.add_argument('--id', type=int, default=0, help='integer ID to assign (default: next available)')
        p.add_argument('--author-name', dest='author_name', type=str, default='', help='author name (default: extract from the script header)')
        p.add_argument('--author-email', dest='author_email', type=str, default='', help='author email (default: extract from the script header)')
        p.add_argument('--title', type=str, default='', help="1-line title (default: the script's 'title' directive)")
        p.add_argument('--db-id', dest='db_id', type=str, default='', help="database ID (default: from the script's fetch task)")
        p.add_argument('--auxiliary-inputs', dest='auxiliary_inputs', type=str, nargs='*', default=[], help='auxiliary input files for the example')
        p.add_argument('--outputs', type=str, nargs='*', default=[], help='output files for the example')
        p = exv.add_parser('update', parents=[contrib], help='update an existing example')
        p.add_argument('example_id', type=int, help='integer ID of the example to update')
        p.add_argument('--name', type=str, default='', help='new short name for the example')
        p.add_argument('--author-name', dest='author_name', type=str, default='')
        p.add_argument('--author-email', dest='author_email', type=str, default='')
        p.add_argument('--title', type=str, default='')
        p.add_argument('--db-id', dest='db_id', type=str, default='')
        p.add_argument('--auxiliary-inputs', dest='auxiliary_inputs', type=str, nargs='*', default=[])
        p.add_argument('--outputs', type=str, nargs='*', default=[])
        p = exv.add_parser('delete', parents=[contrib], help='delete an example by ID')
        p.add_argument('example_id', type=int, help='integer ID of the example to delete')
        p = exv.add_parser('rename', parents=[contrib], help='rename an example')
        p.add_argument('example_id', type=int, help='integer ID of the example to rename')
        p.add_argument('new_name', type=str, help='new name for the example')
        p = exv.add_parser('author', parents=[contrib], help='set an example author')
        p.add_argument('example_id', type=int, help='integer ID of the example')
        p.add_argument('author_name', type=str, help='author name')
        p.add_argument('author_email', type=str, help='author email')

        # ---- pdb-repo ---------------------------------------------------------------
        pr = cats.add_parser('pdb-repo', help='contribute residue coordinates to the PDB repository')
        prv = pr.add_subparsers(dest='verb', required=True, metavar='VERB', help='add-entry')
        p = prv.add_parser('add-entry', parents=[contrib], help='install one or more make-pdb-collection entries into the PDB repository')
        p.add_argument('entry_dir', type=str, metavar='DIR', help='a make-pdb-collection entry directory (named after the residue; contains info.yaml), OR a directory of such entries (e.g. a whole <streamID>/ tree) to install them all at once')
        p.add_argument('--collection', type=str, default=None, metavar='NAME', help='collection/stream tarball to install into (default: each residue segtype; created if absent). With a batch directory, forces all entries into this one collection')
        p.add_argument('--force', action='store_true', help='overwrite an entry already present for this resname in the collection')

        # ---- charmmff ---------------------------------------------------------------
        cf = cats.add_parser('charmmff', help='modify CHARMM force-field content')
        cfv = cf.add_subparsers(dest='verb', required=True, metavar='VERB',
                                help='add-residue | regenerate-segtypes | update-atomselect-macros')
        p = cfv.add_parser('add-residue', parents=[contrib], help='install a custom residue into the force field custom/ directory')
        p.add_argument('file', type=str, metavar='FILE', help='a .str/.rtf/.top file with at least one RESI block')
        p.add_argument('--segtype', type=str, default='ligand', help="segtype to classify the residue's RESI name(s) under (default: ligand)")
        p.add_argument('--force', action='store_true', help='overwrite an existing custom file and permit residue-name collisions with the force field')
        p = cfv.add_parser('regenerate-segtypes', parents=[contrib], help='regenerate the force-field-derived residue->segtype classification (developer use only)')
        p = cfv.add_parser('update-atomselect-macros', parents=[contrib], help='regenerate resources/tcl/macros.tcl from core/labels.py (developer use only)')

        # ---- ledger -----------------------------------------------------------------
        lg = cats.add_parser('ledger', help='inspect or revert recorded modifications')
        lgv = lg.add_subparsers(dest='verb', required=True, metavar='VERB', help='show | revert')
        p = lgv.add_parser('show', help='list recorded modifications')
        p.add_argument('--limit', type=int, default=None, help='show only the most recent N entries')
        p.add_argument('--category', dest='category_filter', type=str, default=None,
                       help='show only entries in this category (example | pdb-repo | charmmff | ledger)')
        p = lgv.add_parser('revert', parents=[contrib],
                           help='reverse a recorded modification by its ledger id (git-revert on a fresh branch)')
        p.add_argument('id', type=int, help='ledger id of the modification to revert (see `ledger show`)')

        return self.parser
