# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
Subcommand to modify the pestifer package.  This subcommand is only exposed if the Pestifer installation includes the entire source tree (i.e., it was installed with pip -e from a git repository).

Modifications are organized under a required *category* positional and a *verb*:

- ``example`` -- manage the example library (``add``/``update``/``delete``/``rename``/``author``).
- ``pdb-repo`` -- contribute residue coordinates to the built-in PDB repository (``add-entry``).
- ``charmmff`` -- modify CHARMM force-field content: contribute a built-in custom residue
  (``add-residue``), or regenerate force-field-derived resources (``regenerate-segtypes``,
  ``update-atomselect-macros``).

For a contribution, the ``--branch`` option (available on the ``pdb-repo`` and ``charmmff``
verbs) folds the git workflow into the command: it verifies your working tree is clean, makes
the change on a new branch, and commits exactly the files it touched.  You then push the branch
and open a pull request for review.
"""
from dataclasses import dataclass

import argparse as ap
import logging
import os

from . import Subcommand

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
        contribute = bool(getattr(args, 'branch', None))
        repo_root = gitutil.package_repo_root()
        touched: list = []
        commit_summary = None

        if contribute:
            if not gitutil.worktree_is_clean(repo_root):
                raise RuntimeError(
                    'working tree is not clean; commit or stash your changes before using --branch '
                    '(the contribution branch must contain only the contribution)')
            if gitutil.branch_exists(repo_root, args.branch):
                raise gitutil.GitError(
                    f"branch '{args.branch}' already exists; choose a different --branch name or delete it first")

        if category == 'example':
            _do_example(RM, args)

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

        if contribute:
            if not touched:
                raise RuntimeError('--branch was given but no package files were modified; nothing to commit')
            if commit_summary is not None and 'resnames' in commit_summary:
                names = ', '.join(commit_summary['resnames'])
                message = (f"contrib(residue): add {names} to built-in custom "
                           f"[segtype: {commit_summary['segtype']}]")
            elif commit_summary is not None and 'entries' in commit_summary:
                entries = commit_summary['entries']
                colls = commit_summary['collections']
                if len(entries) == 1:
                    resname, coll, _kind = entries[0]
                    message = f"contrib(pdb-repo): add {resname} coordinates to the {coll} collection"
                else:
                    message = (f"contrib(pdb-repo): add {len(entries)} entries to "
                               f"{'the ' if len(colls) == 1 else ''}{', '.join(colls)} "
                               f"collection{'s' if len(colls) != 1 else ''}")
            elif verb == 'regenerate-segtypes':
                message = 'contrib(charmmff): regenerate derived segtype classification'
            elif verb == 'update-atomselect-macros':
                message = 'contrib(charmmff): regenerate atomselect macros'
            else:
                message = 'contrib: modify-package change'
            gitutil.create_and_checkout_branch(repo_root, args.branch)
            gitutil.stage_and_commit(repo_root, touched, message)
            print(f"\nCommitted to new branch '{args.branch}':")
            for p in touched:
                print(f"    {os.path.relpath(p, repo_root)}")
            print("\nReview it with `git show`, then push and open a pull request:")
            print(f"    git push -u origin {args.branch}")
            print("    gh pr create      # or open the PR on GitHub")
        return True

    def add_subparser(self, subparsers):
        super().add_subparser(subparsers)

        # shared contribute flag, mixed into the pdb-repo/charmmff verbs
        contrib = ap.ArgumentParser(add_help=False)
        contrib.add_argument('--branch', type=str, default=None, metavar='NAME',
                             help='make the change on a new git branch NAME (off the current HEAD) and commit exactly the files it touches; requires a clean working tree')

        cats = self.parser.add_subparsers(dest='category', required=True, metavar='CATEGORY',
                                          help='category of modification: example | pdb-repo | charmmff')

        # ---- example ----------------------------------------------------------------
        ex = cats.add_parser('example', help='manage the example library')
        exv = ex.add_subparsers(dest='verb', required=True, metavar='VERB',
                                help='add | update | delete | rename | author')
        p = exv.add_parser('add', help='add a new example from a YAML script')
        p.add_argument('scriptname', type=str, help='YAML file of the example')
        p.add_argument('--id', type=int, default=0, help='integer ID to assign (default: next available)')
        p.add_argument('--author-name', dest='author_name', type=str, default='', help='author name (default: extract from the script header)')
        p.add_argument('--author-email', dest='author_email', type=str, default='', help='author email (default: extract from the script header)')
        p.add_argument('--title', type=str, default='', help="1-line title (default: the script's 'title' directive)")
        p.add_argument('--db-id', dest='db_id', type=str, default='', help="database ID (default: from the script's fetch task)")
        p.add_argument('--auxiliary-inputs', dest='auxiliary_inputs', type=str, nargs='*', default=[], help='auxiliary input files for the example')
        p.add_argument('--outputs', type=str, nargs='*', default=[], help='output files for the example')
        p = exv.add_parser('update', help='update an existing example')
        p.add_argument('example_id', type=int, help='integer ID of the example to update')
        p.add_argument('--name', type=str, default='', help='new short name for the example')
        p.add_argument('--author-name', dest='author_name', type=str, default='')
        p.add_argument('--author-email', dest='author_email', type=str, default='')
        p.add_argument('--title', type=str, default='')
        p.add_argument('--db-id', dest='db_id', type=str, default='')
        p.add_argument('--auxiliary-inputs', dest='auxiliary_inputs', type=str, nargs='*', default=[])
        p.add_argument('--outputs', type=str, nargs='*', default=[])
        p = exv.add_parser('delete', help='delete an example by ID')
        p.add_argument('example_id', type=int, help='integer ID of the example to delete')
        p = exv.add_parser('rename', help='rename an example')
        p.add_argument('example_id', type=int, help='integer ID of the example to rename')
        p.add_argument('new_name', type=str, help='new name for the example')
        p = exv.add_parser('author', help='set an example author')
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

        return self.parser
