# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
Subcommand to modify the pestifer package.  This subcommand is only exposed if the Pestifer installation includes the entire source tree (i.e., it was installed with pip -e from a git repository).
The main ways you can modify the package are

1. Managing examples, including adding, updating, and deleting example scripts.  Updating can include name-changes, adding new auxiliary inputs, and modifying expected outputs.
2. Contributing a new built-in *custom* residue: install a CHARMM ``RESI`` definition (from a ``.str``/``.rtf``/``.top`` file) into the force field's ``custom/`` directory and register its segtype.
3. Regenerating atomselect macros based on globals defined in :mod:`core/labels.py <pestifer.core.labels>`.

For a contribution, the ``--branch`` option folds the git workflow into the command: it verifies your working tree is clean, makes the change on a new branch, and commits exactly the files it touched.  You then push the branch and open a pull request for review.
"""
from dataclasses import dataclass

import logging
import os

import argparse as ap

from . import Subcommand

from ..core.resourcemanager import ResourceManager
from ..util import gitutil

logger = logging.getLogger(__name__)


def _add_residue(RM, args, out=print):
    """Install a custom residue; return the list of repository paths it touched."""
    result = RM.add_custom_residue(args.add_residue, segtype=args.segtype, force=args.force)
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
    """Install a make-pdb-collection entry into the PDB repository; return touched paths."""
    result = RM.add_pdb_entry(args.add_pdb_entry, collection=args.collection, force=args.force)
    verb = 'created collection' if result['created_collection'] else 'added to collection'
    out(f"Installed PDB-repo entry for {result['resname']} ({result['nconformers']} conformer"
        f"{'s' if result['nconformers'] != 1 else ''}): {verb} '{result['collection']}'")
    out(f"  tarball: {result['tarball']}")
    out('  resource cache cleared; it will rebuild on the next run.')
    return result['touched_paths'], result

@dataclass
class ModifyPackageSubcommand(Subcommand):
    name: str = 'modify-package'
    short_help: str = "modify the pestifer package"
    long_help: str = "Modify the pestifer package by adding, updating, or deleting examples."

    @staticmethod
    def func(args: ap.Namespace, **kwargs):
        """
        Modify the pestifer source package.  Current modifications allowed involve example addition/update/deletion, and
        updating atomselect macros based on globals defined in :mod:`core/labels.py <pestifer.core.labels>`.

        .. note::
        
            The command ``pestifer modify-package`` can only be invoked if ``pestifer`` is installed as an editable source tree, i.e. the directory containing the ``pestifer`` package.  A simple pip installation of the pestifer package will not allow this command to run, as the pestifer package will not have a full source tree available.  If you want to modify the package, you must clone the repository from GitHub and then execute ``pip install -e .`` from the package root.  Before modifying, it would also be a good idea to create a new branch, rather than modifying in the main branch.
        """
        # First, verify that the user has a full source tree by verifying the existence of the docs directory.  If not, raise an error.
        RM = ResourceManager()
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

        if getattr(args, 'add_residue', None):
            paths, commit_summary = _add_residue(RM, args)
            touched += paths

        if getattr(args, 'add_pdb_entry', None):
            paths, commit_summary = _add_pdb_entry(RM, args)
            touched += paths

        match args.example_action:
            case 'add':
                example_id = args.example_id
                scriptname = args.example_scriptname
                author_name = args.example_author_name
                author_email = args.example_author_email
                title = args.example_title
                db_id = args.example_db_id
                auxiliary_inputs = args.example_auxiliary_inputs
                outputs = args.example_outputs
                if scriptname:
                    RM.add_example(scriptname, example_id=example_id, author_name=author_name, author_email=author_email, title=title, db_id=db_id, auxiliary_inputs=auxiliary_inputs, outputs=outputs)
                else:
                    raise ValueError(f'Invalid parameter for add example action: scriptname={scriptname}. Must be non-empty YAML file name.')
            case 'update':
                example_id = args.example_id
                name = args.example_name
                author_name = args.example_author_name
                author_email = args.example_author_email
                title = args.example_title
                db_id = args.example_db_id
                auxiliary_inputs = args.example_auxiliary_inputs
                outputs = args.example_outputs
                if example_id > 0:
                    RM.update_example(example_id, shortname=name, author_name=author_name, author_email=author_email, title=title, db_id=db_id, auxiliary_inputs=auxiliary_inputs, outputs=outputs)
                else:
                    raise ValueError(f'Invalid parameters for update example action: example_id={example_id}. Must be positive integer.')
            case 'delete':
                example_id = args.example_id
                if example_id > 0:
                    RM.delete_example(example_id)
                else:
                    raise ValueError(f'Invalid parameter for delete example action: example_id={example_id}. Must be positive integer.')
            case 'rename':
                example_id = args.example_id
                new_name = args.example_name
                if example_id > 0 and new_name:
                    RM.rename_example(example_id, new_name=new_name)
                else:
                    raise ValueError(f'Invalid parameters for rename example action: example_id={example_id}, new_name={new_name}. Must be positive integer and non-empty YAML file name.')
            case 'author':
                example_id = args.example_id
                author_name = args.example_author_name
                author_email = args.example_author_email
                if example_id > 0 and author_name and author_email:
                    RM.set_example_author(example_id, author_name, author_email)
                else:
                    raise ValueError(f'Invalid parameters for set-author example action: example_id={example_id}, author_name={author_name}, author_email={author_email}. Must be positive integer and non-empty strings.')
            case _:
                pass
        if args.update_atomselect_macros:
            RM.update_atomselect_macros()

        if getattr(args, 'regenerate_segtypes', False):
            from ..core.labels import _DERIVED_SEGTYPES_PATH
            derived = RM.regenerate_derived_segtypes()
            n = sum(len(v) for v in derived.values())
            print(f"Regenerated force-field-derived segtype classification: {n} residues "
                  f"across {len(derived)} segtypes -> {os.path.relpath(_DERIVED_SEGTYPES_PATH, repo_root)}")
            touched.append(str(_DERIVED_SEGTYPES_PATH))

        if contribute:
            if not touched:
                raise RuntimeError('--branch was given but no package files were modified; nothing to commit')
            if commit_summary is not None and 'resnames' in commit_summary:
                names = ', '.join(commit_summary['resnames'])
                message = (f"contrib(residue): add {names} to built-in custom "
                           f"[segtype: {commit_summary['segtype']}]")
            elif commit_summary is not None and 'collection' in commit_summary:
                message = (f"contrib(pdb-repo): add {commit_summary['resname']} coordinates "
                           f"to the {commit_summary['collection']} collection")
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
        self.parser.add_argument('--add-residue', type=str, default=None, metavar='FILE', help='install FILE (a .str/.rtf/.top file with at least one RESI block) as a pestifer built-in custom residue, copying it into the force field custom/ directory and registering its segtype')
        self.parser.add_argument('--segtype', type=str, default='ligand', help="segtype to classify the added residue's RESI name(s) under (default: ligand)")
        self.parser.add_argument('--add-pdb-entry', type=str, default=None, metavar='DIR', help='install a make-pdb-collection entry directory (named after the residue, containing info.yaml + conformer PDBs) into the built-in PDB repository so make_membrane_system can place the residue')
        self.parser.add_argument('--collection', type=str, default=None, metavar='NAME', help='with --add-pdb-entry: the collection/stream tarball to install into (default: the residue segtype; created if absent)')
        self.parser.add_argument('--force', action='store_true', help='overwrite existing content: with --add-residue, an existing custom file / force-field name collision; with --add-pdb-entry, an entry already present for this resname in the collection')
        self.parser.add_argument('--branch', type=str, default=None, metavar='NAME', help='make the change on a new git branch NAME (created off the current HEAD) and commit exactly the files it touches; requires a clean working tree')
        self.parser.add_argument('--update-atomselect-macros', action='store_true', help='update the resources/tcl/macros.tcl file based on content in core/labels.py; developer use only')
        self.parser.add_argument('--regenerate-segtypes', action='store_true', help='regenerate resources/labels/derived_segtypes.json (the force-field-derived residue->segtype classification) from the installed CHARMM force field(s); developer use only')
        self.parser.add_argument('--example-id', type=int, default=0, help='integer ID of example to modify; developer use only')
        self.parser.add_argument('--example-action', type=str, default=None, choices=[None, 'add', 'update', 'delete', 'rename', 'author'], help='action to perform on the example; choices are [add|update|delete|rename|author]; developer use only')
        self.parser.add_argument('--example-scriptname', type=str, default='', help='yaml file of example; developer use only')
        self.parser.add_argument('--example-name', type=str, default='', help='new name for the example for action \'rename\'; default is basename of the script file (without extension); developer use only')
        self.parser.add_argument('--example-author-name', type=str, default='', help='name of the author; if not given, pestifer attempts to extract it from the script header')
        self.parser.add_argument('--example-author-email', type=str, default='', help='email of the author; if not given, pestifer attempts to extract it from the script header')
        self.parser.add_argument('--example-title', type=str, default='', help='descriptive 1-line title of the example (default: extract from \'title\' directive in the script)')
        self.parser.add_argument('--example-db-id', type=str, default='', help='database ID of the example; if not given, pestifer attempts to extract it from the script\'s "fetch" task')
        self.parser.add_argument('--example-auxiliary-inputs', type=str, nargs='*', default=[], help='list of auxiliary input files for the example')
        self.parser.add_argument('--example-outputs', type=str, nargs='*', default=[], help='list of output files for the example')
        return self.parser