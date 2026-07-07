"""
Some functions for interacting with git
"""
import subprocess

from pathlib import Path


class GitError(RuntimeError):
    """Raised when a git operation fails or a precondition is not met."""


def package_repo_root() -> Path:
    """
    Return the root of the pestifer source repository (the directory containing
    ``.git``).  Only meaningful for an editable/source install.
    """
    return Path(__file__).resolve().parents[2]


def _git(repo_root, *args, check=True):
    """Run ``git -C <repo_root> <args>`` and return its stripped stdout."""
    result = subprocess.run(
        ["git", "-C", str(repo_root), *args],
        capture_output=True, text=True,
    )
    if check and result.returncode != 0:
        raise GitError(f'git {" ".join(args)} failed: {result.stderr.strip() or result.stdout.strip()}')
    return result.stdout.strip()


def worktree_is_clean(repo_root) -> bool:
    """True if the working tree has no staged or unstaged changes."""
    return _git(repo_root, 'status', '--porcelain') == ''


def changed_paths(repo_root) -> list:
    """
    Absolute paths of every file changed, added, or deleted in the working tree, as
    reported by ``git status --porcelain``.  Since callers require a clean tree before
    running an operation, this attributes exactly that operation's changes -- useful when
    an operation touches files it cannot easily enumerate itself (e.g. example management).
    """
    out = _git(repo_root, 'status', '--porcelain')
    paths = []
    for line in out.splitlines():
        entry = line[3:]                       # strip the two-char status + space
        if ' -> ' in entry:                    # rename: "old -> new"
            entry = entry.split(' -> ', 1)[1]
        entry = entry.strip().strip('"')       # git quotes paths with odd characters
        paths.append(str(Path(repo_root) / entry))
    return paths


def current_branch(repo_root) -> str:
    """Name of the currently checked-out branch."""
    return _git(repo_root, 'rev-parse', '--abbrev-ref', 'HEAD')


def branch_exists(repo_root, name: str) -> bool:
    """True if a local branch ``name`` already exists."""
    result = subprocess.run(
        ["git", "-C", str(repo_root), 'rev-parse', '--verify', '--quiet', f'refs/heads/{name}'],
        capture_output=True, text=True,
    )
    return result.returncode == 0


def create_and_checkout_branch(repo_root, name: str, base: str = None):
    """
    Create a new branch ``name`` (off ``base``, or off current HEAD if ``base`` is
    None) and check it out.  Raises :class:`GitError` if the branch already exists.
    """
    if branch_exists(repo_root, name):
        raise GitError(f'branch {name!r} already exists; choose a different --branch name or delete it first')
    args = ['checkout', '-b', name]
    if base:
        args.append(base)
    _git(repo_root, *args)


def stage_and_commit(repo_root, paths, message: str):
    """
    Stage exactly ``paths`` (no more) and create a commit with ``message``.
    Raises :class:`GitError` if there is nothing to commit.
    """
    rel = [str(Path(p).resolve()) for p in paths]
    _git(repo_root, 'add', '--', *rel)
    if _git(repo_root, 'diff', '--cached', '--name-only') == '':
        raise GitError('no staged changes to commit')
    _git(repo_root, 'commit', '-m', message)


def get_git_origin_url():
    """
    Get the URL of the git origin remote.  This is necessary for developers who want to modify the package by adding examples, etc.
    Written by ChatGPT

    Returns
    -------
    str
        The URL of the git origin remote, or None if it cannot be determined.
    """
    try:
        repo_root = Path(__file__).resolve().parent.parent
        result = subprocess.run(
            ["git", "-C", str(repo_root), "remote", "get-url", "origin"],
            capture_output=True,
            text=True,
            check=True
        )
        return result.stdout.strip()
    except Exception:
        return None
