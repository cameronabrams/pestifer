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
    # NB: parse the raw stdout, not _git()'s stripped output -- a leading-space status like
    # " M path" on the first line would otherwise lose its space and shift the path by one char
    result = subprocess.run(
        ["git", "-C", str(repo_root), 'status', '--porcelain'],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise GitError(f'git status --porcelain failed: {result.stderr.strip() or result.stdout.strip()}')
    paths = []
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        entry = line[3:]                       # XY status (2 chars) + a space, then the path
        if ' -> ' in entry:                    # rename: "old -> new"
            entry = entry.split(' -> ', 1)[1]
        entry = entry.strip().strip('"')       # git quotes paths with odd characters
        if entry:
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


def git_identity(repo_root) -> str:
    """Return the configured committer as ``"Name <email>"`` (best effort; may be empty)."""
    name = _git(repo_root, 'config', 'user.name', check=False)
    email = _git(repo_root, 'config', 'user.email', check=False)
    if name and email:
        return f'{name} <{email}>'
    return name or email or ''


def commit_introducing(repo_root, path, needle: str):
    """
    Return the SHA of the most recent commit whose change to ``path`` added or removed a line
    containing ``needle`` (git pickaxe ``-S``), or ``None`` if none is found.  Used to locate
    the commit that recorded a given ledger entry so it can be reverted.
    """
    sha = _git(repo_root, 'log', '-1', '--format=%H', f'-S{needle}', '--', str(path), check=False)
    return sha or None


def revert_into_worktree(repo_root, sha: str):
    """
    Reverse commit ``sha`` into the working tree and index **without committing**
    (``git revert -n``), so the caller can adjust the result (e.g. curate the ledger) before
    committing.  On conflict, restores the pre-revert clean state and raises :class:`GitError`.
    Callers must ensure a clean working tree first.
    """
    result = subprocess.run(
        ["git", "-C", str(repo_root), 'revert', '-n', '--no-edit', sha],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        # clean up the partial revert (the tree was clean before this call)
        subprocess.run(["git", "-C", str(repo_root), 'revert', '--quit'], capture_output=True, text=True)
        subprocess.run(["git", "-C", str(repo_root), 'reset', '--hard', 'HEAD'], capture_output=True, text=True)
        raise GitError(
            f'could not revert commit {sha[:9]} cleanly (likely a conflict with later changes): '
            f'{result.stderr.strip() or result.stdout.strip()}')


def commit_all(repo_root, message: str):
    """Stage every change in the working tree (``git add -A``) and commit it."""
    _git(repo_root, 'add', '-A')
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
