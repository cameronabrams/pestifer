"""
Some functions for interacting with git
"""
import subprocess
from pathlib import Path

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
