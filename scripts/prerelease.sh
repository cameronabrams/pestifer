#!/usr/bin/env bash
# Reversibly stamp the LOCAL editable install with a commit-bound pre-release version, so a
# verification run's logs (and the doc outputs they seed) name precisely the code under test.
#
# Usage:
#   scripts/prerelease.sh stamp <version>   # e.g. 3.9.0rc1  ->  stamps 3.9.0rc1+g<sha>
#   scripts/prerelease.sh unstamp           # restore pyproject.toml and reinstall
#   scripts/prerelease.sh status            # show pyproject / installed / git
#
# The runtime version comes from importlib.metadata (the installed package), NOT pyproject.toml
# directly, so it only changes on reinstall. This stamps <version>+g<short-sha> into
# pyproject.toml and reinstalls, so `pestifer --version` (and every build log) reads it. Requires
# a clean working tree so the embedded commit SHA identifies exactly the code being tested.
#
# It touches ONLY pyproject.toml (left uncommitted) and the editable install. It never edits
# CHANGELOG.md or CITATION.cff, and never commits, tags, or pushes -- i.e. nothing that
# release.sh does. And release.sh refuses to run on a dirty tree, so it cannot collide with a
# live stamp: always `unstamp` before releasing.
#
# Verification flow:
#   git switch main && git pull                # be on the exact commit you intend to release
#   scripts/prerelease.sh stamp 3.9.0rc1       # stamp 3.9.0rc1+g<sha>, reinstall, verify
#   scripts/run_all_examples.sh                # logs now read pestifer v. 3.9.0rc1+g<sha>
#   scripts/prerelease.sh unstamp              # restore clean pyproject + reinstall
#   # only if all examples passed:
#   scripts/release.sh 3.9.0                   # the real, irreversible release

set -euo pipefail
cd "$(dirname "$0")/.."

reinstall() { uv pip install -e . --quiet; }
# read the version the way the logs do: the installed package's metadata, via the same
# `uv run pestifer` entry point the example runner uses.
installed() { uv run pestifer --version 2>/dev/null | awk '{print $2}'; }
pyproj_ver() { grep -m1 '^version = ' pyproject.toml | sed 's/version = "\(.*\)"/\1/'; }
tree_clean() { git diff --quiet && git diff --cached --quiet; }

cmd="${1:-}"
case "$cmd" in
  stamp)
    base="${2:?Usage: scripts/prerelease.sh stamp <version>  (e.g. 3.9.0rc1)}"
    if ! tree_clean; then
      echo "ERROR: working tree is not clean." >&2
      echo "       Commit your changes first (so the stamp's commit SHA identifies exactly the" >&2
      echo "       code under test), or run 'unstamp' if a previous stamp is still applied." >&2
      git status --short >&2
      exit 1
    fi
    sha="$(git rev-parse --short=12 HEAD)"
    stamp="${base}+g${sha}"
    sed -i "s/^version = \".*\"/version = \"${stamp}\"/" pyproject.toml
    reinstall
    got="$(installed)"
    if [ "$got" != "$stamp" ]; then
      echo "ERROR: installed version is '${got}', expected '${stamp}'. Reverting." >&2
      git checkout -- pyproject.toml; reinstall
      exit 1
    fi
    echo "Stamped: ${stamp}   (HEAD ${sha})"
    echo "Every build log will now read:  pestifer v. ${stamp}"
    echo "When the run is done:            scripts/prerelease.sh unstamp"
    ;;
  unstamp)
    git checkout -- pyproject.toml
    reinstall
    echo "Reverted.  pyproject=$(pyproj_ver)   installed=$(installed)"
    ;;
  status)
    echo "pyproject : $(pyproj_ver)"
    echo "installed : $(installed)"
    echo "git HEAD  : $(git rev-parse --short=12 HEAD)   tree_clean=$(tree_clean && echo yes || echo no)"
    ;;
  *)
    echo "Usage: scripts/prerelease.sh {stamp <version>|unstamp|status}" >&2
    exit 1
    ;;
esac
