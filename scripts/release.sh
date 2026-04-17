#!/usr/bin/env bash
# Release pestifer at a given version.
#
# Usage: ./scripts/release.sh <version>
# Example: ./scripts/release.sh 2.2.3
#
# Prerequisites (checked automatically):
#   - Working tree must be clean (no uncommitted changes)
#   - Must be on the main branch
#   - CHANGELOG.md must have a "## [<version>]" entry
#
# What it does:
#   1. Updates the version in pyproject.toml
#   2. Commits the change as "Release v<version>"
#   3. Creates tag v<version>
#   4. Pushes the commit and the tag to origin
#
# The pushed tag triggers the release.yaml CI workflow which builds,
# publishes to PyPI, creates a GitHub Release, and triggers ReadTheDocs.

set -euo pipefail

VERSION="${1:?Usage: scripts/release.sh <version>  (e.g. 2.2.3)}"

# ── Preconditions ────────────────────────────────────────────────────────────

if ! git diff --quiet || ! git diff --cached --quiet; then
    echo "ERROR: working tree has uncommitted changes — commit or stash them first"
    exit 1
fi

BRANCH="$(git branch --show-current)"
if [ "$BRANCH" != "main" ]; then
    echo "ERROR: must be on main branch (currently on '$BRANCH')"
    exit 1
fi

if ! grep -q "^## \[$VERSION\]" CHANGELOG.md; then
    echo "ERROR: no entry for $VERSION found in CHANGELOG.md"
    echo "       Add a '## [$VERSION] - $(date +%Y-%m-%d)' section before releasing"
    exit 1
fi

if git rev-parse "v$VERSION" >/dev/null 2>&1; then
    echo "ERROR: tag v$VERSION already exists"
    exit 1
fi

# ── Version bump ─────────────────────────────────────────────────────────────

echo "Bumping pyproject.toml version to $VERSION"
sed -i "s/^version = \".*\"/version = \"$VERSION\"/" pyproject.toml

# Verify the sed worked
ACTUAL="$(python3 -c "import tomllib; print(tomllib.load(open('pyproject.toml','rb'))['project']['version'])")"
if [ "$ACTUAL" != "$VERSION" ]; then
    echo "ERROR: version in pyproject.toml is '$ACTUAL' after sed — check the file"
    git checkout pyproject.toml
    exit 1
fi

# ── Commit, tag, push ────────────────────────────────────────────────────────

git add pyproject.toml
git commit -m "Release v$VERSION"
git tag "v$VERSION"

echo "Pushing commit and tag v$VERSION to origin..."
git push origin main
git push origin "v$VERSION"

echo ""
echo "Done. The release.yaml workflow will now build and publish v$VERSION."
