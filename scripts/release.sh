#!/usr/bin/env bash
# Release pestifer at a given version.
#
# Usage: ./scripts/release.sh <version>
# Example: ./scripts/release.sh 2.4.4
#
# Prerequisites (checked automatically):
#   - Working tree must be clean (no uncommitted changes)
#   - Must be on the main branch
#   - CHANGELOG.md must have an "## [Unreleased]" section
#   - The release must fit under PyPI's project-size quota (build-free estimate; see below)
#
# What it does:
#   1. Rotates CHANGELOG.md: renames [Unreleased] to [<version>] - <date>
#      and inserts a fresh empty [Unreleased] section above it
#   2. Updates the version in pyproject.toml
#   3. Commits both changes as "Release v<version>"
#   4. Creates tag v<version>
#   5. Pushes the commit and the tag to origin
#
# The pushed tag triggers release.yaml, which runs tests, builds the package,
# publishes to PyPI, creates a GitHub Release with the CHANGELOG notes, and
# triggers a ReadTheDocs build.

set -euo pipefail

VERSION="${1:?Usage: scripts/release.sh <version>  (e.g. 2.4.4)}"
TODAY="$(date +%Y-%m-%d)"

# ── Preconditions ─────────────────────────────────────────────────────────────

if ! git diff --quiet || ! git diff --cached --quiet; then
    echo "ERROR: working tree has uncommitted changes — commit or stash them first"
    exit 1
fi

BRANCH="$(git branch --show-current)"
if [ "$BRANCH" != "main" ]; then
    echo "ERROR: must be on main branch (currently on '$BRANCH')"
    exit 1
fi

if ! grep -q "^## \[Unreleased\]" CHANGELOG.md; then
    echo "ERROR: no '## [Unreleased]' section found in CHANGELOG.md"
    exit 1
fi

if git rev-parse "v$VERSION" >/dev/null 2>&1; then
    echo "ERROR: tag v$VERSION already exists locally"
    exit 1
fi

if git ls-remote --tags origin "refs/tags/v$VERSION" | grep -q .; then
    echo "ERROR: tag v$VERSION already exists on origin"
    exit 1
fi

# ── PyPI project-size preflight ───────────────────────────────────────────────
# PyPI enforces a per-project total-size quota (default 10 GB, summed over every file of
# every release). An upload that would exceed it is rejected with HTTP 400 — but only in CI,
# *after* the tag is pushed and the package is built, leaving a dangling tag. Catch it here
# instead. No build is needed: a release's size is dominated by the static bundled CHARMM
# force field, so the most-recent published release is a reliable size proxy.
#   - Override the assumed limit (e.g. after a granted increase): PESTIFER_PYPI_SIZE_LIMIT_GB=20
#   - Skip the check entirely:                                    SKIP_PYPI_SIZE_CHECK=1
if [ "${SKIP_PYPI_SIZE_CHECK:-0}" != "1" ]; then
    PROJECT="$(grep -m1 '^name = ' pyproject.toml | sed 's/name = "\(.*\)"/\1/')"
    python3 - "$PROJECT" "${PESTIFER_PYPI_SIZE_LIMIT_GB:-10}" <<'PYEOF'
import json, sys, urllib.request
project, limit = sys.argv[1], float(sys.argv[2]) * 1e9
try:
    with urllib.request.urlopen(f"https://pypi.org/pypi/{project}/json", timeout=20) as r:
        data = json.load(r)
except Exception as e:
    print(f"WARNING: could not check PyPI project size ({e}); proceeding without the guard")
    sys.exit(0)
releases = data.get("releases", {})
current = sum(f["size"] for files in releases.values() for f in files)
latest = data.get("info", {}).get("version", "")
est = max(sum(f["size"] for f in releases.get(latest, [])), 40e6) * 1.25  # proxy + floor + margin
projected = current + est
gb = lambda x: x / 1e9
print(f"PyPI size guard: {gb(current):.2f} GB published + ~{est/1e6:.0f} MB new "
      f"= ~{gb(projected):.2f} GB of {gb(limit):.0f} GB "
      f"(free {gb(limit - current):.2f} GB)")
if projected > limit:
    print(f"ERROR: this release would exceed PyPI's {gb(limit):.0f} GB project-size limit for "
          f"'{project}'. Free space by deleting old releases on PyPI, or request an increase "
          f"(https://docs.pypi.org/project-management/storage-limits/), then re-run. "
          f"Adjust the assumed limit with PESTIFER_PYPI_SIZE_LIMIT_GB, or bypass with "
          f"SKIP_PYPI_SIZE_CHECK=1.")
    sys.exit(1)
PYEOF
fi

# ── CHANGELOG rotation ────────────────────────────────────────────────────────

echo "Rotating CHANGELOG.md: [Unreleased] -> [$VERSION] - $TODAY"
sed -i "s/^## \[Unreleased\]/## [$VERSION] - $TODAY/" CHANGELOG.md

# Insert a fresh [Unreleased] section before the new release
sed -i "s/^## \[$VERSION\] - $TODAY/## [Unreleased]\n\n## [$VERSION] - $TODAY/" CHANGELOG.md

# ── Version bump ──────────────────────────────────────────────────────────────

echo "Bumping pyproject.toml version to $VERSION"
sed -i "s/^version = \".*\"/version = \"$VERSION\"/" pyproject.toml

ACTUAL="$(grep '^version = ' pyproject.toml | sed 's/version = \"\(.*\)\"/\1/')"
if [ "$ACTUAL" != "$VERSION" ]; then
    echo "ERROR: version in pyproject.toml is '$ACTUAL' after sed — check the file"
    git checkout pyproject.toml CHANGELOG.md
    exit 1
fi

echo "Bumping CITATION.cff version to $VERSION"
sed -i "s/^version: \".*\"/version: \"$VERSION\"/" CITATION.cff
sed -i "s/^date-released: \".*\"/date-released: \"$TODAY\"/" CITATION.cff

# ── Commit, tag, push ─────────────────────────────────────────────────────────

git add pyproject.toml CHANGELOG.md CITATION.cff
git commit -m "Release v$VERSION"
git tag "v$VERSION"

echo "Pushing commit and tag v$VERSION to origin..."
git push origin main
git push origin "v$VERSION"

echo ""
echo "Done. The release.yaml workflow will now build and publish v$VERSION."
