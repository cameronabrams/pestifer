#!/usr/bin/env bash
# Trim an upstream CHARMM "toppar_c36_<release>.tgz" down to just the files pestifer uses.
#
# The upstream toppar distribution bundles many force fields and methods pestifer never
# loads -- polarizable (drude/, cheq/), implicit solvent (gbsw/, openmm_gbsaobc2/, ace/),
# metal-site parameters (metals/), alternate water models and non-CHARMM ions (non_charmm/),
# silicates/, tamdfff/, larmord/, rush/, charmm_modifications/. These live in subdirectories
# alongside the top-level topology/parameter/stream files and the stream/ subtree, which are
# the only things pestifer parses. Dropping the rest takes the jul24 tarball from 11 MB to
# 3.1 MB (feb26 similarly) and is verified to yield byte-identical build topology.
#
# What is kept: every top-level *.rtf / *.prm / *.str, plus the whole stream/ subtree
# (stream/prot, stream/carb, stream/na, stream/lipid, stream/misc, ...). Everything under
# stream/ is retained because pestifer resolves residue-defining files from there on demand;
# the two files that patches/ targets also live under stream/.
#
# Usage:
#   scripts/trim_toppar.sh <full_toppar_c36_*.tgz> <output_trimmed.tgz>
#
# To restore a dropped subdirectory later, re-run against the full upstream tarball with an
# edited KEEP list, or extract the needed file from the upstream distribution directly.
set -euo pipefail

SRC="${1:?usage: trim_toppar.sh <full_toppar_c36_*.tgz> <output_trimmed.tgz>}"
OUT="${2:?usage: trim_toppar.sh <full_toppar_c36_*.tgz> <output_trimmed.tgz>}"

[ -f "$SRC" ] || { echo "ERROR: source tarball not found: $SRC" >&2; exit 1; }

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

tar xzf "$SRC" -C "$WORK"
ROOT="$WORK/toppar"
[ -d "$ROOT" ] || { echo "ERROR: expected a top-level 'toppar/' directory inside $SRC" >&2; exit 1; }

# Drop every immediate subdirectory of toppar/ except stream/ (top-level files are untouched).
find "$ROOT" -mindepth 1 -maxdepth 1 -type d ! -name stream -exec rm -rf {} +

# Repackage, preserving the toppar/ root the loader expects.
tar czf "$OUT" -C "$WORK" toppar

echo "Trimmed $(du -h "$SRC" | cut -f1) -> $(du -h "$OUT" | cut -f1)  ($OUT)"
