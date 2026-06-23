#!/usr/bin/env bash
#
# Run pestifer build-examples locally (no SLURM), each in its own subdirectory
# so outputs don't collide.  Local companion to run_all_examples.sh (the SLURM
# batch version).
#
# Usage:
#   ./run_all_examples_local.sh                # run every example
#   ./run_all_examples_local.sh 16 17          # run only examples 16 and 17
#   ./run_all_examples_local.sh 1-5 16-17      # ranges are allowed
#
# Environment:
#   PESTIFER   command used to invoke pestifer (default: "pestifer").  With uv,
#              either activate the project venv first (so "pestifer" is on PATH),
#              or set e.g.  PESTIFER="uv run pestifer".
#   OUTROOT    directory under which the per-example subdirs are created
#              (default: the current directory).
#
# pestifer drives VMD and NAMD as subprocesses, so vmd and namd3 (and catdcd,
# for examples that use it) must be on PATH before running.

set -u

PESTIFER="${PESTIFER:-pestifer}"
OUTROOT="${OUTROOT:-.}"

# --- resolve the list of example IDs to run --------------------------------
expand_arg() {            # expand "N" or "N-M" into a sequence of integers
    local a="$1"
    if [[ "$a" =~ ^([0-9]+)-([0-9]+)$ ]]; then
        seq "${BASH_REMATCH[1]}" "${BASH_REMATCH[2]}"
    elif [[ "$a" =~ ^[0-9]+$ ]]; then
        echo "$a"
    else
        echo "ignoring unrecognized argument: $a" >&2
    fi
}

IDS=()
if [ "$#" -gt 0 ]; then
    for a in "$@"; do
        while read -r n; do IDS+=("$n"); done < <(expand_arg "$a")
    done
else
    # auto-detect from the installed example listing (table rows start with the
    # numeric ID followed by an alphanumeric DBID)
    while read -r n; do IDS+=("$n"); done < <(
        $PESTIFER show-resources examples 2>/dev/null \
        | awk '/^[[:space:]]*[0-9]+[[:space:]]+[A-Za-z0-9]/ {print $1}' \
        | sort -n -u)
fi

if [ "${#IDS[@]}" -eq 0 ]; then
    echo "No examples to run (could not detect any; is '$PESTIFER' on PATH?)." >&2
    exit 1
fi

# --- check tooling ---------------------------------------------------------
command -v "${PESTIFER%% *}" >/dev/null 2>&1 \
    || { echo "ERROR: '${PESTIFER%% *}' not found on PATH." >&2; exit 1; }
for tool in vmd namd3; do
    command -v "$tool" >/dev/null 2>&1 \
        || echo "WARNING: '$tool' not on PATH; examples that need it will fail." >&2
done

mkdir -p "$OUTROOT"
echo "Running ${#IDS[@]} example(s): ${IDS[*]}"
echo "pestifer:    $PESTIFER"
echo "output root: $(cd "$OUTROOT" && pwd)"
echo "started at   $(date)"
echo

passed=(); failed=()
overall_start=$SECONDS

for i in "${IDS[@]}"; do
    dir="$OUTROOT/$(printf 'example-%02d' "$i")"
    mkdir -p "$dir"
    echo "=== Example $i: starting at $(date)  ->  $dir ==="
    start=$SECONDS
    ( cd "$dir" && $PESTIFER build-example "$i" ) 2>&1 | tee "$dir/build.log"
    status=${PIPESTATUS[0]}
    elapsed=$((SECONDS - start))
    if [ "$status" -eq 0 ]; then
        echo "=== Example $i: OK in ${elapsed}s ==="
        passed+=("$i")
    else
        echo "=== Example $i: FAILED (exit $status) after ${elapsed}s -- see $dir/build.log ==="
        failed+=("$i")
    fi
    echo
done

total=$((SECONDS - overall_start))
echo "=========================================================="
echo "Finished at $(date)  (total ${total}s)"
echo "Passed (${#passed[@]}): ${passed[*]:-none}"
echo "Failed (${#failed[@]}): ${failed[*]:-none}"
[ "${#failed[@]}" -eq 0 ]
