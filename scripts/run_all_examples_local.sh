#!/usr/bin/env bash
#
# Run pestifer build-examples locally (no SLURM), each in its own subdirectory
# so outputs don't collide.  Local companion to run_all_examples.sh (the SLURM
# batch version).
#
# Usage:
#   ./run_all_examples_local.sh                # run every example, once
#   ./run_all_examples_local.sh 16 17          # run only examples 16 and 17
#   ./run_all_examples_local.sh 1-5 16-17      # ranges are allowed
#   REPLICAS=3 ./run_all_examples_local.sh     # three independent replicas of each
#
# Environment:
#   PESTIFER   command used to invoke pestifer (default: "pestifer").  With uv,
#              either activate the project venv first (so "pestifer" is on PATH),
#              or set e.g.  PESTIFER="uv run pestifer".
#   OUTROOT    directory under which the per-example subdirs are created
#              (default: the current directory).
#   REPLICAS   number of independent replicas of each example (default: 1).
#   SEED_BASE  base NAMD RNG seed (default: 27021972, pestifer's own default).
#              Replica r of every example is built with seed SEED_BASE+r-1.  The
#              seed used is echoed and recorded in each replica's build log.
#              NOTE: re-running with the same seed does NOT reproduce a build
#              exactly.  See "Reproducibility" below.
#
# Replicas vary the MD random-number streams -- velocity assignment and the
# Langevin thermostat.  Model-building seeds (loop closure, membrane packing)
# are NOT varied here; to vary those too, set ligate.ccd.seed and
# make_membrane_system.bilayer.seed in the example YAML.
#
# Reproducibility -- read before quoting any number from a sweep.
#
# Replicas do NOT all start from the same built structure, and the same seed
# does NOT reproduce a build exactly.  NAMD is not bitwise reproducible on
# multiple cores: on identical input with an identical seed it reproduces
# exactly at +p1 and does not at +p24 (measured 2026-08-25).  Every NAMD stage
# inherits this; whether it becomes visible depends on system size and stage
# length, with dynamics amplifying fastest.  Measured over a 27-example
# triplicate sweep:
#
#   - psfgen output IS bit-identical across replicas -- 27 of 27 examples.
#   - Solvated starting structures differ across replicas in 24 of the 25
#     examples that solvate (example 5 is the sole exception).
#   - Where the solvation box is fixed AFTER dynamics, the atom count differs
#     as well -- 11 of 27 examples, spanning 0.76-4.19%.
#
# So the built macromolecule is reproducible, and composition is reproducible
# only for the 16 examples whose box is set before any dynamics runs.  Nothing
# downstream of the first NAMD stage is reproducible bit-for-bit.
#
# With REPLICAS=1 the layout is unchanged: OUTROOT/example-NN.
# With REPLICAS>1 each replica gets its own directory: OUTROOT/example-NN/rep-MM.
#
# pestifer drives VMD and NAMD as subprocesses, so vmd and namd3 (and catdcd,
# for examples that use it) must be on PATH before running.

set -u

PESTIFER="${PESTIFER:-pestifer}"
OUTROOT="${OUTROOT:-.}"
REPLICAS="${REPLICAS:-1}"
SEED_BASE="${SEED_BASE:-27021972}"

case "$REPLICAS" in
    ''|*[!0-9]*) echo "ERROR: REPLICAS must be a positive integer (got '$REPLICAS')" >&2; exit 1 ;;
esac
[ "$REPLICAS" -ge 1 ] || { echo "ERROR: REPLICAS must be >= 1 (got '$REPLICAS')" >&2; exit 1; }
case "$SEED_BASE" in
    ''|*[!0-9]*) echo "ERROR: SEED_BASE must be a non-negative integer (got '$SEED_BASE')" >&2; exit 1 ;;
esac

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
if [ "$REPLICAS" -gt 1 ]; then
    echo "replicas:    $REPLICAS per example, seeds $SEED_BASE..$((SEED_BASE + REPLICAS - 1))"
fi
echo "started at   $(date)"
echo

passed=(); failed=()
overall_start=$SECONDS

for i in "${IDS[@]}"; do
    for r in $(seq 1 "$REPLICAS"); do
        if [ "$REPLICAS" -eq 1 ]; then
            dir="$OUTROOT/$(printf 'example-%02d' "$i")"
            label="Example $i"
            tag="$i"
        else
            dir="$OUTROOT/$(printf 'example-%02d/rep-%02d' "$i" "$r")"
            label="Example $i replica $r"
            tag="$i/$r"
        fi
        seed=$((SEED_BASE + r - 1))
        mkdir -p "$dir"
        echo "=== $label: starting at $(date)  seed $seed  ->  $dir ==="
        start=$SECONDS
        # --seed is passed even for a single run, so the seed that produced these
        # results is always on the record rather than implied by a default.
        ( cd "$dir" && $PESTIFER build-example "$i" --seed "$seed" ) 2>&1 | tee "$dir/build.log"
        status=${PIPESTATUS[0]}
        elapsed=$((SECONDS - start))
        if [ "$status" -eq 0 ]; then
            echo "=== $label: OK in ${elapsed}s ==="
            passed+=("$tag")
        else
            echo "=== $label: FAILED (exit $status) after ${elapsed}s -- see $dir/build.log ==="
            failed+=("$tag")
        fi
        echo
    done
done

total=$((SECONDS - overall_start))
echo "=========================================================="
echo "Finished at $(date)  (total ${total}s)"
echo "Passed (${#passed[@]}): ${passed[*]:-none}"
echo "Failed (${#failed[@]}): ${failed[*]:-none}"
[ "${#failed[@]}" -eq 0 ]
