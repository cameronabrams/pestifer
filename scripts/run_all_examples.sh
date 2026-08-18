#!/bin/bash
#SBATCH --job-name=pestifer-examples
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=pestifer-examples-%j.out
#SBATCH --error=pestifer-examples-%j.err

# Run every pestifer build-example in series, optionally as several independent
# replicas of each.  Each run gets its own subdirectory so outputs don't collide.
# Adjust #SBATCH directives and the environment setup below for your cluster.
#
#   REPLICAS   independent replicas of each example (default: 1)
#   SEED_BASE  base NAMD RNG seed (default: 27021972).  Replica r is built with
#              seed SEED_BASE+r-1, so rerunning with the same SEED_BASE
#              reproduces a replica exactly.  Replicas differ only in their MD
#              random-number streams; model building is not varied.
#
# With REPLICAS=1 the layout is unchanged: example-NN.
# With REPLICAS>1 each replica gets its own directory: example-NN/rep-MM.

# --- environment setup (adapt for your cluster) ---
module purge
module load vmd namd

# Invoke pestifer through uv.  Either activate the project venv beforehand so
# that `pestifer` is on PATH, or point uv at the project here.  Override the
# PESTIFER variable to taste, e.g.:
#   export PESTIFER="uv run --project $HOME/pestifer pestifer"
PESTIFER="${PESTIFER:-uv run pestifer}"
REPLICAS="${REPLICAS:-1}"
SEED_BASE="${SEED_BASE:-27021972}"

# Discover the example count rather than hard-coding it; a stale literal here
# silently skipped every example added since it was written.
N_EXAMPLES=$($PESTIFER show-resources examples 2>/dev/null \
    | awk '/^[[:space:]]*[0-9]+[[:space:]]+[A-Za-z0-9]/ {print $1}' | sort -n | tail -1)
N_EXAMPLES="${N_EXAMPLES:-0}"
if [ "$N_EXAMPLES" -lt 1 ]; then
    echo "ERROR: could not determine the number of examples; is '$PESTIFER' working?" >&2
    exit 1
fi

echo "Starting pestifer example runs at $(date)"
echo "$N_EXAMPLES example(s), $REPLICAS replica(s) each, seeds from $SEED_BASE"

for i in $(seq 1 $N_EXAMPLES); do
    for r in $(seq 1 "$REPLICAS"); do
        if [ "$REPLICAS" -eq 1 ]; then
            dir=$(printf "example-%02d" $i)
            label="Example $i"
        else
            dir=$(printf "example-%02d/rep-%02d" $i $r)
            label="Example $i replica $r"
        fi
        seed=$((SEED_BASE + r - 1))
        mkdir -p "$dir"
        pushd "$dir" > /dev/null
        echo "--- $label: starting at $(date), seed $seed ---"
        $PESTIFER build-example $i --seed "$seed"
        status=$?
        if [ $status -ne 0 ]; then
            echo "--- $label: FAILED (exit $status) at $(date) ---"
        else
            echo "--- $label: completed at $(date) ---"
        fi
        popd > /dev/null
    done
done

echo "All examples finished at $(date)"
