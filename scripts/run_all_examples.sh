#!/bin/bash
#SBATCH --job-name=pestifer-examples
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=pestifer-examples-%j.out
#SBATCH --error=pestifer-examples-%j.err

# Run all 22 pestifer build-examples in series.
# Each example runs in its own subdirectory so outputs don't collide.
# Adjust #SBATCH directives and the environment setup below for your cluster.

# --- environment setup ---
module purge
module load vmd namd
conda activate panacea

PESTIFER=/home/cfa/anaconda3/envs/panacea/bin/pestifer
N_EXAMPLES=22

echo "Starting pestifer example runs at $(date)"

for i in $(seq 1 $N_EXAMPLES); do
    dir=$(printf "example-%02d" $i)
    mkdir -p "$dir"
    pushd "$dir" > /dev/null
    echo "--- Example $i: starting at $(date) ---"
    $PESTIFER build-example $i
    status=$?
    if [ $status -ne 0 ]; then
        echo "--- Example $i: FAILED (exit $status) at $(date) ---"
    else
        echo "--- Example $i: completed at $(date) ---"
    fi
    popd > /dev/null
done

echo "All examples finished at $(date)"
