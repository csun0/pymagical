#!/bin/bash
#SBATCH --mem=256G
#SBATCH --cpus-per-task=16
#SBATCH --output=eval/benchmarks/logs/bench_%A_%a.out
#SBATCH --error=eval/benchmarks/logs/bench_%A_%a.err
#SBATCH --job-name=magical_bench
#SBATCH --array=0-104

# Datasets: 7
DATASETS=("astrocytes" "endothelial" "excitatory_neurons" "inhibitory_neurons" "microglia" "oligodendrocytes" "opcs")

# Implementations: 3
IMPLS=("matlab" "numpy" "numba")

# Iterations: 5
ITERS=(100 500 1000 2000 5000)

# Calculate indices
DATASET_IDX=$((SLURM_ARRAY_TASK_ID / 15))
REMAINING=$((SLURM_ARRAY_TASK_ID % 15))
IMPL_IDX=$((REMAINING / 5))
ITER_IDX=$((REMAINING % 5))

DATASET=${DATASETS[$DATASET_IDX]}
IMPL=${IMPLS[$IMPL_IDX]}
ITER=${ITERS[$ITER_IDX]}

echo "Running Benchmark: Dataset=$DATASET, Impl=$IMPL, Iter=$ITER"

MAIN_DIR="/mnt/home/agebrain/ceph/anderson/snmulti_data/processed/magical/inputs"
mkdir -p eval/benchmarks/logs

if [ "$IMPL" == "matlab" ]; then
    OUTDIR="outputs_bench/matlab"
    mkdir -p "$OUTDIR"
    LABEL="${DATASET}_ml_${ITER}"
    echo "Starting MATLAB..."
    matlab -batch "addpath('eval/benchmarks/scripts'); matlab_runner('$DATASET', $ITER, '$LABEL')"
elif [ "$IMPL" == "numpy" ]; then
    OUTDIR="outputs_bench/numpy"
    mkdir -p "$OUTDIR"
    PREFIX="${DATASET}_py_${ITER}"
    echo "Starting NumPy..."
    uv run pymagical --main-dir "$MAIN_DIR" --cell-dir "$DATASET" --iter "$ITER" --outdir "$OUTDIR" --prefix "$PREFIX"
elif [ "$IMPL" == "numba" ]; then
    OUTDIR="outputs_bench/numba"
    mkdir -p "$OUTDIR"
    PREFIX="${DATASET}_py_${ITER}"
    echo "Starting Numba..."
    uv run pymagical --main-dir "$MAIN_DIR" --cell-dir "$DATASET" --iter "$ITER" --use-numba --outdir "$OUTDIR" --prefix "$PREFIX"
fi

echo "Job $SLURM_ARRAY_TASK_ID completed."
