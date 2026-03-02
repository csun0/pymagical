#!/bin/bash
#SBATCH --mem=512G
#SBATCH --output=eval/benchmarks/logs/py_numpy_5000_%j.out
#SBATCH --job-name=py_magical_numpy_5000

ITERATIONS=5000
OUTDIR=outputs_bench/numpy_5000
PREFIX="astrocytes_numpy"

mkdir -p eval/benchmarks/logs
mkdir -p "$OUTDIR"

echo "Starting pymagical (NumPy) with $ITERATIONS iterations"
uv run pymagical --main-dir /mnt/ceph/users/agebrain/anderson/snmulti_data/pymagical/test_data/ --cell-dir oligodendrocytes --iter "$ITERATIONS" --outdir "$OUTDIR" --prefix "$PREFIX"
echo "Job completed."
