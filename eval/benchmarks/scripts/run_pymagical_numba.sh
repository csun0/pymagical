#!/bin/bash
#SBATCH --mem=512G
#SBATCH --output=eval/benchmarks/logs/py_numba_%j.out
#SBATCH --job-name=py_magical_numba

ITERATIONS=2000
OUTDIR=outputs_bench/numba
PREFIX="astrocytes_numba"

mkdir -p eval/benchmarks/logs
mkdir -p "$OUTDIR"

echo "Starting pymagical (Numba) with $ITERATIONS iterations"
uv run pymagical --iter "$ITERATIONS" --use-numba --outdir "$OUTDIR" --prefix "$PREFIX"
echo "Job completed."
