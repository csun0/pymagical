#!/bin/bash
#SBATCH --mem=512G
#SBATCH --output=eval/benchmarks/logs/py_numpy_%j.out
#SBATCH --job-name=py_magical_numpy

ITERATIONS=2000
OUTDIR=outputs_bench/numpy
PREFIX="astrocytes_numpy"

mkdir -p eval/benchmarks/logs
mkdir -p "$OUTDIR"

echo "Starting pymagical (NumPy) with $ITERATIONS iterations"
uv run pymagical --iter "$ITERATIONS" --outdir "$OUTDIR" --prefix "$PREFIX"
echo "Job completed."
