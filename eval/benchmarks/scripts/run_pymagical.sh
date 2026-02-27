#!/bin/bash
#SBATCH --mem=512G
#SBATCH --output=logs/py_magical_%j.out

# Usage: sbatch run_pymagical.sh <iterations> [outdir] [prefix]

ITERATIONS=${1:-500}
OUTDIR=${2:-outputs}
PREFIX=${3:-"astrocytes"}

mkdir -p logs
mkdir -p "$OUTDIR"

echo "Starting pymagical with $ITERATIONS iterations (prefix: $PREFIX)"
uv run pymagical --iter "$ITERATIONS" --outdir "$OUTDIR" --prefix "$PREFIX"

echo "Job completed."
