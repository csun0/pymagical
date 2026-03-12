#!/bin/bash
#SBATCH --mem=512G
#SBATCH --cpus-per-task=16
#SBATCH --output=eval/benchmarks/logs/py_magical_%j.out
#SBATCH --job-name=py_magical

# ==============================================================================
# pymagical Master Launch Script
# ==============================================================================
# This script is the single entry point for running pymagical on Slurm.
# To alter iterations, datasets, or performance modes, edit the variables below.

# --- 1. Core Configuration ---
ITERATIONS=1000
USE_NUMBA=true         # true/false (set to true for ~28x speedup)
MAIN_DIR="/mnt/home/agebrain/ceph/anderson/snmulti_data/processed/magical/inputs"
OUTDIR_BASE="eval/benchmarks/outputs/python"

# --- 2. Dataset Selection ---
# Array of all available cell types
CELLTYPES=(astrocytes endothelial excitatory_neurons inhibitory_neurons microglia opcs oligodendrocytes)

# OPTION A: Single Cell Type
# Set CELLTYPE manually and do NOT use --array in sbatch
CELLTYPE="astrocytes"

# OPTION B: Slurm Array Job (Uncomment below to enable)
# If using 'sbatch --array=0-6 run_pymagical.sh', CELLTYPE will be set automatically.
if [ ! -z "$SLURM_ARRAY_TASK_ID" ]; then
    CELLTYPE=${CELLTYPES[$SLURM_ARRAY_TASK_ID]}
fi

# --- 3. Output Configuration ---
OUTDIR="$OUTDIR_BASE/$CELLTYPE"
PREFIX="${CELLTYPE}_py"

# ==============================================================================
# Execution Logic (Do not edit below unless necessary)
# ==============================================================================

mkdir -p eval/benchmarks/logs
mkdir -p "$OUTDIR"

NUMBA_FLAG=""
if [ "$USE_NUMBA" = true ]; then
    NUMBA_FLAG="--use-numba"
fi

echo "----------------------------------------------------------------"
echo "Starting pymagical"
echo "Dataset:    $CELLTYPE"
echo "Iterations: $ITERATIONS"
echo "Numba:      $USE_NUMBA"
echo "Output:     $OUTDIR"
echo "----------------------------------------------------------------"

uv run pymagical run \
    --main-dir "$MAIN_DIR" \
    --cell-dir "$CELLTYPE" \
    --iter "$ITERATIONS" \
    $NUMBA_FLAG \
    --outdir "$OUTDIR" \
    --prefix "$PREFIX"

echo "Job completed at $(date)"
