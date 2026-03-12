#!/bin/bash
#SBATCH --mem=512G
#SBATCH --output=eval/benchmarks/logs/ml_magical_%j.out
#SBATCH --job-name=ml_magical

# ==============================================================================
# MATLAB MAGICAL Master Launch Script
# ==============================================================================
# This script is the single entry point for running MATLAB MAGICAL on Slurm.
# To alter iterations, datasets, or labels, edit the variables below.

# --- 1. Core Configuration ---
ITERATIONS=1000
DATASETS=(astrocytes endothelial excitatory_neurons inhibitory_neurons microglia opcs oligodendrocytes)

# --- 2. Dataset Selection ---
# OPTION A: Single Cell Type
# Set CELLTYPE manually and do NOT use --array in sbatch
CELLTYPE="astrocytes"

# OPTION B: Slurm Array Job (Uncomment below to enable)
# If using 'sbatch --array=0-6 run_matlab.sh', CELLTYPE will be set automatically.
if [ ! -z "$SLURM_ARRAY_TASK_ID" ]; then
    CELLTYPE=${DATASETS[$SLURM_ARRAY_TASK_ID]}
fi

# --- 3. Output Label ---
# The label is used to name the output files in eval/benchmarks/outputs/matlab/
LABEL="${CELLTYPE}_ml_${ITERATIONS}"

# ==============================================================================
# Execution Logic (Do not edit below unless necessary)
# ==============================================================================

mkdir -p eval/benchmarks/logs

echo "----------------------------------------------------------------"
echo "Starting MATLAB MAGICAL"
echo "Dataset:    $CELLTYPE"
echo "Iterations: $ITERATIONS"
echo "Label:      $LABEL"
echo "----------------------------------------------------------------"

matlab -batch "addpath('eval/benchmarks/scripts'); matlab_runner('$CELLTYPE', $ITERATIONS, '$LABEL')"

echo "Job completed at $(date)"
