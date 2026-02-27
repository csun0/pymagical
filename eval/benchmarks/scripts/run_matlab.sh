#!/bin/bash
#SBATCH --mem=512G
#SBATCH --output=eval/benchmarks/logs/ml_magical_%j.out

# Usage: sbatch run_matlab.sh <iterations> [output_label]

ITERATIONS=${1:-500}
LABEL=${2:-"astrocytes_ml_$ITERATIONS"}

mkdir -p eval/benchmarks/logs

echo "Starting MATLAB MAGICAL with $ITERATIONS iterations (label: $LABEL)"
matlab -batch "addpath('eval/benchmarks/scripts'); matlab_runner($ITERATIONS, '$LABEL')"

echo "Job completed."
