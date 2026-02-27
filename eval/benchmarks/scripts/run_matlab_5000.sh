#!/bin/bash
#SBATCH --mem=512G
#SBATCH --output=eval/benchmarks/logs/ml_magical_5000_%j.out
#SBATCH --job-name=ml_magical_5000

# Usage: sbatch run_matlab_5000.sh

ITERATIONS=5000
LABEL="astrocytes_ml_5000"

mkdir -p eval/benchmarks/logs

echo "Starting MATLAB MAGICAL with $ITERATIONS iterations (label: $LABEL)"
matlab -batch "addpath('eval/benchmarks/scripts'); matlab_runner($ITERATIONS, '$LABEL')"

echo "Job completed."
