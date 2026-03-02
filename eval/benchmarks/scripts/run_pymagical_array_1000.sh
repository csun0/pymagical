#!/bin/bash
#SBATCH --mem=512G
#SBATCH --output=eval/benchmarks/logs/py_numba_array_1000_%A_%a.out
#SBATCH --job-name=py_magical_1000
#SBATCH --array=0-6
#SBATCH --cpus-per-task=8

CELLTYPES=(astrocytes endothelial excitatory_neurons inhibitory_neurons microglia opcs oligodendrocytes)
CELLTYPE=${CELLTYPES[$SLURM_ARRAY_TASK_ID]}

ITERATIONS=1000
MAIN_DIR="/mnt/home/agebrain/ceph/anderson/snmulti_data/processed/magical/inputs"
OUTDIR="/mnt/home/agebrain/ceph/anderson/snmulti_data/processed/magical/outputs_1000/$CELLTYPE"
PREFIX="${CELLTYPE}_numba"

mkdir -p eval/benchmarks/logs
mkdir -p "$OUTDIR"

echo "Starting pymagical (Numba) for $CELLTYPE with $ITERATIONS iterations"
uv run pymagical --main-dir "$MAIN_DIR" --cell-dir "$CELLTYPE" --iter "$ITERATIONS" --use-numba --outdir "$OUTDIR" --prefix "$PREFIX"
echo "Job completed."
