#!/bin/bash
#SBATCH --mem=512G
#SBATCH --array=1-1

ARRAY_NUM=$SLURM_ARRAY_TASK_ID
TOTAL_ARRAY_NUM=$SLURM_ARRAY_TASK_MAX

time matlab -batch "snmulti_main($ARRAY_NUM, 500)"