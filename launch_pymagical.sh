#!/usr/bin/env bash
set -euo pipefail

# Usage: ./launch_pymagical.sh OUTPUT_DIR [ITERATIONS]
# Relative output paths are resolved from the caller's working directory.
if [[ $# -eq 1 && ( $1 == --help || $1 == -h ) ]]; then
    echo "Usage: $0 OUTPUT_DIR [ITERATIONS (default: 2000)]"
    exit 0
fi
if [[ $# -lt 1 || $# -gt 2 || -z $1 ]]; then
    echo "Usage: $0 OUTPUT_DIR [ITERATIONS (default: 2000)]" >&2
    exit 1
fi

output_dir=$1
iterations=${2:-2000}
if [[ ! $iterations =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: iterations must be a positive integer." >&2
    exit 1
fi
if [[ -e $output_dir || -L $output_dir ]]; then
    echo "Error: output path already exists: $output_dir" >&2
    exit 1
fi
command -v sbatch >/dev/null || { echo "Error: sbatch is required." >&2; exit 1; }
command -v uv >/dev/null || { echo "Error: uv is required." >&2; exit 1; }
repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
celltypes=(astrocytes excitatory_neurons inhibitory_neurons microglia opcs oligodendrocytes)
for ct in "${celltypes[@]}"; do
    if [[ ! -d $repo_dir/inputs/$ct ]]; then
        echo "Error: missing input directory: $repo_dir/inputs/$ct" >&2
        exit 1
    fi
done

# Plain mkdir atomically rejects existing paths, including concurrent launches.
mkdir -- "$output_dir"
output_dir=$(cd -- "$output_dir" && pwd)
for ct in "${celltypes[@]}"; do
    mkdir -- "$output_dir/$ct"
done

# Keep the submitted script and job ID with the results for later inspection.
cat > "$output_dir/run_job.sh" <<'JOB'
#!/usr/bin/env bash
set -euo pipefail
repo_dir=$1
output_dir=$2
iterations=$3
shift 3
celltypes=("$@")
ct=${celltypes[$SLURM_ARRAY_TASK_ID]}
exec > "$output_dir/$ct/run.log" 2>&1
trap 'status=$?; echo "Finished $ct at $(date), exit code $status"' EXIT
cd -- "$repo_dir"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMBA_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "Starting $ct at $(date); iterations=$iterations; job=$SLURM_JOB_ID"
uv run pymagical run \
    --main-dir "$repo_dir/inputs" \
    --cell-dir "$ct" \
    --iter "$iterations" \
    --use-numba \
    --outdir "$output_dir/$ct" \
    --prefix "$ct"
JOB

if ! job_id=$(sbatch --parsable \
    --partition=gen \
    --job-name=pymagical \
    --cpus-per-task=16 \
    --mem=512G \
    --array="0-$((${#celltypes[@]} - 1))" \
    --chdir="$repo_dir" \
    --output="$output_dir/slurm-%A_%a.log" \
    "$output_dir/run_job.sh" "$repo_dir" "$output_dir" "$iterations" "${celltypes[@]}"); then
    echo "Error: SLURM submission failed; launch files remain in $output_dir" >&2
    exit 1
fi
printf '%s\n' "$job_id" > "$output_dir/job_id.txt"
echo "Submitted CPU array job $job_id (${#celltypes[@]} cell types, $iterations iterations each)."
echo "Output: $output_dir"
