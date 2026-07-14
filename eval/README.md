# Evaluation and Benchmarking Suite

This directory contains tools for verifying the statistical fidelity and performance of the `pymagical` Python port compared to the original MATLAB implementation.

## Directory Structure

*   `tests/`: Correctness and fidelity verification.
    *   `test_data_loader.py`: Functional test for the caching data loader.
    *   `compare_results.py`: Main tool for comparing Python and MATLAB output matrices and circuits.
*   `benchmarks/`: Performance profiling and resource scaling.
    *   `profile_run.py`: Breaks down execution time by stage (Loading, Construction, Init, Sampling).
    *   `plot_comparison.py`: Generates comparative bar charts between implementations.
    *   `scripts/`: Slurm (`.sh`) and MATLAB (`.m`) runners for launching cluster jobs.
    *   `logs/`: Centralized directory for Slurm job outputs.

## Usage Examples

### 1. Compare Fidelity and Performance

After running both Python and MATLAB versions for 2000 iterations, use this script to calculate correlations and speedup:

```bash
uv run python eval/tests/compare_results.py \
    --ml-dir path/to/matlab/outputs \
    --py-dir path/to/python/outputs \
    --iter 2000 \
    --ml-prefix astrocytes \
    --py-prefix astrocytes_numba
```

### 2. Profile Python Execution Stages

To see a granular breakdown of where time is spent in the Python pipeline:

```bash
uv run python eval/benchmarks/profile_run.py --iter 100 --output my_profile.png
```

### 3. Launching Cluster Jobs

Submit the Slurm runners. Iterations, cell type, Numba on/off, and paths are set by editing the variables near the top of each script (`ITERATIONS`, `CELLTYPE`, `USE_NUMBA`):

```bash
# Python job (set USE_NUMBA=true/false inside the script)
sbatch eval/benchmarks/scripts/run_pymagical.sh

# Original MATLAB job
sbatch eval/benchmarks/scripts/run_matlab.sh

# Run all seven cell types as a Slurm array
sbatch --array=0-6 eval/benchmarks/scripts/run_pymagical.sh
```

Job output lands in `eval/benchmarks/logs/`. Use `scripts/compare_all.py` to aggregate results across cell types.

## Metrics and Validation

For a detailed history of the validation results and speedups, see:
*   [Performance Report](../docs/performance_report.md)
*   [Optimization Technical Notes](../docs/optimization_tech_notes.md)
*   [Design Decisions Log](../docs/decisions.md)
