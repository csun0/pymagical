# Evaluation and Benchmarking Suite

This directory contains tools for verifying the statistical fidelity and performance of the `pymagical` Python port compared to the original MATLAB implementation.

## Directory Structure

*   `tests/`: Correctness and fidelity verification.
    *   `test_data_loader.py`: Functional test for the caching data loader.
    *   `compare_results.py`: Main tool for comparing Python and MATLAB output matrices and circuits.
*   `benchmarks/`: Performance profiling and resource scaling.
    *   `profile_run.py`: Breaks down execution time by stage (Loading, Construction, Init, Sampling).
    *   `plot_comparison.py`: Generates comparative bar charts between implementations.
    *   `scripts/`: Generic Slurm (`.sh`) and MATLAB (`.m`) runners for launching cluster jobs.

## Usage Examples

### 1. Compare Fidelity and Performance

After running both Python and MATLAB versions for 500 iterations, use this script to calculate correlations and speedup:

```bash
uv run python eval/tests/compare_results.py 
    --ml-dir path/to/matlab/outputs 
    --py-dir path/to/python/outputs 
    --iter 500 
    --prefix astrocytes
```

### 2. Profile Python Execution Stages

To see a granular breakdown of where time is spent in the Python pipeline:

```bash
uv run python eval/benchmarks/profile_run.py --iter 100 --output my_profile.png
```

### 3. Launching Cluster Jobs

Use the provided Slurm wrappers to submit jobs to the queue:

```bash
# Submit Python job for 1000 iterations
sbatch eval/benchmarks/scripts/run_pymagical.sh 1000 /path/to/output astrocytes

# Submit MATLAB job for 1000 iterations
sbatch eval/benchmarks/scripts/run_matlab.sh 1000 astrocytes_ml_1000
```

## Metrics and Validation

For a detailed history of the validation results (Pearson correlations, triad recovery rates, and wall-clock speedups), see [docs/decisions.md](../docs/decisions.md).
