# GEMINI.md - pymagical Project Context

## Project Overview
`pymagical` is a high-performance Python port of the **MAGICAL** (Multiome Accessibility Gene Integration Calling and Looping) algorithm. It is designed to infer functional regulatory circuits—consisting of Transcription Factors (TFs), cis-regulatory elements (Peaks), and target Genes—from paired single-cell RNA-seq and ATAC-seq data using a hierarchical Bayesian Gibbs sampling framework.

The original methodology is documented in [docs/MAGICAL.pdf](docs/MAGICAL.pdf).

### Main Technologies
- **Python 3.10**: Managed via `uv`.
- **Numerical/Scientific**: `numpy`, `scipy`, `pandas`, `statsmodels`, `numba`.
- **Data Optimization**: `pyarrow` for Parquet caching of large genomic matrices.
- **HPC Support**: Standardized Slurm execution scripts with centralized logging.

### Architecture
- `pymagical/`: Core package containing modularized logic:
    - `data_loader.py`: Dual-path IO with Parquet/NPZ caching.
    - `circuits.py`: Candidate circuit construction logic.
    - `initialization.py`: OLS-based model parameter seeding.
    - `estimation.py`: Vectorized and Numba-aware Gibbs sampler.
    - `estimation_kernels.py`: JIT-compiled Numba kernels for high-performance sampling.
    - `magical.py`: High-level execution and output orchestration.
- `eval/`: Comprehensive evaluation suite:
    - `tests/`: Statistical fidelity verification against MATLAB baseline.
    - `benchmarks/`: Performance profiling and resource scaling tools.
    - `benchmarks/logs/`: Centralized Slurm job output directory.
- `docs/`: Technical methodology, design decision records, and performance reports.

## Building and Running

### Environment Setup
Sync the project environment using `uv`:
```bash
uv sync
```

### Running the Pipeline
Run the circuit inference via the installed CLI (use `--use-numba` for 28x speedup):
```bash
uv run pymagical --iter 2000 --use-numba --prefix my_experiment --outdir outputs/
```

### Testing and Evaluation
To compare Python results with a MATLAB baseline:
```bash
uv run python eval/tests/compare_results.py --ml-dir path/to/ml_out --py-dir path/to/py_out --iter 2000
```

## Development Conventions

### Coding Style
- **Indentation**: 4 spaces.
- **Output Naming**: Generated artifacts follow the `{prefix}_{py|ml}_{iter}` pattern.
- **Biological Signage**: Circuit effects are reported as `Overall_Effect [L_dir, B_dir]`.

### Performance Standards
- **Numba Acceleration**: Core sampling loops MUST be implemented in `estimation_kernels.py` using `@njit(parallel=True)` where safe.
- **Memory Contiguity**: Always ensure matrices are `np.ascontiguousarray` (C-order) before passing to Numba kernels.
- **Running Residuals**: The sampler should maintain running residuals to avoid $O(N^3)$ operations, with periodic resets to mitigate numerical drift.

### Statistical Fidelity
- High-iteration runs (2000+) MUST maintain >0.99 Pearson correlation on continuous weights and >95% triad recovery compared to the MATLAB implementation.
