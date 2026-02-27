# GEMINI.md - pymagical Project Context

## Project Overview
`pymagical` is a high-performance Python port of the **MAGICAL** (Multiome Accessibility Gene Integration Calling and Looping) algorithm. It is designed to infer functional regulatory circuits—consisting of Transcription Factors (TFs), cis-regulatory elements (Peaks), and target Genes—from paired single-cell RNA-seq and ATAC-seq data using a hierarchical Bayesian Gibbs sampling framework.

### Main Technologies
- **Python 3.10**: Managed via `uv`.
- **Numerical/Scientific**: `numpy`, `scipy`, `pandas`, `statsmodels`.
- **Data Optimization**: `pyarrow` for Parquet caching of large genomic matrices.
- **HPC Support**: Standardized Slurm execution scripts.

### Architecture
- `pymagical/`: Core package containing modularized logic:
    - `data_loader.py`: Dual-path IO with Parquet/NPZ caching.
    - `circuits.py`: Candidate circuit construction logic.
    - `initialization.py`: OLS-based model parameter seeding.
    - `estimation.py`: Vectorized Gibbs sampler implementation.
    - `magical.py`: High-level execution and output orchestration.
- `eval/`: Comprehensive evaluation suite:
    - `tests/`: Statistical fidelity verification against MATLAB baseline.
    - `benchmarks/`: Performance profiling and resource scaling tools.
- `docs/`: Technical methodology and design decision records.

## Building and Running

### Environment Setup
Sync the project environment using `uv`:
```bash
uv sync
```

### Running the Pipeline
Run the circuit inference via the installed CLI:
```bash
uv run pymagical --iter 500 --prefix my_experiment --outdir outputs/
```
Alternatively, use the module entry point:
```bash
uv run python -m pymagical --iter 500
```

### Testing and Evaluation
To compare Python results with a MATLAB baseline:
```bash
uv run python eval/tests/compare_results.py --ml-dir path/to/ml_out --py-dir path/to/py_out --iter 2000
```

To profile execution stages:
```bash
uv run python eval/benchmarks/profile_run.py --iter 100
```

## Development Conventions

### Coding Style
- **Indentation**: 4 spaces.
- **Output Naming**: All generated artifacts follow the `{prefix}_{py|ml}_{iter}` pattern (e.g., `astrocytes_py_2000.txt`).
- **Biological Signage**: Circuit effects are reported as `Overall_Effect [L_dir, B_dir]`, where `+` is activation and `-` is repression.

### Performance Standards
- **IO Caching**: Always leverage the `data_loader.py` caching mechanism for large datasets.
- **Vectorization**: Favor NumPy vectorization over manual loops in the Gibbs sampler (`estimation.py`) to maintain performance parity with the optimized port.

### Statistical Fidelity
- Convergence and parity with the original MATLAB implementation is the highest priority. High-iteration runs (2000+) should maintain >0.99 Pearson correlation on continuous weights and >95% triad recovery.
