# pymagical

`pymagical` is a high-performance Python port of the **MAGICAL** (Multiome Accessibility Gene Integration Calling and Looping) algorithm. It provides a method for inferring functional regulatory circuits—triads of Transcription Factors (TFs), cis-regulatory elements (Peaks), and target Genes—from single-cell RNA-seq and ATAC-seq data.

The methodology is based on the framework described in:
> **Chen et al., "Mapping disease regulatory circuits at cell-type resolution from single-cell multiomics data," *Nature Computational Science*, 2023.**
> (Available [here](https://www.nature.com/articles/s43588-023-00476-5))

## Key Features

*   **Numba-Accelerated Sampling:** Optional JIT-compiled kernels provide a **~28x speedup** in Gibbs sampling compared to the original MATLAB implementation.
*   **Intelligent IO Caching:** Automatically caches large sparse matrices and genomic metadata into PyArrow-backed Parquet and NumPy formats for near-instant subsequent loads (**15x faster** than MATLAB).
*   **Biological Directionality:** Unlike the original version, `pymagical` automatically classifies inferred circuits as **activators (+)** or **repressors (-)** by analyzing continuous regression weights.

## Installation

The recommended way to install `pymagical` is via [PyPI](https://pypi.org/project/pymagical/):

```bash
pip install pymagical
```

For interactive HTML reports and visualization features, install the `viz` extra:

```bash
pip install "pymagical[viz]"
```

## Quick Start

### 1. Command Line Usage

Once installed, run the circuit inference directly from your terminal. Use `--use-numba` for maximum performance:

```bash
# Run with default data for 500 iterations using Numba
pymagical run --main-dir path/to/data --cell-dir astrocytes --iter 500 --use-numba --outdir results/

# Generate an interactive HTML visualization report (requires [viz] extra)
pymagical viz results/magical_py_500.txt
```

Run `pymagical --help` to see all available flags and subcommands.

### 2. Programmatic Usage

```python
from pymagical import run_magical

run_magical(
    cand_gene_file="genes.txt",
    cand_peak_file="peaks.txt",
    # ... other file paths ...
    iteration_num=2000,
    use_numba=True,
    output_file="my_results.txt"
)
```

## Output Notation

The final triad list includes a biological effect annotation for every identified TF:

`TF_Name (Confidence_Probability, Overall_Effect [L_dir, B_dir])`

*   **Overall Effect:** `+` (Activator) or `-` (Repressor).
*   **L_dir (Looping):** Direction of Peak-to-Gene effect.
*   **B_dir (Binding):** Direction of TF-to-Peak effect.

*Example:* `STAT5B (0.85, + [+,+])` indicates an 85% confident activator that opens a peak which subsequently increases gene expression.

## Documentation

*   [Methodology and Notation Details](docs/methodology.md)
*   [Design Decisions and Benchmarks](docs/decisions.md)
*   [Performance Optimization Technical Notes](docs/optimization_tech_notes.md)
*   [Executive Performance Report](docs/performance_report.md)

## Evaluation and Comparison

The `eval/` directory contains tools for verifying the Python implementation against the original MATLAB baseline.

*   `eval/tests/compare_results.py`: Compare statistical fidelity and triad overlap.
*   `eval/benchmarks/profile_run.py`: Profile stage-by-stage execution runtime.
*   `eval/benchmarks/scripts/`: Slurm submission scripts for large-scale benchmarks.
