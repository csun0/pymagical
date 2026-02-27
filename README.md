# pymagical

`pymagical` is a high-performance Python port of the **MAGICAL** (Multiome Accessibility Gene Integration Calling and Looping) algorithm. It provides an automated pipeline for inferring functional regulatory circuits—triads of Transcription Factors (TFs), cis-regulatory elements (Peaks), and target Genes—from paired single-cell RNA-seq and ATAC-seq data.

## Key Features

*   **Fast & Optimized:** Significantly faster than the original MATLAB implementation (~1.9x speedup in sampling, up to 15x speedup in data loading).
*   **Intelligent IO Caching:** Automatically caches large sparse matrices and genomic metadata into PyArrow-backed Parquet and NumPy formats for near-instant subsequent loads.
*   **Biological Directionality:** Unlike the original version, `pymagical` automatically classifies inferred circuits as **activators (+)** or **repressors (-)** by analyzing the continuous regression weights calculated during Gibbs sampling.
*   **HPC Ready:** Includes built-in support for high-memory Slurm environments and allows for detailed weight-history dumping for downstream statistical distribution analysis.

## Installation

This project uses `uv` for environment management.

```bash
# Clone the repository
git clone <repo-url>
cd pymagical

# Sync the environment
uv sync
```

## Quick Start

### 1. Command Line Usage

Run the circuit inference directly from your terminal:

```bash
# Run with default demo data (astrocytes) for 500 iterations
uv run pymagical --iter 500 --outdir results/

# Run with custom data and dump weight history
uv run pymagical \
    --iter 1000 \
    --prefix my_sample \
    --rna-counts path/to/rna.txt \
    --atac-counts path/to/atac.txt \
    --dump-weights
```

### 2. Programmatic Usage

Integrate `pymagical` into your own Python pipelines:

```python
from pymagical import run_magical

run_magical(
    cand_gene_file="genes.txt",
    cand_peak_file="peaks.txt",
    # ... other file paths ...
    iteration_num=500,
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

## Evaluation and Comparison

The `eval/` directory contains tools for comparing `pymagical` against the original MATLAB implementation and profiling performance.

*   `eval/tests/compare_results.py`: Compare fidelity and performance across implementations.
*   `eval/benchmarks/profile_run.py`: Profile the runtime of different execution stages.
