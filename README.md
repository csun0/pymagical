# pymagical

`pymagical` is a high-performance Python port of the **MAGICAL** (Multiome Accessibility Gene Integration Calling and Looping) algorithm. It provides a method for inferring functional regulatory circuits—triads of Transcription Factors (TFs), cis-regulatory elements (Peaks), and target Genes—from single-cell RNA-seq and ATAC-seq data.

The methodology is based on the framework described in:
> **Chen et al., "Mapping disease regulatory circuits at cell-type resolution from single-cell multiomics data," *Nature Computational Science*, 2023.**
> (Available [here](https://www.nature.com/articles/s43588-023-00476-5))

## Key Features

*   **Numba-Accelerated Sampling:** Optional JIT-compiled kernels provide a **~30x speedup** in Gibbs sampling compared to the original MATLAB implementation (averaged across large-scale benchmarks).
*   **IO Caching:** Automatically caches large sparse matrices and genomic metadata into PyArrow-backed Parquet and NumPy formats for near-instant subsequent loads (**~15x faster** than MATLAB).
*   **Biological Directionality:** Unlike the original version, `pymagical` automatically classifies inferred circuits as **activators (+)** or **repressors (-)** by analyzing continuous regression weights.

## Documentation

For detailed information on setup, biological methodology, and validation, please refer to the following guides:

*   **[Getting Started Tutorial](TUTORIAL.md)**: A complete walkthrough for installing `pymagical` and running your first inference.
*   **[Statistical Fidelity & Matrix Definitions](docs/statistical_fidelity.md)**: Detailed explanation of $B$ and $L$ matrices and validation against MATLAB.
*   **[Methodology Overview](docs/methodology.md)**: Technical details on the hierarchical Bayesian Gibbs sampling framework.
*   **[Performance Report](docs/performance_report.md)**: Benchmarks comparing NumPy and Numba implementations against MATLAB.

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

## Citation
If you use **MAGICAL** in your research, please cite:

```bibtex
@article{chen_mapping_2023,
	title = {Mapping disease regulatory circuits at cell-type resolution from single-cell multiomics data},
	author = {Chen, Xi and Wang, Yuan and Cappuccio, Antonio and Cheng, Wan-Sze and Zamojski, Frederique Ruf and Nair, Venugopalan D. and Miller, Clare M. and Rubenstein, Aliza B. and Nudelman, German and Tadych, Alicja and Theesfeld, Chandra L. and Vornholt, Alexandria and George, Mary-Catherine and Ruffin, Felicia and Dagher, Michael and Chawla, Daniel G. and Soares-Schanoski, Alessandra and Spurbeck, Rachel R. and Ndhlovu, Lishomwa C. and Sebra, Robert and Kleinstein, Steven H. and Letizia, Andrew G. and Ramos, Irene and Fowler, Vance G. and Woods, Christopher W. and Zaslavsky, Elena and Troyanskaya, Olga G. and Sealfon, Stuart C.},
	journal = {Nature Computational Science},
	year = {2023},
	month = jul,
	doi = {10.1038/s43588-023-00476-5},
	url = {https://www.nature.com/articles/s43588-023-00476-5},
}