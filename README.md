# pymagical

`pymagical` is a high-performance Python port of the **MAGICAL** (Multiome Accessibility Gene Integration Calling and Looping) algorithm. It infers functional regulatory circuits (triads of Transcription Factors, cis-regulatory elements or Peaks, and target Genes) from paired single-cell RNA-seq and ATAC-seq data.

The methodology is based on the framework described in:
> **Chen et al., "Mapping disease regulatory circuits at cell-type resolution from single-cell multiomics data" *Nature Computational Science* 2023.**
> (Available [here](https://www.nature.com/articles/s43588-023-00476-5))

## Key Features

*   **IO Caching:** Caches large sparse matrices and genomic metadata into PyArrow-backed Parquet and NumPy formats for near-instant subsequent loads (**~15x faster** than re-parsing text).
*   **Numba-Accelerated Sampling:** JIT-compiled kernels give **~28x faster** Gibbs sampling than the original MATLAB implementation (astrocytes, 2000 iterations).
*   **Biological Directionality:** Classifies inferred circuits as **activators (+)** or **repressors (-)** from the continuous regression weights.

## Documentation

For detailed information on setup, biological methodology, and validation, please refer to the following guides:

*   **[Getting Started Tutorial](TUTORIAL.md)**: A complete walkthrough for installing `pymagical` and running your first inference.
*   **[Statistical Fidelity & Matrix Definitions](docs/statistical_fidelity.md)**: Detailed explanation of $B$ and $L$ matrices and validation against MATLAB.
*   **[Methodology Overview](docs/methodology.md)**: Technical details on the hierarchical Bayesian Gibbs sampling framework.
*   **[Performance Report](docs/performance_report.md)**: Benchmarks comparing NumPy and Numba implementations against MATLAB.

## Installation

### For users (from PyPI)

Install `pymagical` in one step with `pip` (or `uv pip`):

```bash
pip install pymagical
```

For the interactive HTML report (`pymagical viz`), install the optional `viz` extra:

```bash
pip install "pymagical[viz]"
```

Requires Python ≥ 3.10 and a C compiler (Numba needs one for its JIT kernels).

### For developers (from source)

Clone the repository and sync the environment with [uv](https://docs.astral.sh/uv/); this
creates a `.venv` and installs the package (editable) plus the `dev` dependency
group (pytest, plotting, and the viz libraries):

```bash
git clone https://github.com/csun0/pymagical.git
cd pymagical
uv sync
```

Then run commands inside the environment with `uv run`:

```bash
uv run pymagical --help     # CLI
uv run pytest               # test suite
```

> The package version is derived from git tags via `hatch-vcs`. A full clone (or a
> release tarball) builds fine; if you build from a source tree with no git history
> the version falls back to `0.0.0`.

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

`run_magical` takes individual file paths (all required); see the [tutorial](TUTORIAL.md#6-advanced-usage-programmatic-api) for the full argument list.

```python
from pymagical import run_magical

run_magical(
    cand_gene_file="genes.txt",
    cand_peak_file="peaks.txt",
    # ... remaining RNA/ATAC/motif/TAD/refseq file paths (all required) ...
    iteration_num=2000,
    use_numba=True,
    output_file="my_results.txt",
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