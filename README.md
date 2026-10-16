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

Requires Python ≥ 3.10.

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

```bash
# Setup a new workspace
mkdir pymagical_demo
cd pymagical_demo

# Download demo input files
curl -fLC - --retry 3 -O \
https://github.com/csun0/pymagical/releases/download/data-v1.0/magical_demo_input_files.tar.gz
tar -xzf magical_demo_input_files.tar.gz

# Install pymagical in a virtual environment
python -m venv .venv
source .venv/bin/activate
pip install pymagical

# Run with default data for 500 iterations using Numba
pymagical run --main-dir Demo_input_files --cell-dir Demo_input_files --iter 500 --use-numba --outdir results/

# Create the circuit dashboard for all result sets in a directory
pymagical dashboard --input_dir results/
```

Run `pymagical --help` to see all available flags and subcommands.

The dashboard is included in the installed package; no clone or separate update
script is needed. `--input_dir` is required (`--input-dir` also works). It reads
result `.txt` files directly from that directory, skipping matrix and timing
sidecars, saves `results/dashboard.html`, and opens it in your browser. Run again
after adding or updating results to refresh the saved dashboard.

On a remote machine, use `pymagical dashboard --input_dir results/ --no-open`,
then copy the saved HTML to your computer and open it. Use `--output PATH.html`
to save elsewhere. Charts, tables, and fonts use CDNs, so the browser needs
internet access. No `[viz]` extra is required for the dashboard.

The legacy helper accepts the same required input flag:
`uv run src/pymagical/ui/update_ui.py --input_dir results/ --no-open`.

### 2. Run the Downloadable Demo

From the cloned repository root, run `uv sync` as described above, then download
and extract the demo inputs (about 353 MB compressed, 1.5 GB extracted):

```bash
curl --fail --location --retry 3 --continue-at - \
  --output magical_demo_input_files.tar.gz \
  https://github.com/csun0/pymagical/releases/download/data-v1.0/magical_demo_input_files.tar.gz
tar -xzf magical_demo_input_files.tar.gz
```

The CLI defaults match the demo filenames. All files are in the same folder,
so use `--cell-dir Demo_input_files` with `--main-dir Demo_input_files`. Both
directory options resolve independently from the working directory. The bundled
`hg38_Refseq.txt` is header-free and ready to use. Start with a short smoke run:

```bash
uv run pymagical run \
  --main-dir Demo_input_files --cell-dir Demo_input_files \
  --iter 10 --use-numba \
  --prefix demo_smoke --outdir outputs/demo_smoke
```

After the smoke run succeeds, rerun the same command with `--iter 2000`,
`--prefix demo`, and `--outdir outputs/demo`. The main result will be
`outputs/demo/demo_py_2000.txt`, alongside B/L matrices and timing statistics.
Keep `Demo_input_files/.magical_cache/` to reuse parsed inputs on subsequent
runs. Sampling starts a new chain each time; it does not resume from the smoke run.

### 3. Programmatic Usage

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
