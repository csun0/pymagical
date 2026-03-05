# pymagical: Getting Started Tutorial

`pymagical` is a high-performance Python implementation of the MAGICAL algorithm for inferring functional regulatory circuits (TF-Peak-Gene triads) from paired single-cell RNA-seq and ATAC-seq data.

This tutorial walks you through setting up the environment, preparing your data, running the inference engine, and validating the results.

---

## 1. Prerequisites & Installation

### System Requirements
*   **Python 3.10+**
*   **C Compiler** (Required for Numba JIT acceleration)
*   **uv** (Recommended for fast and reproducible environment management)

### Installation Options

#### Option A: Install from PyPI (Recommended for Users)
You can install `pymagical` directly using `pip` or `uv pip`:

```bash
# Using standard pip
pip install pymagical

# Using uv (much faster)
uv pip install pymagical
```

#### Option B: Clone for Development
If you want to run the benchmarks or contribute to the code:

```bash
git clone https://github.com/your-repo/pymagical.git
cd pymagical
uv sync
```

This will automatically create a virtual environment (`.venv`) with all required dependencies, including `pytest`, `numba`, `numpy`, and `pandas`.

---

## 2. Preparing Your Input Data

`pymagical` expects input files in tab-separated text formats (TSV). 

### Required Files and Formats

| File | Description | Format (Columns) | Header? |
| :--- | :--- | :--- | :--- |
| `sig_cr_genes.txt` | Candidate gene list | `gene_symbol` | No |
| `sig_cr_peaks.txt` | Candidate peak list | `chr`, `start`, `end` | No |
| `motif_info.txt` | TF-to-Motif mapping | `motif_index`, `tf_name` | No |
| `motif_prior.txt` | Motif-Peak binding prior | `peak_index`, `motif_index`, `flag` (binary) | No |
| `tad_regions.txt` | TAD boundaries | `chr`, `left_boundary`, `right_boundary` | No |
| `refseq.txt` | Genomic reference | `chr`, `strand`, `start`, `end`, `gene_name` | No |

### Cell-Type Specific Files (Required for each cell type folder)

| File | Description | Format (Columns) | Header? |
| :--- | :--- | :--- | :--- |
| `rna_counts.txt` | scRNA count matrix | `gene_index`, `cell_index`, `read_count` (COO format) | No |
| `rna_genes.txt` | scRNA gene metadata | `gene_index`, `gene_symbol` | No |
| `rna_meta.txt` | scRNA cell metadata | `cell_index`, `barcode`, `type`, `subject_ID`, `condition` | No |
| `atac_counts.txt` | scATAC count matrix | `peak_index`, `cell_index`, `read_count` (COO format) | No |
| `atac_peaks.txt` | scATAC peak metadata | `peak_index`, `chr`, `start`, `end` | No |
| `atac_meta.txt` | scATAC cell metadata | `cell_index`, `barcode`, `type`, `subject_ID`, `condition` | No |

> **Important:** All indices (`gene_index`, `cell_index`, `peak_index`) should be **1-indexed** (starting from 1) to remain compatible with standard MAGICAL data formats.

---

## 3. Running the Inference Pipeline

### Command Line Interface (CLI)
The most efficient way to run the pipeline is using the `pymagical` command (or `uv run pymagical` if using the development environment).

```bash
# Run with Numba acceleration (Recommended)
pymagical run \
    --main-dir ./data \
    --cell-dir astrocytes \
    --iter 2000 \
    --use-numba \
    --prefix my_experiment \
    --outdir outputs/
```

### Key Arguments:
*   `--iter`: Number of Gibbs sampling iterations (2000+ recommended for publication-quality results).
*   `--use-numba`: Enables JIT-compiled kernels. **Provides a ~30x speedup over MATLAB.**
*   `--prefix`: Prefix for the generated output files.

---

## 4. Understanding the Outputs

Once the run completes, you will find several files in your `outdir`:

1.  **`{prefix}_py_{iter}.txt`**: The primary results file. It contains the inferred circuits:
    *   `Peak_Gene_Prob`: The posterior probability of a peak-gene functional link.
    *   `TFs(prob, effect)`: A list of TFs binding that peak, their probability, and their biological direction (e.g., `+ [+,+]` for an activator).
2.  **`{prefix}_py_{iter}_B_matrix.txt`**: The full continuous weight matrix for TF-Peak binding.
3.  **`{prefix}_py_{iter}_L_matrix.txt`**: The full continuous weight matrix for Peak-Gene looping.

---

## 5. Verification & Testing

### Running Unit Tests
To ensure the mathematical kernels and data loaders are working correctly on your system, run the `pytest` suite:

```bash
# From the cloned repository
uv run pytest tests/
```

### Validating against MATLAB (Fidelity)
If you have results from the original MATLAB implementation, you can verify statistical fidelity using the comparison tool:

```bash
uv run python eval/tests/compare_results.py \
    --ml-dir path/to/matlab_out \
    --py-dir ./outputs \
    --iter 2000 \
    --ml-prefix astrocytes \
    --py-prefix astrocytes
```

---

## 6. Advanced Usage: Programmatic API

You can also integrate `pymagical` directly into your Python scripts:

```python
from pymagical import run_magical

run_magical(
    main_dir="data/",
    cell_dir="astrocytes",
    iteration_num=2000,
    use_numba=True,
    output_file="results/astrocytes_circuits.txt"
)
```
