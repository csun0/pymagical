# pymagical: Getting Started Tutorial

`pymagical` is a high-performance Python implementation of the MAGICAL algorithm for inferring functional regulatory circuits (TF-Peak-Gene triads) from paired single-cell RNA-seq and ATAC-seq data.

This tutorial walks you through setting up the environment, preparing your data, running the inference engine, and validating the results.

---

## 1. Prerequisites & Installation

### System Requirements
*   **Python 3.10+**
*   **uv** (Recommended for extremely fast and reproducible environment management)
*   **C Compiler** (For Numba JIT acceleration)

### Installation
Clone the repository and sync the environment using `uv`:

```bash
git clone https://github.com/your-repo/pymagical.git
cd pymagical
uv sync
```

This will automatically create a virtual environment (`.venv`) with all required dependencies, including `pytest`, `numba`, `numpy`, and `pandas`.

---

## 2. Preparing Your Input Data

`pymagical` requires two types of inputs: cell-type-specific count data and genomic metadata.

### Data Layout
Organize your data into a central directory with subdirectories for each cell type:

```text
data/
├── sig_cr_genes.txt        # List of candidate genes
├── sig_cr_peaks.txt        # List of candidate peaks
├── motif_info.txt          # TF names for motifs
├── motif_prior.txt         # TF-Peak binding prior (COO format)
├── tad_regions.txt         # TAD boundary definitions
└── astrocytes/             # Cell-type specific folder
    ├── rna_counts.txt      # scRNA count matrix (COO)
    ├── rna_genes.txt       # scRNA gene metadata
    ├── rna_meta.txt        # scRNA cell metadata
    ├── atac_counts.txt     # scATAC count matrix (COO)
    ├── atac_peaks.txt      # scATAC peak metadata
    └── atac_meta.txt       # scATAC cell metadata
```

> **Note:** For the first run, `pymagical` will parse these text files and cache them as `.parquet` and `.npz` files in a hidden `.magical_cache` folder. Subsequent runs will load near-instantly.

---

## 3. Running the Inference Pipeline

### Command Line Interface (CLI)
The most efficient way to run the pipeline is using the `uv run` command with the `--use-numba` flag for a ~30x speedup.

```bash
uv run pymagical \
    --main-dir ./data \
    --cell-dir astrocytes \
    --iter 2000 \
    --use-numba \
    --prefix my_experiment \
    --outdir outputs/
```

### Key Arguments:
*   `--iter`: Number of Gibbs sampling iterations (2000+ recommended for publication-quality results).
*   `--use-numba`: Enables JIT-compiled kernels. **Always use this for production runs.**
*   `--prefix`: Prefix for the generated output files.

---

## 4. Understanding the Outputs

Once the run completes, you will find several files in your `outdir`:

1.  **`{prefix}_py_{iter}.txt`**: The primary results file. It contains the inferred circuits:
    *   `Peak_Gene_Prob`: The posterior probability of a peak-gene link.
    *   `TFs(prob, effect)`: A list of TFs binding that peak, their probability, and their biological direction (e.g., `+ [+,+]` for an activator).
2.  **`{prefix}_py_{iter}_B_matrix.txt`**: The full continuous weight matrix for TF-Peak binding.
3.  **`{prefix}_py_{iter}_L_matrix.txt`**: The full continuous weight matrix for Peak-Gene looping.

---

## 5. Verification & Testing

### Running Unit Tests
To ensure the mathematical kernels and data loaders are working correctly on your system, run the `pytest` suite:

```bash
uv run pytest tests/
```

### Validating against MATLAB (Fidelity)
If you have results from the original MATLAB implementation, you can verify statistical fidelity (Pearson correlation) using the comparison tool:

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
