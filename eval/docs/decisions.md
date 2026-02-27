# pymagical Project Overview and Decisions Notebook

## PHASE 0: IO Optimization Assessment

### Current Input File Structure
The test input datasets consist of space/tab-delimited text files. The largest files are `atac_counts.txt` (374MB) and `rna_counts.txt` (114MB), which encode sparse matrices in Coordinate format (COO: row, column, value, 1-indexed since it originates from MATLAB). Other files, like `atac_peaks.txt` and `rna_genes.txt`, contain standard row-based metadata.

### Observations & Recommendations
1. **Clarity & Speed:** Parsing hundreds of megabytes of text-based coordinate lists line-by-line is memory and compute-intensive in Python. 
2. **PyArrow / Parquet Storage:** Switching to PyArrow-backed Parquet files would provide tremendous benefits. Parquet's columnar format easily shrinks 300MB+ `.txt` files down to ~20-50MB. Reading Parquet files is an order of magnitude faster because the parser avoids decoding strings to integers and directly loads memory-mapped data.

### Decision for IO
For `pymagical`, we implemented a dual-path IO handler (`data_loader.py`):
- Accepts standard text formats for user-friendliness and compatibility.
- Upon reading text files for the first time, caches the parsed data into an optimized format (`.parquet` for dataframes, SciPy `.npz` for sparse matrices) within a hidden `./.magical_cache` directory, utilizing MD5 hashes of file modification times to validate cache state.
- **Results:** First load (caching) took ~26 seconds total. Subsequent loads took ~1.14 seconds total.

## PHASE 1: LOGICAL PORT
- Replicated MATLAB logic in pure Python using `numpy`, `scipy.sparse`, `pandas`, and `statsmodels`.
- Validated output format and data shapes: successfully generated the requested edge list `.txt`, `B_matrix.txt`, and `L_matrix.txt`.
- Executed the port locally for 10 iterations. Generated correlation metrics against MATLAB 10-iteration baseline. Result showed ~0.95 correlation on both B and L probability matrices.

## PHASE 3: OPTIMIZATION
- Generated a script (`profile_run.py` & `run_matlab_profile.m`) to calculate runtime statistics for both versions to compare their speeds.
- Using standard `numpy` implementations, the translated Python port operates significantly faster than the MATLAB baseline across almost every stage. 
- **10-Iteration Profiling Results:**
  - **Data Loading:** Python (Cache) 1.2s vs MATLAB 18.2s (15x faster)
  - **Circuit Construction:** Python 3.3s vs MATLAB 1.4s (MATLAB is slightly faster here due to highly optimized native set-intersection operations)
  - **Initialization (OLS):** Python 2.4s vs MATLAB 27.5s (11x faster using `statsmodels`)
  - **Gibbs Sampling:** Python 7.5s vs MATLAB 12.9s (1.7x faster)
- The overall Python execution time is massively reduced, scaling exceptionally well for high-iteration HPC cluster compute tasks.
- A visual breakdown chart was outputted to `runtime_comparison.png`.

## PHASE 4: LARGE-SCALE VALIDATION
- Submitted identical 500-iteration Slurm jobs for both the translated `pymagical` Python port and the original MATLAB codebase.
- High-memory configurations (`--mem=512G`) were supplied.
- **Results (500 iterations):** The 500-iteration convergence metrics showed near-perfect statistical parity between the original algorithm and the Python port:
  - **L matrix (Peak-Gene) Pearson correlation:** 0.9987
  - **B matrix (TF-Peak) Pearson correlation:** 0.9945
  - The mean probabilities computed by both frameworks were functionally identical.
- **Runtime Performance (500 iterations):** The Python Slurm job completed in **4 minutes and 52 seconds**, while the MATLAB job took **6 minutes and 43 seconds**.

- **Results (1000 iterations):**
  - **L matrix Pearson correlation:** 0.9990
  - **B matrix Pearson correlation:** 0.9950
  - **Triad Overlap:** Python recovered **95.0%** of the triads found by MATLAB (550 out of 579).
  - **Jaccard Similarity:** 0.7152
  - **Runtime Performance:** Python completed in **6 minutes and 51 seconds**, while MATLAB took **12 minutes and 13 seconds** (nearly **1.8x faster**).

- **Results (2000 iterations):**
  - **L matrix Pearson correlation:** 0.9993
  - **B matrix Pearson correlation:** 0.9952
  - **Triad Overlap:** Python recovered **96.6%** of the triads found by MATLAB (562 out of 582).
  - **Jaccard Similarity:** 0.7159
  - **Runtime Performance:** Python completed in **18 minutes and 21 seconds**, while MATLAB took **35 minutes and 11 seconds** (**~1.9x faster**).

## PHASE 5: FEATURE ENHANCEMENTS
- **Biological Directionality:** The original implementation only provided interaction probabilities. The Python port was modified to calculate whether a TF acts as an **activator (+)** or a **repressor (-)** for its target gene by tracking the continuous linear regression weights generated during the Gibbs sampling loops. 
- The final output string format was updated to include this context, breaking down the interaction into its component parts (e.g., `STAT5B (0.85, + [+,+])` indicating an activator that opens a peak `+` which increases gene expression `+`).
- **Weight History Analysis:** Added a `--dump-weights` CLI flag. When enabled, the sampler outputs the complete, iteration-by-iteration history of continuous B and L weights as `.npy` matrices, allowing for downstream distribution analysis to assess the robustness of taking the mean weight across runs.

