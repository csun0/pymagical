# Computational Scaling and Feasibility Analysis

This document analyzes the computational complexity of the `pymagical` Gibbs sampler and explores the feasibility of running the pipeline without initial candidate filtering.

## 1. Algorithmic Complexity

The `pymagical` engine performs coordinate-wise Gibbs sampling with "Running Residuals" optimization. The three primary sampling stages scale as follows:

| Sampling Stage | Complexity per Iteration | Scaling Factors |
| :--- | :--- | :--- |
| **TF–Peak Binding ($B$)** | $O(M \times P \times S)$ | $M$: TFs, $P$: Peaks, $S$: Samples |
| **Peak–Gene Looping ($L$)** | $O(P \times G \times S)$ | $P$: Peaks, $G$: Genes, $S$: Samples |
| **TF Activity ($T$)** | $O(M \times S \times P)$ | $M$: TFs, $S$: Samples, $P$: Peaks |

The **Peak–Gene Looping ($L$)** stage is the primary bottleneck ($O(P \times G \times S)$).

---

## 2. Feasibility of "Unfiltered" Input

The standard workflow uses DAS/DEG filtering and TAD boundaries. Below is a comparison of runtime feasibility for different input scenarios.

### Numerical Comparison (2,000 Iterations)

| Metric | Filtered (Standard) | **All Genes, Filtered Peaks** | Unfiltered (All Features) |
| :--- | :--- | :--- | :--- |
| **Peaks ($P$)** | ~500 | ~500 | ~100,000 |
| **Genes ($G$)** | ~400 | ~20,000 | ~20,000 |
| **Samples ($S$)** | 20 | 20 | 20 |
| **Ops per Iter** | $4 \times 10^6$ | $2 \times 10^8$ (50x) | $4 \times 10^{10}$ (10,000x) |
| **Numba Runtime** | **~1.5 minutes** | **~75 minutes** | **~10 days** |
| **MATLAB Runtime** | ~47 minutes | ~39 hours | ~300 days |

---

## 3. Implementation Disparities (Python vs. MATLAB)

Users may notice slight deviations in the initial counts of selected peaks and genes (e.g., Python: 384 genes, MATLAB: 383 genes). These are **expected technical artifacts** caused by differences in selection logic:

1.  **TSS Calculation:** 
    *   **MATLAB:** Selects the extreme boundary for genes with multiple transcripts (minimum start for `+` strand, maximum end for `-` strand).
    *   **Python (`pymagical`):** Selects the TSS from the first matching entry in the RefSeq file.
2.  **TAD Boundary Sensitivity:** Because TAD filtering uses strict inequality (`TSS > left & TSS < right`), a shift of even a few base pairs in the calculated TSS can push a gene across a boundary, adding or removing it from the candidate set. 
3.  **Propagation:** A gene being added/removed from a TAD changes the set of peaks that have at least one looping partner, which in turn shifts the final peak count.

**Impact:** These differences only affect a tiny fraction of genes (~1%) and do not impact the overall statistical convergence or biological interpretation of the results.

---

## 4. Practical Implications

The Numba implementation transformations:
*   **Feasible:** Skipping Gene filtering (DEG) is possible if Peaks (DAS) are still filtered (~1 hour runtime).
*   **Infeasible:** Skipping Peak filtering is not recommended (~9 day runtime).
*   **Benefit:** Allows for relaxed statistical thresholds and high-throughput processing of dozens of cell types.
