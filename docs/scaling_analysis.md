# Computational Scaling and Feasibility Analysis

This document analyzes the computational complexity of the `pymagical` Gibbs sampler and explores the feasibility of running the pipeline without initial candidate filtering.

## 1. Algorithmic Complexity

The `pymagical` engine performs coordinate-wise Gibbs sampling with "Running Residuals" optimization. The three primary sampling stages scale as follows:

| Sampling Stage | Complexity per Iteration | Scaling Factors |
| :--- | :--- | :--- |
| **TF-Peak Binding ($B$)** | $O(M \times P \times S)$ | $M$: TFs, $P$: Peaks, $S$: Samples |
| **Peak-Gene Looping ($L$)** | $O(P \times G \times S)$ | $P$: Peaks, $G$: Genes, $S$: Samples |
| **TF Activity ($T$)** | $O(M \times S \times P)$ | $M$: TFs, $S$: Samples, $P$: Peaks |

The **Peak-Gene Looping ($L$)** stage is the primary bottleneck ($O(P \times G \times S)$).

---

## 2. Feasibility of "Unfiltered" Input

The standard workflow uses DAS/DEG filtering and TAD boundaries. The table below projects runtime feasibility for different input scenarios. Only the filtered column is measured; the other runtimes are order-of-magnitude projections scaled from the op counts, not benchmarks.

### Scenario Comparison (2,000 Iterations)

| Metric | Filtered (Standard) | **All Genes, Filtered Peaks** | Unfiltered (All Features) |
| :--- | :--- | :--- | :--- |
| **Peaks ($P$)** | ~500 | ~500 | ~100,000 |
| **Genes ($G$)** | ~400 | ~20,000 | ~20,000 |
| **Samples ($S$)** | 20 | 20 | 20 |
| **Ops per Iter** | $4 \times 10^6$ | $2 \times 10^8$ (50x) | $4 \times 10^{10}$ (10,000x) |
| **Numba Runtime** | **~1.5 min (measured)** | ~1 hour (est.) | ~10 days (est.) |
| **MATLAB Runtime** | ~35 min (measured) | ~1.5 days (est.) | ~300 days (est.) |

---

## 3. Implementation Notes (Python vs. MATLAB)

`pymagical` reproduces MATLAB's candidate-selection logic, including the strand-specific "extreme boundary" TSS rule: for a gene with multiple transcripts it takes the minimum `start` (`+` strand) or maximum `end` (`-` strand) across all transcripts (`circuits.py`). Initial peak and gene counts therefore closely match the MATLAB reference.

Any residual difference (typically within ~1% of genes) comes from tie-breaking when a gene name maps to multiple chromosomes or strands in RefSeq: `pymagical` dedupes on gene name and keeps the first entry. Because TAD filtering uses strict inequality (`TSS > left & TSS < right`), a few-base-pair shift in a borderline TSS can move a gene across a boundary, which in turn changes the set of peaks that retain a looping partner.

**Impact:** These edge cases affect a tiny fraction of genes and do not change statistical convergence or biological interpretation.

---

## 4. Practical Implications

*   **Feasible:** Skipping gene filtering (DEG) works if peaks (DAS) stay filtered (roughly an hour).
*   **Infeasible:** Skipping peak filtering is not recommended (multi-day runtime).
*   **Benefit:** Numba makes relaxed thresholds and high-throughput processing of many cell types practical.
