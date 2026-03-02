# Computational Scaling and Feasibility Analysis

This document analyzes the computational complexity of the `pymagical` Gibbs sampler and explores the feasibility of running the pipeline without initial candidate filtering (i.e., using all genes and all peaks).

## 1. Algorithmic Complexity

The `pymagical` engine performs coordinate-wise Gibbs sampling. Even with the "Running Residuals" optimization, which reduces the cost of updating a single weight from a full matrix-vector product to a vector addition, the total complexity per iteration remains significant as we sweep through the entire model.

The three primary sampling stages scale as follows:

| Sampling Stage | Complexity per Iteration | Scaling Factors |
| :--- | :--- | :--- |
| **TF–Peak Binding ($B$)** | $O(M \times P \times S)$ | $M$: TFs, $P$: Peaks, $S$: Samples |
| **Peak–Gene Looping ($L$)** | $O(P \times G \times S)$ | $P$: Peaks, $G$: Genes, $S$: Samples |
| **TF Activity ($T$)** | $O(M \times S \times P)$ | $M$: TFs, $S$: Samples, $P$: Peaks |

### The Bottleneck
The **Peak–Gene Looping ($L$)** stage is almost always the computational bottleneck because the number of potential interactions between every peak and every gene is far larger than the number of TFs. Its complexity is **$O(P \times G \times S)$**.

---

## 2. Feasibility of "Unfiltered" Input

The standard MAGICAL workflow uses differential analysis (DAS/DEG) and Topologically Associated Domain (TAD) boundaries to restrict the number of candidate peaks ($P$) and genes ($G$). We compare this to an "Unfiltered" scenario below.

### Numerical Comparison (2,000 Iterations)

| Metric | Filtered (Standard) | **All Genes, Filtered Peaks** | Unfiltered (All Features) |
| :--- | :--- | :--- | :--- |
| **Peaks ($P$)** | ~500 | ~500 | ~100,000 |
| **Genes ($G$)** | ~400 | ~20,000 | ~20,000 |
| **Samples ($S$)** | 20 | 20 | 20 |
| **Ops per Iter** | $4 \times 10^6$ | $2 \times 10^8$ (50x) | $4 \times 10^{10}$ (10,000x) |
| **Numba Runtime** | **~1.3 minutes** | **~65 minutes** | **~9 days** |
| **MATLAB Runtime** | ~36 minutes | ~30 hours | ~250 days |

---

## 3. Practical Implications

While the Numba-accelerated version is **~28x faster** than the original MATLAB implementation, it does not change the fundamental $O(N^2)$ scaling of the algorithm.

### Key Takeaway: The "Hybrid" Option is Feasible
The analysis shows that **skipping Gene filtering (DEG) is feasible** as long as you maintain Peak filtering (DAS). Running "All Genes" against a focused set of candidate peaks takes approximately **1 hour** per cell type. This allows researchers to discover regulatory links to genes that might not have passed strict differential expression cutoffs but are still being modulated by distal enhancers.

### Can I use all peaks?
**No.** Using all peaks (~100k) remains the primary driver of computational cost. Even if you filtered genes heavily, the $O(P \times G \times S)$ product would still be too large because $P$ is the largest dimension in genomic data.

### What the Speedup *Does* Allow:
1.  **Relaxed Thresholds:** You can safely relax your differential analysis thresholds (e.g., using an FDR of 0.1 instead of 0.01) to include more marginal candidate circuits without a significant runtime penalty.
2.  **High-Throughput Processing:** You can process 30+ cell types or sub-clusters in the time it previously took to process one.
3.  **Higher Iteration Counts:** You can run 5,000 or 10,000 iterations to ensure absolute convergence on complex datasets, which was previously too slow.

## Conclusion
The Numba implementation transforms MAGICAL from a "bottleneck" step into a "high-speed" step within a bioinformatics pipeline, but the initial biological filtering (selecting DAS and DEG within the same TAD) remains a critical requirement for genomic-scale feasibility.
