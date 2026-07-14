# Executive Summary: pymagical Performance Optimization

This document outlines the architectural improvements and performance gains achieved during the migration of the MAGICAL regulatory circuit inference framework from MATLAB to Python, and the subsequent hyper-optimization using Numba JIT-acceleration.

## 1. Migration: MATLAB to Python (NumPy)
The initial port to Python focused on logical fidelity and IO efficiency. By transitioning to a modern data science stack, we achieved a baseline speedup of **~1.9x** for the primary sampling engine.

### Key Drivers of Improvement:
*   **Intelligent IO Caching:** Replaced line-by-line text parsing with a dual-path IO handler. The first run caches data into PyArrow-backed Parquet and SciPy NPZ formats, reducing subsequent load times from 18 seconds to 1.2 seconds (**15x faster**).
*   **Vectorized Linear Algebra:** Leveraged highly optimized BLAS/LAPACK routines via NumPy and `statsmodels` for OLS-based parameter initialization, reducing the initialization stage from 27 seconds to under 3 seconds.

## 2. Hyper-Optimization: Python to Numba
To enable large-scale circuit discovery without aggressive gene/peak subsetting, we added an optional acceleration layer using **Numba JIT (Just-In-Time)** compilation. For the sampling stage this gave an additional **~15x** speedup over pure NumPy (2000 iterations: 1101s to 72.8s), bringing the total MATLAB-to-Numba sampling speedup to **~28x**.

### Architectural Optimizations:
*   **JIT-Compiled Kernels:** Heavily nested Gibbs sampling loops, previously limited by Python's interpreter overhead, were moved to specialized kernels in `estimation_kernels.py`. These kernels compile to machine code at runtime.
*   **Memory Layout Alignment:** Optimized the data flow by transposing key matrices (Weights, States, and Samples) to ensure **C-contiguous row-major access**. This maximized CPU cache hits and eliminated `NumbaPerformanceWarning` bottlenecks related to non-contiguous slicing.
*   **Running Residuals:** Refactored the sampler to maintain a "running residual" matrix. By updating only the affected components of the residual during coordinate-wise Gibbs steps, we eliminated redundant $O(M \times P \times S)$ matrix-vector multiplications.
*   **Parallelization:** Utilized `numba.prange` to parallelize the sampling of independent peak and gene states across multiple CPU cores.

## 3. Comparative Benchmarks (Astrocytes, 2000 Iterations)
Measured Gibbs sampling time:

*   **MATLAB (Original):** ~34 minutes (2047s)
*   **pymagical (NumPy):** ~18 minutes (1101s)
*   **pymagical (Numba):** **~1.2 minutes (72.8s)**

See [decisions.md](decisions.md) for the full run log (correlations, triad overlap, per-stage timing).

### Conclusion
Numba-accelerated Python turns MAGICAL sampling from a multi-hour HPC task into a sub-two-minute run (**~28x faster than MATLAB**), letting researchers work at higher resolution or across larger genomic regions within the same budget.
