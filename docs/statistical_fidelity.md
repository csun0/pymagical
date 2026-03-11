# Statistical Fidelity: Python vs. MATLAB Implementation

This document logs the specific mathematical and logic discrepancies identified and resolved to ensure `pymagical` provides results identical to the original MAGICAL publication.

## 1. Transcription Start Site (TSS) Calculation
*   **Issue:** Python initially selected the first transcript entry for each gene, leading to arbitrary TSS positions.
*   **Discrepancy:** This caused +/- 1% differences in the initial selection of candidate genes and peaks due to TAD boundary edge cases.
*   **Fix:** Adopted MATLAB's "Extreme Boundary" logic. We now group all transcripts for a gene and select the minimum `start` (for `+` strand) or maximum `end` (for `-` strand) to ensure biological coverage and deterministic results.

## 2. Binary State Sampling Likelihood
*   **Issue:** The likelihood term for the "Off" state ($state=0$) was missing the posterior mean penalty.
*   **Discrepancy:** The sampler was biased toward keeping circuits active, leading to higher circuit counts and lower Jaccard similarity.
*   **Fix:** Updated the kernels to match the MATLAB formula exactly:
    *   `post_0 = exp(-(0.0 - mean)^2 / (2 * var))`
    *   This correctly penalizes the zero state if the data strongly supports a non-zero weight.

## 3. Threshold Logic Inversion
*   **Issue:** The comparison against the random uniform variable was slightly offset and inverted in one case.
*   **Fix:** Synchronized the threshold checks to match MATLAB's exact branching:
    *   **Keep/Activate:** `if rand < P1`
    *   **Deactivate:** `if rand >= P1`

## 4. Iteration Counting & Initialization
*   **Issue:** MATLAB counts the initial state (post-initialization) as the first iteration ($i=1$), whereas Python started accumulating from the first sampling step ($i=0$ after the loop starts).
*   **Fix:** Initialized the frequency and weight accumulators with the initial state values before entering the loop, ensuring that a 2000-iteration run in Python represents the exact same cumulative evidence as MATLAB.

---

## Verification Results (Target)
By resolving these issues, we expect:
1.  **Parity in Initial Counts:** Identical number of TFs, Peaks, and Genes selected.
2.  **Increased Jaccard Similarity:** Higher overlap in discrete triad sets.
3.  **Near-Perfect Pearson Correlation:** Continuous probability matrices should match with $r > 0.999$.
