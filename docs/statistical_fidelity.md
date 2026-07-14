# Statistical Fidelity: Python vs. MATLAB Implementation

This document logs the specific mathematical and logic discrepancies identified and resolved to ensure `pymagical` provides results statistically equivalent to the original MAGICAL implementation (see the correlation targets below; some behavior deliberately diverges, as noted in section 4).

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

## 4. Iteration Counting, Initialization & Burn-in
*   **History:** An earlier port seeded the frequency/weight accumulators with the
    post-initialization state before the loop, to mirror MATLAB counting that state
    as iteration $i=1$. This is **no longer the case**; see the deliberate divergence below.
*   **Current behavior (deliberate divergence from MATLAB):** `pymagical` discards the
    first `burn_in` iterations (default `iteration_num // 5`). Summary accumulators start
    at **zero** (not the OLS/prior seed) and are divided by `iteration_num - burn_in`.
    This removes the starting-state bias that seeding introduced. `burn_in >= iteration_num`
    raises `ValueError`.
*   **Contrast with MATLAB:** the reference implementation discards nothing, seeds the
    accumulators with the prior, and divides by `iteration_num`. The two are therefore
    *not* expected to be bit-identical on the initial iterations; parity is measured on
    the post-burn-in converged matrices (see targets below). Rationale is recorded in
    `docs/decisions.md` and `CLAUDE.md`.

---

## Verification Results (Target)
By resolving these issues, we expect:
1.  **Near-parity in Initial Counts:** Matching counts of TFs, Peaks, and Genes, up to rare tie-breaking edge cases (see `scaling_analysis.md`).
2.  **Increased Jaccard Similarity:** Higher overlap in discrete triad sets.
3.  **Near-Perfect Pearson Correlation:** Continuous probability matrices match with $r > 0.999$ (L) and $r > 0.995$ (B) at 2000 iterations.
