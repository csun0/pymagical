# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`pymagical` is a Python port of **MAGICAL** (Multiome Accessibility Gene Integration Calling and Looping), a hierarchical Bayesian Gibbs sampler that infers regulatory circuits (TF → Peak → Gene triads) from paired single-cell RNA-seq and ATAC-seq data. The Python port is validated against the original MATLAB implementation, which lives in `src/matlabmagical/` and serves as the reference for correctness.

## Commands

```bash
uv sync                                  # create .venv with all deps (incl. numba, pytest)
uv run pytest tests/                     # run the unit test suite
uv run pytest tests/test_estimation.py   # run one test file
uv run pytest tests/test_estimation.py::test_name   # run one test

# Run inference (directory-based input; --use-numba for ~28x faster sampling)
uv run pymagical run --main-dir ./data --cell-dir astrocytes --iter 2000 --use-numba --prefix my_run --outdir outputs/

# `run` is the default subcommand; `pymagical --main-dir ... --iter ...` also works
# Generate interactive HTML report (needs `viz` extra: uv sync then pip install pymagical[viz])
uv run pymagical viz outputs/my_run_py_2000.txt

# Fidelity check vs MATLAB baseline
uv run python eval/tests/compare_results.py --ml-dir <ml_out> --py-dir <py_out> --iter 2000 --ml-prefix astrocytes --py-prefix astrocytes
```

Input files are 1-indexed TSVs (MATLAB origin). `--main-dir` holds shared files (motif_prior.txt, motif_info.txt, tad_regions.txt, refseq); `--cell-dir` holds cell-type-specific scRNA/scATAC files. Individual `--*` flags override resolved paths. See `TUTORIAL.md` for the full file schema.

## Pipeline architecture

`magical.py::run_magical` orchestrates four sequential stages, each in its own module:

1. **`data_loader.py`**: parses input TSVs. COO count matrices (`atac_counts.txt` is ~374MB) are expensive to parse, so first load caches into a hidden `.magical_cache/` dir next to the source (`.parquet` for DataFrames, `.npz` for sparse matrices). Cache validity is keyed on an MD5 fingerprint of source file path/size/mtime; a `completed` flag guards against partial writes. Subsequent loads drop from ~26s to ~1s.
2. **`circuits.py`**: `construct_candidate_circuits_with_tad` intersects candidate peaks/genes with the count matrices and TAD boundaries to build the TF/Peak/Gene binding (`B`) and looping (`L`) candidate state matrices. Only `tad_flag=1` is implemented.
3. **`initialization.py`**: OLS-based (`statsmodels`) seeding of priors (means/variances/probabilities) for the sampler.
4. **`estimation.py::magical_estimation`**: the Gibbs sampler. Returns final `B`/`L` probability matrices, mean weights, and per-edge sign-consistency probabilities.

`write_outputs` (in `magical.py`) thresholds `B`/`L` probabilities at 0.7 and emits the edge list `{prefix}_py_{iter}.txt`, plus full `_B_matrix.txt` / `_L_matrix.txt` and `_timing_stats.txt`.

## The Gibbs sampler: dual implementation (critical)

`magical_estimation` has **two parallel code paths that must produce statistically equivalent results**:

- **NumPy path** (default): readable per-step functions (`tf_activity_t_sampling`, `tf_peak_binding_b_sampling`, `tf_peak_binary_binding_b_state_sampling`, `peak_gene_looping_l_sampling`, `peak_gene_binary_looping_l_state_sampling`) directly in `estimation.py`.
- **Numba path** (`--use-numba`): JIT kernels in `estimation_kernels.py` (`@njit(parallel=True)`). Weights and states are transposed to row-contiguous layout (`b_T` is (M,P), `l_T` is (G,P)) and made `np.ascontiguousarray` before entering kernels, then transposed back each iteration for RSS/summary computation.

Both paths run the same 7-step loop per iteration: (1) TF activity `T`, (2) TF-peak binding weights `B`, (3) TF-peak binary states, (4) ATAC noise variance (inverse-gamma), (5) peak-gene looping weights `L`, (6) peak-gene binary states, (7) RNA noise variance.

When changing sampling logic, **change both paths and keep them in sync**. The Numba kernels are line-for-line ports of the NumPy functions (comments mark spots where they replicate specific MATLAB behavior). Divergence silently breaks fidelity.

### Fidelity constraints (do not "fix" without checking MATLAB)

Several details intentionally mirror MATLAB rather than being statistically "clean". Verify against `src/matlabmagical/` before changing:
- `np.var` uses `ddof=1` where MATLAB `var` does.
- State-sampling smoothing terms differ **per matrix**: B-state uses `+0.25`, L-state uses `+0.1` (matching `TF_peak_binary_binding_B_state_sampling.m` vs `Peak_gene_binary_looping_L_state_samping.m`). Both are followed by `+1e-6`.
- Inside the binary-state sweeps the running `temp` product is **not** updated after a flip; later variables in the sweep condition on the stale full product, exactly as MATLAB does. The commented-out updates in the Numba kernels are left as documentation of this.
- The T-activity posterior mean divides the residual dot-product by `P` only (no `T_var` factor), unlike the B/L means which carry `B_var/S` and `L_var/S`.
- Random normals are clipped to `[-3, 3]`.

Burn-in: the first `burn_in` iterations are discarded; summary accumulators start at zero and are divided by `iteration_num - burn_in`. `burn_in >= iteration_num` raises `ValueError`. (This departs from the original MATLAB, which discarded nothing and divided by `iteration_num` while seeding accumulators with the prior.)

Numba T-sampling maintains `t_sample` directly and omits MATLAB's per-cell `T_A`/`T_R` resample→re-average step; that step injects noise of order `variance_T / cells_per_sample`, negligible at single-cell scale. This is the one deliberate NumPy/Numba divergence.

Target: >0.99 Pearson correlation on `B`/`L` weight/probability matrices vs MATLAB at 2000+ iterations. `docs/statistical_fidelity.md` and `docs/decisions.md` record the rationale; check them before altering sampler math.

## Repo layout notes

- `src/pymagical/`: the package (src layout; version from git tags via hatch-vcs).
- `src/matlabmagical/`: original MATLAB reference implementation. Consult when porting or debugging fidelity.
- `eval/`: separate evaluation suite (fidelity comparison, Slurm/HPC benchmark runners, profiling). Distinct from `tests/`.
- `src/pymagical/ui/frontend/`: React/Vite/TS source for the interactive explorer; `viz.py` generates standalone HTML reports (Plotly + Jinja2, `viz` extra).
- CLI imports are deferred inside `cli.py` so `--help` and arg parsing stay fast and don't require numba/viz deps.

## Conventions

- 4-space indentation.
- Output artifacts follow `{prefix}_{py|ml}_{iter}` naming.
- Circuit direction is reported as `Overall_Effect [L_dir(consistency), B_dir(consistency)]`, where effect = sign of `B_weight * L_weight` and consistency is the fraction of iterations agreeing with the reported sign.
- Numba kernels: keep matrices C-contiguous (`np.ascontiguousarray`) before passing in; use running residuals with periodic resets rather than recomputing full products.
