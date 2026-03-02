import pandas as pd
import numpy as np
import os
import re
import argparse
from scipy.stats import pearsonr

def parse_circuits(filepath):
    """Extract (Gene, Peak, TF) triads from a results text file."""
    circuits = set()
    if not os.path.exists(filepath):
        return circuits
    try:
        with open(filepath, 'r') as f:
            header = next(f)  # skip header
            for line in f:
                parts = line.strip('\n').split('\t')
                if len(parts) < 8:
                    continue
                
                gene = parts[0]
                peak = f"{parts[3]}_{parts[4]}_{parts[5]}"
                tfs_str = parts[7]
                
                # Extract TF names using regex to handle varying formats like 'TF (prob, sign)'
                tfs = re.findall(r'([A-Za-z0-9_-]+)\s*\(', tfs_str)
                
                for tf in tfs:
                    if tf:
                        circuits.add((gene, peak, tf))
    except Exception as e:
        print(f"Error parsing {filepath}: {e}")
                
    return circuits

def load_matrix(filepath):
    """Load L or B matrix as a DataFrame."""
    if not os.path.exists(filepath):
        return None
    try:
        return pd.read_csv(filepath, sep='\t', index_col=0)
    except Exception as e:
        print(f"Error loading matrix {filepath}: {e}")
        return None

def calculate_jaccard(set1, set2):
    """Calculate Jaccard index between two sets."""
    if not set1 and not set2:
        return 1.0
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union if union > 0 else 0.0

def calculate_matrix_correlation(df1, df2):
    """Calculate Pearson correlation between two matrices (common indices/columns)."""
    if df1 is None or df2 is None:
        return np.nan
    
    common_rows = np.intersect1d(df1.index, df2.index)
    common_cols = np.intersect1d(df1.columns, df2.columns)
    
    if len(common_rows) == 0 or len(common_cols) == 0:
        return np.nan
    
    v1 = df1.loc[common_rows, common_cols].values.flatten()
    v2 = df2.loc[common_rows, common_cols].values.flatten()
    
    if np.all(v1 == v1[0]) or np.all(v2 == v2[0]):
        return 1.0 if np.all(v1 == v2) else 0.0
    
    return pearsonr(v1, v2)[0]

def load_timing(filepath):
    """Parse timing_stats.txt file into a dictionary."""
    timing = {}
    if not os.path.exists(filepath):
        return None
    with open(filepath, 'r') as f:
        for line in f:
            if ':' in line:
                key, val = line.strip().split(':', 1)
                timing[key.strip()] = float(val.strip())
    return timing

def main():
    parser = argparse.ArgumentParser(description="Analyze magical benchmarks.")
    parser.add_argument("--bench-dir", type=str, default="outputs_bench", help="Root directory for benchmark outputs")
    parser.add_argument("--out-prefix", type=str, default="benchmark_", help="Output prefix for CSV files")
    args = parser.parse_args()

    datasets = ["astrocytes", "endothelial", "excitatory_neurons", "inhibitory_neurons", "microglia", "oligodendrocytes", "opcs"]
    impls = ["matlab", "numpy", "numba"]
    iters = [100, 500, 1000, 2000, 5000]

    results = []

    # 1. Collect all runs
    all_runs = {}
    for dataset in datasets:
        for impl in impls:
            for it in iters:
                label = f"{'ml' if impl == 'matlab' else 'py'}_{it}"
                prefix = f"{dataset}_{label}"
                path = os.path.join(args.bench_dir, impl, prefix)
                
                circuit_file = f"{path}.txt"
                b_matrix_file = f"{path}_B_matrix.txt"
                l_matrix_file = f"{path}_L_matrix.txt"
                timing_file = f"{path}_timing_stats.txt"
                
                if os.path.exists(circuit_file):
                    all_runs[(dataset, impl, it)] = {
                        "circuits": parse_circuits(circuit_file),
                        "b_matrix": load_matrix(b_matrix_file),
                        "l_matrix": load_matrix(l_matrix_file),
                        "timing": load_timing(timing_file)
                    }

    # 2. Performance Analysis
    perf_data = []
    for (ds, impl, it), data in all_runs.items():
        timing = data["timing"]
        if timing:
            gibbs_keys = [k for k in timing.keys() if "Gibbs" in k]
            if not gibbs_keys: continue
            gibbs_key = gibbs_keys[0]
            perf_data.append({
                "Dataset": ds,
                "Implementation": impl,
                "Iterations": it,
                "Total_Time": sum(timing.values()),
                "Gibbs_Time": timing[gibbs_key],
                "Time_Per_Iter": timing[gibbs_key] / it
            })
    
    if perf_data:
        df_perf = pd.DataFrame(perf_data)
        print("\n=== Performance Summary ===")
        print(df_perf.groupby(["Implementation", "Iterations"])[["Total_Time", "Gibbs_Time", "Time_Per_Iter"]].mean())
        df_perf.to_csv(f"{args.out_prefix}performance.csv", index=False)

    # 3. Fidelity Analysis (MATLAB as baseline)
    fidelity_data = []
    for ds in datasets:
        for it in iters:
            ml_run = all_runs.get((ds, "matlab", it))
            if not ml_run: continue
            
            for impl in ["numpy", "numba"]:
                py_run = all_runs.get((ds, impl, it))
                if not py_run: continue
                
                jaccard = calculate_jaccard(ml_run["circuits"], py_run["circuits"])
                recovery = len(ml_run["circuits"].intersection(py_run["circuits"])) / len(ml_run["circuits"]) if ml_run["circuits"] else 1.0
                b_corr = calculate_matrix_correlation(ml_run["b_matrix"], py_run["b_matrix"])
                l_corr = calculate_matrix_correlation(ml_run["l_matrix"], py_run["l_matrix"])
                
                fidelity_data.append({
                    "Dataset": ds,
                    "Implementation": impl,
                    "Iterations": it,
                    "Jaccard_to_MATLAB": jaccard,
                    "Recovery_to_MATLAB": recovery,
                    "B_Corr_to_MATLAB": b_corr,
                    "L_Corr_to_MATLAB": l_corr
                })
    
    if fidelity_data:
        df_fidelity = pd.DataFrame(fidelity_data)
        print("\n=== Fidelity Summary (vs MATLAB) ===")
        print(df_fidelity.groupby(["Implementation", "Iterations"])[["Jaccard_to_MATLAB", "Recovery_to_MATLAB", "B_Corr_to_MATLAB", "L_Corr_to_MATLAB"]].mean())
        df_fidelity.to_csv(f"{args.out_prefix}fidelity.csv", index=False)

    # 4. Convergence Analysis (Within-method)
    convergence_data = []
    for ds in datasets:
        for impl in impls:
            for i in range(len(iters) - 1):
                it1 = iters[i]
                it2 = iters[i+1]
                
                run1 = all_runs.get((ds, impl, it1))
                run2 = all_runs.get((ds, impl, it2))
                
                if run1 and run2:
                    jaccard = calculate_jaccard(run1["circuits"], run2["circuits"])
                    b_corr = calculate_matrix_correlation(run1["b_matrix"], run2["b_matrix"])
                    l_corr = calculate_matrix_correlation(run1["l_matrix"], run2["l_matrix"])
                    
                    convergence_data.append({
                        "Dataset": ds,
                        "Implementation": impl,
                        "Iter_Pair": f"{it1}_vs_{it2}",
                        "Jaccard_Converge": jaccard,
                        "B_Corr_Converge": b_corr,
                        "L_Corr_Converge": l_corr
                    })
    
    if convergence_data:
        df_conv = pd.DataFrame(convergence_data)
        print("\n=== Convergence Summary (Within-Method) ===")
        print(df_conv.groupby(["Implementation", "Iter_Pair"])[["Jaccard_Converge", "B_Corr_Converge", "L_Corr_Converge"]].mean())
        df_conv.to_csv(f"{args.out_prefix}convergence.csv", index=False)

    print(f"\nAnalysis saved with prefix: {args.out_prefix}")

if __name__ == "__main__":
    main()
