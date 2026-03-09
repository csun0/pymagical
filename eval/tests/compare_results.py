import pandas as pd
import numpy as np
import re
import os
import argparse

def parse_circuits(filepath):
    """Extract (Gene, Peak, TF) triads from a results text file."""
    circuits = set()
    if not os.path.exists(filepath):
        print(f"Warning: {filepath} not found.")
        return circuits
    with open(filepath, 'r') as f:
        next(f)  # skip header
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
                
    return circuits

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

def compare_fidelity(ml_file, py_file, ml_l_mat, py_l_mat, ml_b_mat, py_b_mat):
    print("\n=== Fidelity Comparison ===")
    
    # 1. Matrix Correlation
    print("\n1. Continuous Matrix Correlations:")
    try:
        ml_l = pd.read_csv(ml_l_mat, sep='\t', index_col=0)
        py_l = pd.read_csv(py_l_mat, sep='\t', index_col=0)
        ml_b = pd.read_csv(ml_b_mat, sep='\t', index_col=0)
        py_b = pd.read_csv(py_b_mat, sep='\t', index_col=0)
        
        common_peaks = np.intersect1d(ml_l.index, py_l.index)
        common_genes = np.intersect1d(ml_l.columns, py_l.columns)
        common_tfs = np.intersect1d(ml_b.columns, py_b.columns)
        
        l_val_ml = ml_l.loc[common_peaks, common_genes].values.flatten()
        l_val_py = py_l.loc[common_peaks, common_genes].values.flatten()
        b_val_ml = ml_b.loc[common_peaks, common_tfs].values.flatten()
        b_val_py = py_b.loc[common_peaks, common_tfs].values.flatten()

        l_corr = np.corrcoef(l_val_ml, l_val_py)[0,1]
        b_corr = np.corrcoef(b_val_ml, b_val_py)[0,1]

        # Calculate non-zero correlation (union of non-zero entries)
        l_nonzero = (l_val_ml != 0) | (l_val_py != 0)
        b_nonzero = (b_val_ml != 0) | (b_val_py != 0)
        
        l_corr_nz = np.corrcoef(l_val_ml[l_nonzero], l_val_py[l_nonzero])[0,1] if np.any(l_nonzero) else 0
        b_corr_nz = np.corrcoef(b_val_ml[b_nonzero], b_val_py[b_nonzero])[0,1] if np.any(b_nonzero) else 0

        # Density reporting
        l_dens_py = np.count_nonzero(l_val_py) / len(l_val_py) * 100
        l_dens_ml = np.count_nonzero(l_val_ml) / len(l_val_ml) * 100
        b_dens_py = np.count_nonzero(b_val_py) / len(b_val_py) * 100
        b_dens_ml = np.count_nonzero(b_val_ml) / len(b_val_ml) * 100

        print(f"  L matrix (Peak-Gene) Pearson correlation: {l_corr:.4f} (non-zero only: {l_corr_nz:.4f})")
        print(f"    - Density: Python {l_dens_py:.2f}%, MATLAB {l_dens_ml:.2f}%")
        print(f"  B matrix (TF-Peak) Pearson correlation:   {b_corr:.4f} (non-zero only: {b_corr_nz:.4f})")
        print(f"    - Density: Python {b_dens_py:.2f}%, MATLAB {b_dens_ml:.2f}%")
    except Exception as e:
        print(f"  Error comparing matrices: {e}")
    
    # 2. Triad Overlap
    print("\n2. Inferred Triad (Gene-Peak-TF) Overlap:")
    ml_circuits = parse_circuits(ml_file)
    py_circuits = parse_circuits(py_file)
    
    def get_overlap_metrics(ml_set, py_set, label):
        intersection = ml_set.intersection(py_set)
        union = ml_set.union(py_set)
        jaccard = len(intersection) / len(union) if union else 0
        recovery = len(intersection) / len(ml_set) if ml_set else 0
        
        print(f"\n  {label} Overlap:")
        print(f"    MATLAB total: {len(ml_set)}")
        print(f"    Python total: {len(py_set)}")
        print(f"    Overlapping:  {len(intersection)}")
        print(f"    Jaccard Sim:  {jaccard:.4f}")
        print(f"    Recovery:     {recovery*100:.1f}% of MATLAB {label.lower()} found in Python")

    # Triad Metrics
    get_overlap_metrics(ml_circuits, py_circuits, "Triad")

    # Gene Metrics
    ml_genes = {c[0] for c in ml_circuits}
    py_genes = {c[0] for c in py_circuits}
    get_overlap_metrics(ml_genes, py_genes, "Gene")

    # TF Metrics
    ml_tfs = {c[2] for c in ml_circuits}
    py_tfs = {c[2] for c in py_circuits}
    get_overlap_metrics(ml_tfs, py_tfs, "TF")

def compare_performance(ml_timing_file, py_timing_file):
    print("\n=== Performance Comparison ===")
    ml_timing = load_timing(ml_timing_file)
    py_timing = load_timing(py_timing_file)
    
    if not ml_timing or not py_timing:
        print("  Timing data not available for both implementations.")
        if not ml_timing: print(f"  Missing MATLAB timing: {ml_timing_file}")
        if not py_timing: print(f"  Missing Python timing: {py_timing_file}")
        return
        
    print(f"{'Stage':<30} | {'MATLAB (s)':<10} | {'Python (s)':<10} | {'Speedup'}")
    print("-" * 75)
    
    # Common keys
    stages = ["Data Loading", "Circuit Construction", "Initialization (OLS)"]
    # Handle the variable iter count key
    ml_gibbs_key = [k for k in ml_timing.keys() if "Gibbs Sampling" in k]
    py_gibbs_key = [k for k in py_timing.keys() if "Gibbs Sampling" in k]
    
    if ml_gibbs_key and py_gibbs_key:
        stages.append((ml_gibbs_key[0], py_gibbs_key[0]))
    
    for stage in stages:
        if isinstance(stage, tuple):
            k_ml, k_py = stage
            display_name = "Gibbs Sampling"
        else:
            k_ml = k_py = stage
            display_name = stage
            
        t_ml = ml_timing.get(k_ml, 0)
        t_py = py_timing.get(k_py, 0)
        speedup = t_ml / t_py if t_py > 0 else 0
        print(f"{display_name:<30} | {t_ml:<10.3f} | {t_py:<10.3f} | {speedup:.2f}x")

def main():
    parser = argparse.ArgumentParser(description="Compare Python vs MATLAB results for statistical fidelity and performance.")
    parser.add_argument("--ml-dir", type=str, required=True, help="Directory containing MATLAB outputs")
    parser.add_argument("--py-dir", type=str, required=True, help="Directory containing Python outputs")
    parser.add_argument("--iter", type=int, default=500, help="Number of iterations used in the run")
    parser.add_argument("--ml-prefix", type=str, default="astrocytes", help="Prefix for MATLAB filenames")
    parser.add_argument("--py-prefix", type=str, default="astrocytes", help="Prefix for Python filenames")
    
    args = parser.parse_args()
    
    ml_base = os.path.join(args.ml_dir, f"{args.ml_prefix}_ml_{args.iter}")
    py_base = os.path.join(args.py_dir, f"{args.py_prefix}_py_{args.iter}")
    
    # Check if we need to search deeper (for outputs_bench structure)
    if not os.path.exists(f"{py_base}.txt"):
        # Try looking for prefix inside subdirectory if it exists
        if os.path.exists(os.path.join(args.py_dir, f"{args.py_prefix}_py_{args.iter}.txt")):
             py_base = os.path.join(args.py_dir, f"{args.py_prefix}_py_{args.iter}")
    
    compare_fidelity(
        ml_file=f"{ml_base}.txt",
        py_file=f"{py_base}.txt",
        ml_l_mat=f"{ml_base}_L_matrix.txt",
        py_l_mat=f"{py_base}_L_matrix.txt",
        ml_b_mat=f"{ml_base}_B_matrix.txt",
        py_b_mat=f"{py_base}_B_matrix.txt"
    )
    
    compare_performance(
        ml_timing_file=f"{ml_base}_timing_stats.txt",
        py_timing_file=f"{py_base}_timing_stats.txt"
    )

if __name__ == "__main__":
    main()
