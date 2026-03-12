import pandas as pd
import re
import os
import numpy as np

def get_triads_with_probs(file_path, is_python=True):
    triads = {}
    if not os.path.exists(file_path):
        return triads
    with open(file_path, 'r') as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            gene = parts[0]
            # Standardize peak format
            peak = f"{parts[3]}:{parts[4]}-{parts[5]}"
            l_prob = float(parts[6])
            tfs_raw = parts[7]
            
            # TFs are comma separated
            tfs = [t.strip() for t in tfs_raw.split(',') if t.strip()]
            for t_item in tfs:
                # MATLAB: TF_NAME (prob)
                # Python: TF_NAME (prob, ...)
                match = re.match(r'([A-Za-z0-9_-]+)\s*\(([\d\.]+)', t_item)
                if match:
                    tf_name = match.group(1)
                    tp_prob = float(match.group(2))
                    triads[(gene, peak, tf_name)] = (l_prob, tp_prob)
    return triads

ml_file = 'eval/benchmarks/outputs/matlab/astrocytes_ml_5000.txt'
py_file = 'eval/benchmarks/outputs/convergence_10k/astrocytes_10k_rep1_py_10000.txt'

ml_data = get_triads_with_probs(ml_file, is_python=False)
py_data = get_triads_with_probs(py_file, is_python=True)

ml_set = set(ml_data.keys())
py_set = set(py_data.keys())

print(f"MATLAB triads: {len(ml_set)}")
print(f"Python triads: {len(py_set)}")

# Intersection
common = ml_set.intersection(py_set)
print(f"\nCommon triads: {len(common)}")
print(f"Recovery of MATLAB triads in Python: {len(common) / len(ml_set):.2%}")

# Exclusives
py_only = py_set - ml_set
ml_only = ml_set - py_set
print(f"Python exclusive triads: {len(py_only)}")
print(f"MATLAB exclusive triads: {len(ml_only)}")

# Prob correlation for common triads
if common:
    ml_l = [ml_data[k][0] for k in common]
    py_l = [py_data[k][0] for k in common]
    ml_tp = [ml_data[k][1] for k in common]
    py_tp = [py_data[k][1] for k in common]
    
    corr_l = np.corrcoef(ml_l, py_l)[0, 1]
    corr_tp = np.corrcoef(ml_tp, py_tp)[0, 1]
    
    print(f"\nLooping Prob Pearson Correlation: {corr_l:.4f}")
    print(f"Binding Prob Pearson Correlation: {corr_tp:.4f}")
