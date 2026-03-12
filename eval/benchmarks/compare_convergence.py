import pandas as pd
import re
import os

def get_triads(file_path):
    triads = set()
    if not os.path.exists(file_path):
        return triads
    with open(file_path, 'r') as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            gene = parts[0]
            peak = f"{parts[3]}:{parts[4]}-{parts[5]}"
            tfs_raw = parts[7]
            # Extract TF name: TF_NAME (prob, ...)
            # TFs are comma separated
            tfs = [t.strip() for t in tfs_raw.split(',') if t.strip()]
            for t_item in tfs:
                match = re.match(r'([A-Za-z0-9_-]+)', t_item)
                if match:
                    tf_name = match.group(1)
                    triads.add((gene, peak, tf_name))
    return triads

res_dir = 'eval/benchmarks/outputs/convergence_10k'
reps = [1, 2, 3]
all_triads = []
for r in reps:
    f = os.path.join(res_dir, f"astrocytes_10k_rep{r}_py_10000.txt")
    t = get_triads(f)
    all_triads.append(t)
    print(f"Rep {r}: {len(t)} triads found.")

# Intersection
common = all_triads[0].intersection(all_triads[1]).intersection(all_triads[2])
print(f"\nCommon triads in all 3 runs: {len(common)}")

# Union
total_unique = all_triads[0].union(all_triads[1]).union(all_triads[2])
print(f"Total unique triads across all runs: {len(total_unique)}")

# Percent recovery (mean Jaccard-like)
if len(total_unique) > 0:
    print(f"Overlap consistency: {len(common) / len(total_unique):.2%}")

# pairwise Jaccard
for i in range(3):
    for j in range(i+1, 3):
        intersection = all_triads[i].intersection(all_triads[j])
        union = all_triads[i].union(all_triads[j])
        print(f"Rep {i+1} vs Rep {j+1} Jaccard: {len(intersection)/len(union):.2%}")
