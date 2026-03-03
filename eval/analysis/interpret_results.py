import os
import pandas as pd
import numpy as np
import re
from pathlib import Path

def parse_triads(filepath):
    data = []
    if not os.path.exists(filepath): return pd.DataFrame()
    with open(filepath, 'r') as f:
        next(f)
        for line in f:
            p = line.strip().split('	')
            if len(p) < 8: continue
            gene, peak = p[0], f"{p[3]}_{p[4]}_{p[5]}"
            prob_pg = float(p[6])
            tfs_raw = p[7].strip().rstrip(',')
            matches = re.findall(r'([A-Za-z0-9_-]+)\s*\(([\d\.]+),\s*([\+\-])\s*\[([\+\-]),([\+\-])\]\)', tfs_raw)
            for tf, tf_prob, effect, l_dir, b_dir in matches:
                data.append({
                    'Gene': gene, 'Peak': peak, 'TF': tf,
                    'Combined_Prob': prob_pg * float(tf_prob),
                    'Effect': 1 if effect == '+' else -1,
                    'L_Dir': l_dir, 'B_Dir': b_dir
                })
    return pd.DataFrame(data)

def main():
    results_dir = "outputs_bench/numba"
    datasets = ["astrocytes", "excitatory_neurons", "inhibitory_neurons", "microglia", "oligodendrocytes", "opcs"]
    
    for ds in datasets:
        print(f"\nProcessing {ds.upper()}...")
        df = parse_triads(os.path.join(results_dir, f"{ds}_py_2000.txt"))
        if df.empty: continue

        # 1 & 2: Load and Mode
        tf_stats = df.groupby('TF').agg({
            'Gene': 'nunique',
            'Combined_Prob': 'mean',
            'Effect': lambda x: (x == 1).sum() / len(x)
        }).rename(columns={'Gene': 'Target_Genes', 'Effect': 'Activator_Ratio'})
        
        top_tfs = tf_stats.sort_values('Target_Genes', ascending=False).head(5)
        print("Top 5 Master Regulators:")
        print(top_tfs)

        # 4: GSEA Ranks
        rank_dir = Path("eval/analysis/gsea_ranks") / ds
        rank_dir.mkdir(parents=True, exist_ok=True)
        for tf in top_tfs.index:
            tf_df = df[df['TF'] == tf]
            # Regulatory Potency Score: sum(Prob * Effect) per gene
            gene_ranks = tf_df.groupby('Gene').apply(lambda x: (x['Combined_Prob'] * x['Effect']).sum())
            gene_ranks.sort_values(ascending=False).to_csv(rank_dir / f"{tf}.rnk", sep='	', header=False)

        # 5: Prep for validation (Top 3 circuits)
        top_circuits = df.sort_values('Combined_Prob', ascending=False).head(3)
        print("Top 3 Circuits for Validation:")
        for _, row in top_circuits.iterrows():
            print(f"  {row['TF']} -> {row['Gene']} (Prob: {row['Combined_Prob']:.3f}, Effect: {row['Effect']})")

if __name__ == "__main__":
    main()
