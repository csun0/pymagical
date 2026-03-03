import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys

# Add helpers path
sys.path.append("/mnt/home/csun1/scripts/global_scripts/imports")

try:
    from helpers.init import *
    import helpers.font as hf
    hf.set_font_family("Google Sans")
except ImportError:
    # Fallback if helpers not found during execution (unlikely given find result)
    print("Warning: helpers not found, using standard matplotlib")
    class DummyFont:
        def __getattr__(self, name): return None
    hf = DummyFont()

def main():
    # Load data
    perf_df = pd.read_csv("benchmark_performance.csv")
    conv_df = pd.read_csv("benchmark_convergence.csv")

    datasets = sorted([ds for ds in perf_df["Dataset"].unique() if ds != "endothelial"])
    impl_map = {"matlab": "MATLAB", "numpy": "NumPy", "numba": "Numba"}
    
    # --- Plot 1: Speed per Celltype ---
    fig, axes = plt.subplots(3, 3, figsize=(12, 10), dpi=300, layout="constrained")
    axes_flat = axes.flatten()
    
    for i, ds in enumerate(datasets):
        ax = axes_flat[i]
        df_ds = perf_df[perf_df["Dataset"] == ds].copy()
        df_ds["Implementation"] = df_ds["Implementation"].map(impl_map)
        
        sns.lineplot(data=df_ds, x="Iterations", y="Gibbs_Time", hue="Implementation", 
                     marker="o", ax=ax, linewidth=1.5)
        
        ax.set_yscale("log")
        ax.set_title(f'{ds.replace("_", " ").title()}', loc='left', fontproperties=hf.lf if hasattr(hf, "lf") else None)
        ax.set_xlabel("Iterations", fontproperties=hf.sf if hasattr(hf, "sf") else None)
        ax.set_ylabel("Gibbs Time (s)", fontproperties=hf.sf if hasattr(hf, "sf") else None)
        
        if i != 0: ax.get_legend().remove()
        
        if hasattr(hf, "sf"):
            [l.set_fontproperties(hf.sf) for l in ax.get_yticklabels()]
            [l.set_fontproperties(hf.sf) for l in ax.get_xticklabels()]

    # Hide empty subplots
    for j in range(len(datasets), len(axes_flat)):
        axes_flat[j].axis('off')

    plt.suptitle("Gibbs Sampling Speed Scaling across Celltypes", 
                 fontproperties=hf.mf if hasattr(hf, "mf") else None, x=0.05, ha='left')
    sns.despine()
    plt.savefig("speed_benchmark_all.png")
    print("Saved speed_benchmark_all.png")

    # --- Plot 2: Convergence per Celltype ---
    iter_pair_map = {"100_vs_500": 500, "500_vs_1000": 1000, "1000_vs_2000": 2000, "2000_vs_5000": 5000}
    fig, axes = plt.subplots(3, 3, figsize=(12, 10), dpi=300, layout="constrained")
    axes_flat = axes.flatten()

    for i, ds in enumerate(datasets):
        ax = axes_flat[i]
        df_ds = conv_df[conv_df["Dataset"] == ds].copy()
        df_ds["Higher_Iter"] = df_ds["Iter_Pair"].map(iter_pair_map)
        df_ds["Implementation"] = df_ds["Implementation"].map(impl_map)
        
        sns.lineplot(data=df_ds, x="Higher_Iter", y="Jaccard_Converge", hue="Implementation", 
                     marker="s", linestyle="--", ax=ax)
        
        ax.set_title(f'{ds.replace("_", " ").title()}', loc='left', fontproperties=hf.lf if hasattr(hf, "lf") else None)
        ax.set_xlabel("Iterations", fontproperties=hf.sf if hasattr(hf, "sf") else None)
        ax.set_ylabel("Jaccard Similarity", fontproperties=hf.sf if hasattr(hf, "sf") else None)
        ax.set_ylim(0.7, 1.05)
        
        if i != 0: ax.get_legend().remove()
        
        if hasattr(hf, "sf"):
            [l.set_fontproperties(hf.sf) for l in ax.get_yticklabels()]
            [l.set_fontproperties(hf.sf) for l in ax.get_xticklabels()]

    for j in range(len(datasets), len(axes_flat)):
        axes_flat[j].axis('off')

    plt.suptitle("Circuit Convergence (Jaccard) across Celltypes", 
                 fontproperties=hf.mf if hasattr(hf, "mf") else None, x=0.05, ha='left')
    sns.despine()
    plt.savefig("convergence_benchmark_all.png")
    print("Saved convergence_benchmark_all.png")

if __name__ == "__main__":
    main()
