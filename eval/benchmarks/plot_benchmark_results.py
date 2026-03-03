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
    # Load data from the new organized directory
    data_dir = "eval/benchmarks/data"
    plot_dir = "eval/benchmarks/plots"
    os.makedirs(plot_dir, exist_ok=True)

    perf_df = pd.read_csv(os.path.join(data_dir, "benchmark_performance.csv"))
    conv_df = pd.read_csv(os.path.join(data_dir, "benchmark_convergence.csv"))
    fid_df = pd.read_csv(os.path.join(data_dir, "benchmark_fidelity.csv"))

    datasets = sorted([ds for ds in perf_df["Dataset"].unique() if ds != "endothelial"])
    impl_map = {"matlab": "MATLAB", "numpy": "NumPy", "numba": "Numba"}
    impl_palette = ['grey', 'cornflowerblue', 'dodgerblue']
    
    # --- Plot 1: Speed per Celltype ---
    fig1, axes1 = plt.subplots(3, 2, figsize=(4, 4.5), dpi=300, layout="tight")
    axes_flat1 = axes1.flatten()

    for i, ds in enumerate(datasets):
        ax = axes_flat1[i]
        df_ds = perf_df.loc[perf_df["Dataset"] == ds].copy()
        df_ds["Implementation"] = df_ds["Implementation"].map(impl_map)
        
        sns.lineplot(data=df_ds, x="Iterations", y="Gibbs_Time", hue="Implementation", 
                    marker="s", markersize=4, ax=ax, linewidth=1.5, palette=impl_palette)
        
        ax.set_yscale("log")
        ax.set_title(f'{ds.replace("_", " ").title()}', loc='left', fontproperties=hf.sf)
        ax.set_xlabel("Iterations", fontproperties=hf.sf)
        ax.set_ylabel("Time (s)", fontproperties=hf.sf)
        
        # Always remove the local legend
        ax.get_legend().remove()
        
        [l.set_fontproperties(hf.sf) for l in ax.get_yticklabels()]
        [l.set_fontproperties(hf.sf) for l in ax.get_xticklabels()]

    for j in range(len(datasets), len(axes_flat1)):
        axes_flat1[j].axis('off')

    # Extract handles and labels from the first active axis for the global legend
    handles, labels = axes_flat1[0].get_legend_handles_labels()

    # Create the figure-level legend
    fig1.legend(handles, labels, 
            loc='upper right', 
            bbox_to_anchor=(1,1), 
            ncol=len(labels),            
            prop=hf.sf,
            borderaxespad=0)

    fig1.suptitle("Optimized MAGICAL Runtimes", 
                fontproperties=hf.mf, x=0, y=1, ha='left', va='top')

    sns.despine()
    fig1.savefig(os.path.join(plot_dir, "speed_benchmark_all.png"))
    print(f"Saved {os.path.join(plot_dir, 'speed_benchmark_all.png')}")

    # --- Convergence Plots Helper ---
    iter_pair_map = {"100_vs_500": 500, "500_vs_1000": 1000, "1000_vs_2000": 2000, "2000_vs_5000": 5000}
    
    def plot_metric(df, col_name, x_col, title, filename, y_label, y_lim, palette, target_handles=None, target_labels=None):
        fig, axes = plt.subplots(3, 2, figsize=(4, 4.5), dpi=300, layout="tight")
        axes_flat = axes.flatten()

        for i, ds in enumerate(datasets):
            ax = axes_flat[i]
            df_ds = df[df["Dataset"] == ds].copy()
            if x_col == "Higher_Iter":
                df_ds["Higher_Iter"] = df_ds["Iter_Pair"].map(iter_pair_map)
            
            df_ds["Implementation"] = df_ds["Implementation"].map(impl_map)
            
            sns.lineplot(data=df_ds, x=x_col, y=col_name, hue="Implementation", 
                         marker="s", markersize=4, ax=ax, linewidth=1.5, palette=palette, legend=False)
            
            ax.set_title(f'{ds.replace("_", " ").title()}', loc='left', fontproperties=hf.sf)
            ax.set_xlabel("Iterations", fontproperties=hf.sf)
            ax.set_ylabel(y_label, fontproperties=hf.sf)
            if y_lim: ax.set_ylim(y_lim)
            
            [l.set_fontproperties(hf.sf) for l in ax.get_yticklabels()]
            [l.set_fontproperties(hf.sf) for l in ax.get_xticklabels()]

        for j in range(len(datasets), len(axes_flat)):
            axes_flat[j].axis('off')

        h = target_handles if target_handles else handles
        l = target_labels if target_labels else labels

        fig.legend(h, l, loc='upper right', bbox_to_anchor=(1,1), ncol=len(l), prop=hf.sf, borderaxespad=0)
        fig.suptitle(title, fontproperties=hf.mf, x=0, y=1, ha='left', va='top')
        sns.despine()
        fig.savefig(os.path.join(plot_dir, filename))
        print(f"Saved {os.path.join(plot_dir, filename)}")

    # 1. Within-Method Convergence
    plot_metric(conv_df, "Jaccard_Converge", "Higher_Iter", "Circuit Stability (N vs N-1)", 
                "convergence_jaccard_all.png", "Jaccard", (0.6, 1.05), impl_palette)
    plot_metric(conv_df, "B_Corr_Converge", "Higher_Iter", "Binding Matrix Stability (N vs N-1)", 
                "convergence_pearson_b_all.png", "Pearson r", (0.99, 1.001), impl_palette)
    plot_metric(conv_df, "L_Corr_Converge", "Higher_Iter", "Looping Matrix Stability (N vs N-1)", 
                "convergence_pearson_l_all.png", "Pearson r", (0.99, 1.001), impl_palette)

    # 2. Inter-Method Fidelity (vs MATLAB)
    fid_palette = ['cornflowerblue', 'dodgerblue']
    fid_handles, fid_labels = handles[1:], labels[1:]
    
    plot_metric(fid_df, "Jaccard_to_MATLAB", "Iterations", "Fidelity: Circuit Overlap vs MATLAB", 
                "fidelity_jaccard_all.png", "Jaccard", (0.6, 1.05), fid_palette, fid_handles, fid_labels)
    plot_metric(fid_df, "Recovery_to_MATLAB", "Iterations", "Fidelity: Circuit Recovery vs MATLAB", 
                "fidelity_recovery_all.png", "Recovery %", (0.6, 1.05), fid_palette, fid_handles, fid_labels)
    plot_metric(fid_df, "B_Corr_to_MATLAB", "Iterations", "Fidelity: Binding Matrix vs MATLAB", 
                "fidelity_pearson_b_all.png", "Pearson r", (0.985, 1.001), fid_palette, fid_handles, fid_labels)
    plot_metric(fid_df, "L_Corr_to_MATLAB", "Iterations", "Fidelity: Looping Matrix vs MATLAB", 
                "fidelity_pearson_l_all.png", "Pearson r", (0.985, 1.001), fid_palette, fid_handles, fid_labels)

if __name__ == "__main__":
    main()
