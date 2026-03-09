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
    # Fallback if helpers not found during execution
    print("Warning: helpers not found, using standard matplotlib")
    class DummyFont:
        def __getattr__(self, name): return None
    hf = DummyFont()

def main():
    # Load data from the organized directory
    data_dir = "eval/benchmarks/data"
    plot_dir = "eval/benchmarks/plots"
    os.makedirs(plot_dir, exist_ok=True)

    perf_df = pd.read_csv(os.path.join(data_dir, "benchmark_performance.csv"))
    conv_df = pd.read_csv(os.path.join(data_dir, "benchmark_convergence.csv"))
    fid_df = pd.read_csv(os.path.join(data_dir, "benchmark_fidelity.csv"))

    # Define color palette for datasets (cell types)
    husl_12 = sns.color_palette("husl", 12)
    color_dict = {
        'Excitatory Neurons': husl_12[-4],
        'Inhibitory Neurons': husl_12[0],
        'Astrocytes': husl_12[4],
        'Microglia': husl_12[1],
        'Oligodendrocytes': husl_12[-3],
        'OPCs': husl_12[-2],
        'Endothelial': "goldenrod",
        'Unassigned': 'gray'
    }

    # Filter for Numba only (numpy and numba have same RNG, so we only need one)
    # Exclude endothelial dataset as requested
    fid_numba = fid_df[(fid_df["Implementation"] == "numba") & (fid_df["Dataset"] != "endothelial")].copy()
    
    # Precise mapping to match the provided color_dict keys exactly
    dataset_to_label = {
        "astrocytes": "Astrocytes",
        "excitatory_neurons": "Excitatory Neurons",
        "inhibitory_neurons": "Inhibitory Neurons",
        "microglia": "Microglia",
        "oligodendrocytes": "Oligodendrocytes",
        "opcs": "OPCs"
    }
    fid_numba["Dataset_Label"] = fid_numba["Dataset"].map(dataset_to_label)
    
    def plot_fidelity_combined(df, col_name, title, filename, y_label):
        fig, ax = plt.subplots(1, 1, figsize=(4, 3), dpi=300, layout="constrained")
        
        # Use lineplot with manual color mapping
        sns.lineplot(data=df, x="Iterations", y=col_name, hue="Dataset_Label", 
                     marker="o", markersize=4, ax=ax, linewidth=1.2, palette=color_dict)
        
        ax.set_title(title, loc='left', fontproperties=hf.mf, pad=15)
        ax.set_xlabel("Iterations", fontproperties=hf.sf)
        ax.set_ylabel(y_label, fontproperties=hf.sf)
        
        [l.set_fontproperties(hf.sf) for l in ax.get_yticklabels()]
        [l.set_fontproperties(hf.sf) for l in ax.get_xticklabels()]
        
        # Legend styling
        ax.legend(prop=hf.sf, bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
        
        sns.despine()
        fig.savefig(os.path.join(plot_dir, filename), bbox_inches='tight')
        print(f"Saved {os.path.join(plot_dir, filename)}")
        plt.close(fig)

    # Fidelity: Numba vs MATLAB
    plot_fidelity_combined(fid_numba, "Jaccard_to_MATLAB", "Fidelity: Circuit Overlap (Jaccard)\nNumba vs MATLAB", 
                           "fidelity_jaccard_combined.png", "Jaccard")
    
    plot_fidelity_combined(fid_numba, "Recovery_to_MATLAB", "Fidelity: Circuit Recovery\nNumba vs MATLAB", 
                           "fidelity_recovery_combined.png", "Recovery %")
    
    plot_fidelity_combined(fid_numba, "B_Corr_to_MATLAB", "Fidelity: Binding Matrix (B)\nNumba vs MATLAB", 
                           "fidelity_pearson_b_combined.png", "Pearson r")
    
    plot_fidelity_combined(fid_numba, "L_Corr_to_MATLAB", "Fidelity: Looping Matrix (L)\nNumba vs MATLAB", 
                           "fidelity_pearson_l_combined.png", "Pearson r")

    # New Fidelity Plots: Non-Zero only
    plot_fidelity_combined(fid_numba, "B_Corr_NZ_to_MATLAB", "Fidelity: Binding Matrix (B) [Non-Zero]\nNumba vs MATLAB", 
                           "fidelity_pearson_b_nz_combined.png", "Pearson r (NZ)")
    
    plot_fidelity_combined(fid_numba, "L_Corr_NZ_to_MATLAB", "Fidelity: Looping Matrix (L) [Non-Zero]\nNumba vs MATLAB", 
                           "fidelity_pearson_l_nz_combined.png", "Pearson r (NZ)")

    # --- Keep the original speed plot structure but maybe update to reflect user preference if needed ---
    # The user didn't explicitly ask to change the speed plots, but usually combined is better.
    # For now, I will keep the original speed plots as they were 3x2 grids unless asked otherwise.
    # However, I should update them to follow the same "only numba/matlab" if desired.
    # User said: "for the benchmarking fidelity visualizations... make these changes" 
    # So I will only change fidelity for now.

if __name__ == "__main__":
    main()
