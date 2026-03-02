import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def main():
    # Load data
    perf_df = pd.read_csv("benchmark_performance.csv")
    conv_df = pd.read_csv("benchmark_convergence.csv")

    # Filter for astrocytes (though the summary should be similar across others if aggregated)
    # The user asked specifically for astrocytes
    perf_df = perf_df[perf_df["Dataset"] == "astrocytes"]
    conv_df = conv_df[conv_df["Dataset"] == "astrocytes"]

    # --- Plot 1: Speed vs Iterations ---
    plt.figure(figsize=(10, 6))
    sns.set_style("whitegrid")
    
    # Use Gibbs_Time for pure sampling performance
    # Map implementation names for cleaner legend
    name_map = {"matlab": "MATLAB (Original)", "numpy": "Python (NumPy)", "numba": "Python (Numba)"}
    perf_plot_df = perf_df.copy()
    perf_plot_df["Implementation"] = perf_plot_df["Implementation"].map(name_map)
    
    sns.lineplot(data=perf_plot_df, x="Iterations", y="Gibbs_Time", hue="Implementation", marker="o", linewidth=2.5)
    
    plt.yscale("log")
    plt.title("Gibbs Sampling Performance: Runtime Scaling", fontsize=16, fontweight='bold')
    plt.xlabel("Number of Iterations", fontsize=12)
    plt.ylabel("Gibbs Sampling Time (seconds, log scale)", fontsize=12)
    plt.legend(title="Implementation", frameon=True)
    
    plt.tight_layout()
    plt.savefig("speed_benchmark.png", dpi=300)
    print("Saved speed_benchmark.png")

    # --- Plot 2: Convergence vs Iterations ---
    # We want to show how similarity between N and N+1 increases
    # Map '100_vs_500' to the higher iteration count for X-axis positioning
    iter_pair_map = {
        "100_vs_500": 500,
        "500_vs_1000": 1000,
        "1000_vs_2000": 2000,
        "2000_vs_5000": 5000
    }
    conv_plot_df = conv_df.copy()
    conv_plot_df["Higher_Iter"] = conv_plot_df["Iter_Pair"].map(iter_pair_map)
    conv_plot_df["Implementation"] = conv_plot_df["Implementation"].map(name_map)
    
    fig, ax1 = plt.subplots(figsize=(10, 6))
    
    # Jaccard Similarity (Sets)
    sns.lineplot(data=conv_plot_df, x="Higher_Iter", y="Jaccard_Converge", hue="Implementation", 
                 marker="s", linestyle="--", ax=ax1, legend=False)
    ax1.set_ylabel("Circuit Jaccard Similarity (Within-Method)", fontsize=12, color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')
    ax1.set_xlabel("Iterations (compared to previous level)", fontsize=12)
    
    # Pearson Correlation (Matrices)
    ax2 = ax1.twinx()
    sns.lineplot(data=conv_plot_df, x="Higher_Iter", y="L_Corr_Converge", hue="Implementation", 
                 marker="o", ax=ax2)
    ax2.set_ylabel("Matrix Pearson Correlation (L matrix)", fontsize=12, color='green')
    ax2.tick_params(axis='y', labelcolor='green')
    
    plt.title("Algorithm Convergence: Stability across Iterations", fontsize=16, fontweight='bold')
    
    # Custom legend handling
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='blue', lw=2, linestyle='--'),
                    Line2D([0], [0], color='green', lw=2)]
    ax2.legend(custom_lines, ['Jaccard (Circuit Overlap)', 'Pearson (Matrix Similarity)'], loc='lower right')
    
    plt.tight_layout()
    plt.savefig("convergence_benchmark.png", dpi=300)
    print("Saved convergence_benchmark.png")

if __name__ == "__main__":
    main()
