import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import argparse

# Import custom plotting helpers
sys.path.append('/mnt/home/csun1/scripts/global_scripts')
import helpers.font as hf
hf.set_font_family("Google Sans")

def plot_comparison(py_times, ml_times, py_numba_times, iter_num, out_img):
    stages = ['Data Load', 'Circuit Build', 'Init (OLS)', f'Gibbs ({iter_num} iters)']
    
    # Validation
    if len(py_times) != 4 or len(ml_times) != 4 or len(py_numba_times) != 4:
        print("Error: Expected 4 stages of timing data for all implementations.")
        return

    x = np.arange(len(stages))
    width = 0.25
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 6), dpi=300, layout="constrained")
    
    rects1 = ax.bar(x - width, ml_times, width, label='MATLAB (Original)', color='#FBBC05')
    rects2 = ax.bar(x, py_times, width, label='Python (NumPy)', color='#4285F4')
    rects3 = ax.bar(x + width, py_numba_times, width, label='Python (Numba)', color='#34A853')
    
    ax.set_title(f'Runtime Comparison: MATLAB vs Python vs Numba\n({iter_num} Iterations)', loc='left', fontproperties=hf.mf, pad=15)
    ax.set_ylabel('Time (seconds)', fontproperties=hf.sf)
    
    ax.set_xticks(x)
    ax.set_xticklabels(stages)
    
    [l.set_fontproperties(hf.sf) for l in ax.get_yticklabels()]
    [l.set_fontproperties(hf.sf) for l in ax.get_xticklabels()]
    plt.setp(ax.get_xticklabels(), rotation=15, ha='right')
    
    ax.legend(prop=hf.sf)
    
    def autolabel(rects):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for rect in rects:
            height = rect.get_height()
            if height > 0:
                ax.annotate(f'{height:.1f}s',
                            xy=(rect.get_x() + rect.get_width() / 2, height),
                            xytext=(0, 3),  # 3 points vertical offset
                            textcoords="offset points",
                            ha='center', va='bottom',
                            fontsize=8,
                            fontproperties=hf.sf)
                        
    autolabel(rects1)
    autolabel(rects2)
    autolabel(rects3)
    
    import seaborn as sns
    sns.despine(ax=ax)
    
    plt.savefig(out_img)
    print(f"Saved comparison plot to {out_img}")

def main():
    parser = argparse.ArgumentParser(description="Plot MATLAB vs Python vs Numba performance comparison.")
    parser.add_argument("--iter", type=int, default=10, help="Number of iterations used in profiling")
    parser.add_argument("--ml-times", type=float, nargs=4, default=[18.2, 1.4, 27.5, 12.9], help="MATLAB times: Load, Circuit, Init, Gibbs")
    parser.add_argument("--py-times", type=float, nargs=4, default=[1.2, 3.3, 2.4, 7.5], help="Python NumPy times: Load, Circuit, Init, Gibbs")
    parser.add_argument("--numba-times", type=float, nargs=4, default=[1.2, 3.3, 2.4, 1.5], help="Python Numba times: Load, Circuit, Init, Gibbs")
    parser.add_argument("--output", type=str, default="runtime_comparison.png", help="Output plot filename")
    
    args = parser.parse_args()
    plot_comparison(args.py_times, args.ml_times, args.numba_times, args.iter, args.output)

if __name__ == "__main__":
    main()
