import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import argparse

# Import custom plotting helpers
sys.path.append('/mnt/home/csun1/scripts/global_scripts')
import helpers.font as hf
hf.set_font_family("Google Sans")

def plot_comparison(py_times, ml_times, iter_num, out_img):
    stages = ['Data Load', 'Circuit Build', 'Init (OLS)', f'Gibbs ({iter_num} iters)']
    
    # Validation
    if len(py_times) != 4 or len(ml_times) != 4:
        print("Error: Expected 4 stages of timing data.")
        return

    x = np.arange(len(stages))
    width = 0.35
    
    fig, ax = plt.subplots(1, 1, figsize=(8, 5), dpi=300, layout="constrained")
    
    rects1 = ax.bar(x - width/2, py_times, width, label='Python (pymagical)', color='#4285F4')
    rects2 = ax.bar(x + width/2, ml_times, width, label='MATLAB (Original)', color='#FBBC05')
    
    ax.set_title(f'Runtime Comparison: Python vs MATLAB\n({iter_num} Iterations)', loc='left', fontproperties=hf.mf, pad=15)
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
            ax.annotate(f'{height:.1f}s',
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom',
                        fontproperties=hf.sf)
                        
    autolabel(rects1)
    autolabel(rects2)
    
    import seaborn as sns
    sns.despine(ax=ax)
    
    plt.savefig(out_img)
    print(f"Saved comparison plot to {out_img}")

def main():
    parser = argparse.ArgumentParser(description="Plot Python vs MATLAB performance comparison.")
    parser.add_argument("--iter", type=int, default=10, help="Number of iterations used in profiling")
    parser.add_argument("--py-times", type=float, nargs=4, default=[1.2, 3.3, 2.4, 7.5], help="Python times: Load, Circuit, Init, Gibbs")
    parser.add_argument("--ml-times", type=float, nargs=4, default=[18.2, 1.4, 27.5, 12.9], help="MATLAB times: Load, Circuit, Init, Gibbs")
    parser.add_argument("--output", type=str, default="runtime_comparison.png", help="Output plot filename")
    
    args = parser.parse_args()
    plot_comparison(args.py_times, args.ml_times, args.iter, args.output)

if __name__ == "__main__":
    main()
