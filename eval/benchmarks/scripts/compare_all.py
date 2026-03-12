import os
import subprocess
import pandas as pd
import re

CELLTYPES = ["astrocytes", "endothelial", "excitatory_neurons", "inhibitory_neurons", "microglia", "opcs", "oligodendrocytes"]
ITERATIONS = 1000

PY_OUTDIR = "eval/benchmarks/outputs/python"
ML_OUTDIR = "eval/benchmarks/outputs/matlab"

def run_comparison(celltype):
    py_prefix = f"{celltype}_py"
    # pymagical CLI adds the iteration to the filename
    py_base = os.path.join(PY_OUTDIR, celltype, f"{py_prefix}_py_{ITERATIONS}")
    ml_base = os.path.join(ML_OUTDIR, f"{celltype}_ml_{ITERATIONS}")
    
    cmd = [
        "uv", "run", "python", "eval/tests/compare_results.py",
        "--ml-dir", ML_OUTDIR,
        "--py-dir", os.path.join(PY_OUTDIR, celltype),
        "--iter", str(ITERATIONS),
        "--ml-prefix", celltype,
        "--py-prefix", py_prefix
    ]
    
    print(f"Comparing {celltype}...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.stdout

def parse_metrics(output, celltype):
    metrics = {"Celltype": celltype}
    
    # Extract correlations
    l_corr = re.search(r"L matrix.*Pearson correlation:\s*([\d.]+)", output)
    b_corr = re.search(r"B matrix.*Pearson correlation:\s*([\d.]+)", output)
    if l_corr: metrics["L_Corr"] = float(l_corr.group(1))
    if b_corr: metrics["B_Corr"] = float(b_corr.group(1))
    
    # Extract Triad Overlap
    ml_triads = re.search(r"Triad Overlap:.*?MATLAB total: (\d+)", output, re.DOTALL)
    py_triads = re.search(r"Triad Overlap:.*?Python total: (\d+)", output, re.DOTALL)
    jaccard = re.search(r"Triad Overlap:.*?Jaccard Sim:\s+([\d.]+)", output, re.DOTALL)
    recovery = re.search(r"Triad Overlap:.*?Recovery:\s+([\d.]+)%", output, re.DOTALL)
    
    if ml_triads: metrics["ML_Triads"] = int(ml_triads.group(1))
    if py_triads: metrics["PY_Triads"] = int(py_triads.group(1))
    if jaccard: metrics["Jaccard"] = float(jaccard.group(1))
    if recovery: metrics["Recovery_%"] = float(recovery.group(1))
    
    # Extract Speeds
    ml_speed = re.search(r"Gibbs Sampling\s+\|\s+([\d.]+)", output)
    py_speed = re.search(r"Gibbs Sampling\s+\|\s+[\d.]+\s+\|\s+([\d.]+)", output)
    speedup = re.search(r"Gibbs Sampling\s+\|\s+[\d.]+\s+\|\s+[\d.]+\s+\|\s+([\d.]+)x", output)
    
    if ml_speed: metrics["ML_Time(s)"] = float(ml_speed.group(1))
    if py_speed: metrics["PY_Time(s)"] = float(py_speed.group(1))
    if speedup: metrics["Speedup"] = float(speedup.group(1))
    
    return metrics

def main():
    all_metrics = []
    for ct in CELLTYPES:
        output = run_comparison(ct)
        metrics = parse_metrics(output, ct)
        all_metrics.append(metrics)
    
    df = pd.DataFrame(all_metrics)
    print("\n" + "="*80)
    print(f"SUMMARY REPORT ({ITERATIONS} iterations)")
    print("="*80)
    print(df.to_string(index=False))
    print("="*80)
    
    df.to_csv("eval/benchmarks/fidelity_summary.csv", index=False)
    print("\nSummary saved to eval/benchmarks/fidelity_summary.csv")

if __name__ == "__main__":
    main()
