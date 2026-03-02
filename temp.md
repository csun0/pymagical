 Benchmark Execution Plan

   1. Implementation:
       * Modified matlab_runner.m to be dataset-agnostic.
       * Created eval/benchmarks/scripts/run_comprehensive_benchmark.sh to
         launch a Slurm array (105 jobs) covering all datasets,
         implementations, and iteration levels.
       * Developed eval/benchmarks/analyze_benchmark.py to aggregate results
         and compute fidelity/convergence metrics.

   2. Launch the Benchmark:
      Run the following command to submit the jobs to the Slurm cluster:
   1     sbatch eval/benchmarks/scripts/run_comprehensive_benchmark.sh

   3. Analyze Results:
      Once the jobs complete, run the analysis script to generate summary CSVs
  and console reports:
   1     uv run python eval/benchmarks/analyze_benchmark.py --bench-dir
     outputs_bench

   4. Metrics Captured:
       * Speed: Runtime per iteration and total stage timing.
       * Fidelity: Jaccard similarity and Recovery percentage of (Gene, Peak,
         TF) triads vs. MATLAB.
       * Convergence: Pearson correlation of B (TF-Peak) and L (Peak-Gene)
         matrices between successive iteration levels (e.g., 2000 vs 5000) to
         identify the point of diminishing returns.