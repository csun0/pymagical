import argparse
import os
import sys
from .magical import run_magical
from .viz import generate_report

def main():
    parser = argparse.ArgumentParser(description="pymagical - Python port of MAGICAL regulatory circuit inference.")
    subparsers = parser.add_subparsers(dest="command", help="Commands")

    # --- Run Command ---
    run_parser = subparsers.add_parser("run", help="Run MAGICAL inference (default if no command given)")
    run_parser.add_argument("--iter", type=int, default=500, help="Number of Gibbs sampling iterations (default: 500)")
    run_parser.add_argument("--outdir", type=str, default="outputs", help="Output directory for results")
    run_parser.add_argument("--prefix", type=str, default="oligodendrocytes", help="Prefix for output filenames (default: oligodendrocytes)")
    run_parser.add_argument("--dump-weights", action="store_true", help="Dump history of continuous B and L weights as .npy files")
    run_parser.add_argument("--use-numba", action="store_true", help="Enable Numba JIT optimization for Gibbs sampling (faster, requires numba)")
    
    # Input File Overrides
    test_dir = "/mnt/ceph/users/agebrain/anderson/snmulti_data/pymagical/test_data"
    astrocytes_dir = os.path.join(test_dir, "oligodendrocytes")
    
    run_parser.add_argument("--cand-genes", type=str, default=os.path.join(astrocytes_dir, "sig_cr_genes.txt"))
    run_parser.add_argument("--cand-peaks", type=str, default=os.path.join(astrocytes_dir, "sig_cr_peaks.txt"))
    run_parser.add_argument("--rna-counts", type=str, default=os.path.join(astrocytes_dir, "rna_counts.txt"))
    run_parser.add_argument("--rna-genes", type=str, default=os.path.join(astrocytes_dir, "rna_genes.txt"))
    run_parser.add_argument("--rna-meta", type=str, default=os.path.join(astrocytes_dir, "rna_meta.txt"))
    run_parser.add_argument("--atac-counts", type=str, default=os.path.join(astrocytes_dir, "atac_counts.txt"))
    run_parser.add_argument("--atac-peaks", type=str, default=os.path.join(astrocytes_dir, "atac_peaks.txt"))
    run_parser.add_argument("--atac-meta", type=str, default=os.path.join(astrocytes_dir, "atac_meta.txt"))
    run_parser.add_argument("--motif-mapping", type=str, default=os.path.join(test_dir, "motif_prior.txt"))
    run_parser.add_argument("--motif-info", type=str, default=os.path.join(test_dir, "motif_info.txt"))
    run_parser.add_argument("--tad-file", type=str, default=os.path.join(test_dir, "tad_regions.txt"))
    run_parser.add_argument("--refseq-file", type=str, default=os.path.join(test_dir, "rhemac10_refseq.txt"))

    # --- Viz Command ---
    viz_parser = subparsers.add_parser("viz", help="Generate interactive HTML report from results")
    viz_parser.add_argument("input", type=str, help="Input result text file (e.g., astrocytes_py_2000.txt)")
    viz_parser.add_argument("--output", type=str, default=None, help="Output HTML file path (defaults to input path with .html extension)")

    # Legacy support / default to run
    if len(sys.argv) > 1 and sys.argv[1] not in ["run", "viz", "-h", "--help"]:
        sys.argv.insert(1, "run")
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    if args.command == "run":
        os.makedirs(args.outdir, exist_ok=True)
        out_file = os.path.join(args.outdir, f"{args.prefix}_py_{args.iter}.txt")

        print(f"Running pymagical for {args.iter} iterations...")
        if args.dump_weights:
            print("Weight history dump enabled.")
        
        run_magical(
            cand_gene_file=args.cand_genes,
            cand_peak_file=args.cand_peaks,
            rna_counts_file=args.rna_counts,
            rna_genes_file=args.rna_genes,
            rna_meta_file=args.rna_meta,
            atac_counts_file=args.atac_counts,
            atac_peaks_file=args.atac_peaks,
            atac_meta_file=args.atac_meta,
            motif_mapping_file=args.motif_mapping,
            motif_name_file=args.motif_info,
            tad_flag=1,
            tad_file=args.tad_file,
            refseq_file=args.refseq_file,
            output_file=out_file,
            iteration_num=args.iter,
            dump_weight_history=args.dump_weights,
            use_numba=args.use_numba
        )
        
        # Optional: Auto-generate report after run
        # report_file = out_file.replace('.txt', '.html')
        # generate_report(out_file, report_file)

    elif args.command == "viz":
        output = args.output if args.output else args.input.replace(".txt", ".html")
        generate_report(args.input, output)

if __name__ == "__main__":
    main()
