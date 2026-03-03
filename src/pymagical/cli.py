import argparse
import os
import sys

def main():
    parser = argparse.ArgumentParser(description="pymagical - Python port of MAGICAL regulatory circuit inference.")
    subparsers = parser.add_subparsers(dest="command", help="Commands")

    # --- Run Command ---
    run_parser = subparsers.add_parser("run", help="Run MAGICAL inference (default if no command given)")
    run_parser.add_argument("--iter", type=int, default=500, help="Number of Gibbs sampling iterations (default: 500)")
    run_parser.add_argument("--outdir", type=str, default="outputs", help="Output directory for results")
    run_parser.add_argument("--prefix", type=str, default="magical", help="Prefix for output filenames (default: magical)")
    run_parser.add_argument("--dump-weights", action="store_true", help="Dump history of continuous B and L weights as .npy files")
    run_parser.add_argument("--use-numba", action="store_true", help="Enable Numba JIT optimization for Gibbs sampling (faster, requires numba)")
    
    # Directory-based input
    run_parser.add_argument("--main-dir", type=str, help="Main data folder containing motifs, tad, and refseq info")
    run_parser.add_argument("--cell-dir", type=str, help="Celltype specific folder containing scRNA and scATAC data")

    # Input File Overrides (defaults will be set if directories are provided)
    run_parser.add_argument("--cand-genes", type=str, help="Override: sig_cr_genes.txt path")
    run_parser.add_argument("--cand-peaks", type=str, help="Override: sig_cr_peaks.txt path")
    run_parser.add_argument("--rna-counts", type=str, help="Override: rna_counts.txt path")
    run_parser.add_argument("--rna-genes", type=str, help="Override: rna_genes.txt path")
    run_parser.add_argument("--rna-meta", type=str, help="Override: rna_meta.txt path")
    run_parser.add_argument("--atac-counts", type=str, help="Override: atac_counts.txt path")
    run_parser.add_argument("--atac-peaks", type=str, help="Override: atac_peaks.txt path")
    run_parser.add_argument("--atac-meta", type=str, help="Override: atac_meta.txt path")
    run_parser.add_argument("--motif-mapping", type=str, help="Override: motif_prior.txt path")
    run_parser.add_argument("--motif-info", type=str, help="Override: motif_info.txt path")
    run_parser.add_argument("--tad-file", type=str, help="Override: tad_regions.txt path")
    run_parser.add_argument("--refseq-file", type=str, help="Override: rhemac10_refseq.txt path")

    # --- Viz Command ---
    viz_parser = subparsers.add_parser("viz", help="Generate interactive HTML report from results")
    viz_parser.add_argument("input", type=str, help="Input result text file (e.g., astrocytes_py_2000.txt)")
    viz_parser.add_argument("--output", type=str, default=None, help="Output HTML file path (defaults to input path with .html extension)")

    # Legacy support / default to run
    # If -h/--help is the only argument, we want to show help for 'run' too
    # since it is the primary command and users expect to see its flags.
    if len(sys.argv) == 2 and sys.argv[1] in ["-h", "--help"]:
        # Print top-level help, then explicitly print 'run' help
        parser.print_help()
        print("\n'run' command options (default):")
        run_parser.print_help()
        sys.exit(0)

    if len(sys.argv) > 1 and sys.argv[1] not in ["run", "viz", "-h", "--help"]:
        sys.argv.insert(1, "run")
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        print("\nUse 'pymagical run --help' to see all inference flags.")
        sys.exit(1)

    args = parser.parse_args()

    if args.command == "run":
        # Defer imports
        from .magical import run_magical

        # Resolve paths
        if args.main_dir and args.cell_dir:
            # If cell_dir is provided as a relative path, assume it is under main_dir
            if not os.path.isabs(args.cell_dir):
                args.cell_dir = os.path.join(args.main_dir, args.cell_dir)

        def resolve(val, folder, filename):
            if val is not None:
                return val
            if folder is not None:
                return os.path.join(folder, filename)
            return None

        cand_genes = resolve(args.cand_genes, args.cell_dir, "sig_cr_genes.txt")
        cand_peaks = resolve(args.cand_peaks, args.cell_dir, "sig_cr_peaks.txt")
        rna_counts = resolve(args.rna_counts, args.cell_dir, "rna_counts.txt")
        rna_genes = resolve(args.rna_genes, args.cell_dir, "rna_genes.txt")
        rna_meta = resolve(args.rna_meta, args.cell_dir, "rna_meta.txt")
        atac_counts = resolve(args.atac_counts, args.cell_dir, "atac_counts.txt")
        atac_peaks = resolve(args.atac_peaks, args.cell_dir, "atac_peaks.txt")
        atac_meta = resolve(args.atac_meta, args.cell_dir, "atac_meta.txt")
        
        motif_mapping = resolve(args.motif_mapping, args.main_dir, "motif_prior.txt")
        motif_info = resolve(args.motif_info, args.main_dir, "motif_info.txt")
        tad_file = resolve(args.tad_file, args.main_dir, "tad_regions.txt")
        refseq_file = resolve(args.refseq_file, args.main_dir, "rhemac10_refseq.txt")

        # Validate that we have all required files
        required_files = {
            "Candidate Genes": cand_genes,
            "Candidate Peaks": cand_peaks,
            "RNA Counts": rna_counts,
            "RNA Genes": rna_genes,
            "RNA Meta": rna_meta,
            "ATAC Counts": atac_counts,
            "ATAC Peaks": atac_peaks,
            "ATAC Meta": atac_meta,
            "Motif Mapping": motif_mapping,
            "Motif Info": motif_info,
            "TAD File": tad_file,
            "RefSeq File": refseq_file
        }

        missing = [name for name, path in required_files.items() if path is None or not os.path.exists(path)]
        if missing:
            print(f"Error: Missing required input files or directories:\n  " + "\n  ".join(missing))
            print("\nPlease provide --main-dir and --cell-dir, or individual file overrides.")
            sys.exit(1)

        os.makedirs(args.outdir, exist_ok=True)
        out_file = os.path.join(args.outdir, f"{args.prefix}_py_{args.iter}.txt")

        print(f"Running pymagical for {args.iter} iterations...")
        if args.dump_weights:
            print("Weight history dump enabled.")
        
        run_magical(
            cand_gene_file=cand_genes,
            cand_peak_file=cand_peaks,
            rna_counts_file=rna_counts,
            rna_genes_file=rna_genes,
            rna_meta_file=rna_meta,
            atac_counts_file=atac_counts,
            atac_peaks_file=atac_peaks,
            atac_meta_file=atac_meta,
            motif_mapping_file=motif_mapping,
            motif_name_file=motif_info,
            tad_flag=1,
            tad_file=tad_file,
            refseq_file=refseq_file,
            output_file=out_file,
            iteration_num=args.iter,
            dump_weight_history=args.dump_weights,
            use_numba=args.use_numba
        )
        
    elif args.command == "viz":
        try:
            # Defer imports
            from .viz import generate_report
            output = args.output if args.output else args.input.replace(".txt", ".html")
            generate_report(args.input, output)
        except ImportError:
            print("Error: The 'viz' dependencies (plotly, jinja2) are not installed.")
            print("Please install them using: pip install pymagical[viz]")
            sys.exit(1)

if __name__ == "__main__":
    main()
