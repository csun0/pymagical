import argparse
import os
from .magical import run_magical

def main():
    parser = argparse.ArgumentParser(description="pymagical - Python port of MAGICAL regulatory circuit inference.")
    
    # Core Parameters
    parser.add_argument("--iter", type=int, default=500, help="Number of Gibbs sampling iterations (default: 500)")
    parser.add_argument("--outdir", type=str, default="outputs", help="Output directory for results")
    parser.add_argument("--prefix", type=str, default="astrocytes", help="Prefix for output filenames (default: astrocytes)")
    parser.add_argument("--dump-weights", action="store_true", help="Dump history of continuous B and L weights as .npy files")
    
    # Input File Overrides (Defaults to astrocytes demo data)
    test_dir = "/mnt/ceph/users/agebrain/anderson/snmulti_data/pymagical/test_data"
    astrocytes_dir = os.path.join(test_dir, "astrocytes")
    
    parser.add_argument("--cand-genes", type=str, default=os.path.join(astrocytes_dir, "sig_cr_genes.txt"))
    parser.add_argument("--cand-peaks", type=str, default=os.path.join(astrocytes_dir, "sig_cr_peaks.txt"))
    parser.add_argument("--rna-counts", type=str, default=os.path.join(astrocytes_dir, "rna_counts.txt"))
    parser.add_argument("--rna-genes", type=str, default=os.path.join(astrocytes_dir, "rna_genes.txt"))
    parser.add_argument("--rna-meta", type=str, default=os.path.join(astrocytes_dir, "rna_meta.txt"))
    parser.add_argument("--atac-counts", type=str, default=os.path.join(astrocytes_dir, "atac_counts.txt"))
    parser.add_argument("--atac-peaks", type=str, default=os.path.join(astrocytes_dir, "atac_peaks.txt"))
    parser.add_argument("--atac-meta", type=str, default=os.path.join(astrocytes_dir, "atac_meta.txt"))
    parser.add_argument("--motif-mapping", type=str, default=os.path.join(test_dir, "motif_prior.txt"))
    parser.add_argument("--motif-info", type=str, default=os.path.join(test_dir, "motif_info.txt"))
    parser.add_argument("--tad-file", type=str, default=os.path.join(test_dir, "tad_regions.txt"))
    parser.add_argument("--refseq-file", type=str, default=os.path.join(test_dir, "rhemac10_refseq.txt"))
    
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    # Use the prefix and iteration count to name the final result file, coherent with MATLAB naming
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
        dump_weight_history=args.dump_weights
    )

if __name__ == "__main__":
    main()
