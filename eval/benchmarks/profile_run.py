import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
import argparse

# Add package root to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))

from pymagical.data_loader import (
    load_candidate_genes, load_candidate_peaks, load_scrna_data, 
    load_scatac_data, load_motif_prior, load_tad_regions, load_refseq
)
from pymagical.circuits import construct_candidate_circuits_with_tad
from pymagical.initialization import initialize_magical
from pymagical.estimation import magical_estimation

# Import custom plotting helpers
sys.path.append('/mnt/home/csun1/scripts/global_scripts')
import helpers.font as hf
hf.set_font_family("Google Sans")

def profile_magical(iter_num=10, out_img='runtime_statistics.png'):
    test_dir = "/mnt/ceph/users/agebrain/anderson/snmulti_data/pymagical/test_data"
    astrocytes_dir = os.path.join(test_dir, "astrocytes")
    
    times = {}
    
    print(f"Profiling Python execution for {iter_num} iterations...")
    
    t0 = time.time()
    cand_genes = load_candidate_genes(os.path.join(astrocytes_dir, "sig_cr_genes.txt"))
    cand_peaks = load_candidate_peaks(os.path.join(astrocytes_dir, "sig_cr_peaks.txt"))
    rna_genes, rna_cells, rna_counts = load_scrna_data(
        os.path.join(astrocytes_dir, "rna_counts.txt"),
        os.path.join(astrocytes_dir, "rna_genes.txt"),
        os.path.join(astrocytes_dir, "rna_meta.txt")
    )
    atac_peaks, atac_cells, atac_counts = load_scatac_data(
        os.path.join(astrocytes_dir, "atac_counts.txt"),
        os.path.join(astrocytes_dir, "atac_peaks.txt"),
        os.path.join(astrocytes_dir, "atac_meta.txt")
    )
    motifs, motif_prior = load_motif_prior(
        os.path.join(test_dir, "motif_info.txt"),
        os.path.join(test_dir, "motif_prior.txt"),
        len(atac_peaks)
    )
    refseq = load_refseq(os.path.join(test_dir, "rhemac10_refseq.txt"))
    tads = load_tad_regions(os.path.join(test_dir, "tad_regions.txt"))
    times['Data Loading'] = time.time() - t0
    
    t0 = time.time()
    common_samples = np.intersect1d(rna_cells['subject_ID'].unique(), atac_cells['subject_ID'].unique())
    circuit_res = construct_candidate_circuits_with_tad(
        common_samples, cand_genes, cand_peaks,
        rna_genes, rna_cells, rna_counts,
        atac_peaks, atac_cells, atac_counts,
        motifs, motif_prior, refseq, tads
    )
    times['Circuit Construction'] = time.time() - t0
    
    (cand_tfs, cand_tf_log2count, curr_cand_peaks, cand_peak_log2count,
     curr_cand_genes_dict, cand_gene_log2count, curr_cand_tf_binding, curr_cand_peak_gene_looping,
     atac_cell_vector, scatac_read_count_matrix, rna_cell_vector, scrna_read_count_matrix) = circuit_res
     
    M = len(cand_tfs)
    S = len(common_samples)
    P = len(curr_cand_peaks)
    G = len(curr_cand_genes_dict['symbols'])
    
    t0 = time.time()
    init_res = initialize_magical(
        cand_tf_log2count, cand_peak_log2count, cand_gene_log2count,
        curr_cand_tf_binding, curr_cand_peak_gene_looping,
        M, S, atac_cell_vector, rna_cell_vector
    )
    times['Initialization (OLS)'] = time.time() - t0
    
    (t_a_prior, t_r_prior, t_sample_mean, t_sample_var,
     b_prior, b_mean, b_var, b_prob,
     l_prior, l_mean, l_var, l_prob) = init_res
     
    t0 = time.time()
    magical_estimation(
        atac_cell_vector, cand_peak_log2count, rna_cell_vector, cand_gene_log2count,
        curr_cand_tf_binding, curr_cand_peak_gene_looping,
        t_a_prior, t_r_prior, t_sample_mean, t_sample_var,
        b_prior, b_mean, b_var, b_prob,
        l_prior, l_mean, l_var, l_prob,
        M, S, P, G, iter_num
    )
    times[f'Gibbs Sampling ({iter_num} iters)'] = time.time() - t0
    
    print("Runtime Statistics:")
    for k, v in times.items():
        print(f"  {k}: {v:.3f} seconds")
        
    labels = list(times.keys())
    values = list(times.values())
    
    fig, axes = plt.subplots(1, 1, figsize=(6, 4), dpi=300, layout="constrained")
    axes = np.ravel(axes)
    ax = axes[0]
    
    ax.bar(labels, values, color=['#4285F4', '#EA4335', '#FBBC05', '#34A853'])
    
    ax.set_ylabel('Time (seconds)', fontproperties=hf.sf)
    ax.set_title(f'Python Execution Breakdown\n({iter_num} iterations)', loc='left', fontproperties=hf.mf, pad=15)
    
    [l.set_fontproperties(hf.sf) for l in ax.get_yticklabels()]
    [l.set_fontproperties(hf.sf) for l in ax.get_xticklabels()]
    plt.setp(ax.get_xticklabels(), rotation=15, ha='right')
    
    import seaborn as sns
    sns.despine(trim=True, ax=ax)
    
    plt.savefig(out_img)
    print(f"Saved plot to {out_img}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Profile Python implementation stages.")
    parser.add_argument("--iter", type=int, default=10, help="Number of iterations")
    parser.add_argument("--output", type=str, default="runtime_statistics.png", help="Output plot filename")
    args = parser.parse_args()
    
    profile_magical(args.iter, args.output)
