import os
import time
import numpy as np

from .data_loader import (
    load_candidate_genes, load_candidate_peaks, load_scrna_data, 
    load_scatac_data, load_motif_prior, load_tad_regions, load_refseq
)
from .circuits import construct_candidate_circuits_with_tad
from .initialization import initialize_magical
from .estimation import magical_estimation

def run_magical(
    cand_gene_file, cand_peak_file,
    rna_counts_file, rna_genes_file, rna_meta_file,
    atac_counts_file, atac_peaks_file, atac_meta_file,
    motif_mapping_file, motif_name_file,
    tad_flag, tad_file, refseq_file,
    output_file, iteration_num, burn_in=None, dump_weight_history=False,
    use_numba=False
):
    timing = {}
    
    # Default burn-in: 20% of iterations
    if burn_in is None:
        burn_in = iteration_num // 5

    t_start = time.time()
    print('loading all input data ...\n')
    cand_genes = load_candidate_genes(cand_gene_file)
    cand_peaks = load_candidate_peaks(cand_peak_file)
    
    print('Loading scRNAseq data ...\n')
    rna_genes, rna_cells, rna_counts = load_scrna_data(rna_counts_file, rna_genes_file, rna_meta_file)
    
    print('Loading scATACseq data ...\n')
    atac_peaks, atac_cells, atac_counts = load_scatac_data(atac_counts_file, atac_peaks_file, atac_meta_file)
    
    print('Loading TF motif prior ...\n')
    motifs, motif_prior = load_motif_prior(motif_name_file, motif_mapping_file, len(atac_peaks))
    
    common_samples = np.intersect1d(rna_cells['subject_ID'].unique(), atac_cells['subject_ID'].unique())
    print(f"There are paired data for {len(common_samples)} samples/subjects.\n")
    
    refseq = load_refseq(refseq_file)
    timing['Data Loading'] = time.time() - t_start
    
    t_stage = time.time()
    if tad_flag == 1:
        tads = load_tad_regions(tad_file)
        
        circuit_res = construct_candidate_circuits_with_tad(
            common_samples, cand_genes, cand_peaks,
            rna_genes, rna_cells, rna_counts,
            atac_peaks, atac_cells, atac_counts,
            motifs, motif_prior, refseq, tads
        )
    else:
        raise NotImplementedError("TAD_flag=0 not implemented")
        
    if circuit_res is None:
        return
    timing['Circuit Construction'] = time.time() - t_stage
        
    (cand_tfs, cand_tf_log2count, curr_cand_peaks, cand_peak_log2count,
     curr_cand_genes_dict, cand_gene_log2count, curr_cand_tf_binding, curr_cand_peak_gene_looping,
     atac_cell_vector, scatac_read_count_matrix, rna_cell_vector, scrna_read_count_matrix) = circuit_res
     
    M = len(cand_tfs)
    S = len(common_samples)
    P = len(curr_cand_peaks)
    G = len(curr_cand_genes_dict['symbols'])
    
    t_stage = time.time()
    init_res = initialize_magical(
        cand_tf_log2count, cand_peak_log2count, cand_gene_log2count,
        curr_cand_tf_binding, curr_cand_peak_gene_looping,
        M, S, atac_cell_vector, rna_cell_vector
    )
    timing['Initialization (OLS)'] = time.time() - t_stage
    
    (t_a_prior, t_r_prior, t_sample_mean, t_sample_var,
     b_prior, b_mean, b_var, b_prob,
     l_prior, l_mean, l_var, l_prob) = init_res
     
    t_stage = time.time()
    (b_prob_final, l_prob_final, b_weight_final, l_weight_final, 
     b_history, l_history, b_pos_prob, l_pos_prob) = magical_estimation(
        atac_cell_vector, cand_peak_log2count, rna_cell_vector, cand_gene_log2count,
        curr_cand_tf_binding, curr_cand_peak_gene_looping,
        t_a_prior, t_r_prior, t_sample_mean, t_sample_var,
        b_prior, b_mean, b_var, b_prob,
        l_prior, l_mean, l_var, l_prob,
        M, S, P, G, iteration_num, burn_in=burn_in,
        dump_weight_history=dump_weight_history,
        use_numba=use_numba
    )
    timing[f'Gibbs Sampling ({iteration_num} iters)'] = time.time() - t_stage
    
    print('\nWriting output matrices...')
    write_outputs(output_file, b_prob_final, l_prob_final, b_weight_final, l_weight_final, 
                  b_pos_prob, l_pos_prob, curr_cand_peaks, curr_cand_genes_dict, cand_tfs, timing)
    
    if dump_weight_history and b_history is not None and l_history is not None:
        print('Writing weight histories...')
        base_name = os.path.splitext(output_file)[0]
        np.save(f"{base_name}_B_history.npy", b_history)
        np.save(f"{base_name}_L_history.npy", l_history)
        print(f"Histories saved to {base_name}_B_history.npy and {base_name}_L_history.npy")
        
    print("Done!")
    
def write_outputs(output_file, b_prob_final, l_prob_final, b_weight_final, l_weight_final, 
                  b_pos_prob, l_pos_prob, peaks, genes, tfs, timing=None):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    prob_threshold_tf_peak = 0.7
    prob_threshold_peak_gene = 0.7
    
    with open(output_file, 'w') as f:
        f.write('Gene_symbol\tGene_chr\tGene_TSS\tPeak_chr\tPeak_start\tPeak_end\tPeak_Gene_Prob\tTFs(prob, effect [L(consistency), B(consistency)])\n')
        
        xx, yy = np.where(l_prob_final > prob_threshold_peak_gene)
        for i in range(len(xx)):
            p = xx[i]
            g = yy[i]
            
            tf_probs = b_prob_final[p, :]
            tf_indices = np.argsort(tf_probs)[::-1]
            
            valid_tfs = []
            for t_idx in tf_indices:
                if tf_probs[t_idx] > prob_threshold_tf_peak:
                    # Biological effect: TF -> Peak weight * Peak -> Gene weight
                    b_w = b_weight_final[p, t_idx]
                    l_w = l_weight_final[p, g]
                    effect_val = b_w * l_w
                    
                    overall_dir = "+" if effect_val >= 0 else "-"
                    l_dir = "+" if l_w >= 0 else "-"
                    b_dir = "+" if b_w >= 0 else "-"
                    
                    # Consistency is fraction of iterations agreeing with reported sign
                    l_cons = l_pos_prob[p, g]
                    if l_dir == "-": l_cons = 1.0 - l_cons
                    
                    b_cons = b_pos_prob[p, t_idx]
                    if b_dir == "-": b_cons = 1.0 - b_cons
                    
                    valid_tfs.append(f"{tfs[t_idx]} ({tf_probs[t_idx]:.2f}, {overall_dir} [{l_dir}({l_cons:.2f}),{b_dir}({b_cons:.2f})])")
            
            if valid_tfs:
                tf_str = ", ".join(valid_tfs) + ", "
                f.write(f"{genes['symbols'][g]}\tchr{genes['tss'][g, 0]}\t{genes['tss'][g, 1]}\t"
                        f"chr{peaks['chr_num'].iloc[p]}\t{peaks['point1'].iloc[p]}\t{peaks['point2'].iloc[p]}\t"
                        f"{l_prob_final[p, g]:g}\t{tf_str}\n")
                        
    base_name = os.path.splitext(output_file)[0]
    b_file = f"{base_name}_B_matrix.txt"
    l_file = f"{base_name}_L_matrix.txt"
    timing_file = f"{base_name}_timing_stats.txt"
    
    with open(b_file, 'w') as f:
        f.write("Peak_Coordinates\t" + "\t".join(tfs) + "\n")
        for p in range(len(peaks)):
            peak_str = f"chr{peaks['chr_num'].iloc[p]}_{peaks['point1'].iloc[p]}_{peaks['point2'].iloc[p]}"
            probs_str = "\t".join([f"{val:g}" for val in b_prob_final[p, :]])
            f.write(f"{peak_str}\t{probs_str}\n")
            
    with open(l_file, 'w') as f:
        f.write("Peak_Coordinates\t" + "\t".join(genes['symbols']) + "\n")
        for p in range(len(peaks)):
            peak_str = f"chr{peaks['chr_num'].iloc[p]}_{peaks['point1'].iloc[p]}_{peaks['point2'].iloc[p]}"
            probs_str = "\t".join([f"{val:g}" for val in l_prob_final[p, :]])
            f.write(f"{peak_str}\t{probs_str}\n")

    if timing:
        with open(timing_file, 'w') as f:
            for stage, duration in timing.items():
                f.write(f"{stage}: {duration:.6f}\n")
