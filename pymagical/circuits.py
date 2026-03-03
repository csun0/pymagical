import numpy as np
import pandas as pd
from scipy import sparse

def construct_candidate_circuits_with_tad(
    common_samples,
    cand_genes, cand_peaks, 
    rna_genes, rna_cells, rna_count_matrix,
    atac_peaks, atac_cells, atac_count_matrix,
    motifs, tf_peak_binding_matrix,
    refseq, tad_regions
):
    print("Starting candidate circuits construction...")
    
    # --- 1. TF-peak binding ---
    # Intersect candidate peaks with scATAC peaks
    # Merge on chr_num, point1, point2
    atac_peaks_merge = atac_peaks[['peak_index', 'chr_num', 'point1', 'point2']].copy()
    cand_peaks_merge = cand_peaks[['chr_num', 'point1', 'point2']].copy()
    cand_peaks_merge['cand_idx'] = np.arange(len(cand_peaks))
    
    merged_peaks = pd.merge(atac_peaks_merge, cand_peaks_merge, on=['chr_num', 'point1', 'point2'], how='inner')
    # Using 'inner' retains order if we sort by cand_idx, but to match MATLAB 'stable' we keep first occurrence
    merged_peaks = merged_peaks.drop_duplicates(subset=['chr_num', 'point1', 'point2']).sort_values('cand_idx')
    
    bidex = merged_peaks['peak_index'].values - 1 # 0-indexed ATAC peak indices
    cidex = merged_peaks['cand_idx'].values # 0-indexed candidate peak indices
    
    curr_cand_peaks = cand_peaks.iloc[cidex].copy().reset_index(drop=True)
    curr_cand_tf_binding = tf_peak_binding_matrix[bidex, :]
    
    num_peaks = curr_cand_tf_binding.shape[0]
    total_atac_peaks = len(atac_peaks)
    
    # TF filtering
    tf_num = np.array(curr_cand_tf_binding.sum(axis=0)).flatten()
    tf_pct = tf_num / num_peaks
    
    # TF enrichment: (tf_num / cand_peaks) / (total_tf_bindings / total_atac_peaks)
    total_tf_bindings = np.array(tf_peak_binding_matrix.sum(axis=0)).flatten()
    # Avoid division by zero
    bg_freq = total_tf_bindings / total_atac_peaks
    bg_freq[bg_freq == 0] = 1.0 # prevent nan
    tf_enrichment_fc = (tf_num / num_peaks) / bg_freq
    
    pct_threshold = 0.05
    num_threshold = 30
    enrichment_threshold = 0.8
    
    tf_index = np.where((tf_pct > pct_threshold) & (tf_num > num_threshold) & (tf_enrichment_fc > enrichment_threshold))[0]
    
    if len(tf_index) > 0:
        cand_tfs = motifs.iloc[tf_index]['name'].values
        curr_cand_tf_binding = curr_cand_tf_binding[:, tf_index]
    else:
        print("Too few peaks with TF binding sites. MAGICAL not applicable to this cell type!")
        return None
    
    # Filter peaks with >0 TF binding
    peak_sums = np.array(curr_cand_tf_binding.sum(axis=1)).flatten()
    peak_idx = np.where(peak_sums > 0)[0]
    
    curr_cand_peaks = curr_cand_peaks.iloc[peak_idx].reset_index(drop=True)
    curr_cand_tf_binding = curr_cand_tf_binding[peak_idx, :]
    
    # --- 2. Peak-gene looping ---
    # Intersect candidate genes with Refseq
    # Use MATLAB-style logic: for each gene name, use min(start) for + strand or max(end) for - strand
    # across all possible transcripts to capture the outermost TSS.
    refseq_extreme = refseq.groupby(['gene_name', 'strand', 'chr_num']).agg({
        'start': 'min',
        'end': 'max'
    }).reset_index()
    
    # In case a gene has multiple entries with different strands/chrs, we follow MATLAB's 
    # 'pick first chr/strand but aggregate coordinates' behavior by deduping on gene_name.
    refseq_extreme = refseq_extreme.drop_duplicates(subset=['gene_name'])
    
    cand_genes_df = pd.DataFrame({'gene_symbol': cand_genes})
    merged_genes = pd.merge(cand_genes_df, refseq_extreme, left_on='gene_symbol', right_on='gene_name', how='inner')
    
    # Calculate TSS based on the strand-specific extreme boundary
    merged_genes['gene_TSS'] = np.where(merged_genes['strand'] == '+', merged_genes['start'], merged_genes['end'])
    curr_cand_genes = merged_genes['gene_symbol'].values
    gene_tss = merged_genes[['chr_num', 'gene_TSS']].values
    
    num_cand_peaks = len(curr_cand_peaks)
    num_cand_genes = len(curr_cand_genes)
    
    peak_gene_looping_tad = np.zeros((num_cand_peaks, num_cand_genes), dtype=int)
    
    peak_chr = curr_cand_peaks['chr_num'].values
    peak_center = (curr_cand_peaks['point1'].values + curr_cand_peaks['point2'].values) / 2
    peak_p1 = curr_cand_peaks['point1'].values
    peak_p2 = curr_cand_peaks['point2'].values
    
    # Build TAD looping mask
    for _, tad in tad_regions.iterrows():
        tad_chr = tad['chr_num']
        t_left = tad['left_boundary']
        t_right = tad['right_boundary']
        
        # Peaks in TAD
        p_idx = np.where((peak_chr == tad_chr) & (peak_p1 > t_left) & (peak_p2 < t_right))[0]
        # Genes in TAD
        g_idx = np.where((gene_tss[:, 0] == tad_chr) & (gene_tss[:, 1] > t_left) & (gene_tss[:, 1] < t_right))[0]
        
        if len(p_idx) > 0 and len(g_idx) > 0:
            for p in p_idx:
                peak_gene_looping_tad[p, g_idx] = 1
                
    # Build Distance looping mask (< 1e6)
    peak_gene_looping_dist = np.zeros((num_cand_peaks, num_cand_genes), dtype=int)
    for g in range(num_cand_genes):
        g_chr = gene_tss[g, 0]
        g_tss = gene_tss[g, 1]
        dist_mask = (peak_chr == g_chr) & (np.abs(peak_center - g_tss) < 1e6)
        peak_gene_looping_dist[dist_mask, g] = 1
        
    curr_cand_peak_gene_looping = peak_gene_looping_tad * peak_gene_looping_dist
    
    # Filter peaks and genes with >0 looping
    p_sums = curr_cand_peak_gene_looping.sum(axis=1)
    g_sums = curr_cand_peak_gene_looping.sum(axis=0)
    
    p_idx = np.where(p_sums > 0)[0]
    g_idx = np.where(g_sums > 0)[0]
    
    curr_cand_peak_gene_looping = curr_cand_peak_gene_looping[p_idx, :][:, g_idx]
    curr_cand_peaks = curr_cand_peaks.iloc[p_idx].reset_index(drop=True)
    curr_cand_genes = curr_cand_genes[g_idx]
    gene_tss = gene_tss[g_idx, :]
    curr_cand_tf_binding = curr_cand_tf_binding[p_idx, :]
    
    # Filter TFs
    tf_sums = np.array(curr_cand_tf_binding.sum(axis=0)).flatten()
    tf_idx = np.where(tf_sums > 0)[0]
    cand_tfs = cand_tfs[tf_idx]
    curr_cand_tf_binding = curr_cand_tf_binding[:, tf_idx]
    
    # --- 3. Pseudo-bulk calculation ---
    num_samples = len(common_samples)
    atac_cell_vector = np.zeros(len(atac_cells), dtype=int)
    rna_cell_vector = np.zeros(len(rna_cells), dtype=int)
    
    atac_counts = np.zeros((atac_count_matrix.shape[0], num_samples))
    rna_counts = np.zeros((rna_count_matrix.shape[0], num_samples))
    
    for s_idx, sample in enumerate(common_samples):
        # ATAC
        a_mask = (atac_cells['subject_ID'] == sample).values
        atac_cell_vector[a_mask] = s_idx + 1 # 1-based internally for output tracking? Let's use 0-based
        if np.any(a_mask):
            atac_counts[:, s_idx] = np.array(atac_count_matrix[:, a_mask].sum(axis=1)).flatten()
            
        # RNA
        r_mask = (rna_cells['subject_ID'] == sample).values
        rna_cell_vector[r_mask] = s_idx + 1
        if np.any(r_mask):
            rna_counts[:, s_idx] = np.array(rna_count_matrix[:, r_mask].sum(axis=1)).flatten()
            
    # Use 0-based for python logic
    atac_cell_vector = atac_cell_vector - 1
    rna_cell_vector = rna_cell_vector - 1
    
    total_atac_reads = 5e6
    atac_raw_count_sum = atac_counts.sum(axis=0) + 1
    for s in range(num_samples):
        atac_counts[:, s] = atac_counts[:, s] / atac_raw_count_sum[s] * total_atac_reads
    atac_log2 = np.log2(atac_counts + 1)
    
    total_rna_reads = 5e6
    rna_raw_count_sum = rna_counts.sum(axis=0) + 1
    for s in range(num_samples):
        rna_counts[:, s] = rna_counts[:, s] / rna_raw_count_sum[s] * total_rna_reads
    rna_log2 = np.log2(rna_counts + 1)
    
    # --- 4. Select actively accessible peaks ---
    # The peaks in curr_cand_peaks need their original index in the atac matrix
    # Re-merge to find original atac peak indices
    atac_peaks_merge = atac_peaks[['peak_index', 'chr_num', 'point1', 'point2']].copy()
    curr_cand_peaks_merge = curr_cand_peaks[['chr_num', 'point1', 'point2']].copy()
    curr_cand_peaks_merge['cand_idx'] = np.arange(len(curr_cand_peaks))
    
    merged = pd.merge(curr_cand_peaks_merge, atac_peaks_merge, on=['chr_num', 'point1', 'point2'], how='inner')
    merged = merged.drop_duplicates(subset=['chr_num', 'point1', 'point2']).sort_values('cand_idx')
    
    bidex = merged['cand_idx'].values
    cidex = merged['peak_index'].values - 1
    
    curr_cand_peaks = curr_cand_peaks.iloc[bidex].reset_index(drop=True)
    curr_cand_tf_binding = curr_cand_tf_binding[bidex, :]
    curr_cand_peak_gene_looping = curr_cand_peak_gene_looping[bidex, :]
    cand_peak_log2count = atac_log2[cidex, :]
    scatac_read_count_matrix = atac_count_matrix[cidex, :]
    
    active_mask = cand_peak_log2count.sum(axis=1) > 0
    curr_cand_peaks = curr_cand_peaks[active_mask].reset_index(drop=True)
    curr_cand_tf_binding = curr_cand_tf_binding[active_mask, :]
    curr_cand_peak_gene_looping = curr_cand_peak_gene_looping[active_mask, :]
    cand_peak_log2count = cand_peak_log2count[active_mask, :]
    scatac_read_count_matrix = scatac_read_count_matrix[active_mask, :]
    
    # mean center
    cand_peak_log2count = cand_peak_log2count - cand_peak_log2count.mean(axis=1, keepdims=True)
    
    # --- 5. Select actively expressed genes ---
    rna_genes_df = rna_genes.copy()
    rna_genes_df['rna_idx'] = np.arange(len(rna_genes_df))
    curr_cand_genes_df = pd.DataFrame({'gene_symbol': curr_cand_genes, 'cand_idx': np.arange(len(curr_cand_genes))})
    
    merged = pd.merge(curr_cand_genes_df, rna_genes_df, on='gene_symbol', how='inner').sort_values('cand_idx')
    bidex = merged['cand_idx'].values
    cidex = merged['rna_idx'].values
    
    curr_cand_genes = curr_cand_genes[bidex]
    gene_tss = gene_tss[bidex, :]
    curr_cand_peak_gene_looping = curr_cand_peak_gene_looping[:, bidex]
    cand_gene_log2count = rna_log2[cidex, :]
    scrna_read_count_matrix = rna_count_matrix[cidex, :]
    
    active_mask = cand_gene_log2count.sum(axis=1) > 0
    curr_cand_genes = curr_cand_genes[active_mask]
    gene_tss = gene_tss[active_mask, :]
    curr_cand_peak_gene_looping = curr_cand_peak_gene_looping[:, active_mask]
    cand_gene_log2count = cand_gene_log2count[active_mask, :]
    scrna_read_count_matrix = scrna_read_count_matrix[active_mask, :]
    
    cand_gene_log2count = cand_gene_log2count - cand_gene_log2count.mean(axis=1, keepdims=True)
    
    # --- 6. Select actively expressed TFs ---
    cand_tfs_df = pd.DataFrame({'gene_symbol': cand_tfs, 'tf_idx': np.arange(len(cand_tfs))})
    merged = pd.merge(cand_tfs_df, rna_genes_df, on='gene_symbol', how='inner').sort_values('tf_idx')
    
    bidex = merged['tf_idx'].values
    cidex = merged['rna_idx'].values
    
    cand_tfs = cand_tfs[bidex]
    curr_cand_tf_binding = curr_cand_tf_binding[:, bidex]
    cand_tf_log2count = rna_log2[cidex, :]
    
    active_mask = cand_tf_log2count.sum(axis=1) > 0
    cand_tfs = cand_tfs[active_mask]
    curr_cand_tf_binding = curr_cand_tf_binding[:, active_mask]
    cand_tf_log2count = cand_tf_log2count[active_mask, :]
    
    cand_tf_log2count = cand_tf_log2count - cand_tf_log2count.mean(axis=1, keepdims=True)
    
    # --- 7. Final filter ---
    p_mask = (curr_cand_peak_gene_looping.sum(axis=1) > 0) & (np.array(curr_cand_tf_binding.sum(axis=1)).flatten() > 0)
    curr_cand_peaks = curr_cand_peaks[p_mask].reset_index(drop=True)
    curr_cand_peak_gene_looping = curr_cand_peak_gene_looping[p_mask, :]
    curr_cand_tf_binding = curr_cand_tf_binding[p_mask, :]
    cand_peak_log2count = cand_peak_log2count[p_mask, :]
    scatac_read_count_matrix = scatac_read_count_matrix[p_mask, :]
    
    g_mask = curr_cand_peak_gene_looping.sum(axis=0) > 0
    curr_cand_genes = curr_cand_genes[g_mask]
    gene_tss = gene_tss[g_mask, :]
    curr_cand_peak_gene_looping = curr_cand_peak_gene_looping[:, g_mask]
    cand_gene_log2count = cand_gene_log2count[g_mask, :]
    scrna_read_count_matrix = scrna_read_count_matrix[g_mask, :]
    
    tf_mask = np.array(curr_cand_tf_binding.sum(axis=0)).flatten() > 0
    cand_tfs = cand_tfs[tf_mask]
    curr_cand_tf_binding = curr_cand_tf_binding[:, tf_mask]
    cand_tf_log2count = cand_tf_log2count[tf_mask, :]
    
    print(f"MAGICAL initially selected {len(cand_tfs)} TFs, {len(curr_cand_peaks)} peaks, and {len(curr_cand_genes)} genes for circuit inference.")
    
    # Add TSS back to genes for returning
    genes_out = {'symbols': curr_cand_genes, 'tss': gene_tss}
    
    return (
        cand_tfs, cand_tf_log2count,
        curr_cand_peaks, cand_peak_log2count,
        genes_out, cand_gene_log2count,
        curr_cand_tf_binding, curr_cand_peak_gene_looping,
        atac_cell_vector, scatac_read_count_matrix,
        rna_cell_vector, scrna_read_count_matrix
    )
