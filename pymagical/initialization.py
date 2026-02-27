import numpy as np
import statsmodels.api as sm

def initialize_magical(
    cand_tf_log2count, cand_peak_log2count, cand_gene_log2count,
    cand_tf_peak_binding, cand_peak_gene_looping,
    M, S, atac_cell_vector, rna_cell_vector
):
    print("Initializing MAGICAL model parameters...")
    
    # TF activity prior: TF RNA expression
    t_sample_mean = cand_tf_log2count.copy()
    # MATLAB var(X, 0, 2) is N-1 degrees of freedom
    t_sample_var = np.var(cand_tf_log2count, axis=1, ddof=1)
    
    t_a_prior = np.zeros((M, len(atac_cell_vector)))
    t_r_prior = np.zeros((M, len(rna_cell_vector)))
    
    for m in range(M):
        for s in range(S):
            a_index = np.where(atac_cell_vector == s)[0]
            t_a_prior[m, a_index] = np.random.randn(len(a_index)) * np.sqrt(t_sample_var[m]) + t_sample_mean[m, s]
            
            r_index = np.where(rna_cell_vector == s)[0]
            t_r_prior[m, r_index] = np.random.randn(len(r_index)) * np.sqrt(t_sample_var[m]) + t_sample_mean[m, s]
            
    # TF-peak binding prior: regression weight
    # If cand_tf_peak_binding is sparse, convert to dense
    if hasattr(cand_tf_peak_binding, 'toarray'):
        b_prior = cand_tf_peak_binding.toarray().astype(float)
    else:
        b_prior = cand_tf_peak_binding.astype(float).copy()
        
    b_prob = b_prior.copy()
    b_var = np.ones(M)
    
    xx, yy = np.where(b_prior > 0)
    for i in range(len(xx)):
        p_idx = xx[i]
        t_idx = yy[i]
        
        X = sm.add_constant(cand_tf_log2count[t_idx, :])
        y = cand_peak_log2count[p_idx, :]
        mdl = sm.OLS(y, X).fit()
        
        # mdl.params[1] is slope, mdl.pvalues[1] is p-value
        if len(mdl.params) > 1:
            b_prior[p_idx, t_idx] = mdl.params[1]
            b_prob[p_idx, t_idx] = 1 - mdl.pvalues[1]
        else:
            b_prior[p_idx, t_idx] = 0.0
            b_prob[p_idx, t_idx] = 0.0
            
    b_mean = b_prior.copy()
    for m in range(M):
        mask = b_prior[:, m] > 0
        if np.any(mask):
            b_var[m] = np.var(b_prior[mask, m], ddof=1)
        if np.isnan(b_var[m]) or b_var[m] == 0:
            b_var[m] = 1.0 # fallback
            
    # Peak-Gene looping prior
    if hasattr(cand_peak_gene_looping, 'toarray'):
        l_prior = cand_peak_gene_looping.toarray().astype(float)
    else:
        l_prior = cand_peak_gene_looping.astype(float).copy()
        
    l_prob = l_prior.copy()
    
    xx, yy = np.where(l_prior > 0)
    for i in range(len(xx)):
        p_idx = xx[i]
        g_idx = yy[i]
        
        X = sm.add_constant(cand_peak_log2count[p_idx, :])
        y = cand_gene_log2count[g_idx, :]
        mdl = sm.OLS(y, X).fit()
        
        if len(mdl.params) > 1:
            l_prior[p_idx, g_idx] = mdl.params[1]
            l_prob[p_idx, g_idx] = 1 - mdl.pvalues[1]
        else:
            l_prior[p_idx, g_idx] = 0.0
            l_prob[p_idx, g_idx] = 0.0
            
    l_mean = l_prior.copy()
    nonzero_l = l_prior[l_prior != 0]
    if len(nonzero_l) > 1:
        l_var = np.var(nonzero_l, ddof=1)
    else:
        l_var = 1.0
        
    return (
        t_a_prior, t_r_prior, t_sample_mean, t_sample_var,
        b_prior, b_mean, b_var, b_prob,
        l_prior, l_mean, l_var, l_prob
    )
