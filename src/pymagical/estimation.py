import numpy as np
import scipy.stats as stats

try:
    from .estimation_kernels import (
        sample_b_state_kernel, sample_l_state_kernel, 
        sample_b_weight_kernel, sample_l_weight_kernel, 
        sample_t_kernel, update_ta_tr_kernel
    )
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False

def tf_activity_t_sampling(atac_cell_vector, a_sample, rna_cell_vector, r_sample, b, t_a, t_r, t_mean, t_var, sigma_a_noise, M, S, P, G):
    t_sample = np.zeros((M, S))
    for s in range(S):
        t_sample[:, s] = np.mean(t_a[:, atac_cell_vector == s], axis=1)
        
    tf_index = np.random.permutation(M)
    for m in tf_index:
        bm = b[:, m] 
        bm_dot_bm = np.dot(bm, bm)
        temp_var = (bm_dot_bm * t_var[m] / P) + sigma_a_noise
        
        resid = a_sample - b @ t_sample + np.outer(bm, t_sample[m, :])
        mean_t = (bm @ resid / P + t_mean[m, :] * sigma_a_noise) / temp_var
        variance_t = t_var[m] * sigma_a_noise / temp_var
        
        aa = np.random.randn(S)
        aa = np.clip(aa, -3, 3)
        t_sample[m, :] = aa * np.sqrt(variance_t) + mean_t
        
        for s in range(S):
            a_idx = np.where(atac_cell_vector == s)[0]
            t_a[m, a_idx] = np.random.randn(len(a_idx)) * np.sqrt(variance_t) + t_sample[m, s]
            
            r_idx = np.where(rna_cell_vector == s)[0]
            t_r[m, r_idx] = np.random.randn(len(r_idx)) * np.sqrt(variance_t) + t_sample[m, s]
            
    return t_a, t_r, t_sample

def tf_peak_binding_b_sampling(atac_cell_vector, a_sample, b, t_a, b_state, b_mean, b_var, sigma_a_noise, M, S, P):
    t_sample = np.zeros((M, S))
    for s in range(S):
        t_sample[:, s] = np.mean(t_a[:, atac_cell_vector == s], axis=1)
        
    tf_index = np.random.permutation(M)
    for m in tf_index:
        tm = t_sample[m, :] 
        tm_dot_tm = np.dot(tm, tm)
        temp_var = tm_dot_tm * b_var[m] / S + sigma_a_noise
        
        resid = a_sample - b @ t_sample + np.outer(b[:, m], tm)
        mean_b = (resid @ tm * b_var[m] / S + b_mean[:, m] * sigma_a_noise) / temp_var
        variance_b = b_var[m] * sigma_a_noise / temp_var
        
        bb = np.random.randn(P)
        bb = np.clip(bb, -3, 3)
        b[:, m] = (bb * np.sqrt(variance_b) + mean_b) * b_state[:, m]
        
    return b

def tf_peak_binary_binding_b_state_sampling(atac_cell_vector, a_sample, b, t_a, b_state, b_mean, b_var, b_prob, sigma_a_noise, M, S, P):
    t_sample = np.zeros((M, S))
    for s in range(S):
        t_sample[:, s] = np.mean(t_a[:, atac_cell_vector == s], axis=1)
        
    tf_index = np.random.permutation(M)
    for p in range(P):
        temp = b[p, :] @ t_sample
        for m in tf_index:
            tm = t_sample[m, :]
            tm_dot_tm = np.dot(tm, tm)
            temp_var = (tm_dot_tm * b_var[m] / S) + sigma_a_noise
            
            if b_state[p, m] == 1:
                mean_b = ((a_sample[p, :] - temp + b[p, m]*tm) @ tm * b_var[m] / S + b_mean[p, m]*sigma_a_noise) / temp_var
                variance_b = b_var[m] * sigma_a_noise / temp_var
                
                post_b1 = np.exp(-(b[p, m]*1 - mean_b)**2 / (2*variance_b)) * (b_prob[p, m] + 0.25) + 1e-6
                post_b0 = np.exp(-(b[p, m]*0 - mean_b)**2 / (2*variance_b)) * (1 - b_prob[p, m] + 0.25) + 1e-6
                
                p1 = post_b1 / (post_b1 + post_b0)
                if p1 < np.random.rand():
                    temp = temp - b[p, m]*tm
                    b[p, m] = 0
                    b_state[p, m] = 0
                    
            elif b_state[p, m] == 0 and b_mean[p, m] != 0:
                mean_b = ((a_sample[p, :] - temp) @ tm * b_var[m] / S + b_mean[p, m]*sigma_a_noise) / temp_var
                variance_b = b_var[m] * sigma_a_noise / temp_var
                
                bb = np.random.randn()
                bb = np.clip(bb, -3, 3)
                b_temp = bb * np.sqrt(variance_b) + mean_b
                
                post_b1 = np.exp(-(b_temp*1 - mean_b)**2 / (2*variance_b)) * (b_prob[p, m] + 0.25) + 1e-6
                post_b0 = np.exp(-(b_temp*0 - mean_b)**2 / (2*variance_b)) * (1 - b_prob[p, m] + 0.25) + 1e-6
                
                p1 = post_b1 / (post_b1 + post_b0)
                if p1 >= np.random.rand():
                    b[p, m] = b_temp
                    b_state[p, m] = 1
                    temp = temp + b[p, m]*tm
            else:
                b[p, m] = 0
                b_state[p, m] = 0
                
    return b_state, b

def peak_gene_looping_l_sampling(rna_cell_vector, r_sample, l, b, t_r, l_state, l_mean, l_var, sigma_r_noise, M, S, P, G):
    t_sample = np.zeros((M, S))
    for s in range(S):
        t_sample[:, s] = np.mean(t_r[:, rna_cell_vector == s], axis=1)
        
    a_estimate = b @ t_sample 
    
    p_index = np.random.permutation(P)
    for p in p_index:
        ap = a_estimate[p, :] 
        ap_dot_ap = np.dot(ap, ap)
        temp_var = ap_dot_ap * l_var / S + sigma_r_noise
        
        resid = r_sample - l.T @ a_estimate + np.outer(l[p, :], ap) 
        
        mean_l = (resid @ ap * l_var / S + l_mean[p, :] * sigma_r_noise) / temp_var 
        variance_l = l_var * sigma_r_noise / temp_var
        
        ll = np.random.randn(G)
        ll = np.clip(ll, -3, 3)
        l[p, :] = (ll * np.sqrt(variance_l) + mean_l) * l_state[p, :]
        
    return l

def peak_gene_binary_looping_l_state_sampling(rna_cell_vector, r_sample, l, b, t_r, l_state, l_mean, l_var, l_prob, sigma_r_noise, M, S, P, G):
    t_sample = np.zeros((M, S))
    for s in range(S):
        t_sample[:, s] = np.mean(t_r[:, rna_cell_vector == s], axis=1)
        
    a_estimate = b @ t_sample
    ap_dot_ap_arr = np.sum(a_estimate**2, axis=1)
    
    for g in range(G):
        temp = l[:, g] @ a_estimate 
        for p in range(P):
            ap = a_estimate[p, :]
            temp_var = ap_dot_ap_arr[p] * l_var / S + sigma_r_noise
            
            if l_state[p, g] == 1:
                mean_l = ((r_sample[g, :] - temp + l[p, g]*ap) @ ap * l_var / S + l_mean[p, g]*sigma_r_noise) / temp_var
                variance_l = l_var * sigma_r_noise / temp_var
                
                post_l1 = np.exp(-(l[p, g]*1 - mean_l)**2 / (2*variance_l)) * (l_prob[p, g] + 0.1) + 1e-6
                post_l0 = np.exp(-(l[p, g]*0 - mean_l)**2 / (2*variance_l)) * (1 - l_prob[p, g] + 0.1) + 1e-6
                
                p1 = post_l1 / (post_l1 + post_l0)
                if p1 < np.random.rand():
                    temp = temp - l[p, g]*ap
                    l[p, g] = 0
                    l_state[p, g] = 0
            elif l_state[p, g] == 0 and l_mean[p, g] != 0:
                mean_l = ((r_sample[g, :] - temp) @ ap * l_var / S + l_mean[p, g]*sigma_r_noise) / temp_var
                variance_l = l_var * sigma_r_noise / temp_var
                
                ll = np.random.randn()
                ll = np.clip(ll, -3, 3)
                l_temp = ll * np.sqrt(variance_l) + mean_l
                
                post_l1 = np.exp(-(l_temp*1 - mean_l)**2 / (2*variance_l)) * (l_prob[p, g] + 0.1) + 1e-6
                post_l0 = np.exp(-(l_temp*0 - mean_l)**2 / (2*variance_l)) * (1 - l_prob[p, g] + 0.1) + 1e-6
                
                p1 = post_l1 / (post_l1 + post_l0)
                if p1 >= np.random.rand():
                    l[p, g] = l_temp
                    l_state[p, g] = 1
                    temp = temp + l[p, g]*ap
            else:
                l[p, g] = 0
                l_state[p, g] = 0
                
    return l_state, l

def magical_estimation(
    atac_cell_vector, cand_peak_log2count, rna_cell_vector, cand_gene_log2count,
    cand_tf_peak_binding, cand_peak_gene_looping,
    t_a_prior, t_r_prior, t_prior_mean, t_prior_var,
    b_prior, b_mean, b_var, b_prob,
    l_prior, l_mean, l_var, l_prob,
    M, S, P, G, iteration_num,
    dump_weight_history=False,
    use_numba=False
):
    if use_numba and not HAS_NUMBA:
        print("Warning: Numba not found. Falling back to NumPy implementation.")
        use_numba = False

    if use_numba:
        print("Using Numba-accelerated kernels.")
    else:
        print("Using standard NumPy implementation.")

    print("MAGICAL work starts...")
    
    a_sample = cand_peak_log2count
    r_sample = cand_gene_log2count
    
    b = b_prior.copy()
    l = l_prior.copy()
    t_a = t_a_prior.copy()
    t_r = t_r_prior.copy()
    
    # Pre-calculate t_sample if using Numba
    if use_numba:
        t_sample = np.zeros((M, S))
        for s in range(S):
            t_sample[:, s] = np.mean(t_a[:, atac_cell_vector == s], axis=1)
        
        # Optimize for memory access:
        a_sample_opt = np.ascontiguousarray(a_sample.T) # (S, P)
        r_sample_opt = np.ascontiguousarray(r_sample.T) # (S, G)
        
        # Transposed weights for row-contiguous TF/Gene access
        b_T = np.ascontiguousarray(b.T) # (M, P)
        l_T = np.ascontiguousarray(l.T) # (G, P)
        
        # Transposed priors
        b_mean_T = np.ascontiguousarray(b_mean.T)
        b_prob_T = np.ascontiguousarray(b_prob.T)
        l_mean_T = np.ascontiguousarray(l_mean.T)
        l_prob_T = np.ascontiguousarray(l_prob.T)
        
        t_sample = np.ascontiguousarray(t_sample)
        
        # Sign consistency accumulators
        b_pos_count = np.zeros((M, P))
        l_pos_count = np.zeros((G, P))
    else:
        a_sample_opt = a_sample
        r_sample_opt = r_sample

    # States
    if hasattr(cand_tf_peak_binding, 'toarray'):
        b_state = cand_tf_peak_binding.toarray().astype(float)
    else:
        b_state = cand_tf_peak_binding.astype(float).copy()
        
    if hasattr(cand_peak_gene_looping, 'toarray'):
        l_state = cand_peak_gene_looping.toarray().astype(float)
    else:
        l_state = cand_peak_gene_looping.astype(float).copy()
    
    if use_numba:
        b_state_T = np.ascontiguousarray(b_state.T) # (M, P)
        l_state_T = np.ascontiguousarray(l_state.T) # (G, P)
        
    b_state_frq = np.zeros_like(b_state)
    l_state_frq = np.zeros_like(l_state)
    
    b_weight_sum = np.zeros_like(b_prior)
    l_weight_sum = np.zeros_like(l_prior)
    
    if dump_weight_history:
        b_history = np.zeros((iteration_num, P, M), dtype=np.float32)
        l_history = np.zeros((iteration_num, P, G), dtype=np.float32)
    else:
        b_history = None
        l_history = None
    
    alpha_a = 1.0
    beta_a = 1.0
    sigma_a_noise = np.var(a_sample - b_prior @ t_prior_mean, ddof=1)
    
    alpha_r = 1.0
    beta_r = 1.0
    sigma_r_noise = np.var(r_sample - l_prior.T @ (b_prior @ t_prior_mean), ddof=1)
    
    iteration_seg = max(1, iteration_num // 10)
    
    for i in range(iteration_num):
        if use_numba:
            # Step 1: TF activity
            tf_index = np.random.permutation(M)
            t_sample = sample_t_kernel(a_sample_opt, b_T, t_sample, t_prior_mean, t_prior_var, sigma_a_noise, tf_index, M, S, P)
            t_a, t_r = update_ta_tr_kernel(t_a, t_r, atac_cell_vector, rna_cell_vector, t_sample, t_prior_var[0], M, S)
            
            # Step 2: TF-peak binding weights
            tf_index = np.random.permutation(M)
            b_T, b_pos_count = sample_b_weight_kernel(a_sample_opt, b_T, t_sample, b_state_T, b_mean_T, b_var, sigma_a_noise, tf_index, M, S, P, b_pos_count)
            
            # Step 3: TF-peak binary states
            tf_index = np.random.permutation(M)
            b_state_T, b_T = sample_b_state_kernel(a_sample_opt, b_T, t_sample, b_state_T, b_mean_T, b_var, b_prob_T, sigma_a_noise, tf_index, P, S, M)
            
            # Update b and b_state from T for RSS and summary
            b = b_T.T
            b_state = b_state_T.T
            
            # Step 4: ATAC variance control
            rss_a = np.sum((a_sample_opt - t_sample.T @ b_T)**2)
            scale_a = 1.0 / (beta_a + rss_a / (2*P*S))
            sigma_a_noise = 1.0 / np.random.gamma(shape=alpha_a + 0.5, scale=scale_a)
            
            # Step 5: Peak-Gene looping weights
            a_estimate = b @ t_sample 
            p_index = np.random.permutation(P)
            l_T, l_pos_count = sample_l_weight_kernel(r_sample_opt, l_T, a_estimate, l_state_T, l_mean_T, l_var, sigma_r_noise, p_index, M, S, P, G, l_pos_count)
            
            # Step 6: Peak-Gene looping binary state
            ap_dot_ap_arr = np.sum(a_estimate**2, axis=1)
            l_state_T, l_T = sample_l_state_kernel(r_sample_opt, l_T, a_estimate, l_state_T, l_mean_T, l_var, l_prob_T, sigma_r_noise, P, S, G, ap_dot_ap_arr)
            
            # Update l and l_state
            l = l_T.T
            l_state = l_state_T.T
            
            # Step 7: RNA variance control
            rss_r = np.sum((r_sample_opt - a_estimate.T @ l_T.T)**2)
            scale_r = 1.0 / (beta_r + rss_r / (2*G*S))
            sigma_r_noise = 1.0 / np.random.gamma(shape=alpha_r + 0.5, scale=scale_r)
        else:
            # Step 1: TF activity
            t_a, t_r, t_sample = tf_activity_t_sampling(
                atac_cell_vector, a_sample, rna_cell_vector, r_sample,
                b, t_a, t_r, t_prior_mean, t_prior_var, sigma_a_noise, M, S, P, G
            )
            
            # Step 2: TF-peak binding
            b = tf_peak_binding_b_sampling(
                atac_cell_vector, a_sample, b, t_a, b_state, b_mean, b_var, sigma_a_noise, M, S, P
            )
            
            # Step 3: TF-peak binary state
            b_state, b = tf_peak_binary_binding_b_state_sampling(
                atac_cell_vector, a_sample, b, t_a, b_state, b_mean, b_var, b_prob, sigma_a_noise, M, S, P
            )
            
            # Step 4: ATAC variance control
            rss_a = np.sum((a_sample - b @ t_sample)**2)
            scale_a = 1.0 / (beta_a + rss_a / (2*P*S))
            sigma_a_noise = 1.0 / np.random.gamma(shape=alpha_a + 0.5, scale=scale_a)
            
            # Step 5: Peak-Gene looping
            l = peak_gene_looping_l_sampling(
                rna_cell_vector, r_sample, l, b, t_r, l_state, l_mean, l_var, sigma_r_noise, M, S, P, G
            )
            
            # Step 6: Peak-Gene looping binary state
            l_state, l = peak_gene_binary_looping_l_state_sampling(
                rna_cell_vector, r_sample, l, b, t_r, l_state, l_mean, l_var, l_prob, sigma_r_noise, M, S, P, G
            )
            
            # Step 7: RNA variance control
            rss_r = np.sum((r_sample - l.T @ (b @ t_sample))**2)
            scale_r = 1.0 / (beta_r + rss_r / (2*G*S))
            sigma_r_noise = 1.0 / np.random.gamma(shape=alpha_r + 0.5, scale=scale_r)
        
        # Summary
        b_state_frq += b_state
        l_state_frq += l_state
        
        b_weight_sum += b
        l_weight_sum += l
        
        if dump_weight_history:
            b_history[i, :, :] = b
            l_history[i, :, :] = l
        
        if (i + 1) % iteration_seg == 0:
            print(f"MAGICAL finished {int(100 * (i + 1) / iteration_num)} percent")
            
    cand_tf_peak_binding_prob = b_state_frq / iteration_num
    cand_peak_gene_looping_prob = l_state_frq / iteration_num
    
    cand_tf_peak_binding_weight = b_weight_sum / iteration_num
    cand_peak_gene_looping_weight = l_weight_sum / iteration_num
    
    if use_numba:
        b_pos_prob = b_pos_count.T / iteration_num
        l_pos_prob = l_pos_count.T / iteration_num
    else:
        # For numpy, return zeros for now
        b_pos_prob = np.zeros_like(cand_tf_peak_binding_prob)
        l_pos_prob = np.zeros_like(cand_peak_gene_looping_prob)
    
    return cand_tf_peak_binding_prob, cand_peak_gene_looping_prob, cand_tf_peak_binding_weight, cand_peak_gene_looping_weight, b_history, l_history, b_pos_prob, l_pos_prob
