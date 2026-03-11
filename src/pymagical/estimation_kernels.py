import numpy as np
from numba import njit, prange

@njit(fastmath=True)
def _clip(x, low, high):
    if x < low: return low
    if x > high: return high
    return x

@njit(fastmath=True, parallel=True)
def sample_b_state_kernel(
    a_sample_T, b_T, t_sample, b_state_T, b_mean_T, b_var, b_prob_T, 
    sigma_a_noise, tf_index, P, S, M
):
    for p in prange(P):
        temp = np.zeros(S)
        for m in range(M):
            if b_T[m, p] != 0.0:
                for s in range(S):
                    temp[s] += t_sample[m, s] * b_T[m, p]
                    
        for m in tf_index:
            tm = t_sample[m, :] 
            tm_dot_tm = np.dot(tm, tm)
            temp_var = (tm_dot_tm * b_var[m] / S) + sigma_a_noise
            
            if b_state_T[m, p] == 1.0:
                dot_val = 0.0
                for s in range(S):
                    dot_val += (a_sample_T[s, p] - temp[s] + b_T[m, p] * tm[s]) * tm[s]
                
                mean_b = (dot_val * b_var[m] / S + b_mean_T[m, p] * sigma_a_noise) / temp_var
                variance_b = b_var[m] * sigma_a_noise / temp_var
                
                # Match MATLAB exactly: likelihood of current weight vs likelihood of 0
                post_b1 = np.exp(-(b_T[m, p] - mean_b)**2 / (2 * variance_b)) * (b_prob_T[m, p] + 0.25) + 1e-6
                post_b0 = np.exp(-(0.0 - mean_b)**2 / (2 * variance_b)) * (1.0 - b_prob_T[m, p] + 0.25) + 1e-6
                
                prob_1 = post_b1 / (post_b1 + post_b0)
                if np.random.random() >= prob_1:
                    for s in range(S):
                        temp[s] -= b_T[m, p] * tm[s]
                    b_T[m, p] = 0.0
                    b_state_T[m, p] = 0.0
                    
            elif b_state_T[m, p] == 0.0 and b_mean_T[m, p] != 0.0:
                dot_val = 0.0
                for s in range(S):
                    dot_val += (a_sample_T[s, p] - temp[s]) * tm[s]
                
                mean_b = (dot_val * b_var[m] / S + b_mean_T[m, p] * sigma_a_noise) / temp_var
                variance_b = b_var[m] * sigma_a_noise / temp_var
                
                b_temp = _clip(np.random.standard_normal(), -3.0, 3.0) * np.sqrt(variance_b) + mean_b
                
                post_b1 = np.exp(-(b_temp - mean_b)**2 / (2 * variance_b)) * (b_prob_T[m, p] + 0.25) + 1e-6
                post_b0 = np.exp(-(0.0 - mean_b)**2 / (2 * variance_b)) * (1.0 - b_prob_T[m, p] + 0.25) + 1e-6
                
                prob_1 = post_b1 / (post_b1 + post_b0)
                if np.random.random() < prob_1:
                    b_T[m, p] = b_temp
                    b_state_T[m, p] = 1.0
                    for s in range(S):
                        temp[s] += b_T[m, p] * tm[s]
    return b_state_T, b_T

@njit(fastmath=True, parallel=True)
def sample_l_state_kernel(
    r_sample_T, l_T, a_estimate, l_state_T, l_mean_T, l_var, l_prob_T,
    sigma_r_noise, P, S, G, ap_dot_ap_arr
):
    for g in prange(G):
        temp = np.zeros(S)
        for p in range(P):
            if l_T[g, p] != 0.0:
                for s in range(S):
                    temp[s] += a_estimate[p, s] * l_T[g, p]
                    
        for p in range(P):
            ap = a_estimate[p, :] 
            temp_var = ap_dot_ap_arr[p] * l_var / S + sigma_r_noise
            
            if l_state_T[g, p] == 1.0:
                dot_val = 0.0
                for s in range(S):
                    dot_val += (r_sample_T[s, g] - temp[s] + l_T[g, p] * ap[s]) * ap[s]
                
                mean_l = (dot_val * l_var / S + l_mean_T[g, p] * sigma_r_noise) / temp_var
                variance_l = l_var * sigma_r_noise / temp_var
                
                post_l1 = np.exp(-(l_T[g, p] - mean_l)**2 / (2 * variance_l)) * (l_prob_T[g, p] + 0.1) + 1e-6
                post_l0 = np.exp(-(0.0 - mean_l)**2 / (2 * variance_l)) * (1.0 - l_prob_T[g, p] + 0.1) + 1e-6
                
                prob_1 = post_l1 / (post_l1 + post_l0)
                if np.random.random() >= prob_1:
                    for s in range(S):
                        temp[s] -= l_T[g, p] * ap[s]
                    l_T[g, p] = 0.0
                    l_state_T[g, p] = 0.0
            elif l_state_T[g, p] == 0.0 and l_mean_T[g, p] != 0.0:
                dot_val = 0.0
                for s in range(S):
                    dot_val += (r_sample_T[s, g] - temp[s]) * ap[s]
                
                mean_l = (dot_val * l_var / S + l_mean_T[g, p] * sigma_r_noise) / temp_var
                variance_l = l_var * sigma_r_noise / temp_var
                
                l_temp = _clip(np.random.standard_normal(), -3.0, 3.0) * np.sqrt(variance_l) + mean_l
                
                post_l1 = np.exp(-(l_temp - mean_l)**2 / (2 * variance_l)) * (l_prob_T[g, p] + 0.1) + 1e-6
                post_l0 = np.exp(-(0.0 - mean_l)**2 / (2 * variance_l)) * (1.0 - l_prob_T[g, p] + 0.1) + 1e-6
                
                prob_1 = post_l1 / (post_l1 + post_l0)
                if np.random.random() < prob_1:
                    l_T[g, p] = l_temp
                    l_state_T[g, p] = 1.0
                    for s in range(S):
                        temp[s] += l_T[g, p] * ap[s]
    return l_state_T, l_T

@njit(fastmath=True)
def sample_b_weight_kernel(
    a_sample_T, b_T, t_sample, b_state_T, b_mean_T, b_var, sigma_a_noise, tf_index, M, S, P,
    b_pos_count
):
    a_resid_T = a_sample_T - t_sample.T @ b_T
    for m in tf_index:
        tm = t_sample[m, :]
        tm_dot_tm = np.dot(tm, tm)
        temp_var = tm_dot_tm * b_var[m] / S + sigma_a_noise
        variance_b = b_var[m] * sigma_a_noise / temp_var
        sqrt_var_b = np.sqrt(variance_b)
        
        for p in range(P):
            if b_state_T[m, p] != 0.0:
                old_bm_p = b_T[m, p]
                dot_val = 0.0
                for s in range(S):
                    dot_val += (a_resid_T[s, p] + old_bm_p * tm[s]) * tm[s]
                
                mean_b = (dot_val * b_var[m] / S + b_mean_T[m, p] * sigma_a_noise) / temp_var
                new_bm_p = _clip(np.random.standard_normal(), -3.0, 3.0) * sqrt_var_b + mean_b
                b_T[m, p] = new_bm_p
                
                if b_T[m, p] > 0:
                    b_pos_count[m, p] += 1.0
                
                diff = old_bm_p - new_bm_p
                if diff != 0.0:
                    for s in range(S):
                        a_resid_T[s, p] += diff * tm[s]
            else:
                b_T[m, p] = 0.0
    return b_T, b_pos_count

@njit(fastmath=True)
def sample_l_weight_kernel(
    r_sample_T, l_T, a_estimate, l_state_T, l_mean_T, l_var, sigma_r_noise, p_index, M, S, P, G,
    l_pos_count
):
    r_resid_T = r_sample_T - a_estimate.T @ l_T.T
    for p in p_index:
        ap = a_estimate[p, :] 
        ap_dot_ap = np.dot(ap, ap)
        temp_var = ap_dot_ap * l_var / S + sigma_r_noise
        curr_var_l = l_var * sigma_r_noise / temp_var
        sqrt_var_l = np.sqrt(curr_var_l)
        
        for g in range(G):
            if l_state_T[g, p] != 0.0:
                old_lp_g = l_T[g, p]
                dot_val = 0.0
                for s in range(S):
                    dot_val += (r_resid_T[s, g] + old_lp_g * ap[s]) * ap[s]
                
                mean_l = (dot_val * l_var / S + l_mean_T[g, p] * sigma_r_noise) / temp_var
                new_lp_g = _clip(np.random.standard_normal(), -3.0, 3.0) * sqrt_var_l + mean_l
                l_T[g, p] = new_lp_g
                
                if l_T[g, p] > 0:
                    l_pos_count[g, p] += 1.0
                
                diff = old_lp_g - new_lp_g
                if diff != 0.0:
                    for s in range(S):
                        r_resid_T[s, g] += diff * ap[s]
            else:
                l_T[g, p] = 0.0
    return l_T, l_pos_count

@njit(fastmath=True)
def sample_t_kernel(
    a_sample_T, b_T, t_sample, t_mean, t_var, sigma_a_noise, tf_index, M, S, P
):
    a_resid_T = a_sample_T - t_sample.T @ b_T
    for m in tf_index:
        bm = b_T[m, :] 
        bm_dot_bm = np.dot(bm, bm)
        temp_var = (bm_dot_bm * t_var[m] / P) + sigma_a_noise
        curr_var_t = t_var[m] * sigma_a_noise / temp_var
        sqrt_var_t = np.sqrt(curr_var_t)
        
        for s in range(S):
            old_tm_s = t_sample[m, s]
            dot_val = 0.0
            for p in range(P):
                dot_val += (a_resid_T[s, p] + bm[p] * old_tm_s) * bm[p]
                
            mean_val = (dot_val / P + t_mean[m, s] * sigma_a_noise) / temp_var
            new_tm_s = _clip(np.random.standard_normal(), -3.0, 3.0) * sqrt_var_t + mean_val
            t_sample[m, s] = new_tm_s
            
            diff = old_tm_s - new_tm_s
            if diff != 0.0:
                for p in range(P):
                    a_resid_T[s, p] += diff * bm[p]
    return t_sample

@njit(fastmath=True, parallel=True)
def update_ta_tr_kernel(t_a, t_r, atac_cell_vector, rna_cell_vector, t_sample, variance_t, M, S):
    v_t = np.sqrt(variance_t)
    for m in prange(M):
        for i in range(atac_cell_vector.shape[0]):
            s = atac_cell_vector[i]
            t_a[m, i] = np.random.standard_normal() * v_t + t_sample[m, s]
        for i in range(rna_cell_vector.shape[0]):
            s = rna_cell_vector[i]
            t_r[m, i] = np.random.standard_normal() * v_t + t_sample[m, s]
    return t_a, t_r
