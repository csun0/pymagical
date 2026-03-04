import pytest
import numpy as np
from pymagical.estimation import (
    peak_gene_looping_l_sampling,
    tf_peak_binding_b_sampling,
    tf_activity_t_sampling
)

@pytest.fixture
def sampling_inputs(small_matrices):
    m = small_matrices
    N_peaks, N_genes = m['gene_peak_prior'].shape
    _, N_tfs = m['motif_prior'].shape
    N_cells = m['rna'].shape[1]
    S = 2 # 2 pseudo-samples
    
    # Setup vectors and samples
    atac_cell_vector = np.repeat(np.arange(S), N_cells // S)
    rna_cell_vector = np.repeat(np.arange(S), N_cells // S)
    
    # Priors
    t_mean = np.zeros((N_tfs, S))
    t_var = np.ones(N_tfs)
    b_state = m['motif_prior'].astype(float)
    b_mean = np.zeros((N_peaks, N_tfs))
    b_var = np.ones(N_tfs)
    l_state = m['gene_peak_prior'].astype(float)
    l_mean = np.zeros((N_peaks, N_genes))
    l_var = 1.0
    
    # Weights
    b = np.zeros((N_peaks, N_tfs))
    l = np.zeros((N_peaks, N_genes))
    t_a = np.random.randn(N_tfs, N_cells).astype(np.float32)
    t_r = np.random.randn(N_tfs, N_cells).astype(np.float32)
    
    return {
        'atac_cell_vector': atac_cell_vector,
        'rna_cell_vector': rna_cell_vector,
        'a_sample': m['atac'][:, :S], # use first S columns as pseudo-bulk
        'r_sample': m['rna'][:, :S],
        'b': b, 'l': l, 't_a': t_a, 't_r': t_r,
        'b_state': b_state, 'b_mean': b_mean, 'b_var': b_var,
        'l_state': l_state, 'l_mean': l_mean, 'l_var': l_var,
        't_mean': t_mean, 't_var': t_var,
        'M': N_tfs, 'S': S, 'P': N_peaks, 'G': N_genes,
        'sigma_a_noise': 1.0, 'sigma_r_noise': 1.0
    }

def test_peak_gene_looping_l_sampling(sampling_inputs):
    i = sampling_inputs
    l_new = peak_gene_looping_l_sampling(
        i['rna_cell_vector'], i['r_sample'], i['l'], i['b'], i['t_r'],
        i['l_state'], i['l_mean'], i['l_var'], i['sigma_r_noise'],
        i['M'], i['S'], i['P'], i['G']
    )
    assert l_new.shape == (i['P'], i['G'])

def test_tf_peak_binding_b_sampling(sampling_inputs):
    i = sampling_inputs
    b_new = tf_peak_binding_b_sampling(
        i['atac_cell_vector'], i['a_sample'], i['b'], i['t_a'],
        i['b_state'], i['b_mean'], i['b_var'], i['sigma_a_noise'],
        i['M'], i['S'], i['P']
    )
    assert b_new.shape == (i['P'], i['M'])

def test_tf_activity_t_sampling(sampling_inputs):
    i = sampling_inputs
    ta, tr, t_sample = tf_activity_t_sampling(
        i['atac_cell_vector'], i['a_sample'], i['rna_cell_vector'], i['r_sample'],
        i['b'], i['t_a'], i['t_r'], i['t_mean'], i['t_var'], i['sigma_a_noise'],
        i['M'], i['S'], i['P'], i['G']
    )
    assert t_sample.shape == (i['M'], i['S'])
