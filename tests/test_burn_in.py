import pytest
import numpy as np
from pymagical.estimation import magical_estimation

@pytest.fixture
def mock_estimation_inputs():
    M, S, P, G = 2, 2, 5, 3
    return {
        'atac_cell_vector': np.array([0, 1]),
        'cand_peak_log2count': np.random.randn(P, S),
        'rna_cell_vector': np.array([0, 1]),
        'cand_gene_log2count': np.random.randn(G, S),
        'cand_tf_peak_binding': np.ones((P, M)),
        'cand_peak_gene_looping': np.ones((P, G)),
        't_a_prior': np.random.randn(M, S),
        't_r_prior': np.random.randn(M, S),
        't_prior_mean': np.zeros((M, S)),
        't_prior_var': np.ones(M),
        'b_prior': np.zeros((P, M)),
        'b_mean': np.zeros((P, M)),
        'b_var': np.ones(M),
        'b_prob': np.full((P, M), 0.5),
        'l_prior': np.zeros((P, G)),
        'l_mean': np.zeros((P, G)),
        'l_var': 1.0,
        'l_prob': np.full((P, G), 0.5),
        'M': M, 'S': S, 'P': P, 'G': G
    }

def test_burn_in_logic(mock_estimation_inputs):
    i = mock_estimation_inputs
    iteration_num = 10
    burn_in = 5
    
    # We want to check if the outputs are averaged correctly (over iteration_num - burn_in)
    # We also want to check if it doesn't crash
    res = magical_estimation(
        **i, iteration_num=iteration_num, burn_in=burn_in, use_numba=False
    )
    
    b_prob, l_prob, b_weight, l_weight, b_history, l_history, b_pos_prob, l_pos_prob = res
    
    # Probabilities should be between 0 and 1
    assert np.all(b_prob >= 0) and np.all(b_prob <= 1)
    assert np.all(l_prob >= 0) and np.all(l_prob <= 1)
    
    # History should have all iterations
    if b_history is not None:
        assert b_history.shape[0] == iteration_num

def test_zero_burn_in(mock_estimation_inputs):
    i = mock_estimation_inputs
    res = magical_estimation(
        **i, iteration_num=5, burn_in=0, use_numba=False
    )
    b_prob, l_prob, b_weight, l_weight, b_history, l_history, b_pos_prob, l_pos_prob = res
    assert b_prob.shape == (i['P'], i['M'])

def test_invalid_burn_in(mock_estimation_inputs):
    i = mock_estimation_inputs
    with pytest.raises(ValueError, match="Burn-in period must be less than total iterations"):
        magical_estimation(
            **i, iteration_num=5, burn_in=5, use_numba=False
        )

def test_numba_burn_in(mock_estimation_inputs):
    # Check if numba is available
    try:
        import numba
        HAS_NUMBA = True
    except ImportError:
        HAS_NUMBA = False
        
    if not HAS_NUMBA:
        pytest.skip("Numba not installed")
        
    i = mock_estimation_inputs
    iteration_num = 10
    burn_in = 5
    
    # Run with Numba
    res = magical_estimation(
        **i, iteration_num=iteration_num, burn_in=burn_in, use_numba=True
    )
    
    b_prob, l_prob, b_weight, l_weight, b_history, l_history, b_pos_prob, l_pos_prob = res
    
    assert b_prob.shape == (i['P'], i['M'])
    assert np.all(b_pos_prob >= 0) and np.all(b_pos_prob <= 1)
    assert np.all(l_pos_prob >= 0) and np.all(l_pos_prob <= 1)
