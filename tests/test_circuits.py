import pytest
import pandas as pd
import numpy as np
from pymagical.circuits import (
    construct_candidate_circuits_with_tad, construct_candidate_circuits_without_tad
)

def test_construct_candidate_circuits_with_tad_basic():
    # This function is complex and requires many inputs.
    # We will test if we can call it with minimal mock data.
    
    # Mock data setup
    common_samples = ["Sample1", "Sample2"]
    cand_genes = ["GeneA", "GeneB"]
    
    rna_genes = pd.DataFrame({'gene_symbol': ["GeneA", "GeneB"]})
    rna_cells = pd.DataFrame({
        'subject_ID': ["Sample1", "Sample2"],
        'cell_ID': ["C1", "C2"]
    })
    rna_count_matrix = np.random.randint(0, 10, size=(2, 2))
    
    atac_cells = pd.DataFrame({
        'subject_ID': ["Sample1", "Sample2"],
        'cell_ID': ["C1", "C2"]
    })
    
    motifs = pd.DataFrame({'name': ["TF1", "TF2"]})
    
    # Actually, let's just mock enough peaks to pass internal filters.
    N_test_peaks = 40
    cand_peaks = pd.DataFrame({
        'chr_num': [1] * N_test_peaks,
        'point1': np.arange(N_test_peaks) * 100,
        'point2': np.arange(N_test_peaks) * 100 + 50
    })
    atac_peaks = cand_peaks.copy()
    atac_peaks['peak_index'] = np.arange(N_test_peaks) + 1
    
    tf_peak_binding_matrix = np.ones((N_test_peaks, 2)) # All TFs bind all peaks
    atac_count_matrix = np.random.randint(1, 10, size=(N_test_peaks, 2))
    
    refseq = pd.DataFrame({
        'gene_name': ["GeneA", "GeneB"],
        'strand': ["+", "-"],
        'chr_num': [1, 1],
        'start': [1000, 5000],
        'end': [1100, 5100]
    })
    
    tad_regions = pd.DataFrame({
        'chr_num': [1],
        'left_boundary': [0],
        'right_boundary': [100000] # Large enough to cover our mock coordinates
    })
    
    results = construct_candidate_circuits_with_tad(
        common_samples,
        cand_genes, cand_peaks, 
        rna_genes, rna_cells, rna_count_matrix,
        atac_peaks, atac_cells, atac_count_matrix,
        motifs, tf_peak_binding_matrix,
        refseq, tad_regions
    )
    
    # If successful, it returns a tuple of many items
    assert results is not None
    assert len(results) == 12

def test_construct_candidate_circuits_without_tad_uses_distance_only():
    # Peaks every 100 bp on chr1 (centers 25..3925). GeneA TSS = 1000 (+ strand start),
    # GeneB TSS = 5100 (- strand end). With a 1 kb cutoff only GeneA has peaks in range.
    np.random.seed(0)
    common_samples = ["Sample1", "Sample2"]
    n_peaks = 40
    cand_peaks = pd.DataFrame({
        'chr_num': [1] * n_peaks,
        'point1': np.arange(n_peaks) * 100,
        'point2': np.arange(n_peaks) * 100 + 50
    })
    atac_peaks = cand_peaks.copy()
    atac_peaks['peak_index'] = np.arange(n_peaks) + 1
    cells = pd.DataFrame({'subject_ID': ["Sample1", "Sample2"]})
    refseq = pd.DataFrame({
        'gene_name': ["GeneA", "GeneB"], 'strand': ["+", "-"], 'chr_num': [1, 1],
        'start': [1000, 5000], 'end': [1100, 5100]
    })

    results = construct_candidate_circuits_without_tad(
        common_samples,
        ["GeneA", "GeneB"], cand_peaks,
        pd.DataFrame({'gene_symbol': ["GeneA", "GeneB"]}), cells, np.random.randint(1, 10, size=(2, 2)),
        atac_peaks, cells, np.random.randint(1, 10, size=(n_peaks, 2)),
        pd.DataFrame({'name': ["GeneA", "GeneB"]}), np.ones((n_peaks, 2)),
        refseq, distance_control=1000
    )

    peaks, genes, looping = results[2], results[4], results[7]
    assert list(genes['symbols']) == ["GeneA"]
    centers = (peaks['point1'] + peaks['point2']).values / 2
    assert np.all(np.abs(centers - 1000) < 1000)
    assert len(peaks) == 20  # centers 25..1925
    assert looping.shape == (20, 1) and looping.all()
