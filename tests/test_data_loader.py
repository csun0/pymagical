import os
import pytest
import numpy as np
import pandas as pd
from pymagical.data_loader import (
    load_candidate_genes, load_candidate_peaks, load_scrna_data, 
    load_scatac_data, load_motif_prior, load_tad_regions, load_refseq
)

def test_load_candidate_genes(mock_data_dir):
    gene_file = os.path.join(mock_data_dir, "sig_cr_genes.txt")
    genes = load_candidate_genes(gene_file)
    assert len(genes) == 3
    assert "GeneA" in genes

def test_load_candidate_peaks(mock_data_dir):
    peak_file = os.path.join(mock_data_dir, "sig_cr_peaks.txt")
    peaks_df = load_candidate_peaks(peak_file)
    assert len(peaks_df) == 3
    # The loader returns a DataFrame with 'chr', 'point1', 'point2'
    assert peaks_df.iloc[0]['chr'] == "chr1"

def test_load_scrna_data(mock_data_dir):
    counts_file = os.path.join(mock_data_dir, "rna_counts.txt")
    genes_file = os.path.join(mock_data_dir, "rna_genes.txt")
    meta_file = os.path.join(mock_data_dir, "rna_meta.txt")
    
    genes, cells, counts = load_scrna_data(counts_file, genes_file, meta_file)
    assert len(genes) == 3
    assert len(cells) == 5
    assert counts.shape == (3, 5)

def test_load_scatac_data(mock_data_dir):
    counts_file = os.path.join(mock_data_dir, "atac_counts.txt")
    peaks_file = os.path.join(mock_data_dir, "atac_peaks.txt")
    meta_file = os.path.join(mock_data_dir, "atac_meta.txt")
    
    peaks, cells, counts = load_scatac_data(counts_file, peaks_file, meta_file)
    assert len(peaks) == 3
    assert len(cells) == 5
    assert counts.shape == (3, 5)

def test_load_motif_prior(mock_data_dir):
    info_file = os.path.join(mock_data_dir, "motif_info.txt")
    prior_file = os.path.join(mock_data_dir, "motif_prior.txt")
    
    # We created 3 peaks in mock_data_dir
    motifs, prior = load_motif_prior(info_file, prior_file, num_peaks=3)
    assert len(motifs) == 2
    assert prior.shape == (3, 2)

def test_load_tad_regions(mock_data_dir):
    tad_file = os.path.join(mock_data_dir, "tad_regions.txt")
    tads = load_tad_regions(tad_file)
    assert len(tads) == 1
    # DataFrame returns rows as dictionaries with .iloc[0] or similar
    assert tads.iloc[0]['chr'] == "chr1"

def test_load_refseq(mock_data_dir):
    refseq_file = os.path.join(mock_data_dir, "refseq.txt")
    refseq = load_refseq(refseq_file)
    assert len(refseq) == 2
    assert "GeneA" in refseq['gene_name'].values
