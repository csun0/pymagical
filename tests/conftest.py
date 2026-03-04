import pytest
import numpy as np
import pandas as pd
import tempfile
import os

@pytest.fixture
def mock_data_dir():
    with tempfile.TemporaryDirectory() as tmp_dir:
        # Create small synthetic files for testing data_loader
        
        # 1. Candidate genes (No header)
        genes = ["GeneA", "GeneB", "GeneC"]
        with open(os.path.join(tmp_dir, "sig_cr_genes.txt"), "w") as f:
            f.write("\n".join(genes))
            
        # 2. Candidate peaks (4 columns, no header: peak_idx, chr, p1, p2)
        peaks = ["1\tchr1\t100\t200", "2\tchr1\t300\t400", "3\tchr1\t500\t600"]
        with open(os.path.join(tmp_dir, "sig_cr_peaks.txt"), "w") as f:
            f.write("\n".join(peaks))
            
        # 3. scRNA data
        # rna_genes: gene_index, gene_symbol
        with open(os.path.join(tmp_dir, "rna_genes.txt"), "w") as f:
            for i, g in enumerate(genes):
                f.write(f"{i+1}\t{g}\n")
        
        # rna_meta: cell_index, cell_barcode, cell_type, subject_ID, condition
        cells = ["C1", "C2", "C3", "C4", "C5"]
        with open(os.path.join(tmp_dir, "rna_meta.txt"), "w") as f:
            for i, c in enumerate(cells):
                f.write(f"{i+1}\t{c}\tType1\tS1\tCtrl\n")
        
        # rna_counts: gene_index, cell_index, readcount (COOR format)
        with open(os.path.join(tmp_dir, "rna_counts.txt"), "w") as f:
            for g_idx in range(1, 4):
                for c_idx in range(1, 6):
                    f.write(f"{g_idx}\t{c_idx}\t{np.random.randint(1, 10)}\n")
                    
        # 4. scATAC data
        # atac_peaks: peak_index, chr, point1, point2
        with open(os.path.join(tmp_dir, "atac_peaks.txt"), "w") as f:
            f.write("1\tchr1\t100\t200\n")
            f.write("2\tchr1\t300\t400\n")
            f.write("3\tchr1\t500\t600\n")
            
        # atac_meta: cell_index, cell_barcode, cell_type, subject_ID, condition
        with open(os.path.join(tmp_dir, "atac_meta.txt"), "w") as f:
            for i, c in enumerate(cells):
                f.write(f"{i+1}\t{c}\tType1\tS1\tCtrl\n")
                
        # atac_counts: peak_index, cell_index, readcount
        with open(os.path.join(tmp_dir, "atac_counts.txt"), "w") as f:
            for p_idx in range(1, 4):
                for c_idx in range(1, 6):
                    f.write(f"{p_idx}\t{c_idx}\t{np.random.randint(1, 10)}\n")
                
        # 5. Motif data
        # motif_info: motif_index, name
        tfs = ["TF1", "TF2"]
        with open(os.path.join(tmp_dir, "motif_info.txt"), "w") as f:
            f.write("1\tTF1\n2\tTF2\n")
        
        # motif_prior: peak_index, motif_index, flag
        with open(os.path.join(tmp_dir, "motif_prior.txt"), "w") as f:
            f.write("1\t1\t1\n1\t2\t0\n2\t1\t0\n2\t2\t1\n3\t1\t1\n3\t2\t1\n")
        
        # 6. TAD regions: chr, left, right
        with open(os.path.join(tmp_dir, "tad_regions.txt"), "w") as f:
            f.write("chr1\t0\t1000\n")
            
        # 7. RefSeq: chr, strand, start, end, gene_name
        with open(os.path.join(tmp_dir, "refseq.txt"), "w") as f:
            f.write("chr1\t+\t100\t200\tGeneA\n")
            f.write("chr1\t-\t300\t400\tGeneB\n")

        yield tmp_dir

@pytest.fixture
def small_matrices():
    """Generates small matrices for testing estimation logic."""
    np.random.seed(42)
    N_peaks = 10
    N_genes = 5
    N_tfs = 3
    N_cells = 20
    
    return {
        'rna': np.random.randn(N_genes, N_cells).astype(np.float32),
        'atac': np.random.randn(N_peaks, N_cells).astype(np.float32),
        'motif_prior': np.random.randint(0, 2, size=(N_peaks, N_tfs)).astype(np.float32),
        'gene_peak_prior': np.random.randint(0, 2, size=(N_peaks, N_genes)).astype(np.float32)
    }
