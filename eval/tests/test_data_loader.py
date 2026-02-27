import os
import time
from pymagical.data_loader import (
    load_candidate_genes, load_candidate_peaks, load_scrna_data, 
    load_scatac_data, load_motif_prior, load_tad_regions, load_refseq
)

def test_data_loader():
    test_dir = "/mnt/ceph/users/agebrain/anderson/snmulti_data/pymagical/test_data"
    astrocytes_dir = os.path.join(test_dir, "astrocytes")
    
    print("Loading candidate genes...")
    cand_genes = load_candidate_genes(os.path.join(astrocytes_dir, "sig_cr_genes.txt"))
    print(f"Loaded {len(cand_genes)} candidate genes.")
    
    print("Loading candidate peaks...")
    cand_peaks = load_candidate_peaks(os.path.join(astrocytes_dir, "sig_cr_peaks.txt"))
    print(f"Loaded {len(cand_peaks)} candidate peaks.")
    
    print("Loading scRNA data...")
    t0 = time.time()
    rna_genes, rna_cells, rna_counts = load_scrna_data(
        os.path.join(astrocytes_dir, "rna_counts.txt"),
        os.path.join(astrocytes_dir, "rna_genes.txt"),
        os.path.join(astrocytes_dir, "rna_meta.txt")
    )
    print(f"Loaded scRNA data in {time.time()-t0:.2f}s. Matrix shape: {rna_counts.shape}")
    
    print("Loading scATAC data...")
    t0 = time.time()
    atac_peaks, atac_cells, atac_counts = load_scatac_data(
        os.path.join(astrocytes_dir, "atac_counts.txt"),
        os.path.join(astrocytes_dir, "atac_peaks.txt"),
        os.path.join(astrocytes_dir, "atac_meta.txt")
    )
    print(f"Loaded scATAC data in {time.time()-t0:.2f}s. Matrix shape: {atac_counts.shape}")
    
    print("Loading Motif data...")
    t0 = time.time()
    motifs, motif_prior = load_motif_prior(
        os.path.join(test_dir, "motif_info.txt"),
        os.path.join(test_dir, "motif_prior.txt"),
        len(atac_peaks)
    )
    print(f"Loaded Motif prior in {time.time()-t0:.2f}s. Matrix shape: {motif_prior.shape}")
    
    print("Loading TAD regions...")
    tads = load_tad_regions(os.path.join(test_dir, "tad_regions.txt"))
    print(f"Loaded {len(tads)} TAD regions.")
    
    print("Loading RefSeq...")
    refseq = load_refseq(os.path.join(test_dir, "rhemac10_refseq.txt"))
    print(f"Loaded {len(refseq)} RefSeq genes.")
    
if __name__ == "__main__":
    test_data_loader()
