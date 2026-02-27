import os
import hashlib
import pandas as pd
import numpy as np
from scipy import sparse

def _get_cache_prefix(source_files, prefix=""):
    """Generate cache paths based on the hash of the source file paths."""
    data_dir = os.path.dirname(os.path.abspath(source_files[0]))
    cache_dir = os.path.join(data_dir, ".magical_cache")
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir, exist_ok=True)
        
    mtimes = [os.path.getmtime(f) for f in source_files]
    max_mtime = max(mtimes)
    
    path_str = "|".join(source_files)
    path_hash = hashlib.md5(path_str.encode('utf-8')).hexdigest()
    
    return os.path.join(cache_dir, f"{prefix}_{path_hash}"), max_mtime

def _is_cache_valid(cache_file, source_mtime):
    if os.path.exists(cache_file):
        return os.path.getmtime(cache_file) >= source_mtime
    return False

def parse_chr_to_num(chr_series):
    """Convert chr string (e.g. 'chr1') to number. Exclude X/Y by returning 0."""
    chr_num = np.zeros(len(chr_series), dtype=int)
    for i in range(1, 23):  # up to 22 for Rhesus Macaque/Human autosomes
        chr_num[chr_series == f'chr{i}'] = i
    return chr_num

def load_candidate_genes(filepath):
    """Load candidate gene symbols."""
    cache_prefix, src_mtime = _get_cache_prefix([filepath], "cand_genes")
    cache_path = f"{cache_prefix}.parquet"
    
    if _is_cache_valid(cache_path, src_mtime):
        print(f"Using cached file for candidate genes: {os.path.basename(cache_path)}")
        return pd.read_parquet(cache_path)['gene_symbol'].values
        
    df = pd.read_csv(filepath, header=None, names=['gene_symbol'])
    df.to_parquet(cache_path, engine='pyarrow', index=False)
    return df['gene_symbol'].values

def load_candidate_peaks(filepath):
    """Load candidate peaks."""
    cache_prefix, src_mtime = _get_cache_prefix([filepath], "cand_peaks")
    cache_path = f"{cache_prefix}.parquet"
    
    if _is_cache_valid(cache_path, src_mtime):
        print(f"Using cached file for candidate peaks: {os.path.basename(cache_path)}")
        return pd.read_parquet(cache_path)
        
    df = pd.read_csv(filepath, sep=r'\s+', header=None, names=['chr', 'point1', 'point2'])
    df['chr_num'] = parse_chr_to_num(df['chr'])
    df.to_parquet(cache_path, engine='pyarrow', index=False)
    return df

def load_scrna_data(counts_file, genes_file, meta_file):
    """Load scRNA genes, cells, and sparse count matrix."""
    cache_prefix, src_mtime = _get_cache_prefix([counts_file, genes_file, meta_file], "scrna")
    genes_cache = f"{cache_prefix}_genes.parquet"
    cells_cache = f"{cache_prefix}_cells.parquet"
    counts_cache = f"{cache_prefix}_counts.npz"
    
    if _is_cache_valid(genes_cache, src_mtime) and _is_cache_valid(counts_cache, src_mtime):
        print(f"Using cached files for scRNA data: {os.path.basename(cache_prefix)}_*")
        genes_df = pd.read_parquet(genes_cache)
        cells_df = pd.read_parquet(cells_cache)
        count_matrix = sparse.load_npz(counts_cache)
        return genes_df, cells_df, count_matrix
        
    genes_df = pd.read_csv(genes_file, sep=r'\s+', header=None, names=['gene_index', 'gene_symbol'])
    cells_df = pd.read_csv(meta_file, sep='\t', header=None, 
                           names=['cell_index', 'cell_barcode', 'cell_type', 'subject_ID', 'condition'])
    
    counts_df = pd.read_csv(counts_file, sep=r'\s+', header=None, names=['gene_index', 'cell_index', 'readcount'])
    row = counts_df['gene_index'].values - 1
    col = counts_df['cell_index'].values - 1
    data = counts_df['readcount'].values
    
    count_matrix = sparse.coo_matrix((data, (row, col)), shape=(len(genes_df), len(cells_df))).tocsr()
    
    genes_df.to_parquet(genes_cache, engine='pyarrow', index=False)
    cells_df.to_parquet(cells_cache, engine='pyarrow', index=False)
    sparse.save_npz(counts_cache, count_matrix)
    
    return genes_df, cells_df, count_matrix

def load_scatac_data(counts_file, peaks_file, meta_file):
    """Load scATAC peaks, cells, and sparse count matrix."""
    cache_prefix, src_mtime = _get_cache_prefix([counts_file, peaks_file, meta_file], "scatac")
    peaks_cache = f"{cache_prefix}_peaks.parquet"
    cells_cache = f"{cache_prefix}_cells.parquet"
    counts_cache = f"{cache_prefix}_counts.npz"
    
    if _is_cache_valid(peaks_cache, src_mtime) and _is_cache_valid(counts_cache, src_mtime):
        print(f"Using cached files for scATAC data: {os.path.basename(cache_prefix)}_*")
        peaks_df = pd.read_parquet(peaks_cache)
        cells_df = pd.read_parquet(cells_cache)
        count_matrix = sparse.load_npz(counts_cache)
        return peaks_df, cells_df, count_matrix

    peaks_df = pd.read_csv(peaks_file, sep=r'\s+', header=None, 
                           names=['peak_index', 'chr', 'point1', 'point2'])
    peaks_df['chr_num'] = parse_chr_to_num(peaks_df['chr'])
    
    cells_df = pd.read_csv(meta_file, sep='\t', header=None, 
                           names=['cell_index', 'cell_barcode', 'cell_type', 'subject_ID', 'condition'])
    
    counts_df = pd.read_csv(counts_file, sep=r'\s+', header=None, names=['peak_index', 'cell_index', 'readcount'])
    row = counts_df['peak_index'].values - 1
    col = counts_df['cell_index'].values - 1
    data = counts_df['readcount'].values
    
    count_matrix = sparse.coo_matrix((data, (row, col)), shape=(len(peaks_df), len(cells_df))).tocsr()
    
    peaks_df.to_parquet(peaks_cache, engine='pyarrow', index=False)
    cells_df.to_parquet(cells_cache, engine='pyarrow', index=False)
    sparse.save_npz(counts_cache, count_matrix)
    
    return peaks_df, cells_df, count_matrix

def load_motif_prior(name_file, mapping_file, num_peaks):
    """Load motif names and TF-peak binding prior sparse matrix."""
    cache_prefix, src_mtime = _get_cache_prefix([name_file, mapping_file], "motif")
    motifs_cache = f"{cache_prefix}_motifs.parquet"
    mapping_cache = f"{cache_prefix}_mapping.npz"
    
    if _is_cache_valid(motifs_cache, src_mtime) and _is_cache_valid(mapping_cache, src_mtime):
        print(f"Using cached files for motif prior: {os.path.basename(cache_prefix)}_*")
        motifs_df = pd.read_parquet(motifs_cache)
        tf_peak_binding_matrix = sparse.load_npz(mapping_cache)
        return motifs_df, tf_peak_binding_matrix
        
    motifs_df = pd.read_csv(name_file, sep=r'\s+', header=None, names=['motif_index', 'name'])
    mapping_df = pd.read_csv(mapping_file, sep=r'\s+', header=None, names=['peak_index', 'motif_index', 'flag'])
    
    row = mapping_df['peak_index'].values - 1
    col = mapping_df['motif_index'].values - 1
    data = mapping_df['flag'].values
    
    tf_peak_binding_matrix = sparse.coo_matrix((data, (row, col)), shape=(num_peaks, len(motifs_df))).tocsr()
    
    motifs_df.to_parquet(motifs_cache, engine='pyarrow', index=False)
    sparse.save_npz(mapping_cache, tf_peak_binding_matrix)
    
    return motifs_df, tf_peak_binding_matrix

def load_tad_regions(filepath):
    """Load TAD regions."""
    cache_prefix, src_mtime = _get_cache_prefix([filepath], "tad")
    cache_path = f"{cache_prefix}.parquet"
    
    if _is_cache_valid(cache_path, src_mtime):
        print(f"Using cached file for TAD regions: {os.path.basename(cache_path)}")
        return pd.read_parquet(cache_path)
        
    df = pd.read_csv(filepath, sep=r'\s+', header=None, names=['chr', 'left_boundary', 'right_boundary'])
    df['chr_num'] = parse_chr_to_num(df['chr'])
    df.to_parquet(cache_path, engine='pyarrow', index=False)
    return df

def load_refseq(filepath):
    """Load Refseq info."""
    cache_prefix, src_mtime = _get_cache_prefix([filepath], "refseq")
    cache_path = f"{cache_prefix}.parquet"
    
    if _is_cache_valid(cache_path, src_mtime):
        print(f"Using cached file for Refseq: {os.path.basename(cache_path)}")
        return pd.read_parquet(cache_path)
        
    df = pd.read_csv(filepath, sep=r'\s+', header=0, 
                     names=['chr', 'strand', 'start', 'end', 'gene_name'])
    df['chr_num'] = parse_chr_to_num(df['chr'])
    df.to_parquet(cache_path, engine='pyarrow', index=False)
    return df
