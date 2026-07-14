import os
import hashlib
import json
import tempfile
import pandas as pd
import numpy as np
from scipy import sparse

def _get_source_fingerprint(source_files):
    """Generate a fingerprint based on file paths, sizes, and modification times."""
    fingerprint_parts = []
    for f in sorted(source_files):
        abs_path = os.path.abspath(f)
        stat = os.stat(abs_path)
        fingerprint_parts.append(f"{abs_path}|{stat.st_size}|{stat.st_mtime}")
    
    fingerprint_str = "||".join(fingerprint_parts)
    return hashlib.md5(fingerprint_str.encode('utf-8')).hexdigest()

def _get_cache_dir(source_file):
    """Get or create the .magical_cache directory."""
    data_dir = os.path.dirname(os.path.abspath(source_file))
    cache_dir = os.path.join(data_dir, ".magical_cache")
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir, exist_ok=True)
    return cache_dir

def _check_cache_integrity(cache_files, metadata_file, source_fingerprint):
    """Check if all cache files exist, have non-zero size, and metadata matches."""
    if not os.path.exists(metadata_file):
        return False
    
    for f in cache_files:
        if not os.path.exists(f) or os.path.getsize(f) == 0:
            return False
            
    try:
        with open(metadata_file, 'r') as f:
            meta = json.load(f)
        return meta.get('source_fingerprint') == source_fingerprint and meta.get('completed', False)
    except (json.JSONDecodeError, IOError, ValueError):
        return False

def _atomic_save_metadata(metadata_file, source_fingerprint):
    """Save metadata file atomically."""
    meta = {
        'source_fingerprint': source_fingerprint,
        'completed': True,
        'timestamp': pd.Timestamp.now().isoformat()
    }
    # Use a safer temp file approach
    target_dir = os.path.dirname(metadata_file)
    with tempfile.NamedTemporaryFile(mode='w', dir=target_dir, delete=False, suffix=".tmp") as f:
        json.dump(meta, f)
        tmp_path = f.name
    
    try:
        os.replace(tmp_path, metadata_file)
    except Exception as e:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise e

def _atomic_save_parquet(df, cache_path, **kwargs):
    """Save a DataFrame to Parquet atomically."""
    target_dir = os.path.dirname(cache_path)
    with tempfile.NamedTemporaryFile(dir=target_dir, delete=False, suffix=".tmp.parquet") as f:
        tmp_path = f.name
    
    try:
        df.to_parquet(tmp_path, engine='pyarrow', index=False, **kwargs)
        if os.path.getsize(tmp_path) == 0:
            raise IOError(f"Failed to write data to {tmp_path} (size 0)")
        os.replace(tmp_path, cache_path)
    except Exception as e:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise e

def _atomic_save_npz(matrix, cache_path):
    """Save a sparse matrix to NPZ atomically."""
    target_dir = os.path.dirname(cache_path)
    # scipy.save_npz adds .npz if not present, so we handle suffix carefully
    with tempfile.NamedTemporaryFile(dir=target_dir, delete=False, suffix=".tmp.npz") as f:
        tmp_path = f.name
    
    try:
        sparse.save_npz(tmp_path, matrix)
        if os.path.getsize(tmp_path) == 0:
            raise IOError(f"Failed to write data to {tmp_path} (size 0)")
        os.replace(tmp_path, cache_path)
    except Exception as e:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise e

def parse_chr_to_num(chr_series):
    """Convert chr string (e.g. 'chr1') to number. Exclude X/Y by returning 0."""
    chr_num = np.zeros(len(chr_series), dtype=int)
    for i in range(1, 23):  # up to 22 for Rhesus Macaque/Human autosomes
        chr_num[chr_series == f'chr{i}'] = i
    return chr_num

def load_candidate_genes(filepath):
    """Load candidate gene symbols."""
    fingerprint = _get_source_fingerprint([filepath])
    cache_dir = _get_cache_dir(filepath)
    cache_path = os.path.join(cache_dir, f"cand_genes_{fingerprint}.parquet")
    meta_path = cache_path + ".meta"
    
    if _check_cache_integrity([cache_path], meta_path, fingerprint):
        print(f"Using cached file for candidate genes: {os.path.basename(cache_path)}")
        return pd.read_parquet(cache_path)['gene_symbol'].values
        
    print(f"Loading candidate genes from {os.path.basename(filepath)} ...")
    df = pd.read_csv(filepath, header=None, names=['gene_symbol'])
    _atomic_save_parquet(df, cache_path)
    _atomic_save_metadata(meta_path, fingerprint)
    return df['gene_symbol'].values

def load_candidate_peaks(filepath):
    """Load candidate peaks."""
    fingerprint = _get_source_fingerprint([filepath])
    cache_dir = _get_cache_dir(filepath)
    cache_path = os.path.join(cache_dir, f"cand_peaks_{fingerprint}.parquet")
    meta_path = cache_path + ".meta"
    
    if _check_cache_integrity([cache_path], meta_path, fingerprint):
        print(f"Using cached file for candidate peaks: {os.path.basename(cache_path)}")
        return pd.read_parquet(cache_path)
        
    print(f"Loading candidate peaks from {os.path.basename(filepath)} ...")
    df = pd.read_csv(filepath, sep=r'\s+', header=None, names=['chr', 'point1', 'point2'])
    df['chr_num'] = parse_chr_to_num(df['chr'])
    _atomic_save_parquet(df, cache_path)
    _atomic_save_metadata(meta_path, fingerprint)
    return df

def load_scrna_data(counts_file, genes_file, meta_file):
    """Load scRNA genes, cells, and sparse count matrix."""
    sources = [counts_file, genes_file, meta_file]
    fingerprint = _get_source_fingerprint(sources)
    cache_dir = _get_cache_dir(counts_file)
    
    genes_cache = os.path.join(cache_dir, f"scrna_{fingerprint}_genes.parquet")
    cells_cache = os.path.join(cache_dir, f"scrna_{fingerprint}_cells.parquet")
    counts_cache = os.path.join(cache_dir, f"scrna_{fingerprint}_counts.npz")
    meta_path = os.path.join(cache_dir, f"scrna_{fingerprint}.meta")
    
    cache_files = [genes_cache, cells_cache, counts_cache]
    
    if _check_cache_integrity(cache_files, meta_path, fingerprint):
        print(f"Using cached files for scRNA data: {os.path.basename(meta_path)}")
        genes_df = pd.read_parquet(genes_cache)
        cells_df = pd.read_parquet(cells_cache)
        count_matrix = sparse.load_npz(counts_cache)
        return genes_df, cells_df, count_matrix
        
    print(f"Loading scRNA genes from {os.path.basename(genes_file)} ...")
    genes_df = pd.read_csv(genes_file, sep='\t', header=None, names=['gene_index', 'gene_symbol'])
    print(f"Loading scRNA meta from {os.path.basename(meta_file)} ...")
    cells_df = pd.read_csv(meta_file, sep='\t', header=None, 
                           names=['cell_index', 'cell_barcode', 'cell_type', 'subject_ID', 'condition'])
    
    print(f"Loading scRNA counts from {os.path.basename(counts_file)} ...")
    counts_df = pd.read_csv(counts_file, sep='\t', header=None, names=['gene_index', 'cell_index', 'readcount'])
    row = counts_df['gene_index'].values - 1
    col = counts_df['cell_index'].values - 1
    data = counts_df['readcount'].values

    # Inputs are 1-indexed (MATLAB origin); an out-of-range index would silently
    # grow the sparse matrix instead of erroring, corrupting downstream shapes.
    if len(row) and (row.min() < 0 or row.max() >= len(genes_df)
                     or col.min() < 0 or col.max() >= len(cells_df)):
        raise ValueError(
            f"scRNA counts index out of bounds: gene_index in "
            f"[{row.min() + 1}, {row.max() + 1}] vs {len(genes_df)} genes; "
            f"cell_index in [{col.min() + 1}, {col.max() + 1}] vs {len(cells_df)} cells."
        )

    count_matrix = sparse.coo_matrix((data, (row, col)), shape=(len(genes_df), len(cells_df))).tocsr()
    
    _atomic_save_parquet(genes_df, genes_cache)
    _atomic_save_parquet(cells_df, cells_cache)
    _atomic_save_npz(count_matrix, counts_cache)
    _atomic_save_metadata(meta_path, fingerprint)
    
    return genes_df, cells_df, count_matrix

def load_scatac_data(counts_file, peaks_file, meta_file):
    """Load scATAC peaks, cells, and sparse count matrix."""
    sources = [counts_file, peaks_file, meta_file]
    fingerprint = _get_source_fingerprint(sources)
    cache_dir = _get_cache_dir(counts_file)
    
    peaks_cache = os.path.join(cache_dir, f"scatac_{fingerprint}_peaks.parquet")
    cells_cache = os.path.join(cache_dir, f"scatac_{fingerprint}_cells.parquet")
    counts_cache = os.path.join(cache_dir, f"scatac_{fingerprint}_counts.npz")
    meta_path = os.path.join(cache_dir, f"scatac_{fingerprint}.meta")
    
    cache_files = [peaks_cache, cells_cache, counts_cache]
    
    if _check_cache_integrity(cache_files, meta_path, fingerprint):
        print(f"Using cached files for scATAC data: {os.path.basename(meta_path)}")
        peaks_df = pd.read_parquet(peaks_cache)
        cells_df = pd.read_parquet(cells_cache)
        count_matrix = sparse.load_npz(counts_cache)
        return peaks_df, cells_df, count_matrix

    print(f"Loading scATAC peaks from {os.path.basename(peaks_file)} ...")
    peaks_df = pd.read_csv(peaks_file, sep='\t', header=None, 
                           names=['peak_index', 'chr', 'point1', 'point2'])
    peaks_df['chr_num'] = parse_chr_to_num(peaks_df['chr'])
    
    print(f"Loading scATAC meta from {os.path.basename(meta_file)} ...")
    cells_df = pd.read_csv(meta_file, sep='\t', header=None, 
                           names=['cell_index', 'cell_barcode', 'cell_type', 'subject_ID', 'condition'])
    
    print(f"Loading scATAC counts from {os.path.basename(counts_file)} ...")
    counts_df = pd.read_csv(counts_file, sep='\t', header=None, names=['peak_index', 'cell_index', 'readcount'])
    row = counts_df['peak_index'].values - 1
    col = counts_df['cell_index'].values - 1
    data = counts_df['readcount'].values

    # Inputs are 1-indexed (MATLAB origin); an out-of-range index would silently
    # grow the sparse matrix instead of erroring, corrupting downstream shapes.
    if len(row) and (row.min() < 0 or row.max() >= len(peaks_df)
                     or col.min() < 0 or col.max() >= len(cells_df)):
        raise ValueError(
            f"scATAC counts index out of bounds: peak_index in "
            f"[{row.min() + 1}, {row.max() + 1}] vs {len(peaks_df)} peaks; "
            f"cell_index in [{col.min() + 1}, {col.max() + 1}] vs {len(cells_df)} cells."
        )

    count_matrix = sparse.coo_matrix((data, (row, col)), shape=(len(peaks_df), len(cells_df))).tocsr()
    
    _atomic_save_parquet(peaks_df, peaks_cache)
    _atomic_save_parquet(cells_df, cells_cache)
    _atomic_save_npz(count_matrix, counts_cache)
    _atomic_save_metadata(meta_path, fingerprint)
    
    return peaks_df, cells_df, count_matrix

def load_motif_prior(name_file, mapping_file, num_peaks):
    """Load motif names and TF-peak binding prior sparse matrix."""
    sources = [name_file, mapping_file]
    fingerprint = _get_source_fingerprint(sources)
    cache_dir = _get_cache_dir(name_file)
    
    motifs_cache = os.path.join(cache_dir, f"motif_{fingerprint}_motifs.parquet")
    mapping_cache = os.path.join(cache_dir, f"motif_{fingerprint}_mapping.npz")
    meta_path = os.path.join(cache_dir, f"motif_{fingerprint}.meta")
    
    cache_files = [motifs_cache, mapping_cache]
    
    if _check_cache_integrity(cache_files, meta_path, fingerprint):
        print(f"Using cached files for motif prior: {os.path.basename(meta_path)}")
        motifs_df = pd.read_parquet(motifs_cache)
        tf_peak_binding_matrix = sparse.load_npz(mapping_cache)
        return motifs_df, tf_peak_binding_matrix
        
    print(f"Loading motif names from {os.path.basename(name_file)} ...")
    motifs_df = pd.read_csv(name_file, sep=r'\s+', header=None, names=['motif_index', 'name'])
    print(f"Loading motif mapping from {os.path.basename(mapping_file)} ...")
    mapping_df = pd.read_csv(mapping_file, sep=r'\s+', header=None, names=['peak_index', 'motif_index', 'flag'])
    
    row = mapping_df['peak_index'].values - 1
    col = mapping_df['motif_index'].values - 1
    data = mapping_df['flag'].values
    
    tf_peak_binding_matrix = sparse.coo_matrix((data, (row, col)), shape=(num_peaks, len(motifs_df))).tocsr()
    
    _atomic_save_parquet(motifs_df, motifs_cache)
    _atomic_save_npz(tf_peak_binding_matrix, mapping_cache)
    _atomic_save_metadata(meta_path, fingerprint)
    
    return motifs_df, tf_peak_binding_matrix

def load_tad_regions(filepath):
    """Load TAD regions."""
    fingerprint = _get_source_fingerprint([filepath])
    cache_dir = _get_cache_dir(filepath)
    cache_path = os.path.join(cache_dir, f"tad_{fingerprint}.parquet")
    meta_path = cache_path + ".meta"
    
    if _check_cache_integrity([cache_path], meta_path, fingerprint):
        print(f"Using cached file for TAD regions: {os.path.basename(cache_path)}")
        return pd.read_parquet(cache_path)
        
    print(f"Loading TAD regions from {os.path.basename(filepath)} ...")
    df = pd.read_csv(filepath, sep=r'\s+', header=None, names=['chr', 'left_boundary', 'right_boundary'])
    df['chr_num'] = parse_chr_to_num(df['chr'])
    _atomic_save_parquet(df, cache_path)
    _atomic_save_metadata(meta_path, fingerprint)
    return df

def load_refseq(filepath):
    """Load Refseq info."""
    fingerprint = _get_source_fingerprint([filepath])
    cache_dir = _get_cache_dir(filepath)
    cache_path = os.path.join(cache_dir, f"refseq_{fingerprint}.parquet")
    meta_path = cache_path + ".meta"
    
    if _check_cache_integrity([cache_path], meta_path, fingerprint):
        print(f"Using cached file for Refseq: {os.path.basename(cache_path)}")
        return pd.read_parquet(cache_path)
        
    print(f"Loading Refseq info from {os.path.basename(filepath)} ...")
    df = pd.read_csv(filepath, sep=r'\s+', header=None, 
                     names=['chr', 'strand', 'start', 'end', 'gene_name'])
    df['chr_num'] = parse_chr_to_num(df['chr'])
    _atomic_save_parquet(df, cache_path)
    _atomic_save_metadata(meta_path, fingerprint)
    return df

