function snmulti_main(celltype_idx, iteration_num)

    PRO_DIR = "/mnt/home/agebrain/ceph/anderson/snmulti_data/processed/";
    INPUT_DIR = fullfile(PRO_DIR, "magical/inputs");
    OUTPUT_DIR = fullfile(PRO_DIR, "magical/outputs");
    if ~isfolder(OUTPUT_DIR)
        mkdir(OUTPUT_DIR)
    end

    short_cell_types = fullfile(PRO_DIR, "metadata/names_ct_lower.tsv");
    
    % Use 'FileType', 'text' to force MATLAB to read the .tsv file
    cellTypes = readcell(short_cell_types, 'FileType', 'text', 'Delimiter', '\t');
    
    celltype = cellTypes{celltype_idx};
    fprintf("cell type: %s\n", celltype);

    CT_DIR = fullfile(INPUT_DIR, celltype);

    scRNA_readcount_file_path = fullfile(CT_DIR, "rna_counts.txt");
    scRNA_gene_file_path = fullfile(CT_DIR, "rna_genes.txt");
    scRNA_cellmeta_file_path = fullfile(CT_DIR, "rna_meta.txt");

    scATAC_readcount_file_path = fullfile(CT_DIR, "atac_counts.txt");
    scATAC_peak_file_path = fullfile(CT_DIR, "atac_peaks.txt");
    scATAC_cellmeta_file_path = fullfile(CT_DIR, "atac_meta.txt");

    Motif_mapping_file_path = fullfile(INPUT_DIR, "motif_prior.txt");
    Motif_name_file_path = fullfile(INPUT_DIR, "motif_info.txt");
    TAD_flag = 1; 
    TAD_file_path = fullfile(INPUT_DIR, "tad_regions.txt");
    Ref_seq_file_path = fullfile(INPUT_DIR, "rhemac10_refseq.txt");

    Candidate_gene_file_path = fullfile(CT_DIR, "sig_cr_genes.txt");
    Candidate_peak_file_path = fullfile(CT_DIR, "sig_cr_peaks.txt");

    Output_file_path = fullfile(OUTPUT_DIR, sprintf("%s_%d.txt", celltype, iteration_num));
    
    % Using fprintf for cleaner console output
    fprintf("output file: %s\n", Output_file_path);
    fprintf("iteration number: %d\n", iteration_num);

    MAGICAL(Candidate_gene_file_path, Candidate_peak_file_path,...
    scRNA_readcount_file_path, scRNA_gene_file_path, scRNA_cellmeta_file_path,...
    scATAC_readcount_file_path, scATAC_peak_file_path, scATAC_cellmeta_file_path,...
    Motif_mapping_file_path, Motif_name_file_path, TAD_flag, TAD_file_path, Ref_seq_file_path, ...
    Output_file_path, iteration_num);
    
end