function matlab_runner(dataset_name, iteration_num, output_label)
    if nargin < 1; dataset_name = 'astrocytes'; end
    if nargin < 2; iteration_num = 500; end
    if nargin < 3; output_label = sprintf('%s_ml_%d', dataset_name, iteration_num); end

    test_dir = '/mnt/home/agebrain/ceph/anderson/snmulti_data/processed/magical/inputs';
    ct_dir = fullfile(test_dir, dataset_name);
    output_dir = 'eval/benchmarks/outputs/matlab';
    if ~isfolder(output_dir); mkdir(output_dir); end

    Candidate_gene_file_path = fullfile(ct_dir, "sig_cr_genes.txt");
    Candidate_peak_file_path = fullfile(ct_dir, "sig_cr_peaks.txt");
    scRNA_readcount_file_path = fullfile(ct_dir, "rna_counts.txt");
    scRNA_gene_file_path = fullfile(ct_dir, "rna_genes.txt");
    scRNA_cellmeta_file_path = fullfile(ct_dir, "rna_meta.txt");
    scATAC_readcount_file_path = fullfile(ct_dir, "atac_counts.txt");
    scATAC_peak_file_path = fullfile(ct_dir, "atac_peaks.txt");
    scATAC_cellmeta_file_path = fullfile(ct_dir, "atac_meta.txt");
    Motif_mapping_file_path = fullfile(test_dir, "motif_prior.txt");
    Motif_name_file_path = fullfile(test_dir, "motif_info.txt");
    TAD_flag = 1; 
    TAD_file_path = fullfile(test_dir, "tad_regions.txt");
    Ref_seq_file_path = fullfile(test_dir, "rhemac10_refseq.txt");
    Output_file_path = fullfile(output_dir, [output_label, '.txt']);
    timing_file_path = fullfile(output_dir, [output_label, '_timing_stats.txt']);

    addpath(fullfile(pwd, 'src/matlabmagical'));

    fprintf('Starting MATLAB MAGICAL with %d iterations...\n', iteration_num);
    rng('shuffle');

    % --- Data Loading ---
    tic;
    fprintf('loading all input data ...\n');
    Candidate_genes.gene_symbols=textread(Candidate_gene_file_path, '%s');
    [Candidate_peaks.chr, Candidate_peaks.point1, Candidate_peaks.point2]=textread(Candidate_peak_file_path, '%s %d %d');
    Candidate_peaks.chr_num=zeros(length(Candidate_peaks.chr), 1);
    for i=1:22; Candidate_peaks.chr_num(strcmp(Candidate_peaks.chr, ['chr', num2str(i)])>0,1)=i; end
    [scRNA_genes.gene_index, scRNA_genes.gene_symbols]=textread(scRNA_gene_file_path, '%d %s');
    [scRNA_cells.cell_index, scRNA_cells.cell_barcode, scRNA_cells.cell_type, scRNA_cells.subject_ID, scRNA_cells.condition]=textread(scRNA_cellmeta_file_path, '%d %s %s %s %s', 'delimiter', '\t');
    scATAC_assay_temp = readtable(scATAC_peak_file_path);
    scATAC_peaks.peak_index=scATAC_assay_temp{:,1};
    scATAC_peaks.chr=scATAC_assay_temp{:,2};
    scATAC_peaks.point1=scATAC_assay_temp{:,3};
    scATAC_peaks.point2=scATAC_assay_temp{:,4};
    scATAC_peaks.chr_num=zeros(length(scATAC_peaks.chr), 1);
    for i=1:22; scATAC_peaks.chr_num(strcmp(scATAC_peaks.chr, ['chr', num2str(i)])>0,1)=i; end
    [scATAC_cells.cell_index, scATAC_cells.cell_barcode, scATAC_cells.cell_type, scATAC_cells.subject_ID, scATAC_cells.condition]=textread(scATAC_cellmeta_file_path, '%d %s %s %s %s', 'delimiter', '\t');
    [scRNA_read_count_table.gene_index, scRNA_read_count_table.cell_index, scRNA_read_count_table.readcount]=textread(scRNA_readcount_file_path, '%d %d %d');
    scRNA_read_count_matrix=sparse(scRNA_read_count_table.gene_index,scRNA_read_count_table.cell_index,scRNA_read_count_table.readcount);
    if size(scRNA_read_count_matrix,1)<length(scRNA_genes.gene_symbols); scRNA_read_count_matrix(end+1:length(scRNA_genes.gene_symbols),:)=0; end
    if size(scRNA_read_count_matrix,2)<length(scRNA_cells.cell_barcode); scRNA_read_count_matrix(:,end+1:length(scRNA_cells.cell_barcode))=0; end
    [scATAC_read_count_table.peak_index, scATAC_read_count_table.cell_index, scATAC_read_count_table.readcount]=textread(scATAC_readcount_file_path, '%d %d %d');
    scATAC_read_count_matrix=sparse(scATAC_read_count_table.peak_index,scATAC_read_count_table.cell_index,scATAC_read_count_table.readcount);
    if size(scATAC_read_count_matrix,1)<length(scATAC_peaks.chr); scATAC_read_count_matrix(end+1:length(scATAC_peaks.chr),:)=0; end
    if size(scATAC_read_count_matrix,2)<length(scATAC_cells.cell_barcode); scATAC_read_count_matrix(:,end+1:length(scATAC_cells.cell_barcode))=0; end
    [Motifs.motif_index, Motifs.name]=textread(Motif_name_file_path, '%d %s');
    [Motif_mapping_table.peak_index, Motif_mapping_table.motif_index, Motif_mapping_table.flag]=textread(Motif_mapping_file_path, '%d %d %d');
    TF_peak_binding_matrix=sparse(Motif_mapping_table.peak_index, Motif_mapping_table.motif_index, Motif_mapping_table.flag);
    if size(TF_peak_binding_matrix,1)<length(scATAC_peaks.chr); TF_peak_binding_matrix(end+1:length(scATAC_peaks.chr),:)=0; end
    if size(TF_peak_binding_matrix,2)<length(Motifs.name); TF_peak_binding_matrix(:,end+1:length(Motifs.name))=0; end
    [TAD.chr, TAD.left_boundary, TAD.right_boundary]=textread(TAD_file_path, '%s %d %d');
    TAD.chr_num=zeros(length(TAD.chr), 1);
    for i=1:22; TAD.chr_num(strcmp(TAD.chr, ['chr', num2str(i)])>0,1)=i; end
    [Refseq.chr, Refseq.strand, Refseq.start, Refseq.end, Refseq.gene_name]=textread(Ref_seq_file_path, '%s %s %d %d %s', 'headerlines', 1);
    Refseq.chr_num=zeros(length(Refseq.chr), 1);
    for i=1:22; Refseq.chr_num(strcmp(Refseq.chr, ['chr', num2str(i)])>0,1)=i; end
    Common_samples=intersect(unique(scRNA_cells.subject_ID), unique(scATAC_cells.subject_ID));
    time_load = toc;

    % --- Circuit Construction ---
    tic;
    [Candidate_TFs, Candidate_TF_log2Count,...
        Candidate_peaks, Candidate_Peak_log2Count,...
        Candidate_genes, Candidate_Gene_log2Count,...
        Candidate_TF_Peak_Binding, Candidate_Peak_Gene_looping,...
        ATAC_cell_vector, scATAC_read_count_matrix, RNA_cell_vector, scRNA_read_count_matrix]=...
        Candidate_circuits_construction_with_TAD(Common_samples, Candidate_genes, Candidate_peaks,...
                                      scRNA_genes, scRNA_cells, scRNA_read_count_matrix, ...
                                      scATAC_peaks, scATAC_cells, scATAC_read_count_matrix,...
                                      Motifs, TF_peak_binding_matrix,...
                                      Refseq, TAD);
    time_circuit = toc;

    % --- Initialization ---
    tic;
    S=length(Common_samples); M=length(Candidate_TFs); P=length(Candidate_peaks.peak_index); G=length(Candidate_genes.gene_symbols);
    [T_A_initial, T_R_initial, T_sample_mean, T_sample_var, B_initial, B_prior_mean, B_prior_var, B_prob, L_initial, L_prior_mean, L_prior_var, L_prob]=...
        MAGICAL_initialization(Candidate_TF_log2Count, Candidate_Peak_log2Count, Candidate_Gene_log2Count, Candidate_TF_Peak_Binding, Candidate_Peak_Gene_looping, M, S, ATAC_cell_vector, RNA_cell_vector);
    time_init = toc;

    % --- Gibbs Sampling ---
    tic;
    [Candidate_TF_Peak_Binding_prob, Candidate_Peak_Gene_Looping_prob]=...
        MAGICAL_estimation(scATAC_read_count_matrix, ATAC_cell_vector, Candidate_Peak_log2Count, scRNA_read_count_matrix, RNA_cell_vector, Candidate_Gene_log2Count, Candidate_TF_Peak_Binding, Candidate_Peak_Gene_looping,...
        T_A_initial, T_R_initial, T_sample_mean, T_sample_var, B_initial, B_prior_mean, B_prior_var, B_prob, L_initial, L_prior_mean, L_prior_var, L_prob, M, S, P, G, iteration_num);
    time_gibbs = toc;

    % --- Output ---
    prob_threshold_TF_peak_binding=0.7;
    prob_threshold_peak_gene_looping=0.7;

    fid=fopen(Output_file_path, 'w');
    fprintf(fid, 'Gene_symbol\tGene_chr\tGene_TSS\tPeak_chr\tPeak_start\tPeak_end\tLooping_prob\tTFs(binding prob)\n');

    [xx,yy]=find(Candidate_Peak_Gene_Looping_prob>prob_threshold_peak_gene_looping);
    circuit_flag=zeros(length(xx),1);
    TF_vector=zeros(length(Candidate_TFs), 1);
    for i=1:length(xx)
        [TF_prob, TF_index]=sort(full(Candidate_TF_Peak_Binding_prob(xx(i),:)), 'descend');
        index=find(TF_prob>prob_threshold_TF_peak_binding);
        if ~isempty(index)
            circuit_flag(i)=1;
            fprintf(fid, '%s\tchr%d\t%d\tchr%d\t%d\t%d\t%G\t', ...
                Candidate_genes.gene_symbols{yy(i)}, Candidate_genes.gene_TSS(yy(i),:),...
                Candidate_peaks.chr_num(xx(i)), Candidate_peaks.point1(xx(i)),  Candidate_peaks.point2(xx(i)), full(Candidate_Peak_Gene_Looping_prob(xx(i),yy(i))));
            
            for t=1:length(index)
                fprintf(fid, '%s (%.2f), ', Candidate_TFs{TF_index(t)}, TF_prob(t));
                TF_vector(TF_index(t))=TF_vector(TF_index(t))+1;
            end
            fprintf(fid, '\n');
        end
    end
    fclose(fid);

    % Write matrices
    [fpath, fname, fext] = fileparts(Output_file_path);
    B_matrix_file = fullfile(fpath, sprintf('%s_B_matrix.txt', fname));
    L_matrix_file = fullfile(fpath, sprintf('%s_L_matrix.txt', fname));

    fid_B = fopen(B_matrix_file, 'w');
    fprintf(fid_B, 'Peak_Coordinates');
    for j = 1:length(Candidate_TFs); fprintf(fid_B, '\t%s', Candidate_TFs{j}); end
    fprintf(fid_B, '\n');
    for i = 1:P
        fprintf(fid_B, 'chr%d_%d_%d', Candidate_peaks.chr_num(i), Candidate_peaks.point1(i), Candidate_peaks.point2(i));
        row_data = full(Candidate_TF_Peak_Binding_prob(i, :));
        fprintf(fid_B, '\t%G', row_data);
        fprintf(fid_B, '\n');
    end
    fclose(fid_B);

    fid_L = fopen(L_matrix_file, 'w');
    fprintf(fid_L, 'Peak_Coordinates');
    for j = 1:G; fprintf(fid_L, '\t%s', Candidate_genes.gene_symbols{j}); end
    fprintf(fid_L, '\n');
    for i = 1:P
        fprintf(fid_L, 'chr%d_%d_%d', Candidate_peaks.chr_num(i), Candidate_peaks.point1(i), Candidate_peaks.point2(i));
        row_data = full(Candidate_Peak_Gene_Looping_prob(i, :));
        fprintf(fid_L, '\t%G', row_data);
        fprintf(fid_L, '\n');
    end
    fclose(fid_L);

    % Write timing
    fid = fopen(timing_file_path, 'w');
    fprintf(fid, 'Data Loading: %f\n', time_load);
    fprintf(fid, 'Circuit Construction: %f\n', time_circuit);
    fprintf(fid, 'Initialization (OLS): %f\n', time_init);
    fprintf(fid, 'Gibbs Sampling (%d iters): %f\n', iteration_num, time_gibbs);
    fclose(fid);
end
