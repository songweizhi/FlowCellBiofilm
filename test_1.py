
################################################## get_mutated_genes ###################################################

SNV_matrix_cdc_list = [SNV_210_mono_uniq_file, SNV_210_co_uniq_file, SNV_210_concurrence_file, SNV_D2_mono_uniq_file, SNV_D2_co_uniq_file, SNV_D2_concurrence_file]

SNV_to_gene_dict = {}
SNV_to_ko_id_dict = {}
SNV_to_ko_desc_dict = {}
for SNV_matrix_cdc in SNV_matrix_cdc_list:

    # output files
    SNV_matrix_cdc_path, SNV_matrix_cdc_basename, SNV_matrix_cdc_ext = sep_path_basename_ext(SNV_matrix_cdc)

    output_mutated_genes_cate =         '%s/%s_mutated_genes_category.txt'     % (output_dir, SNV_matrix_cdc_basename)
    output_mutated_genes_cate_fun =     '%s/%s_mutated_genes_cate_fun.txt'     % (output_dir, SNV_matrix_cdc_basename)
    output_summary =                    '%s/%s_summary.txt'                    % (output_dir, SNV_matrix_cdc_basename)
    effect_file =                       '%s/%s_mutation_effect.txt'            % (output_dir, SNV_matrix_cdc_basename)
    output_seq_nc =                     '%s/%s_mutated_genes_nc.fasta'         % (output_dir, SNV_matrix_cdc_basename)
    output_seq_aa =                     '%s/%s_mutated_genes_aa.fasta'         % (output_dir, SNV_matrix_cdc_basename)

    # store sequences in dict
    ref_seq_dict = {}
    for each_seq in SeqIO.parse(combined_ref_fasta, 'fasta'):
        ref_seq_dict[each_seq.id] = str(each_seq.seq)

    # get seq id list
    seq_id_list = []
    for cds in open(combined_ref_gff):
        if (cds.startswith('2.10')) or (cds.startswith('D2')):
            seq_ID = cds.strip().split('\t')[0]
            if seq_ID not in seq_id_list:
                seq_id_list.append(seq_ID)

    # get ending positions for all genes
    ORF_ending_pos_dict = {}
    ORF_seq_id_dict = {}
    ORF_strand_dict = {}
    for cds2 in open(combined_ref_gff):
        if (cds2.startswith('2.10')) or (cds2.startswith('D2')):
            cds2_split = cds2.strip().split('\t')
            seq_id2 = cds2_split[0]
            gene_id2 = cds2_split[8].split(';')[0].split('=')[1]
            start_pos2 = int(cds2_split[3])
            end_pos2 = int(cds2_split[4])
            strand = cds2_split[6]
            ORF_ending_pos_dict[gene_id2] = [start_pos2, end_pos2]
            ORF_seq_id_dict[gene_id2] = seq_id2
            ORF_strand_dict[gene_id2] = strand

    # initialize coding region dict
    coding_region_dict = {}
    for each_id in seq_id_list:
        coding_region_dict[each_id] = []

    # add coding pos to coding_region_dict
    for each_cds in open(combined_ref_gff):
        if (each_cds.startswith('2.10')) or (each_cds.startswith('D2')):
            each_cds_split = each_cds.strip().split('\t')
            seq_id = each_cds_split[0]
            start_pos = int(each_cds_split[3])
            end_pos = int(each_cds_split[4])
            pos_list = list(range(start_pos, end_pos + 1))
            for each_pos in pos_list:
                coding_region_dict[seq_id].append(each_pos)

    output_handle = open(output_summary, 'w')
    for each_snv in open(SNV_matrix_cdc):
        if not each_snv.startswith('\t'):
            each_snv_id =  each_snv.strip().split('\t')[0]
            each_snv_seq = each_snv.strip().split(',')[0].split('|')[0]
            each_snv_pos = int(each_snv.strip().split(',')[0].split('|')[1])
            each_snv_pos_wt = each_snv.strip().split(',')[0].split('|')[2]
            each_snv_pos_v = each_snv.strip().split(',')[0].split('|')[3]

            # get SNV location
            location = ''
            if each_snv_pos in coding_region_dict[each_snv_seq]:
                location = 'Coding'
            else:
                location = 'Intergenic'
                output_handle.write('%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\n' % (each_snv_id,
                                                                              SNV_parallel_dict_combined[each_snv.strip()],
                                                                              location))
                SNV_to_gene_dict[each_snv_id]    = 'Intergenic'
                SNV_to_ko_id_dict[each_snv_id]   = 'NA'
                SNV_to_ko_desc_dict[each_snv_id] = 'NA'

            # get current gene's COG ID, description, mutation type and mutation effect
            if location == 'Coding':
                for each_gene in ORF_ending_pos_dict:
                    if (each_snv_seq == ORF_seq_id_dict[each_gene]) and (ORF_ending_pos_dict[each_gene][0] <= each_snv_pos <= ORF_ending_pos_dict[each_gene][1]):

                        start_pos_raw = ORF_ending_pos_dict[each_gene][0]
                        end_pos_raw = ORF_ending_pos_dict[each_gene][1]
                        snv_pos_raw = each_snv_pos
                        start_pos_rescaled = ORF_ending_pos_dict[each_gene][0] - (ORF_ending_pos_dict[each_gene][0] - 1)
                        end_pos_rescaled = ORF_ending_pos_dict[each_gene][1] - (ORF_ending_pos_dict[each_gene][0] - 1)
                        snv_pos_rescaled = each_snv_pos - (ORF_ending_pos_dict[each_gene][0] - 1)
                        mutation_type = ''

                        # get current gene's COG ID and description
                        current_COG_id = ''
                        current_COG_function = ''
                        if each_gene in gene_KEGG_id_dict:
                            current_COG_id = gene_KEGG_id_dict[each_gene]
                            current_COG_function = gene_KEGG_function_dict[each_gene]
                        else:
                            current_COG_id = 'NA'
                            current_COG_function = 'NA'

                        # get current gene's mutation type
                        if each_snv_pos_v == '-':
                            mutation_type = 'Frameshift'
                            output_handle.write('%s\t%s\t%s\t%s\t%s\tNA\t%s\t%s\t%s\n' % (each_snv_id,
                                                                                          SNV_parallel_dict_combined[each_snv.strip()],
                                                                                          location,
                                                                                          ORF_strand_dict[each_gene],
                                                                                          each_gene,
                                                                                          mutation_type,
                                                                                          gene_KEGG_id_dict[each_gene],
                                                                                          gene_KEGG_function_dict[each_gene]))
                            SNV_to_gene_dict[each_snv_id]    = each_gene
                            SNV_to_ko_id_dict[each_snv_id]   = gene_KEGG_id_dict[each_gene]
                            SNV_to_ko_desc_dict[each_snv_id] = gene_KEGG_function_dict[each_gene]

                        elif each_snv_pos_v in ['A', 'T', 'C', 'G']:

                            # get all reading frame
                            rf_pos_list = get_rf_pos_list(snv_pos_rescaled)

                            # get the mutated reading frame with rescaled position
                            mutated_rf_rescaled = []
                            for each_rf in rf_pos_list:
                                if snv_pos_rescaled in each_rf:
                                    mutated_rf_rescaled = each_rf

                            # get the mutated reading frame with raw position
                            mutated_rf_raw = []
                            for each_rescaled_pos in mutated_rf_rescaled:
                                raw_pos = each_rescaled_pos + ORF_ending_pos_dict[each_gene][0] - 1
                                mutated_rf_raw.append(raw_pos)

                            # get sequence of raw and mutated reading frame
                            rf_seq_raw = ref_seq_dict[each_snv_seq][(mutated_rf_raw[0] - 1):(mutated_rf_raw[0] + 2)]
                            rf_seq_mutated = ''
                            for each_bp in mutated_rf_raw:
                                current_bp = ''
                                if each_bp == each_snv_pos:
                                    current_bp = each_snv_pos_v
                                else:
                                    current_bp = ref_seq_dict[each_snv_seq][each_bp - 1]
                                rf_seq_mutated += current_bp

                            if ORF_strand_dict[each_gene] == '-':
                                rf_seq_raw = str(Seq(rf_seq_raw, generic_dna).reverse_complement())
                                rf_seq_mutated = str(Seq(rf_seq_mutated, generic_dna).reverse_complement())

                            rf_seq_raw_aa = str(SeqRecord(Seq(rf_seq_raw)).seq.translate(table=transl_table))
                            rf_seq_mutated_aa = str(SeqRecord(Seq(rf_seq_mutated)).seq.translate(table=transl_table))

                            mutation_type_term = ''
                            if rf_seq_raw_aa == rf_seq_mutated_aa:
                                mutation_type_term = 'Silent'
                            elif rf_seq_mutated_aa == '*':
                                mutation_type_term = 'Nonsense'
                            else:
                                mutation_type_term = 'Missense'
                            aa_mutation = '%s(%s)->%s(%s)' % (rf_seq_raw, rf_seq_raw_aa, rf_seq_mutated, rf_seq_mutated_aa)

                            # print out
                            for_write = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (each_snv_id,
                                                                                  SNV_parallel_dict_combined[each_snv.strip()],
                                                                                  location,
                                                                                  ORF_strand_dict[each_gene],
                                                                                  each_gene,
                                                                                  aa_mutation,
                                                                                  mutation_type_term,
                                                                                  gene_KEGG_id_dict[each_gene],
                                                                                  gene_KEGG_function_dict[each_gene])
                            SNV_to_gene_dict[each_snv_id]    = each_gene
                            SNV_to_ko_id_dict[each_snv_id]   = gene_KEGG_id_dict[each_gene]
                            SNV_to_ko_desc_dict[each_snv_id] = gene_KEGG_function_dict[each_gene]

                            output_handle.write(for_write)
    output_handle.close()

    # get mutation_cate_summary
    mutation_cate_dict = get_mutation_cate_summary(output_summary, output_mutated_genes_cate)
    os.system('rm %s' % output_mutated_genes_cate)



