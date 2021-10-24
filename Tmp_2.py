

def blast_results_to_dict(blastn_results, iden_cutoff, query_cov_cutoff):

    query_to_subject_list_dict = {}

    for blast_hit in open(blastn_results):
        blast_hit_split = blast_hit.strip().split('\t')
        query = blast_hit_split[0]
        subject = blast_hit_split[1]
        subject_with_prefix = 'GenomicSeq__%s' % subject
        iden = float(blast_hit_split[2])
        align_len = int(blast_hit_split[3])
        query_len = int(blast_hit_split[12])
        coverage_q = float(align_len) * 100 / float(query_len)
        if (iden >= iden_cutoff) and (coverage_q >= query_cov_cutoff):
            if query not in query_to_subject_list_dict:
                query_to_subject_list_dict[query] = [subject_with_prefix]
            else:
                query_to_subject_list_dict[query].append(subject_with_prefix)

    return query_to_subject_list_dict

