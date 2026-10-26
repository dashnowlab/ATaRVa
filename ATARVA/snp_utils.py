import numpy as np

def haplocluster_reads(cooper, locus_key):
    """
    cluster reads into haplotypes using SNP allele information.

    :param cooper:      cooper object
    :param locus_key:   locus key
    :return:            minimum SNP position used for clustering, or -1 if clustering fails
    """

    locus = cooper.cooper_loci_info[locus_key]
    locus_data = cooper.cooper_loci_data[locus_key]
    locus_cov  = locus_data.depth

    relevant_snps = set()    # SNPs relevant to the locus based on proximity
    for rindex in locus_data.reads:
        relevant_snps |= (cooper.cooper_read_data[rindex].snps)
    relevant_snps = sorted(list(filter(lambda x: (x in cooper.cooper_snp_data) and (cooper.cooper_snp_data[x].cov >= 3) and
                                                 (locus.start - cooper.args.snp_dist < x < locus.end + cooper.args.snp_dist),
                            relevant_snps)))

    relevant_snp_data = {}
    alt_snp_cov = {}
    for pos in relevant_snps:
        qual_check = False
        snp_data = cooper.cooper_snp_data[pos]
        max_coverage = 0

        passed_alts = set()
        for alt in snp_data.sub:
            alt_reads = snp_data.sub[alt].intersection(set(locus_data.reads))
            alt_cov   = len(alt_reads)
            if alt_cov == 0: continue
            if alt_cov > max_coverage: max_coverage = alt_cov

            avg_alt_qual = sum([snp_data.qual[idx] for idx in alt_reads])/alt_cov

            if avg_alt_qual <= cooper.args.snp_qual:
                qual_check = True
                break
            passed_alts.add(alt)
        if max_coverage == 0 or qual_check:
            continue
        alt_snp_cov[pos] = max_coverage

        relevant_snp_data[pos] = { 'cov': 0, 'alleles': {}, 'qual': {} }
        for alt in passed_alts:
            alt_reads = snp_data.sub[alt].intersection(set(locus_data.reads))
            relevant_snp_data[pos]['cov']           += len(alt_reads)
            relevant_snp_data[pos]['alleles'][alt]   = alt_reads
            relevant_snp_data[pos]['qual'][alt]      = snp_data.qual
        ref_reads = snp_data.ref.intersection(set(locus_data.reads))
        if len(ref_reads) > 0:
            relevant_snp_data[pos]['alleles']['r'] = ref_reads
            relevant_snp_data[pos]['cov']         += len(ref_reads)

    del_positions = list(filter(lambda x: relevant_snp_data[x]['cov'] < 5, relevant_snp_data.keys()))

    for pos in del_positions:
        del relevant_snp_data[pos]
    ordered_snp_on_cov = sorted(relevant_snp_data.keys(), key = lambda item : relevant_snp_data[item]['cov'] + alt_snp_cov[item], reverse = True)

    sig_snp_data    = {}
    ordered_sig_snps = []

    for pos in ordered_snp_on_cov:
        snp_data  = relevant_snp_data[pos]
        pos_cov   = snp_data['cov']
        if pos_cov < 0.6 * locus_cov: break
        ordered_sig_snps.append(pos)
        sig_snp_data[pos] = { 'cov': snp_data['cov'], 'alleles': snp_data['alleles'], 'qual': snp_data['qual'] }

    min_snp_pos       = qvalue_phasing(cooper, locus, locus_data, sig_snp_data, ordered_sig_snps)

    THRESHOLD_RANGES  = [(0.3, 0.7), (0.25, 0.75), (0.2, 0.8)]
    if not locus_data.is_genotyped:
        for tier_idx, (lower_thresh, upper_thresh) in enumerate(THRESHOLD_RANGES):

            sig_snp_data    = {}
            ordered_sig_snps = []

            for pos in ordered_snp_on_cov:
                snp_data  = relevant_snp_data[pos]
                pos_cov   = snp_data['cov']
                if pos_cov < 0.6 * locus_cov: break

                # count alleles where read fraction is within threshold bounds
                balanced_alleles = sum(
                    lower_thresh * pos_cov <= len(reads) <= upper_thresh * pos_cov
                    for reads in snp_data['alleles'].values()
                )

                if balanced_alleles >= 2:
                    ordered_sig_snps.append(pos)
                    sig_snp_data[pos] = { 'cov': snp_data['cov'], 'alleles': snp_data['alleles'], 'qual': snp_data['qual'] }

            if not ordered_sig_snps:
                if tier_idx < 2: continue
                return -1

            min_snp_pos = merge_snpreadsets(cooper, locus, locus_data, sig_snp_data, ordered_sig_snps)

            if locus_data.is_genotyped or tier_idx == 2:
                break

    return min_snp_pos


def merge_snpreadsets(cooper, locus, locus_data, sig_snp_data, ordered_sig_snps):
    """
    merge SNP read sets to phase reads into two haplotype clusters.

    :param cooper:              cooper object
    :param locus_data:          locus data object
    :param sig_snp_data:        significant SNP position data
    :param ordered_sig_snps:    SNP positions ordered by significance
    :return:                    [haplotypes, success, min_snp, skip_point, snp_quals, phased_reads, snp_count]
    """
    CLUSTER_JOIN_HIGH     = 0.7    # min intersection to join a cluster
    CLUSTER_JOIN_LOW      = 0.05   # max intersection to confirm non-membership

    del_positions = []
    for pos in ordered_sig_snps:
        remove_alts = []
        for alt, reads in sig_snp_data[pos]['alleles'].items():
            if len(reads) < 0.2 * sig_snp_data[pos]['cov']:
                remove_alts.append(alt)
        for alt in remove_alts:
            del sig_snp_data[pos]['alleles'][alt]
            if alt in sig_snp_data[pos]['qual']:
                del sig_snp_data[pos]['qual'][alt]
        if len(sig_snp_data[pos]['alleles']) < 2:
            del_positions.append(pos)
    for pos in del_positions:
        del sig_snp_data[pos]
        ordered_sig_snps.remove(pos)

    snps     = ordered_sig_snps[:cooper.args.snp_count]

    # --- compute quality values for top SNPs ---
    snp_quals = [max(list([max(list(x.values())) for x in sig_snp_data[pos]['qual'].values()])) for pos in ordered_sig_snps]
    snp_quals = ','.join(str(int(q)) for q in snp_quals)

    # --- compute pairwise mismatch scores between SNPs ---
    mismatch_scores = {}                  # {pos_a: {pos_b: mismatch_score}}

    for i, pos_a in enumerate(snps):
        if mismatch_scores and i == len(snps) - 1:
            break
        mismatch_scores[pos_a] = {}
        readsets_a = list(sig_snp_data[pos_a]['alleles'].values())

        for pos_b in snps[i + 1:]:
            score = sum(
                min(len(reads_b & reads_a) for reads_a in readsets_a)
                for reads_b in sig_snp_data[pos_b]['alleles'].values()
            )
            mismatch_scores[pos_a][pos_b] = score

    # --- select best SNPs for phasing ---
    # prefer positions with at least 2 zero-mismatch neighbours
    sig_snps = []
    for pos, neighbours in mismatch_scores.items():
        if sum(1 for s in neighbours.values() if s == 0) >= 2:
            sig_snps.append(pos)
            sig_snps.extend(sorted(neighbours, key=lambda p: neighbours[p]))
            break

    if cooper.args.amplicon and len(sig_snps) < 3:
        return -1

    # fallback — pick position with lowest sum of two best scores
    if not sig_snps:
        best_pos = min(
            mismatch_scores,
            key=lambda p: sum(sorted(mismatch_scores[p].values())[:2])
        )
        sig_snps.append(best_pos)
        sig_snps.extend(
            sorted(mismatch_scores[best_pos], key=lambda p: mismatch_scores[best_pos][p])
        )

    # if len(sig_snps) < 2: return -1

    # --- build final ordered SNP dict ---
    final_snp_dict = {
        pos: sig_snp_data[pos]['alleles']
        for pos in sig_snps
        if pos in sig_snp_data
    }
    min_snp_pos = min(final_snp_dict)

    # --- cluster reads into two haplotypes ---
    cluster1: set = set()
    cluster2: set = set()

    for alleles in final_snp_dict.values():
        for read_set in alleles.values():
            if not cluster1:
                cluster1 |= read_set
            elif not cluster2:
                cluster2 |= read_set
            elif (len(read_set & cluster1) > CLUSTER_JOIN_HIGH * len(read_set) and
                  len(read_set & cluster2) < CLUSTER_JOIN_LOW  * len(read_set)):
                cluster1 |= read_set
            elif (len(read_set & cluster2) > CLUSTER_JOIN_HIGH * len(read_set) and
                  len(read_set & cluster1) < CLUSTER_JOIN_LOW  * len(read_set)):
                cluster2 |= read_set

    # --- validate phasing coverage ---
    total_phased = len(cluster1) + len(cluster2)

    if cooper.args.amplicon and (locus_data.depth - total_phased > 10 or total_phased < 0.9 * locus_data.depth):
        locus_data.is_genotyped    = False
        locus_data.skip_code     = 1
        return min_snp_pos

    if total_phased >= cooper.args.phased_frac * locus_data.depth and \
        len(cluster1) >= cooper.args.cluster_frac*locus_data.depth and len(cluster2) >= cooper.args.cluster_frac*locus_data.depth:
        locus_data.hap_read_sets      = (list(cluster1), list(cluster2))
        locus_data.hap_alen_sets      = ([locus_data.read_alens[ridx][0] for ridx in cluster1], [locus_data.read_alens[ridx][0] for ridx in cluster2])
        locus_data.is_genotyped       = True
        locus_data.phase_mode         = 'snp'
        locus_data.hap_category       = 0
        locus_data.n_phasing_snps     = len(final_snp_dict)
        locus_data.phasing_snp_quals  = snp_quals

        return min_snp_pos

    locus_data.is_genotyped    = False
    locus_data.skip_code     = 1
    return min_snp_pos


def _cluster_reads(score, num_reads):
    """
    Greedy 2-way clustering over read pairs ordered by most negative score.
    Works on matrix indices; returns (cluster1_indices, cluster2_indices).
    
    :param score:      pairwise score matrix (numpy array)
    :param num_reads:  number of reads (size of score matrix)
    :return:           (cluster1_indices, cluster2_indices)
    """
    flat = score.ravel()
    neg = np.flatnonzero(flat < 0)
    order = neg[np.argsort(flat[neg], kind='stable')]
    rows = (order // num_reads).tolist()
    cols = (order % num_reads).tolist()

    cluster1, cluster2, phased = [], [], set()
    sum1 = np.zeros(num_reads)
    sum2 = np.zeros(num_reads)

    def place(j):
        s1 = sum1[j] / len(cluster1) if cluster1 else 0.0
        s2 = sum2[j] / len(cluster2) if cluster2 else 0.0
        if s1 < s2:
            cluster2.append(j); sum2.__iadd__(score[j])
        elif s1 > s2:
            cluster1.append(j); sum1.__iadd__(score[j])
        phased.add(j)

    for itr, (ia, ib) in enumerate(zip(rows, cols)):
        if ia in phased and ib in phased:
            continue
        if itr == 0:
            cluster1.append(ia); sum1 += score[ia]
            cluster2.append(ib); sum2 += score[ib]
            phased.update((ia, ib))
            continue
        if ib not in phased: place(ib)
        if ia not in phased: place(ia)
        if len(phased) == num_reads:
            break

    for j in range(num_reads):
        if j not in phased:
            place(j)

    return cluster1, cluster2


def qvalue_phasing(cooper, locus, locus_data, sig_snp_data, ordered_sig_snps):
    """
    Phase reads into two haplotype clusters using SNP allele information and quality values.
    
    :param cooper:          Cooper instance
    :param locus:           Locus instance
    :param locus_data:      Locus data
    :param sig_snp_data:    Significant SNP data
    :param ordered_sig_snps: Ordered significant SNPs
    :return:                -1 if phasing fails, otherwise the number of phased reads
    """
    args = cooper.args

    # --- drop low-support alleles, then SNPs left with <2 alleles ---
    del_positions = set()
    for pos in ordered_sig_snps:
        d = sig_snp_data[pos]
        alleles, qual = d['alleles'], d['qual']
        thr = 0.2 * d['cov']
        for alt in [a for a, reads in alleles.items() if len(reads) < thr]:
            del alleles[alt]
            qual.pop(alt, None)
        if len(alleles) < 2:
            del_positions.add(pos)
    if del_positions:
        for pos in del_positions:
            del sig_snp_data[pos]
        ordered_sig_snps[:] = [p for p in ordered_sig_snps if p not in del_positions]

    if not ordered_sig_snps: return -1
    if args.amplicon and len(ordered_sig_snps) < 3: return -1

    depth = locus_data.depth
    if len(ordered_sig_snps) == 1:
        max_alt_cov = max(len(r) for r in sig_snp_data[ordered_sig_snps[0]]['alleles'].values())
        if max_alt_cov > 0.7 * depth or max_alt_cov < 0.3 * depth:
            return -1

    ordered_sig_snps = ordered_sig_snps[:args.snp_count]
    min_snp_pos  = min(ordered_sig_snps)
    sorted_reads = sorted(locus_data.reads)
    num_reads    = len(sorted_reads)

    snp_quals = ','.join(
        str(int(max(max(x.values()) for x in sig_snp_data[pos]['qual'].values())))
        for pos in ordered_sig_snps
    )

    # --- pairwise score matrix (vectorized) ---
    read_idx   = {r: i for i, r in enumerate(sorted_reads)}
    mean_quals = None
    score = np.zeros((num_reads, num_reads), dtype=np.float64)

    for snp_pos in ordered_sig_snps:
        d = sig_snp_data[snp_pos]
        codes = np.full(num_reads, -1, dtype=np.int32)
        quals = np.zeros(num_reads, dtype=np.float64)
        for code, (allele, reads) in enumerate(d['alleles'].items()):
            if allele == 'r':
                if mean_quals is None:
                    mean_quals = np.array([cooper.cooper_read_data[r].mean_qual for r in sorted_reads], dtype=np.float64)
                for r in reads:
                    i = read_idx.get(r)
                    if i is not None:
                        codes[i] = code; quals[i] = mean_quals[i]
            else:
                aq = d['qual'][allele]
                for r in reads:
                    i = read_idx.get(r)
                    if i is not None:
                        codes[i] = code; quals[i] = aq[r]

        valid = np.flatnonzero((codes >= 0) & (quals >= args.snp_qual))
        if valid.size < 2: continue
        cv, qv = codes[valid], quals[valid]
        score[np.ix_(valid, valid)] -= (qv[:, None] + qv[None, :]) * (cv[:, None] != cv[None, :])

    # --- clustering (replaces the old loops and score_cluster_membership) ---
    c1_idx, c2_idx = _cluster_reads(score, num_reads)
    cluster1 = [sorted_reads[i] for i in c1_idx]
    cluster2 = [sorted_reads[i] for i in c2_idx]
    n1, n2 = len(cluster1), len(cluster2)
    total_phased = n1 + n2

    # --- validate phasing coverage ---
    if args.amplicon and (depth - total_phased > 10 or total_phased < 0.9 * depth):
        locus_data.is_genotyped = False
        locus_data.skip_code    = 1
        return min_snp_pos

    if total_phased >= args.phased_frac * depth and \
       n1 >= args.cluster_frac * depth and n2 >= args.cluster_frac * depth:
        alens = locus_data.read_alens
        locus_data.hap_read_sets     = (cluster1, cluster2)
        locus_data.hap_alen_sets     = ([alens[r][0] for r in cluster1], [alens[r][0] for r in cluster2])
        locus_data.is_genotyped      = True
        locus_data.phase_mode        = 'snp'
        locus_data.hap_category      = 0
        locus_data.n_phasing_snps    = len(sig_snp_data)
        locus_data.phasing_snp_quals = snp_quals
        return min_snp_pos

    locus_data.is_genotyped = False
    locus_data.skip_code    = 1
    return min_snp_pos
