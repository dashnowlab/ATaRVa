import pyabpoa as pa
import statistics as stats
from collections import Counter
from itertools import chain

def seq_organiser(seqs):
    """
    organises the sequences in descending order of their lengths and returns a flat list of sequences.

    param seqs: list of sequences to be organised
    return: flat list of sequences organised in descending order of their lengths
    """

    seq_lens = {}
    for seq in seqs:
        slen = len(seq)
        if slen in seq_lens:
            seq_lens[slen].append(seq)
        else:
            seq_lens[slen] = [seq]
    sorted_mode_seqs = sorted(seq_lens.values(), key=lambda x : len(x), reverse=True)
    flat_list = list(chain.from_iterable(sorted_mode_seqs))
    del seq_lens, sorted_mode_seqs
    return flat_list


def seq_alter_median(seqs):
    """
    Organises the sequences in descending order of their lengths and returns a flat list of sequences with
    the median length sequence in the middle and the rest of the sequences alternatively arranged on either side.
    
    param seqs: list of sequences to be organised
    return: flat list of sequences with median length sequence in the middle and the rest of the sequences alternatively arranged on either side
    """
    
    odd         = len(seqs) %2 == 1
    sorted_seqs = sorted(seqs)
    median_idx  = len(sorted_seqs)//2
    upper_batch = sorted_seqs[median_idx:] # 1st half
    lower_batch = sorted_seqs[:median_idx][::-1] # reversed 2nd half
    insert_idx  = []
    increment   = 1
    for x in list(range(len(upper_batch))): # adjusted indices to insert the lower batch seqs alternatively
        insert_idx.append(x + increment)
        increment += 1
    low_idx   = 0
    break_idx = insert_idx[-1] # Have to break the follwoing loop for even numbers of the total reads, because the lower batch will have one seq lesser than upper batch
    for idx in insert_idx:
        if odd and (break_idx == idx): break
        upper_batch.insert(idx, lower_batch[low_idx]) # Direct alternative insertion of lower batch seqs into upper batch
        low_idx += 1
    del sorted_seqs, lower_batch
    return upper_batch


def consensus_seq_poa(seqs, amplicon):
    """
    Builds a consensus sequence from a list of sequences using Partial Order Alignment (POA) algorithm.

    param seqs: list of sequences to build consensus from
    param amplicon: boolean indicating if the sequences are from an amplicon
    return: consensus sequence
    """

    if len(seqs)<7:
        cons_algrm='MF'
    else:
        cons_algrm='HB'

    median_len = len(stats.median_high(seqs))
    if median_len > 10000:
        fseqs = [seq for seq in seqs if len(seq) == median_len]
        return ''.join(Counter(bases).most_common(1)[0][0] for bases in zip(*fseqs))

    if not amplicon:
        sorted_seqs = seq_organiser(seqs)
    else:
        sorted_seqs = seq_alter_median(seqs)
    abpoa  = pa.msa_aligner(cons_algrm=cons_algrm)
    result = abpoa.msa(sorted_seqs, out_cons=True, out_msa=False)
    return result.cons_seq[0]
