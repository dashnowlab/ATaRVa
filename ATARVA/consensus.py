import pyabpoa as pa
import statistics as stats
from collections import Counter
from itertools import chain

import random

MAX_POA_LEN   = 10000    # POA memory grows ~quadratically with length; tune on your data
MAX_POA_READS = 50

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
    sorted_seqs = sorted(seq_lens.values(), key=lambda x : len(x), reverse=True)
    flat_list   = list(chain.from_iterable(sorted_seqs))
    del seq_lens, sorted_seqs
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
    cons_algrm = 'MF' if len(seqs) < 7 else 'HB'

    median_seq = stats.median_high(seqs)             # was: len(median_high(seqs))
    median_len = len(median_seq)
    lens       = [len(s) for s in seqs]
    max_len    = max(lens)

    if median_len > MAX_POA_LEN or max_len - median_len > MAX_POA_LEN * 0.5:
        # cheap path: column-wise majority over reads of the modal-ish length
        fseqs = [s for s, l in zip(seqs, lens) if l == median_len]
        return ''.join(Counter(bases).most_common(1)[0][0] for bases in zip(*fseqs))

    if len(seqs) > MAX_POA_READS:                    # cap depth
        seqs = random.Random(0).sample(seqs, MAX_POA_READS)   # fixed seed keeps runs reproducible

    sorted_seqs = seq_organiser(seqs) if not amplicon else seq_alter_median(seqs)
    result = pa.msa_aligner(cons_algrm=cons_algrm).msa(sorted_seqs, out_cons=True, out_msa=False)
    return result.cons_seq[0]
