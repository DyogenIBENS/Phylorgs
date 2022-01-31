#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Get compositional frequencies of an alignment: GC content, gaps, 'N's.

For the gap count to be meaningful, you should remove columns containing only gaps
(use `ungap.py`)
"""


from sys import stdin
import numpy as np
import argparse
from Bio import AlignIO, SeqIO
from collections import Counter
from seqtools.symbols import GAPS, NUCLEOTIDES, NUCL_UNKNOWN, NUCL_AMBIGUOUS, CODONS_STOP, AA, AA_UNKNOWN, AA_AMBIGUOUS

import logging
logger = logging.getLogger(__name__)


def weighted_std(values, weights, axis=None):
    average = np.average(values, axis, weights=weights)
    variance = np.average((values-average)**2, axis, weights=weights)
    return np.sqrt(variance)


def get_seq_counts(alignment, is_aa=False):
    """Counts for each sequence"""
    alphabet = AA if is_aa else NUCLEOTIDES
    alphab_len = len(alphabet)
    unknown = AA_UNKNOWN if is_aa else NUCL_UNKNOWN
    ambiguous = AA_AMBIGUOUS if is_aa else NUCL_AMBIGUOUS
    seq_resid = np.zeros((alphab_len, len(alignment)))  # counts of A, C, G, T or amino-acids
    seq_gap   = np.zeros(len(alignment))
    seq_N     = np.zeros(len(alignment))
    if is_aa:
        seq_CpG, seq_stops = None, None
    else:
        seq_CpG   = np.zeros(len(alignment))
        seq_stops = np.zeros(len(alignment))

    try:
        length = alignment.get_alignment_length()
    except AttributeError:
        length = [len(record) for record in alignment]

    for i, record in enumerate(alignment):

        counts = Counter(record.seq.upper())

        seq_resid[:, i] = [counts[n] for n in alphabet]
        seq_gap[i]  = sum(counts[g] for g in GAPS)
        seq_N[i]    = counts[unknown]
        for a, possible in ambiguous.items():
            seq_N[i] += counts[a] * (1. - len(possible)/alphab_len)
            for p in possible:
                seq_resid[alphabet.index(p), i] += counts[a] / alphab_len

        if seq_resid[:, i].sum() + seq_gap[i] + seq_N[i] != len(record):
            raise ValueError('Residue counts do not add up: at #%d %s: len=%d VS %d resid + %d gaps %d unknown' %(i, record.name, len(record), seq_resid[:,i].sum(), seq_gap[i], seq_N[i]))
        if not is_aa:
            seq_stops[i] = sum(str(record.seq[k:k+3]) in CODONS_STOP for k in range(0,len(record),3))
            seq_CpG[i]  = record.seq.upper().count('CG')

    return length, seq_resid, seq_gap, seq_N, seq_CpG, seq_stops


def get_seq_freqs(length, seq_nucl, seq_gap, seq_N, seq_CpG=None, seq_stops=None):
    """composition frequencies **per sequence**"""
    seq_nucltot = seq_nucl.sum(axis=0)
    seq_lengths = seq_nucltot + seq_N
    
    seq_nucl_freq = seq_nucl.astype(float) / seq_nucltot
    seq_N_freq    = seq_N.astype(float)    / seq_lengths
    seq_gap_freq  = seq_gap.astype(float)  / length
    freqs = (seq_lengths, seq_gap_freq, seq_N_freq) + tuple(seq_nucl_freq)
    try:
        seq_stop_freq = seq_stops * 3. / seq_lengths
        seq_CpG_freq  = seq_CpG * 2. / seq_nucltot  # Is the *2 necessary?
        GC_index = [NUCLEOTIDES.index('C'), NUCLEOTIDES.index('G')]
        seq_GC_freq   = seq_nucl_freq[GC_index].sum(axis=0)
        freqs += (seq_GC_freq, seq_CpG_freq, seq_stop_freq)
    except TypeError:
        pass

    return freqs


def get_al_compo_summary(length, seq_counts, seq_freqs):
    """Stats over the whole alignment
    
    Return: for each summary type, a list of values for each compositional element.
    """
    #seq_nucl, seq_gap, seq_N, seq_CpG = seq_counts
    seq_nucl, seq_gap, seq_N, seq_CpG, seq_stops = seq_counts
    
    seq_lengths = seq_freqs[0]

    nucl_count = seq_nucl.sum(axis=1)  # total number of each nucleotide (summed over sequences)
    nucl_tot   = nucl_count.sum()
    gap_count  = seq_gap.sum()  # seq_gap.mean() --> / length
    N_count    = seq_N.sum()

    global_stats = (np.sum(length),
                    float(gap_count) / np.array(length),
                    float(N_count)   / (nucl_tot+N_count)) \
                   + tuple(nucl_count.astype(float) / nucl_tot)
    try:
        CpG_count  = seq_CpG.sum()
        stop_count = seq_stops.sum()
        GC_count   = nucl_count[NUCLEOTIDES.index('C')] + nucl_count[NUCLEOTIDES.index('G')]
        global_stats += (float(GC_count)      / nucl_tot,
                         float(CpG_count) * 2 / nucl_tot,
                         float(stop_count)* 3 / nucl_tot)
    except AttributeError as err:
        if 'NoneType' in err.args[0]:
            CpG_count = stop_count = GC_count = np.NaN
        else:
            raise

    seq_means   = [s.mean()     for s in seq_freqs]
    seq_medians = [np.median(s) for s in seq_freqs]
    seq_stds    = [s.std()      for s in seq_freqs]

    # Same as glob in all cases but "len" and "gaps".
    # Both variables don't need to be weighted by the number of non-gaps.
    seq_w_means = [np.average(s, weights=seq_lengths)   for s in seq_freqs]

    seq_w_stds  = [weighted_std(s, weights=seq_lengths) for s in seq_freqs]

    # TODO: return a np.array
    return global_stats, seq_means, seq_medians, seq_stds, seq_w_means, seq_w_stds


def make_al_compo(alignment, byseq=False, is_aa=False):
    length, *seq_counts = get_seq_counts(alignment, is_aa)

    seq_freqs = get_seq_freqs(length, *seq_counts)
    
    if byseq:
        stats = np.array(seq_freqs).T
        stat_names = [r.name for r in alignment]  # 1 seq per row
    else:
        stats = get_al_compo_summary(length, seq_counts, seq_freqs)
        stat_names = ['glob', 'mean', 'med', 'std', 'w_mean', 'w_std']  # summary stats
    return stat_names, stats


def print_al_compo(alignment_file, format='fasta', byseq=False, is_aa=False):
    sequences = list(SeqIO.parse(alignment_file, format)) if byseq else AlignIO.read(alignment_file, format=format)
    stat_names, stats = make_al_compo(sequences, byseq, is_aa)

    n_stats = len(stats[0])
    
    # column names
    header = ('', 'len', 'gaps')
    if is_aa:
        header += ('X',) + tuple(AA)
    else:
        header += ('N',) + tuple(NUCLEOTIDES) + ('GC', 'CpG', 'stops')
    print('\t'.join(['%8s']*(n_stats+1)) % header)

    row_template = '\t'.join(['%8s'] + ['%8g']*n_stats)

    for desc, stat in zip(stat_names, stats):
        print(row_template % tuple([desc] + list(stat)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('alignment_file', nargs='?',
                        type=argparse.FileType('r'), default=stdin)
    parser.add_argument('-f', '--format', default='fasta', help='[%(default)s]')
    parser.add_argument('-s', '--byseq', action='store_true',
                        help='Show stats for each sequence. [Default: by alignment site]')
    parser.add_argument('-a', '--aa', action='store_true', dest='is_aa',
                        help='Interpret input sequences as amino-acids')

    args = parser.parse_args()
    print_al_compo(**vars(args))


if __name__ == '__main__':
    main()
