#!/usr/bin/env python3


"""
Perform Chi² test of compositional heterogeneity of biological sequences:
    a sequence composition VS the mean composition
"""


from sys import stdin, stderr, exit
import warnings
import bz2
import argparse as ap
import numpy as np
from collections import Counter
from scipy.stats import chisquare
from seqtools.symbols import GAPS, AA, AA_UNKNOWN, AA_AMBIGUOUS, NUCLEOTIDES, NUCL_UNKNOWN, NUCL_AMBIGUOUS
from seqtools.arrayal import AA2INT, NUCL2INT


def iter_multifasta(seqfile):
    f = stdin if seqfile == '-' else bz2.open(seqfile) if seqfile.endswith('.bz2') else open(seqfile)
    line = next(f)
    lineno = 0
    while not line.startswith('>'):
        line = next(f)
        lineno += 1
    if lineno > 0:
        print('WARNING: lines 1-%d are not a sequence label: SKIP.' % lineno, file=stderr)

    header = line.rstrip()[1:]

    seq = ''
    for line in f:
        if line.startswith('>'):
            yield header, seq
            header = line.rstrip()[1:]
            seq = ''
        elif line.rstrip() == r'\\':
            # Sometimes used to separate multiple sequence alignments inside a single fasta file.
            # See seqtools.seqconv.split_multidata
            #TODO
            continue
        else:
            seq += line.strip()

    if seq:
        yield header, seq


def get_residue_numbers(seqfile, as_freq=False):
    """Parse proportions of residues per sequence. Assume fasta input.
    If as_freq=True, convert to frequencies instead of absolute counts"""
    seqlengths = {}
    seqcounts = {}
    seen_residues = set()
    seen_labels = []  # Mostly to store the order of sequences
    for header, seq in iter_multifasta(seqfile):
        label, description = header.split(maxsplit=1)
        if label in seqlengths:
            raise ValueError('Duplicated sequence label at %r' % header)
        seen_labels.append(label)
        seqcounts[label] = Counter(seq.upper())
        seen_residues.update(seqcounts[label])
        Ngaps = sum(seqcounts[label].pop(g, 0) for g in GAPS)
        seqlengths[label] = len(seq) - Ngaps

    # Autodetect nucleotide or AA
    seen_residues.difference_update(GAPS)
    if seen_residues <= set(NUCLEOTIDES + NUCL_UNKNOWN).union(NUCL_AMBIGUOUS):
        converter = NUCL2INT
    elif seen_residues <= set(AA + AA_UNKNOWN).union(AA_AMBIGUOUS):
        converter = AA2INT
    else:
        raise ValueError('sequence characters not allowed: ' + ' '.join(seen_residues.difference(AA + AA_UNKNOWN + AA_AMBIGUOUS)))

    max_resid_id = max(converter.values()) + 1  # corresponds to an unknown residue (N for DNA, X for peptide)
    #print('DEBUG: max_resid_id=%d' % max_resid_id, file=stderr)
    residnumbers = np.zeros((len(seen_labels), max_resid_id), dtype=(float if as_freq else np.uint))
    for i, label in enumerate(seen_labels):
        for resid, count in seqcounts[label].items():
            j = converter.get(resid, max_resid_id) - 1  # -1 because we won't store gaps
            residnumbers[i, j] = count
            if as_freq:
                residnumbers[i, j] /= seqlengths[label]
    return seen_labels, seqlengths, residnumbers


def compo_test_each_row(mat, average_without_self=False, weight=False):
    # Drop residues never found (columns of zeros)
    if (mat > 1.).any() and (mat[0 < mat] < 1.).any():
        raise ValueError('Input numbers must be either frequencies or absolute counts, not a mix of both')
    elif (mat[0 < mat] < 1.).any():
        raise NotImplementedError('Frequencies given, use counts instead.')

    # Weight values by the sequence lengths
    totals = mat.sum(axis=1)
    weights = totals if weight else np.ones(mat.shape[0], dtype=np.uint8)

    mat = mat[:,(mat>0).any(axis=0)]
    #print('DEBUG: mat.shape = %s,%s' % mat.shape, file=stderr)

    rowfreqs = mat / totals[:,np.newaxis]

    if average_without_self:
        rowidx = list(range(mat.shape[0]))
        def ref_compo(i):
            kept = rowidx[:i] + rowidx[i+1:]
            return np.average(rowfreqs[kept,:], axis=0, weights=weights[kept]) * totals[i]
    else:
        average_freq = np.average(rowfreqs, axis=0, weights=weights)
        def ref_compo(i): return average_freq * totals[i]

    test = chisquare
    # Verify that there are no frequencies too low (would be better with original counts)
    lowerlimit = 5
    underlimit = (mat < lowerlimit).sum(axis=1)
    if underlimit.any():
        warnings.warn('%d sequences have low frequency residues (< 5 counts). Maybe use an exact test.' % (underlimit>0).sum())
        # Auto switch to exact test?

    pvalues = []
    for i, row in enumerate(mat):
        #print('DEBUG: # %s    row size %s VS ref_compo size %s\nrow=%s sum=%s\nref_compo=%s sum=%s' % (
        #      i, row.shape, ref_compo(i).shape,
        #      ', '.join('%g' % x for x in row), row.sum(),
        #      ', '.join('%g' % x for x in ref_compo(i)), ref_compo(i).sum()), file=stderr)
        pvalues.append(test(row, ref_compo(i)).pvalue)
        #if i > 4: break

    return np.array(pvalues)


def compo_hetero(seqfile, average_without_self, weight_length=False, correct_multi=None, exact=False, plot=False):
        labels, lengths, resid_nb = get_residue_numbers(seqfile, exact)
        #print('DEBUG:%d sequences x %d residues' % resid_nb.shape, file=stderr)
        #print('AA = %s (n=%d)' % (AA, len(AA)), file=stderr)
        #print('NUCL = %s' % NUCLEOTIDES, file=stderr)

        if plot:
            import matplotlib as mpl
            mpl.use('TkAgg')
            import matplotlib.pyplot as plt
            im = plt.pcolormesh(resid_nb.T, cmap='RdPu_r', vmin=0)
            ax = plt.gca()
            plt.colorbar(im)
            ax.set_ylabel('Residue')
            ax.set_yticks(np.arange(resid_nb.shape[1]) + 0.5)
            residues = NUCLEOTIDES + 'N' if resid_nb.shape[1] < 6 else AA + 'X'
            ax.set_yticklabels(list(residues))
            ax.invert_yaxis()  # Letters ordered from top to bottom
            ax.set_xticks(np.arange(len(labels)) + 0.5)
            ax.set_xticklabels(labels, fontdict=dict(fontsize='x-small', rotation=45, va='top', ha='right'))
            ax.figure.set_size_inches(16,8)
            plt.tight_layout()
            plt.show(block=True)

        with warnings.catch_warnings(record=True) as warned:
            warnings.simplefilter('always')  # Always warn
            pvalues = compo_test_each_row(resid_nb, average_without_self, weight_length)
            for w in warned:
                print('WARNING:%s: %s' % (seqfile, w.message), file=stderr)

        if correct_multi:
            from statsmodels.stats.multitest import multipletests
            reject, adj_pvalues, _, _ = multipletests(pvalues, method=correct_multi)
            pvalues = adj_pvalues

        return labels, pvalues


def main():
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('seqfiles', nargs='*', help="Fasta files. '-' for stdin, or bzipped files are accepted. Multiple files are analysed independently.")
    parser.add_argument('-c', '--count', action='store_true', help='Output the count of rejected sequences [default: list sequences]')
    parser.add_argument('-w', '--weight-length', action='store_true',
                        help='The average composition is weighted by the sequence lengths')
    parser.add_argument('-L', '--leave-one-out-average', action='store_true',
                        help='Leave the sequence being tested out of the average composition')
    parser.add_argument('-m', '--multi-test-correction',
                        help=('method name (see statsmodels.stats.multitest): ',
                            'bonferroni, holm, fdr_bh (Benjamini-Hochberg)'))
    parser.add_argument('-s', '--signif-cutoff', default=0.05, type=float)
    parser.add_argument('-e', '--exact', '--exact-test', action='store_true')
    parser.add_argument('-p', '--plot', '--plot-freqs', action='store_true')
    # TODO: average weighted based on sequence proximity in a tree?
    args = parser.parse_args()

    seqfiles = ['-'] if not args.seqfiles else args.seqfiles

    failed = 0  # Number of alignments containing at least one rejected sequence
    for n, seqfile in enumerate(seqfiles):
        labels, pvalues = compo_hetero(seqfile, args.leave_one_out_average, args.weight_length, args.multi_test_correction, args.exact, args.plot)
        nonhomogeneous = np.flatnonzero(pvalues<args.signif_cutoff)

        prefix = seqfile + ': ' if len(seqfiles)>1 else ''
        if args.count:
            print(prefix + str(len(nonhomogeneous)))
        elif len(nonhomogeneous):
            print('# %s%d non-homogeneous sequences:' % (prefix, len(nonhomogeneous)))
            for i in nonhomogeneous:
                print('%s\t%g' % (labels[i], pvalues[i]))
            failed += 1
        else:
            print('# %s0 non-homogeneous sequence.' % prefix)
    #return failed


if __name__ == '__main__':
    exit(main())

