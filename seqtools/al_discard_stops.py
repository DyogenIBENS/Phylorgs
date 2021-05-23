#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from sys import stdin, stderr, stdout
import argparse as ap
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq


stops = set(('TAA', 'TAG', 'TGA',))


def discard_stops(al, to_gaps=False, no_ambiguous=False):
    labels = [record.name for record in al]

    new_al = al[:,:0]
    prev_start = 0
    n_stop_cols = 0

    if no_ambiguous:
        def filter_codon(codon):
            return 'N' in codon or codon in stops
    else:
        def filter_codon(codon):
            return codon in stops

    for i in range(0, al.get_alignment_length(), 3):
        column = list(al[:,i:(i+3)]) # list() to be able to assign.

        is_stop = [filter_codon(str(c.seq)) for c in column]
        if any(is_stop):
            n_stop_cols += 1
            for stop_j, check in enumerate(is_stop):
                if check:
                    print('Codon %d, seq %d %r (%s)' %(i//3, stop_j, labels[stop_j], column[stop_j].seq),
                            file=stderr)
                    column[stop_j].seq = Seq('---')

            new_al += al[:,prev_start:i]
            if to_gaps:
                new_al += MultipleSeqAlignment(column)
            prev_start = i+3

    if not any(is_stop):
        # Only if the last column is not a stop:
        # Include the remaining columns, up to the last one.
        i+=3
        new_al += al[:,prev_start:i]

    L = new_al.get_alignment_length()
    if to_gaps:
        assert L == al.get_alignment_length(), 'new length %d ' % L
    else:
        assert L + 3*n_stop_cols == al.get_alignment_length(), 'new length = %d; n_stop_cols = %d' % (L, n_stop_cols)
    print('new length = %d; n_stop_cols = %d' % (L, n_stop_cols), file=stderr)

    return new_al


def main():
    p = ap.ArgumentParser(description=__doc__)
    p.add_argument('infile', nargs='?', default=stdin, help='[stdin]')
    p.add_argument('-f', '--fmt', default='fasta')
    p.add_argument('-g', '--to-gaps', action='store_true')
    p.add_argument('-a', '--no-ambiguous', action='store_true')

    args = p.parse_args()

    al = AlignIO.read(args.infile, args.fmt)

    new_al = discard_stops(al, args.to_gaps, args.no_ambiguous)

    AlignIO.write(new_al, stdout, args.fmt)


def test_discard_stops():
    from io import StringIO
    al = AlignIO.read(StringIO('>s1\nCCCTAA\n>s2\nGGGGGG\n'), 'fasta')
    print('# TEST last is stop, to_gaps=True', file=stderr)
    print(discard_stops(al, to_gaps=True).format('fasta'))
    print('# TEST last is stop, to_gaps=False', file=stderr)
    print(discard_stops(al, to_gaps=False).format('fasta'))
    al = AlignIO.read(StringIO('>s1\nCCCTAATTT\n>s2\nGGGGGGCCC\n'), 'fasta')
    print('# TEST stop in the middle, to_gaps=True', file=stderr)
    print(discard_stops(al, to_gaps=True).format('fasta'))
    print('# TEST stop in the middle, to_gaps=False', file=stderr)
    print(discard_stops(al, to_gaps=False).format('fasta'))

    al = AlignIO.read(StringIO('>s1\nCCCTTNTTT\n>s2\nGGGGGGCCC\n'), 'fasta')
    print('# TEST ambiguous in the middle, to_gaps=True', file=stderr)
    print(discard_stops(al, to_gaps=True, no_ambiguous=True).format('fasta'))

if __name__ == '__main__':
    main()
