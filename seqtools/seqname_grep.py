#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Select sequences from file (e.g. fasta) based on their name.
"""


from sys import stdout, stdin
import re
import argparse
from Bio import SeqIO, Seq, AlignIO, Align


def seqrecords_grep(records_iterable, pattern, negate=False):
    reg = re.compile(pattern)
    
    for record in records_iterable:
        match = reg.search(record.name)
        if (match and not negate) or (not match and negate):
            yield record


def algrep(alignment, pattern, negate=False):
    out_records = seqrecords_grep(alignment, pattern, negate)
    return Align.MultipleSeqAlignment(out_records)


def main(pattern, filename, negate=False, fromfile=False, fixed=False, word=False, label=False, seq_format='fasta', outfile=stdout):
    #al = AlignIO.read(alignment_file, format=seq_format)

    if label:
        begin, end = r'^', r'$'
    elif word:
        begin, end = r'\b', r'\b'
    else:
        begin, end = '', ''
    if fixed:
        def wrap_pattern_elem(pat):
            return begin + re.escape(pat) + end
    else:
        def wrap_pattern_elem(pat):
            return begin + pat + end

    if fromfile:
        if pattern == '-':
            pattern = '|'.join(wrap_pattern_elem(line.rstrip()) for line in stdin)
        else:
            with open(pattern) as f:
                pattern = '|'.join(wrap_pattern_elem(line.rstrip()) for line in f)
    sequences = SeqIO.parse(filename, format=seq_format)
    SeqIO.write(seqrecords_grep(sequences, pattern, negate), outfile, seq_format)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pattern')
    parser.add_argument('filename', nargs='?', type=argparse.FileType('r'),
                        default=stdin)
    parser.add_argument('-o', '--outfile', nargs='?', type=argparse.FileType('w'),
                        default=stdout)
    parser.add_argument('-f', '--fromfile', action='store_true',
                        help='Read pattern from file (one alternative per line).')
    parser.add_argument('-F', '--fixed', action='store_true', help='Disable regex')
    parser.add_argument('-w', '--word', action='store_true', help='Only match whole words')
    parser.add_argument('-x', '--label', action='store_true', help='Only match whole label')
    parser.add_argument('-v', '--negate', action='store_true')
    parser.add_argument('-s', '--seq-format', default='fasta')

    args = parser.parse_args()
    main(**vars(args))
