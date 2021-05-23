#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from sys import stdin
import argparse as ap
from Bio import AlignIO


def main():
    p = ap.ArgumentParser(description=__doc__)
    p.add_argument('infiles', nargs='*', help='[stdin]')
    p.add_argument('-f', '--fmt', default='fasta', help='[%(default)s]')
    p.add_argument('-l', '--from-list', action='store_true',
                   help='Each line of the input files designate an alignment file')

    args = p.parse_args()

    if args.from_list:
        if not args.infiles:
            infiles = [line.rstrip() for line in stdin if not line.startswith('#')]
        else:
            infiles = []
            for infile in args.infiles:
                with open(infile) as f:
                    for line in f:
                        if not line.startswith('#'):
                            infiles.append(line.rstrip())
    else:
        infiles = args.infiles if args.infiles else [stdin]

    for infile in infiles:
        lengths = [str(al.get_alignment_length()) for al in AlignIO.parse(infile, args.fmt)]
        print(('-' if (infile is stdin) else str(infile)) + '\t' + ' '.join(lengths))


if __name__ == '__main__':
    main()

