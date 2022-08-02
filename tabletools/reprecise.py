#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from sys import stdin
import argparse as ap
from collections import OrderedDict
from itertools import zip_longest
import bz2
import logging
logger = logging.getLogger(__name__)


def main():
    logging.basicConfig(format=logging.BASIC_FORMAT)
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('infile', help="'-' for stdin")
    parser.add_argument('-F', '--formats', default='',
                        help='space-separated list of instructions column-range:format')
    parser.add_argument('-H', '--header', action='store_true', help='First row as header')
    parser.add_argument('-s', '--sep', default='\t', help="['\\t']")
    parser.add_argument('-d', '--default', default='%s', help="Default format [%(default)s]")
    args = parser.parse_args()
    sep = args.sep

    formatters = OrderedDict()
    lastcol = 0
    for cmd in args.formats.split():
        cols, fmt = cmd.split(':')
        if fmt in formatters:
            logger.warning('Duplicate format given: %r. Columns might be updated.', fmt)
        else:
            formatters[fmt] = []
        for colrange_part in cols.split(','):
            try:
                start, end = colrange_part.split('-', maxsplit=1)
                start = 0 if not start else int(start)
                end = None if not end else int(end)
                if end is not None and end > lastcol:
                    lastcol = end
                elif end is None and start >= lastcol:
                    lastcol = start + 1
                formatters[fmt].append(slice(start, end))
            except ValueError:
                col = int(colrange_part)
                if col > lastcol:
                    lastcol = col
                formatters[fmt].append(col)

    row_fmts = [args.default] * (lastcol)
    tail_fmt = args.default
    for fmt, locations in formatters.items():
        for colrange in locations:
            if isinstance(colrange, slice):
                if colrange.stop is None:
                    colrange = slice(colrange.start, lastcol)
                    tail_fmt = fmt
                for col in range(colrange.start, colrange.stop):
                    row_fmts[col] = fmt
    logger.debug('row_fmts = %s, tail_fmt = %s', ' '.join(row_fmts), tail_fmt)

    f = stdin if args.infile == '-' else open(args.infile)
    if args.header:
        print(next(f), end='')
    for line in f:
        if line.lstrip().startswith('#'):
            print(line, end='')
        else:
            row = line.rstrip().split(sep)
            for i in range(len(row)):
                try:
                    row[i] = int(row[i])
                except ValueError:
                    try:
                        row[i] = float(row[i])
                    except ValueError:
                        pass

            if len(row) > len(row_fmts):
                row_fmts += [tail_fmt] * (len(row) - len(row_fmts))
            #logger.debug('len(row) = %d ; len(row_fmts) = %d', len(row), len(row_fmts))
            #logger.debug('row: %s', row)

            print(sep.join(fmt % val for fmt,val in zip(row_fmts, row)))


if __name__ == '__main__':
    main()
