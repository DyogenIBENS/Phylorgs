#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from sys import stdin, stdout, stderr
import argparse
import re
from io import StringIO
from Bio import SeqIO, AlignIO


EXT2FMT = {'fa': 'fasta', 'fasta': 'fasta', 'mfa': 'fasta',
           'phylip': 'phylip-relaxed',
           'nx': 'nexus',
           'paml': 'evolver'}


PHYLIP_HEAD_REG = re.compile(r'^\s*\d+\s\d+$')


def guess_format(filename, default='fasta'):
    if filename in (stdin, stdout, '-'):
        return default
    ext = filename.rsplit('.', 1)[1].split(':')[0]
    # specifications for the number of alignments are specified after colon
    return EXT2FMT[ext]


def split_evolver(f, start=0, end=None):
    al_lines = []
    n_al = 0
    for line in f:
        if line.rstrip():
            if PHYLIP_HEAD_REG.match(line.rstrip()):
                #print('Phylip start (len(al_lines) = %d)' % len(al_lines))
                if al_lines:
                    assert len(al_lines) == nseq+1

                    if n_al>=start:
                        #yield SeqIO.parse(StringIO(''.join(al_lines)), 'phylip-sequential')
                        yield StringIO(''.join(al_lines))

                    al_lines = []
                    n_al += 1
                    if end is not None and n_al > end:
                        print('Parsed %d alignments.' % (n_al-start), file=stderr)
                        return
                nseq, length = [int(x) for x in line.split()]
                al_lines = [line]
            else:
                #print('+ %r...' % line[:20])
                al_lines.append(line)
    assert len(al_lines) == nseq+1
    #print(''.join(line[:20].rstrip() + '\n' for line in al_lines))

    yield StringIO(''.join(al_lines))
    print('Parsed %d alignments.' % (n_al+1), file=stderr)


def write_phylip_sequential_relaxed(al, outfile):
    """Bio implements:
    - phylip (strict, interleaved),
    - phylip-sequential (strict, sequential),
    - phylip-relaxed (relaxed, interleaved),

    but not a relaxed sequential (as needed for PhyloBayes).
    
    Write sequences on a single line.
    """
    
    if outfile is not stdout:
        outfile = open(outfile)
    try:
        outfile.write('%d\t%d\n' % (len(al), al.get_alignment_length()))
        for record in al:
            outfile.write('%s  %s\n' % (record.name, record.seq))
    finally:
        if outfile is not stdout:
            outfile.close()


def convert_one_dataset(infile, out, fro='fasta', to='phylip-relaxed'):
    if to == 'phylip-sequential-relaxed':
        write_phylip_sequential_relaxed(AlignIO.read(infile, fro), outfile)
    else:
        AlignIO.convert(infile, fro, out, to)


def main(infile, outfile, fro=None, to=None):
    if fro is None:
        fro = guess_format(infile)
    if to is None:
        to = guess_format(outfile, 'phylip-relaxed')

    if infile == '-':
        infile = stdin
    if outfile not in (stdout, '-'):
        out = open(outfile, 'w')
    else:
        out = stdout

    if fro == 'evolver':
        ext = infile.rsplit('.', 1)[1]
        if ':' in ext:
            infile, selection = infile.rsplit(':', 1)
            start, end = re.match(r'(\d+)(?:-(\d*))?$', selection).groups()
            start = int(start)
            if end is None:
                end = start
            elif end == '':
                end = None  # read the entire file
            else:
                end = int(end)
        else:
            start, end = 0, None
        with open(infile) as f:
            for aldata in split_evolver(f, start, end):
                convert_one_dataset(aldata, out, 'phylip-sequential', to)
                out.write('\n\n')
    else:
        convert_one_dataset(infile, out, fro, to)

    if outfile not in (stdout, '-'):
        out.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile', nargs='?', #type=argparse.FileType('r'),
                        default=stdin)
    parser.add_argument('outfile', nargs='?', #type=argparse.FileType('w'),
                        default=stdout)
    parser.add_argument('-f', '--from', dest='fro',
                        help='[guess from extension, or fasta if stdin]')
    parser.add_argument('-t', '--to',
                        help='[guess from extension, or phylip-relaxed if stdout]')
    
    args = parser.parse_args()
    main(**vars(args))

