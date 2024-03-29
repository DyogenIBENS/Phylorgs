#!/usr/bin/env python3
# -*- coding: utf-8 -*-

_EPILOG="""
For multidata files (formats "evolver"/"multifasta"), a range can be specified
at the end of the input filename:

filename:0    (the first dataset)
filename:5-10 (the 5 to 10 excluded)
filename:5-   (the 5th to the end)
"""

from sys import stdin, stdout, stderr
import os.path as op
import argparse
import re
from io import StringIO
from Bio import SeqIO, AlignIO

import logging
logger = logging.getLogger(__name__)


EXT2FMT = {'fa': 'fasta',
           'mfa': 'multifasta',
           'phy': 'phylip-sequential-relaxed', # 'phylip': 'phylip-relaxed',
           'nx': 'nexus', 'nex': 'nexus',
           'paml': 'evolver', 'mphy': 'evolver',
           'sto': 'stockholm' # input not implemented yet
           }


PHYLIP_HEAD_REG = re.compile(r'^\s*\d+\s+\d*$')  # nseq, nsites


def guess_format(filename, default='fasta'):
    if filename in (stdin, stdout, '-'):
        return default
    ext = op.splitext(filename)[1].split(':')[0].lstrip('.')
    if not ext:
        return default
    # specifications for the number of alignments are specified after colon
    return EXT2FMT.get(ext, ext)


def parse_filename_range(infile):
    try:
        infile, selection = infile.rsplit(':', 1)
        # ValueError if there is no ':'
        start, end = re.match(r'(\d+)(?:-(\d*))?$', selection).groups()
        start = int(start)
        if end is None:
            end = start
        elif end == '':
            end = None  # read the entire file
        else:
            end = int(end)
    except ValueError:
        start, end = 0, None
    return infile, start, end


def split_evolver(f, start=0, end=None):
    al_lines = []
    n_al = 0
    nseq = -1
    for line in f:
        if line.rstrip():
            if PHYLIP_HEAD_REG.match(line.rstrip()):
                #print('Phylip start (len(al_lines) = %d)' % len(al_lines))
                if al_lines:
                    assert len(al_lines) == nseq+1

                    if n_al>=start:
                        #yield SeqIO.parse(StringIO(''.join(al_lines)), 'phylip-sequential')
                        logger.debug('Data #%d: %d input lines.', n_al, len(al_lines))
                        yield StringIO(''.join(al_lines))

                    al_lines = []
                    n_al += 1
                    if end is not None and n_al > end:
                        logger.info('Parsed %d alignments.' % (n_al-start))
                        return
                nseq, length = [int(x) for x in line.split()]
                al_lines = [line]
            else:
                #print('+ %r...' % line[:20])
                al_lines.append(line)
    assert len(al_lines) == nseq+1, "nseq=%s, %d al_lines at data #%d" % (nseq, len(al_lines), n_al)
    #print(''.join(line[:20].rstrip() + '\n' for line in al_lines))

    if n_al >= start:
        logger.debug('Data #%d: %d input lines.', n_al, len(al_lines))
        yield StringIO(''.join(al_lines))
        n_al += 1
    logger.info('Parsed %d alignments.' % n_al)


def split_multidata(f, sep='', start=0, end=None):
    al_lines = []
    n_al = 0
    for line in f:
        if line.rstrip() != sep:
            al_lines.append(line)
        else:
            if al_lines:
                if n_al>=start:
                    yield StringIO(''.join(al_lines))
                n_al += 1

                if end is not None and n_al > end:
                    logger.info('Parsed %d alignments.' % (n_al-start))
                    return
                al_lines = []
    if al_lines and n_al >= start:
        yield StringIO(''.join(al_lines))
        n_al += 1
    logger.info('Parsed %d alignments.' % (n_al-start))
    #print(''.join(line[:20].rstrip() + '\n' for line in al_lines))


def write_phylip_sequential_relaxed(al, outfile, stockholm_data=None):
    """Bio implements:
    - phylip (strict, interleaved),
    - phylip-sequential (strict, sequential),
    - phylip-relaxed (relaxed, interleaved),

    but not a relaxed sequential (as needed for PhyloBayes).
    
    Write sequences on a single line.
    """
    
    try:
        outfile = open(outfile, 'w')
        was_closed = True
    except TypeError:
        was_closed = False
        pass  # Already an opened file object.
    try:
        if stockholm_data is not None:
            # Stockholm formatted metadata (comments)
            outfile.write('# STOCKHOLM 1.0\n')
            outfile.write('\n'.join('#=%s %s' % item for item in stockholm_data) + '\n')
        else:
            outfile.write('%d\t%d\n' % (len(al), al.get_alignment_length()))
        for record in al:
            outfile.write('%s  %s\n' % (record.name, record.seq))
    finally:
        if was_closed:
            outfile.close()


def convert_one_dataset(infile, out, fro='fasta', to='phylip-sequential-relaxed'):
    if to in ('phylip-sequential-relaxed', 'stockholm'):
        metadata = None
        if to == 'stockholm' and isinstance(infile, (str, bytes)):
            metadata = [('GF ID', op.basename(op.splitext(infile)[0]))]
        write_phylip_sequential_relaxed(AlignIO.read(infile, fro), out, metadata)
    else:
        AlignIO.convert(infile, fro, out, to)


def main(infiles, outfile=None, fro=None, to=None):
    if to is None:
        to = guess_format(outfile, 'phylip-sequential-relaxed')

    n_al = 0  # Counter of the input alignments.

    if outfile is None:
        joint_out = True
        out = stdout
        if to == 'evolver': to = 'phylip-sequential-relaxed'
    elif to[:7] in ('evolver', 'multifa'):
        # FIXME: the above check is unnecessary: just specify a uniq outfile, and -to=phy
        to = 'phylip-sequential' if to == 'evolver' else 'fasta'
        joint_out = True
        out = open(outfile, 'w')
    #FIXME: output to Nexus requires an Alphabet to be specified for Biopython...
    else:
        try:
            outfile_n = outfile % n_al
            joint_out = False
        except TypeError:
            out = open(outfile, 'w')
            joint_out = True

    for infile in infiles:
        if fro is None:
            fro = guess_format(infile)
            if fro == 'phylip-sequential-relaxed':
                fro = 'phylip-relaxed'
        if fro.startswith('evolver') or fro.startswith('multifa'):
            iterdata = split_evolver if fro.startswith('evolver') else split_multidata
            srcfmt = 'phylip-relaxed' if fro.startswith('evolver') else 'fasta'
            # NOTE: PAML's Phylip parser allows up to 30 characters in sequence names.
            # In addition, it separates the name and sequence by 2 spaces.
            infile, start, end = parse_filename_range(infile)

            f = stdin if infile == '-' else open(infile)
            try:
                for aldata in iterdata(f, start=start, end=end):
                    if not joint_out:
                        out = open(outfile % n_al, 'w')
                    logger.debug('Data #%d: out=%s %s -> %s (joint=%s)', n_al, out.name, srcfmt, to, joint_out)
                    n_al += 1
                    convert_one_dataset(aldata, out, srcfmt, to)
                    if joint_out:
                        out.write('//\n\n' if to=='stockholm' else '\n\n')
                    else:
                        out.close()
            finally:
                if infile != '-': f.close()
        else:
            if infile == '-': infile = stdin
            if not joint_out:
                out = open(outfile % n_al, 'w')
            logger.debug('Data #%d: out=%s %s -> %s (joint=%s)', n_al, out.name, fro, to, joint_out)
            n_al += 1
            convert_one_dataset(infile, out, fro, to)
            if joint_out:
                out.write('//\n\n' if to=='stockholm' else '\n\n')
            else:
                out.close()

    if outfile is not None:
        out.close()


if __name__ == '__main__':
    logging.basicConfig(format='%(levelname)s:%(funcName)-20s:%(message)s')

    parser = argparse.ArgumentParser(description=__doc__, epilog=_EPILOG,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('infiles', nargs='*', #type=argparse.FileType('r'),
                        default=['-'], help='[stdin]')
    parser.add_argument('-o', '--outfile',
            help='Useful to specify splitted outfiles, e.g. out_%%02d.fa [stdout]')
    parser.add_argument('-f', '--from', dest='fro',
                        help='[guess from extension, or fasta if stdin]')
    parser.add_argument('-t', '--to',
                        help='[guess from extension, or phylip-sequential-relaxed if stdout]')
    parser.add_argument('-v', '--verbose', action='count', default=0)

    args = parser.parse_args()

    logger.setLevel(30 - 10*args.verbose)
    delattr(args, 'verbose')

    main(**vars(args))

