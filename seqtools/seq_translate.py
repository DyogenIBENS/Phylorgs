#!/usr/bin/env python3
# -*- coding: utf-8 -*-


_epilog = ("Note: if input format is \"evolver\"/\"multifasta\" and output format requires "
"one file per alignment (e.g."
"regular fasta or phylip), the outfile argument should include a formatter string:\n\n"
"Example: ./seq_translate.py -f evolver mc.paml output_prot_%02d.phy")


import sys
import argparse
from Bio import SeqIO


try:
    from seqtools.seqconv import parse_filename_range, split_multidata, split_evolver
except ImportError:
    print('seqtools.seqconv.split_evolver not available for PAML input format.', file=sys.stderr)


def iter_translate(records, format="fasta"):
    for record in records:
        tr_record = record.translate(gap='-', id=True, name=True, description=True)
        assert len(tr_record.seq) > 0
        yield tr_record


def translate(inputfile, outfile, input_format="fasta", output_format='fasta'):
    if inputfile == '-':
        inputfile = sys.stdin
    SeqIO.write(iter_translate(SeqIO.parse(inputfile, input_format)), outfile, output_format)


def main():
    parser = argparse.ArgumentParser(description=__doc__, epilog=_epilog)
    parser.add_argument('inputfile', nargs='?', default=sys.stdin)
    parser.add_argument('outfile', nargs='?', default=sys.stdout)
    parser.add_argument('-if', '--input-format', default="fasta", 
                        help=('Formats for Biopython (fasta, phylip, phylip-'
                        'sequential, phylip-relaxed, etc.) or "evolver" for '
                        'the output of PAML evolver, and "multifasta" for '
                        'concatenated fastas [%(default)s]'))
    parser.add_argument('-of', '--output-format', default="fasta",
                        help='[%(default)s]')

    args = parser.parse_args()
    if args.input_format == 'evolver' or args.input_format.startswith('multifa'):
        # This is just a concatenated Phylip/fasta, with extra blank lines.
        if args.input_format == 'evolver':
            iterdata = split_evolver
            datumformat = 'phylip-sequential'
        else:
            iterdata = split_multidata
            datumformat = 'fasta'
        try:
            if args.outfile == sys.stdout:
                out = sys.stdout
                joint_out = True
            elif args.output_format == 'evolver':
                out = open(args.outfile, 'w')
                joint_out = True
            else:
                joint_out = False

            infile, start, end = parse_filename_range(args.inputfile)
            if infile == '-' or infile == sys.stdin:
                f = sys.stdin
            else:
                f = open(infile)
            for n_al, aldata in enumerate(iterdata(f, start=start, end=end)):
                if not joint_out:
                    out = open(args.outfile % n_al, 'w')
                else:
                    out.write('\n\n')
                #print(''.join(line[:20].rstrip() + '\n' for line in al_lines))
                translate(aldata, out, datumformat, args.output_format)
                if joint_out:
                    out.write('\n\n')
                else:
                    out.close()
        except BaseException as err:
            if args.outfile == sys.stdout:
                f.close()
                out.close()
            raise

    else:
        translate(**vars(args))


if __name__=='__main__':
    main()
