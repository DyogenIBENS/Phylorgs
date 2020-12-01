#!/usr/bin/env python3
# -*- coding: utf-8 -*-


_epilog = ("Note: if input format is \"evolver\" and output format requires "
"one file per alignment (e.g."
"regular fasta or phylip), the outfile argument should include a formatter string:\n\n"
"Example: ./seq_translate.py -f evolver mc.paml output_prot_%02d.phy")


import sys
import argparse
from Bio import SeqIO
from Bio.Alphabet import IUPAC

try:
    from seqtools.seqconv import split_evolver
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
    SeqIO.write(iter_translate(SeqIO.parse(inputfile, input_format, alphabet=IUPAC.ambiguous_dna)), outfile, output_format)


def main():
    parser = argparse.ArgumentParser(description=__doc__, epilog=_epilog)
    parser.add_argument('inputfile', nargs='?', default=sys.stdin)
    parser.add_argument('outfile', nargs='?', default=sys.stdout)
    parser.add_argument('-if', '--input-format', default="fasta", 
                        help=('Formats for Biopython (fasta, phylip, phylip-'
                        'sequential, phylip-relaxed, etc.) or "evolver" for '
                        'the output of PAML evolver. [%(default)s]'))
    parser.add_argument('-of', '--output-format', default="fasta",
                        help='[%(default)s]')

    args = parser.parse_args()
    if args.input_format == 'evolver':
        # This is just a concatenated Phylip, with extra blank lines.
        try:
            if args.outfile == sys.stdout:
                out = sys.stdout
                joint_out = True
            elif args.output_format == 'evolver':
                out = open(args.outfile, 'w')
                joint_out = True
            else:
                joint_out = False

            if args.inputfile == '-' or args.inputfile == sys.stdin:
                f = sys.stdin
            else:
                f = open(args.inputfile)
            for aldata in split_evolver(f):
                if not joint_out:
                    out = open(args.outfile % n_al, 'w')
                else:
                    out.write('\n\n')
                #print(''.join(line[:20].rstrip() + '\n' for line in al_lines))
                translate(aldata, out, 'phylip-sequential', args.output_format)
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
