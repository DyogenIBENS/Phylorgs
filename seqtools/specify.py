#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Append the species name to a gene name in a sequence alignment.
Can convert ensembl gene names to species, as well as assembly names
(ex: loxAfr3).

CAVEATS:

Because of Biopython Bio.SeqIO parser, spaces are not allowed in sequence
labels. To include the whole fasta identifier line, use record.description.

"""

EXAMPLES="""
EXAMPLES

# For short species names + shorter numeric ids:
  -t 'gene=s/^ENS(G|[A-Z]{3}G)0{5}|^MGP_(SPRET|CAROLI|Pahari)EiJ_G0|^WBGene00/G/'

# For parsing the PhylTree extended newick format (long labels):

  -i '^[.*]?(?P<sp>[A-Z][a-z _.-]+)(?:\|.*)?$'
"""

from sys import stdin, stdout
import os.path as op
import re
import argparse
from Bio import SeqIO
from genomicustools.identify import ultimate_seq2sp, convert_gene2species, SP2GENEID
from UItools.sed import parse_simple_subst


ENSEMBL_VERSION = 87

KEYS = ('gene', 'sp', 'shortsp', 'ENSG_PREFIX')  # Sorted by dependency order
DEFAULT_IN_FMT = r'(?P<gene>.*)'  # Use the named group feature.
# Here we capture the whole string as the gene name.
DEFAULT_OUT_FMT = TREEBEST_OUT_FMT = '{gene}_{sp}'
ALE_OUT_FMT = '{sp}_{gene}'



#TODO: move to genomicustools.identify
special_short2long = {
       'Scere': 'Saccharomyces cerevisiae',
       'Celeg': 'Caenorhabditis elegans',
       'Dmela': 'Drosophila melanogaster',
       'Csav' : 'Ciona savignyi',
       'Hsa'  : 'Homo sapiens',
       'Mspre': 'Mus spretus',
       'Mcaro': 'Mus caroli',
       'Mpaha': 'Mus pahari'
       }
special_long2short = {
       'Saccharomyces cerevisiae': 'Scere',
       'Caenorhabditis elegans'  : 'Celeg',
       'Drosophila melanogaster' : 'Dmela',
       'Ciona savignyi'          : 'Csav',
       'Homo sapiens'            : 'Hsa',
       'Mus spretus': 'Mspre',
       'Mus caroli' : 'Mcaro',
       'Mus pahari' : 'Mpaha',
       'Mus musculus': 'Mmusc'  # Because there is the 'Mus' genus!
       }

# The conversions follow a logical directed graph. It would be possible in practice to
# convert back in the other direction, but currently this requires to change the graph
# (ie. change the input key of each converter function)

# label parsing:
# sp_dot    -> sp
# sp_under  ->

# Then it will attempt conversions in this order: 'sp', 'shortsp', 'ENSG_PREFIX'

# gene -> sp -> shortsp
#            -> ENSG_PREFIX
# shortsp -> sp (if no 'gene')


def identify_sp(infos, ensembl_version=ENSEMBL_VERSION):
    try:
        return ultimate_seq2sp(infos['gene'], ensembl_version)
    except KeyError:
        try:
            return special_short2long[infos['shortsp']]
        except KeyError:
            return convert_gene2species('ENS'+infos['shortsp'].upper()+'G', ensembl_version)

def identify_shortsp(infos, ensembl_version):
    sp = infos['sp']
    return special_long2short.get(sp,
                                  SP2GENEID[ensembl_version][sp][3:6].capitalize())

def identify_ENSG_PREFIX(infos, ensembl_version):
    sp = infos['sp']
    return SP2GENEID[ensembl_version][sp]

identifiers = dict(sp=identify_sp, shortsp=identify_shortsp,
                   ENSG_PREFIX=identify_ENSG_PREFIX)


def parse_label(label, input_fmt=DEFAULT_IN_FMT, ensembl_version=ENSEMBL_VERSION):
    # NOTE: this is impossible to ask for an output gene label from a species label.
    try:
        infos = re.compile(input_fmt).match(label).groupdict()
    except AttributeError as err:
        err.args += ('label=%r' % label,)
        raise
    if 'sp' not in infos and 'sp_under' in infos:
        infos['sp'] = infos.pop('sp_under').replace('_', ' ')
    elif 'sp' not in infos and 'sp_dot' in infos:
        infos['sp'] = infos.pop('sp_dot').replace('.', ' ')
    for key in KEYS[1:]:
        if key not in infos:
            try:
                infos[key] = identifiers[key](infos, ensembl_version)
            except KeyError as err:
                err.args += ('%s is unidentifiable from label %r using pattern %r'
                             % (key, label, input_fmt),)
                raise
    return infos


DEFAULT_TRANSFORMS = dict(sp=lambda string: string.replace(' ', '.'))
# NOTE that CAFE v5 (Hahn et al) does not accept dots in species names...


# ~~> genomicustools.identify
def makelabel(elements, get_label, input_fmt=DEFAULT_IN_FMT, label_fmt=DEFAULT_OUT_FMT,
            transforms=DEFAULT_TRANSFORMS,
            ensembl_version=ENSEMBL_VERSION):  # uniq=True
    seen_labels = {}
    for elem in elements:
        infos = parse_label(get_label(elem), input_fmt, ensembl_version)
        for key, transform in transforms.items():
            infos[key] = transform(infos[key])
        label = label_fmt.format(**infos)

        if label in seen_labels:
            seen_labels[label] += 1
            label += '.%d' % (seen_labels[label] - 1)
        else:
            seen_labels[label] = 1
        yield elem, label, infos


def seq_specify(inputfile, file_fmt="fasta",
                input_fmt=DEFAULT_IN_FMT, label_fmt=DEFAULT_OUT_FMT,
                transforms=DEFAULT_TRANSFORMS,
                ensembl_version=ENSEMBL_VERSION):

    sequences = SeqIO.parse(inputfile, file_fmt)

    for record, label, _ in makelabel(sequences, lambda record: record.id,
                                    input_fmt, label_fmt, transforms,
                                    ensembl_version):
        record.id = label
        record.description = ''
        yield record


def tree_specify(tree, input_fmt=DEFAULT_IN_FMT, label_fmt=DEFAULT_OUT_FMT,
                 transforms=DEFAULT_TRANSFORMS,
                 ensembl_version=ENSEMBL_VERSION, inplace=True):
    if not inplace:
        tree = tree.copy()
    for leaf, label, infos in makelabel(tree.iter_leaves(), lambda leaf: leaf.name,
            input_fmt, label_fmt, transforms, ensembl_version):
        leaf.name = label
        leaf.add_feature('S', infos['sp'])

    if not inplace:
        return tree


EXT2FMT = {'.fa': 'fasta',
           '.phy': 'phylip-relaxed',
           '.nw': 'nwk',
           '.nwk': 'nwk',
           '.newick': 'nwk',
           '.tree': 'nwk'}


def specify(inputfile, outfile, input_fmt=DEFAULT_IN_FMT,
            label_fmt=DEFAULT_OUT_FMT, transform=None, file_fmt=None,
            ensembl_version=ENSEMBL_VERSION):

    if file_fmt is None:
        ext = op.splitext(inputfile)[1]
        file_fmt = EXT2FMT[ext]
    transforms = DEFAULT_TRANSFORMS
    if transform is not None:
        for arg in transform:
            key, val = arg.split('=', 1)
            transforms[key] = parse_simple_subst(val)

    if file_fmt == 'nwk':
        import ete3
        try:
            tree = ete3.Tree(stdin.read() if inputfile=='-' else inputfile, format=1)
        except ete3.parser.newick.NewickError as err:
            err.args = (err.args[0] + ' ERROR with treefile %s%s' % (
                inputfile[:70],
                '...'+inputfile[-10:] if len(inputfile)>80 else inputfile[70:]),)
            raise
        tree_specify(tree, input_fmt, label_fmt, transforms, ensembl_version)
        if outfile is stdout:
            print(tree.write(format=1))
        else:
            tree.write(outfile=outfile, format=1, format_root_node=True)
    else:
        if inputfile == '-':
            inputfile = stdin
        SeqIO.write(seq_specify(inputfile, file_fmt, input_fmt, label_fmt, transforms,
                                ensembl_version),
                    outfile, file_fmt)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, epilog=EXAMPLES,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('inputfile', help='"-" for stdin (requires -f).')
    parser.add_argument('outfile', nargs='?', default=stdout)
    parser.add_argument('-i', '--input-fmt', default=DEFAULT_IN_FMT,
                        help=('sequence label regular expression, with named '
                            'groups [%(default)r]. Valid keys: gene sp shortsp'
                            ' sp_dot sp_under.'))
    parser.add_argument('-l', '--label-fmt', default=DEFAULT_OUT_FMT,
            help='Valid extra keys: %s [%%(default)r]' % ' '.join(identifiers.keys()))
    parser.add_argument('-t', '--transform', action='append',
                        help='substitution commands applied to input variables [%(default)s]')
    parser.add_argument('-f', '--file-fmt',
                        help='[automatic detection of nwk/fasta/phylip unless using stdin]')
    parser.add_argument('-e', '--ensembl-version', type=int,
                        default=ENSEMBL_VERSION, help='[%(default)r]')

    args = parser.parse_args()
    specify(**vars(args))
