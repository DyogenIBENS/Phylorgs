#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Convert between miscellaneous phylogenetic tree file formats.

Details:
--------
Phylogenetic trees can be represented as string in various formats:
- Newick (many variants) -> Ete3 is a good tool for this.
- Nexus
- PhyloXML

Also, the embedding of attributes in the parenthesized representations may vary.
- NHX in newick trees: '[&&NHX:key=value:key2=value2]'
- another format in nexus trees such as produced by BEAST/LSD: [&key=value,key2=value2].

The support values may be placed inconsistently between formats/programs:
- after the node label / after the branch length.

Notes:
The Bio.Phylo.convert doesn't properly convert the comments from BEAST Nexus files.

"""


from sys import stdin, stdout
from functools import partial
from Bio import Phylo
from Bio.Phylo import NewickIO
from Bio.Nexus import Nexus
from dendro.parsers import beast_comment_parser, nhx_comment_parser, read_newick2nexus
from dendro.converters import BioNexusTrees_to_BioPhylo




def comment_formatter(variables: dict, float_fmt='%g', sep=':', begin='&&NHX:'):
    """Serialize a dict object into a string, designed for phylogenetic tree comments.

    Values that are tuple/list are assumed to contain numbers that will be collapsed
    using '..' separator. This is supposed to occur when representing ranges.
    """
    if not variables:
        return ''
    text = begin  # Bio.Phylo does not want '[' or ']' at the start/end.
    for var, value in variables.items():
        text += var + '='
        if isinstance(value, (list, tuple)):
            text += '{' + '..'.join(float_fmt % x for x in value) + '}'
        else:
            try:
                text += float_fmt % value
            except TypeError:
                text += str(value)
        text += sep

    return text.rstrip(sep)


NHX_comment_formatter = comment_formatter
nexus_comment_formatter = partial(comment_formatter, sep=',', begin='&')


def apply_node_formats(tree, float_fmt='%g', comment_parser=beast_comment_parser, comment_formatter=NHX_comment_formatter):
    """Rename nodes (=Clade) of a Bio.Phylo tree using the original newick labels.

    Bio.Phylo [version 1.73] parses node numbers as Support values.
    Use this function when it is not intended by the format (e.g. Beast/Treeannotator output tree)
    """
    for node in tree.get_terminals() + tree.get_nonterminals():
        if not node.is_terminal(): # For lsd2/Treeannotator output
            try:
                node.name = '%d' % node.confidence
            except TypeError as err:
                if node.confidence:
                    # confidence is not a integer, put it into the name
                    node.name = node.confidence
            node.confidence = None
        node.comment = comment_formatter(comment_parser(node.comment), float_fmt)


def apply_nexus_node_formats(tree, float_fmt='%g', comment_parser=beast_comment_parser, comment_formatter=nexus_comment_formatter):
    """Rename nodes of a Bio.Nexus tree using the original newick labels

    Bio.Phylo [version 1.73] parses node numbers as Support values.
    Use this function when it is not intended by the format (e.g. Beast/Treeannotator output tree)
    """
    for nodeid, node in tree.chain.items():
        data = node.get_data()
        #if tree.is_internal(nodeid):
            #if data.support is not None:
        data.comment = comment_formatter(comment_parser(data.comment), float_fmt)


def nexus2nhx(infile, outfile, float_fmt='%g'):
    """Designed to work with Beast/treeannotator output"""
    try:
        if outfile.writable():
            out = outfile
        else:
            raise ValueError('File %r not writable' % outfile)
    except AttributeError:
        out = open(outfile, 'w')
    try:
        for tree in Phylo.parse(infile, 'nexus'):
            apply_node_formats(tree, float_fmt, beast_comment_parser, NHX_comment_formatter)
            Phylo.write(tree, out, 'newick', format_branch_length=float_fmt)
    finally:
        if isinstance(outfile, (str, bytes)):
            out.close()

## TODO: tests:
# 1. Check if the resulting newick can be parsed by Ete3/Biopython/Newick-utils.

def index_nexus_treenodes(nx):
    # See:
    #   nxtree.all_ids()
    #   nxtree.get_taxa()
    #   nxtree.id
    translate = {}
    taxa = {}
    next_id = 1
    for ntree in nx.trees:
        #FIXME: index the leaves first.
        for _, node in ntree.chain.items():
            nodename = node.data.taxon
            if nodename and nodename not in taxa:
                taxa[nodename] = next_id
                translate[next_id] = nodename
                next_id += 1
    nx.translate = translate  # Overwrite current translate


# I use this block only in the bayestraits function
NEX_TAXA_BLOCK = """Begin Taxa;
 Dimensions NTax=%(count)d;
 TaxLabels
\t%(labels)s;
End;
"""


def write_nexus_trees(nx, handle, float_fmt='%g', **kwargs):
    """
    Write Nexus trees to file in Nexus format, using the NewickIO writer as a backend.
    Translate and Numerical Labels are used if defined.
    
    Note: wrote this function because Biopython did not have any way of
    writing a list of trees into a Nexus file, including the translate block.
    See:
    - Bio.Phylo.NexusIO.write: requires a list of Phylo.Trees.Tree as input, lacks translate block.
    - Nexus.Nexus.write_nexus_data: does not output trees.
    """
    handle.write("#NEXUS\n")
    if nx.taxlabels:
        handle.write(NEX_TAXA_BLOCK % {
            "count": len(nx.taxlabels),
            "labels": "\n\t".join(nx.taxlabels)})

    trees = BioNexusTrees_to_BioPhylo(nx.trees, nx.translate)
    writer = NewickIO.Writer(trees)
    
    # Build the Trees block
    nex_lines = ["Begin Trees;"]

    # Optionally add a Translate command
    if nx.translate:
        translate = ["%4d %s" % id_name for id_name in sorted(nx.translate.items(), key=lambda item: item[0])]
        nex_lines += ["\tTranslate\n\t" + ",\n\t".join(translate) + ";"]

    # Add the trees
    nex_lines += [
            "\n".join(writer.to_strings(plain_newick=False, format_branch_length=float_fmt, **kwargs)),
            "End;\n"]

    handle.write("\n".join(nex_lines))


def write_nexus_trees_to_bayestraits(nx, handle, float_fmt='%g', **kwargs):
    """Modified from Bio.Phylo.NexusIO.write():

    Compared to write_nexus_tree, it adds a Taxa block.
    """
    taxlabels = set()
    for nt in nx.trees:
        taxlabels.update(nt.get_taxa())
    nx.taxlabels = sorted(taxlabels)

    write_nexus_trees(nx, handle, float_fmt, **kwargs)
    # confidence actually contains the node label, so should not be formatted as float.
    return len(nx.trees)


def index_nexus_trees(nx, namefmt='t%d', start=1, step=1):
    idx = start
    for tree in nx.trees:
        tree.name = namefmt % idx
        idx += step

def nhx2bayestraits(infile, outfile, float_fmt='%g', **nwk_kwargs):
    """
    Phylo.convert(infile, 'newick', outfile, 'nexus') does not reindex labels,
    as required by some programs (BayesTraits).
    """
    nx = read_newick2nexus(infile)
    for tree in nx.trees:
        apply_nexus_node_formats(tree, comment_parser=nhx_comment_parser, comment_formatter=nexus_comment_formatter)
    index_nexus_trees(nx)
    index_nexus_treenodes(nx)
    try:
        if outfile.writable():
            out = outfile
        else:
            raise ValueError('File %r not writable' % outfile)
    except AttributeError:
        out = open(outfile, 'w')
    try:
        write_nexus_trees_to_bayestraits(nx, out, float_fmt, **nwk_kwargs)
    finally:
        if isinstance(outfile, (str, bytes)):
            out.close()


def nexus_rewrite(infile, outfile, float_fmt='%g'):
    nx = Nexus.Nexus(infile)
    if not nx.taxlabels and nx.translate:
        nx.taxlabels = [tax for tax_id, tax in sorted(nx.translate.items(), key=lambda item: item[0])]
    for tree in nx.trees:
        apply_nexus_node_formats(tree, float_fmt, beast_comment_parser, nexus_comment_formatter)
    try:
        if outfile.writable():
            out = outfile
        else:
            raise ValueError('File %r not writable' % outfile)
    except AttributeError:
        out = open(outfile, 'w')
    try:
        write_nexus_trees(nx, out, float_fmt)
    finally:
        if isinstance(outfile, (str, bytes)):
            out.close()

    #nx.write_nexus_data(outfile, float_fmt)

# Checks:
# 1. load the resulting nexus
# nx2 = Nexus.Nexus(nexfile)
# 2. plot the converted tree:
# plottree(nx2.trees[0], bionexus_methods.get_items, bionexus_methods.get_label, nx2.trees[0].node(nx2.trees[0].root))

def main():
    import argparse as ap
    parser = ap.ArgumentParser(__doc__)
    
    common_parser = ap.ArgumentParser(add_help=False)
    common_parser.add_argument('infile', nargs='?', default=stdin, help='[stdin]')
    common_parser.add_argument('-o', '--outfile', default=stdout, help='[stdout]')
    common_parser.add_argument('-F', '--float-fmt', default='%g', help='float format (branch lengths, comments)')

    subp = parser.add_subparsers(dest='command')
    subp.add_parser('nexus2nhx', help='Tested on Beast output', parents=[common_parser])
    subp.add_parser('nhx2bayestraits', parents=[common_parser])
    subp.add_parser('nexus', parents=[common_parser])

    args = parser.parse_args()

    if args.command == 'nexus2nhx':
        nexus2nhx(args.infile, args.outfile, args.float_fmt)
    elif args.command == 'nhx2bayestraits':
        nhx2bayestraits(args.infile, args.outfile, args.float_fmt)
    elif args.command == 'nexus':
        nexus_rewrite(args.infile, args.outfile, args.float_fmt)
    else:
        raise ValueError('Invalid command %r' % args.command)


if __name__ == '__main__':
    main()
