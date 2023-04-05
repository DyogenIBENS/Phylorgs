#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""""""

from dendro.parsers import chooseparser, eval_optiontext, expand_short_options, LONG_KWARGS_ETE3, PARSERSPEC_HELP
from dendro.any import methodchoice
from dendro.bates import dfw_descendants_generalized
import numpy as np

import argparse as ap


def discretize_branchlengths(tree_methods, tree, nbins=5):
    get_items = tree_methods.get_items
    set_items = tree_methods.set_items
    get_root = tree_methods.get_root
    get_dist = tree_methods.get_dist
    root = get_root(tree)
    rootdist = get_dist(tree, root)  # Unused.

    distrib = []
    iter_tree = list(dfw_descendants_generalized(tree, get_items, queue=[(root, rootdist)]))

    for (parent, dist), items in iter_tree:
        distrib.extend((d for ch,d in items))

    bins = np.histogram_bin_edges(distrib, bins=nbins)
    newdists = np.digitize(distrib, bins)
    i = 0

    for (parent, dist), items in iter_tree:
        newitems = [(ch, nd) for (ch, _), nd in zip(items, newdists[i:])]
        set_items(tree, (parent, dist), newitems)
        i+=len(newitems)

    tree_methods.set_dist(tree, root, np.digitize(rootdist, bins))
    return tree

### TODO: new_age = log(10 + age)

def multiply_branchlengths(tree_methods, tree, factor=1):
    get_items = tree_methods.get_items
    set_items = tree_methods.set_items
    get_root = tree_methods.get_root
    get_dist = tree_methods.get_dist
    root = get_root(tree)
    rootdist = get_dist(tree, root)  # Unused.

    iter_tree = list(dfw_descendants_generalized(tree, get_items, queue=[(root, rootdist)]))

    for (parent, dist), items in iter_tree:
        newitems = [(ch, d*factor) for ch, d in items]
        set_items(tree, (parent, dist), newitems)

    tree_methods.set_dist(tree, root, rootdist*factor)
    return tree


def main():
    parser = ap.ArgumentParser(description=__doc__, epilog=PARSERSPEC_HELP,
                               formatter_class=ap.RawTextHelpFormatter)
    #common_parser = ap.ArgumentParser(add_help=False)
    parser.add_argument('treefile')
    parser.add_argument('-p', '--parser', default='ete3:1', help='[%(default)s]')

    subp = parser.add_subparsers(dest='transform')

    pars_discret = subp.add_parser('discretize', aliases=['disc'],
                                   help='Convert branch lengths to a limited set of integers')
    pars_discret.add_argument('-n', '--nbins', type=int, default=5,
                              help='number of bins [%(default)s]')

    pars_multiply = subp.add_parser('multiply', aliases=['mult'], help='')
    pars_multiply.add_argument('factor', type=float, help='')

    #TODO:
    #pars_replace = subp.add_parser('replace', aliases=['repl'], help='')
    #pars_replace.add_argument('lengths_table',
    #                    help=('3 columns tab-separated file. The 2 first are '
    #                    'for the parent and child node labels of the branch. '
    #                    "'-' for stdin.")

    args = parser.parse_args()
    transform = args.transform
    delattr(args, 'transform')

    parse_tree = chooseparser(args.parser, silent=True)
    treetype, _, tree_optiontext = args.parser.partition(':')
    tree_methods = methodchoice[treetype.strip().lower()]

    write_kwargs = eval_optiontext(tree_optiontext)
    expand_short_options(write_kwargs, LONG_KWARGS_ETE3)

    if transform in ('disc', 'discretize'):
        func = discretize_branchlengths
        kwargs = dict(nbins=args.nbins)
    elif transform in ('mult', 'multiply'):
        func = multiply_branchlengths
        kwargs = dict(factor=args.factor)
    else:
        raise ValueError('Invalid transform command %r' % transform)

    for tree in parse_tree(args.treefile):
        func(tree_methods, tree, **kwargs)  # inplace
        tree_methods.print_newick(tree, **write_kwargs)


if __name__=='__main__':
    main()
