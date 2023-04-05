#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Set new branch length from a given table."""


from sys import stdin
import argparse as ap
from copy import deepcopy
from math import nan

from dendro.parsers import chooseparser, eval_optiontext, expand_short_options, LONG_KWARGS_ETE3
from dendro.any import methodchoice
from dendro.bates import dfw_descendants_generalized

import logging
logger = logging.getLogger(__name__)


def set_branch_lengths(tree, lengths, treetype='ete3', multiply=False, add=False, unfound='nan'):
    methods = methodchoice[treetype]
    logger.debug('methods=%s (%r)', methods, treetype)

    get_children = methods.get_children
    set_items = methods.set_items
    get_label = methods.get_label
    get_dist = methods.get_dist

    root = methods.get_root(tree)
    #newtree = deepcopy(tree)

    for node, children in dfw_descendants_generalized(tree, get_children, queue=[root]):
        nodename = get_label(tree, node)
        newitems = []
        for ch in children:
            childname = get_label(tree, ch)
            if (nodename, childname) in lengths:
                new_length = lengths[(nodename, childname)]
            elif childname in lengths:
                new_length = lengths[childname]
            else:
                msg = 'Unfound branch: %r -> %r' % (nodename, childname)
                if unfound == 'raise':
                    raise ValueError(msg)
                logger.warning(msg)
                if unfound == 'nan':
                    new_length = nan
                elif unfound == 'keep':
                    new_length = 0 if add else 1 if multiply else get_dist(tree, ch)
                else:
                    raise ValueError("'unfound' argument must be 'raise|warn|ignore'")

            if multiply:
                new_length *= get_dist(tree, ch)
            elif add:
                new_length += get_dist(tree, ch)
            #logger.debug('set dist %r %s', ch, new_length)
            #set_dist(tree, ch, new_length)
            newitems.append((ch, new_length))
        set_items(tree, (node, nan), newitems)

    #logger.debug('Root items=%s; new items=%s', newtree.items[root], newitems)
    return tree


def main():
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('lengths_table',
                        help=('3 columns tab-separated file. The 2 first are '
                        'for the parent and child node labels of the branch. '
                        'if only 2 columns: node label, new value (do not match'
                        "parent). '-' for stdin."))
    parser.add_argument('treefiles', nargs='+')
    parser.add_argument('-p', '--parser', default='ete3:1', help='[%(default)s]')
    parser.add_argument('-u', '--unfound', default='nan', choices=['raise', 'nan', 'keep'],
                        help='How to update tree edges not found in table [%(default)s]')
    modif_parser = parser.add_mutually_exclusive_group()
    modif_parser.add_argument('-m', '--multiply', action='store_true',
                        help='multiply instead of replacing (useful for rates)')
    modif_parser.add_argument('-a', '--add', action='store_true',
                        help='add instead of replacing')

    args = parser.parse_args()

    f = stdin if args.lengths_table == '-' else open(args.lengths_table)
    try:
        lengths = {}
        n_col = None
        for line in f:
            fields = line.rstrip().split('\t')
            if n_col is None:
                n_col = len(fields)
                logger.debug("table format: %d columns", n_col)
            elif n_col != len(fields):
                raise ValueError('Inconsistent number of columns: %d != %d ' % (len(fields), n_col))
            if n_col == 3:
                lengths[(fields[0], fields[1])] = float(fields[2])
            else:
                lengths[fields[0]] = float(fields[1])
    finally:
        if args.lengths_table == '-':
            f.close()

    treetype, _, tree_optiontext = args.parser.partition(':')
    treetype = treetype.strip().lower()
    parse_tree = chooseparser(args.parser, silent=True)

    print_newick = methodchoice[treetype].print_newick
    write_kwargs = eval_optiontext(tree_optiontext)
    expand_short_options(write_kwargs, LONG_KWARGS_ETE3)

    for treefile in args.treefiles:
        for tree in parse_tree(treefile):
            print_newick(set_branch_lengths(tree, lengths, treetype,
                                            args.multiply, args.add, args.unfound),
                        **write_kwargs)


if __name__ == '__main__':
    logging.basicConfig(format='%(levelname)s:%(funcName)-20s:%(message)s')
    #logger.setLevel(logging.DEBUG)
    main()
