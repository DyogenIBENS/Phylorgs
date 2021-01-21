#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Set new branch length from a given table."""


from sys import stdin
import argparse as ap
from copy import deepcopy
from math import nan

from dendro.parsers import parserchoice
from dendro.any import methodchoice
from dendro.bates import dfw_descendants_generalized

import logging
logger = logging.getLogger(__name__)


def set_branch_lengths(tree, lengths, treetype='ete3', multiply=False):
    methods = methodchoice[treetype]
    logger.debug('methods=%s (%r)', methods, treetype)

    get_children = methods.get_children
    set_items = methods.set_items
    get_label = methods.get_label
    get_dist = methods.get_dist

    root = methods.get_root(tree)
    #newtree = deepcopy(tree)

    for node, children in dfw_descendants_generalized(tree, get_children, [(None, root)]):
        nodename = get_label(tree, node)
        newitems = []
        for ch in children:
            childname = get_label(tree, ch)
            try:
                new_length = lengths[(nodename, childname)]
            except KeyError:
                logger.warning('Unfound branch: %s -> %s', nodename, childname)
                new_length = nan
            if multiply:
                new_length *= get_dist(tree, ch)
            logger.debug('set dist %s %s', ch, new_length)
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
                        "'-' for stdin."))
    parser.add_argument('treefiles', nargs='+')
    parser.add_argument('-p', '--parser', default='ete3_f1',
                        help='[%(default)s]')
    parser.add_argument('-m', '--multiply', action='store_true',
                        help='multiply instead of replacing (useful for rates)')

    args = parser.parse_args()

    f = stdin if args.lengths_table == '-' else open(args.lengths_table)
    try:
        lengths = {}
        for line in f:
            fields = line.rstrip().split('\t')
            lengths[(fields[0], fields[1])] = float(fields[2])
    finally:
        if args.lengths_table == '-':
            f.close()

    print_newick = methodchoice[args.parser].print_newick

    for treefile in args.treefiles:
        for tree in parserchoice[args.parser](treefile):
            print_newick(set_branch_lengths(tree, lengths, treetype=args.parser,
                                            multiply=args.multiply))


if __name__ == '__main__':
    logging.basicConfig(format='%(levelname)s:%(funcName)-20s:%(message)s')
    #logger.setLevel(logging.DEBUG)
    main()


