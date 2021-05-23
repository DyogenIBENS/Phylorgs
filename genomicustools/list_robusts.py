#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse as ap
from LibsDyogen import myProteinTree, myPhylTree
from dendro.reconciled import prottree_list_robusts

import logging


def main():
    logging.basicConfig(format=logging.BASIC_FORMAT)#, level=logging.INFO)
    p = ap.ArgumentParser(description=__doc__)
    p.add_argument('prottreefile')
    p.add_argument('phyltreefile')
    p.add_argument('-a', '--ancestor')
    args = p.parse_args()

    trees = myProteinTree.loadTree(args.prottreefile)
    phyltree = myPhylTree.PhylogeneticTree(args.phyltreefile)
    for taxon, nodename in prottree_list_robusts(trees, phyltree, args.ancestor):
        print(taxon, nodename, sep='\t')

if __name__ == '__main__':
    main()
