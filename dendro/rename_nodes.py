#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from sys import stdin
import argparse as ap
import ete3
from dendro.parsers import parserchoice
from dendro.converters import converterchoice


def load_conversion(conversionfile):
    with open(conversionfile) as stream:
        conversion = {}
        for line in stream:
            if not line.startswith('#'):
                field1, field2, *extra = line.rstrip('\r\n').split('\t')
                if field1 and field2:
                    conversion[field1] = field2
    return conversion


def rename(tree, conversion):
    for node in tree.traverse():
        try:
            node.name = conversion[node.name]
        except KeyError:
            pass


def main(conversionfile, treefile=stdin, parser='ete3_f1'):
    conversion = load_conversion(conversionfile)

    parse = parserchoice[parser]
    convert = converterchoice[parser]['ete3']

    for tree_object in parse(treefile):
        tree = convert(tree_object)
        rename(tree, conversion)
        print(tree.write(format=1, format_root_node=True))


if __name__ == '__main__':
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('conversionfile')
    parser.add_argument('treefile', nargs='?')
    parser.add_argument('-p', '--parser', default='ete3_f1',
                        choices=['ete3_f1', 'ete3', 'phyltree', 'proteintree'],
                        help='How to load the tree. `myPhylTree` allows to ' \
                             'automatically set unique node names [%(default)s].')
    
    args = parser.parse_args()
    main(**vars(args))

