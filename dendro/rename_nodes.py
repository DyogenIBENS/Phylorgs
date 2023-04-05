#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from sys import stdin
import argparse as ap
import ete3
from dendro.parsers import chooseparser, eval_optiontext, expand_short_options, LONG_KWARGS_ETE3
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


def main(conversionfile, treefile=stdin, parser='ete3:1'):
    conversion = load_conversion(conversionfile)

    parse = chooseparser(parser)
    treetype, _, tree_optiontext = parser.partition(':')
    write_kwargs = {'format': 1, 'format_root_node': True,
                    **eval_optiontext(tree_optiontext)}
    expand_short_options(write_kwargs, LONG_KWARGS_ETE3)
    convert = converterchoice[treetype.strip().lower()]['ete3']

    for tree_object in parse(treefile):
        tree = convert(tree_object)
        rename(tree, conversion)
        print(tree.write(**write_kwargs))


if __name__ == '__main__':
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('conversionfile')
    parser.add_argument('treefile', nargs='?')
    parser.add_argument('-p', '--parser', default='ete3:1',
                        help=('How to load the tree (ete3/phyltree/prottree '
                              '+ options). `phyltree` allows to automatically '
                              'set unique node names [%(default)s].'))

    args = parser.parse_args()
    main(**vars(args))

