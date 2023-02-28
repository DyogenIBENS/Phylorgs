#!/usr/bin/env python3


"""
Trim down a tree to N leaves while keeping the N most basal divergences.
Keep the closest leaf under each retained basal node.
"""


from sys import stdin
import argparse as ap
from ete3 import Tree
from dendro.trimmer import get_basal


def main():
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('treefiles', nargs='*', default=['-'], help='[stdin]')
    parser.add_argument('-t', '--tree-format', type=int, default=1, help='Ete3 newick flavor [%(default)s]')
    parser.add_argument('-q', '--quoted-labels', action='store_true', help='Parse quoted node names')
    parser.add_argument('-N', type=int, default=3, help='The final number of leaves [%(default)s]')
    args = parser.parse_args()

    for treefile in args.treefiles:
        newick = stdin.read() if treefile == '-' else treefile
        tree = Tree(newick, format=args.tree_format, quoted_node_names=args.quoted_labels)

        basal_nodes = get_basal(tree, args.N)

        for basal in basal_nodes:
            if not basal.is_leaf():
                closest, dist = basal.get_closest_leaf()
                for child in basal.get_children():
                    child.detach()
                basal.add_child(closest, dist=dist)

        print(tree.write(format=args.tree_format, quoted_node_names=args.quoted_labels), format_root_node=True)


if __name__ == '__main__':
    main()
