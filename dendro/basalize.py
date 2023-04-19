#!/usr/bin/env python3


"""
Trim down a tree to N leaves while keeping the N most basal divergences.
Keep the closest leaf under each retained basal node.
"""


from sys import stdin, stderr
import argparse as ap
from ete3 import Tree
from dendro.trimmer import get_basal


def main():
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('treefiles', nargs='*', default=['-'], help='[stdin]')
    parser.add_argument('-t', '--tree-format', type=int, default=1, help='Ete3 newick flavor [%(default)s]')
    parser.add_argument('-q', '--quoted-labels', action='store_true', help='Parse quoted node names')
    parser.add_argument('-N', type=int, default=3, help='The final number of leaves [%(default)s]')
    parser.add_argument('-v', '--reverse', action='store_true', help='Output the terminal subtrees instead of the parent tree')
    args = parser.parse_args()

    nwk_format = args.tree_format
    quoted_labels = args.quoted_labels
    reverse = args.reverse

    for treefile in args.treefiles:
        newick = stdin.read() if treefile == '-' else treefile
        tree = Tree(newick, format=nwk_format, quoted_node_names=quoted_labels)

        basal_nodes, excluded = get_basal(tree.children, args.N)

        print('retained %d basal nodes VS %d excluded' % (len(basal_nodes), len(excluded)), file=stderr)
        #print('-', '\n- '.join(n.name for n in basal_nodes), file=stderr)

        for basal in basal_nodes:
            if reverse:
                print(basal.write(format=nwk_format, quoted_node_names=quoted_labels, format_root_node=True))
            if not basal.is_leaf():
                closest, dist = basal.get_closest_leaf(topology_only=False)
                for child in basal.get_children():
                    child.detach()
                basal.add_child(closest, dist=dist)
        for node in excluded:
            node.detach()

        if len(tree) > args.N:
            raise RuntimeError('Failed to trim tree: %d leaves' % len(tree))

        if reverse:
            print('')
        else:
            print(tree.write(format=nwk_format, quoted_node_names=quoted_labels, format_root_node=True))


if __name__ == '__main__':
    main()
