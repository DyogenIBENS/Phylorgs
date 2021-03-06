#!/usr/bin/env python3

"""USAGE:
    ./reshape_tree.py [-h] <input_tree> <output_tree>

OPTIONS:
    -h: print this help and exit

DESCRIPTION:
    produces reshaped newick tree compatible with codeml:
    Reshaping means removing the inner node labels, and removing nodes that
    have only one child (codeml fails on such a tree).
    Keep foreground branch labels such as '#1', '#2' or '$1'.
    Requires the python module 'ete3'.
NOTE:
    It is not an equivalent of reshape_tree.sh, since the latter only removes
    nodes with a single child if they are directly connected to the root.
"""

import sys
import re
import ete3


if __name__=='__main__':
    if len(sys.argv) != 3:
        print(__doc__)
        if len(sys.argv) == 2 and sys.argv[1] in ('-h', '--help'):
            sys.exit(0)
        sys.exit(1)

    input_tree = sys.argv[1]
    output_tree = sys.argv[2]

    re_mark = re.compile(r'[#$][0-9]')

    tree = ete3.Tree(input_tree, format=1)
    for node in tree.traverse():
        newroot = node
        if len(node.children) > 1:
            break
    for node in newroot.traverse():
        if len(node.children) == 1:
            node.delete(prevent_nondicotomic=False, preserve_branch_length=True)
        elif not node.is_leaf():
            m = re_mark.search(node.name)
            if m:
                node.name = m.group()
            else:
                node.name = ''
            

    # The format compatible with codeml is format=5
    #newroot.write(outfile=output_tree, format=5)
    newroot.write(outfile=output_tree, format=1)

