#!/usr/bin/env python3


from sys import stdin, stderr
import argparse as ap
from ete3 import Tree
from ete3.parser.newick import NewickError
from itertools import cycle


def iter_multinewick(filename):
    f = stdin if filename == '-' else open(filename)
    try:
        for line in f:
            if not line.startswith('#'):
                yield line.rstrip()
    finally:
        if filename != '-':
            f.close()


def zip_if_aligned_else_onetomany(list1, iterable2):
    """If both iterators have many elements, zip them.
        Otherwise, perform all possible pairwise comparisons"""
    iter1 = cycle(list1) if len(list1) == 1 else list1
    iter2 = iter(iterable2)  # We make it consumable.

    iter_pairs = zip(iter1, iter2)

    for i, (elem1, elem2) in enumerate(iter_pairs):
        yield elem1, elem2

    if i==0 and len(list1)>1:
        # There was a single element in iterable2, we should make all comparison
        # between this element and each element of list1 (recycle elem2)
        for elem1 in list1[1:]:
            yield elem1, elem2

    elif len(list1) != 1 and i != len(list1) - 1 :
        print('Warning: Extra elements in list 1. Ignored', file=stderr)
    else:
        try:
            next(iter2)
        except StopIteration:
            return
        print('Warning: Extra elements in list 2. Ignored', file=stderr)



def main():
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('reftrees', help="'-' for stdin.")
    parser.add_argument('targettrees', help="'-' for stdin.")
    parser.add_argument('-t', '--format', type=int, default=0, help='Ete3 Newick flavor [%(default)s]')
    parser.add_argument('-r', '--rooted', action='store_true', help='Switch to rooted tree comparison [unrooted]')
    parser.add_argument('--correct-poly-size', action='store_true')
    args = parser.parse_args()
    unrooted = not args.rooted

    def load_reftree(newick):
        reftree = Tree(newick, format=args.format, quoted_node_names=True)
        reftree.add_feature('newick', newick[:100])
        if len(newick) > 100:
            reftree.newick += '[...]'
        if args.correct_poly_size and unrooted and len(reftree.children) > 2:
            # Probable bug in Ete3 which considers that the root must be corrected.
            # It should not as it's conventional to have the polytomic root for unrooted trees.
            reftree.set_outgroup(reftree.children[0])
        return reftree

    ref_trees = [load_reftree(newick) for newick in iter_multinewick(args.reftrees)]

    iter_targets = iter_multinewick(args.targettrees)

    count = 0
    for ref_tree, newick in zip_if_aligned_else_onetomany(ref_trees, iter_targets):
        count += 1
        try:
            tree = Tree(newick, format=args.format, quoted_node_names=True)
            if args.correct_poly_size and unrooted and len(tree.children) > 2:
                tree.set_outgroup(tree.children[0])

            RF = ref_tree.robinson_foulds(tree, unrooted_trees=unrooted, correct_by_polytomy_size=args.correct_poly_size)

            if RF[1] == 0:
                if not set(ref_tree.iter_leaves()).intersection(tree.get_leaves()):
                    raise ValueError("Uncomparable trees (no common leaves)")
                raise RuntimeError("Unexpected: RFmax = 0 (RF=%s)" % RF[0])

            # Output: normalized, absolute, max)
            print('%g\t%d\t%d' % (float(RF[0])/RF[1], RF[0], RF[1]))

        except BaseException as err:
            info = ". Comparison %d ref %r VS newick %r" % (count, ref_tree.newick, newick)
            newargs = (err.args[0] + info,) + err.args[1:]
            raise type(err)(*newargs)


if __name__ == '__main__':
    main()
