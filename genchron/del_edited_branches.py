#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Read the tree forest, delete the longest branches stemming from edited nodes,
and output the new forest."""


from sys import stdin, stdout
import argparse
import LibsDyogen.myProteinTree as ProteinTree

from dendro.bates import iter_leaves


#def prottree_getchildren(tree, node):
#    return [child for child, _ in tree.data.get(node, [])]

INFINITE_DIST = 10000
EDITED_NODE_ID = 100000000


def select_farthest_child(nodedata):  #keep=1
    """Return the child node that has the shortest distance""".
    children, dists = zip(*nodedata)
    sorted_children = sorted(enumerate(dists), key=lambda x: x[1])
    else:
        selected = 0
    return children[selected], selected


def select_farthest_leaves(nodedata, tree):
    pass


def select_less_leaves(nodedata, tree):  #keep=1
    """Return the child node whose descendant species are the least numerous."""
    # Select based on the number of species in the leaves
    species_counts = []
    for ch,_ in nodedata:
        child_leaves = list(iter_leaves(tree,
                                lambda tr,n: [c for c,_ in tr.data.get(n, [])],
                                queue=[ch]))
        child_species = set(tree.info[leaf]['taxon_name'] for leaf in child_leaves)
        species_counts.append(len(child_species))

    selected, least_species = min(enumerate(species_counts), key=lambda x: x[1])
    return nodedata[selected][0], selected


def lock_targets(node, tree, edited_node_id=EDITED_NODE_ID, infinite_dist=INFINITE_DIST,
                 fromset=None):
    """Return which children to delete"""
    nodeinfo = tree.info[node]
    nodedata = tree.data[node]

    # Too long branches
    try:
        targets = {child: i for i,(child,dist) in enumerate(nodedata) if dist >= infinite_dist}
    except ValueError as err:
        err.args += (nodedata,)
        raise

    # Gene splits
    targets.update({child: i for i,(child,_) in enumerate(nodedata) \
                    if tree.info[child]['Duplication'] == 10})

    # Gene edition
    if nodeinfo['Duplication'] == 3 or \
            (nodeinfo['Duplication'] > 0 and node >= edited_node_id):
        edited_child, edited = select_farthest(nodedata)
        targets[edited_child] = edited


    if fromset is not None:
        for i, (child, _) in enumerate(nodedata):
            if child in fromset:
                fromset.delete(child)
                targets[child] = i

    return targets


def lock_targets_filter(node, tree, filter):
    """Return which children to delete"""
    nodeinfo = tree.info[node]
    nodedata = tree.data[node]

    children_info = {}
    for child, dist in nodedata:
        children_info[child] = tree.info[child]
        children_info[child]['dist'] = dist

    # Gene splits
    targets = {child: i for i,(child,childinfo) in
                    enumerate(children_info.items())
                    if eval(filter.format(**childinfo))}

    # Gene edition
    return targets


def knock_targets(targets, tree, nodedata, nodeinfo):
    """Detach target nodes (i.e delete subtrees) from tree"""
    leaf_count = 0
    for target, target_index in sorted(tuple(targets.items()),
                                       key=lambda item: -item[1]):
        tree.info.pop(target)
        try:
            tree.data.pop(target)
        except KeyError as err:
            leaf_count += 1

        nodedata.pop(target_index)
    return leaf_count


def filterbranches(tree, node, edited_node_id=EDITED_NODE_ID,
                   infinite_dist=INFINITE_DIST, filter=None, fromset=None):
    """Descend the tree from root to branches. Discard any anomalous node (and its descendants)"""
    del_count = 0
    del_leaf_count = 0

    #if node is None:
    #    node = tree.root

    data = tree.data.get(node)

    if data:
        info = tree.info[node]
        targets = lock_targets(node, tree, edited_node_id, infinite_dist)
        if targets:
            del_count      += len(targets)
            del_leaf_count += knock_targets(targets, tree, data, info)
        
        for child,_ in data:
            ch_dc, ch_dlc = filterbranches(tree, child, edited_node_id, infinite_dist)
            del_count      += ch_dc
            del_leaf_count += ch_dlc

    return del_count, del_leaf_count


def main(treeforestfile, outfile, dryrun=False, edited_node_id=EDITED_NODE_ID,
         infinite_dist=INFINITE_DIST):
    total_deleted = 0
    total_leaves_deleted = 0
    if dryrun and outfile is not stdout: outfile.close()
    if treeforestfile == '-': treeforestfile = stdin
    for tree in ProteinTree.loadTree(treeforestfile):
        del_count, del_leaf_count = filterbranches(tree, tree.root,
                                                   edited_node_id,
                                                   infinite_dist)
        total_deleted        += del_count
        total_leaves_deleted += del_leaf_count
        if not dryrun:
            tree.printTree(outfile)
        #break
    print("Deleted %d branches, of which %d leaves." % (total_deleted, total_leaves_deleted))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('treeforestfile', nargs='?', default='-')
    parser.add_argument('outfile', nargs='?', default=stdout,
                        type=argparse.FileType('w'))
    parser.add_argument('-n', '--dryrun', action='store_true',
                        help='Do not delete or output anything, just print '\
                             'counts.')
    parser.add_argument('-e', '--edited_node_id', type=int, default=EDITED_NODE_ID,
                        help="Use 100000000 [default] or 500000000000000 (5e14).")
    parser.add_argument('-i', '--infinite-dist', type=float, default=INFINITE_DIST,
                        help="branch length considered aberrant [%(default)s].")
    parser.add_argument('-f', '--filter',
                        help='Logical expression to select nodes.')
    
    args = parser.parse_args()
    main(**vars(args))
