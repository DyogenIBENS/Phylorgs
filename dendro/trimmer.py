#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Cut tree branches in various ways"""
# Anc. gr. name: tom

from copy import copy, deepcopy

from dendro.bates import rev_dfw_descendants, iter_distleaves, dfw_pairs_generalized, dfw_descendants_generalized
import logging
logger = logging.getLogger(__name__)


def reset_data(tree, node, data):
    node.children = []
    for child, dist in data:
        node.add_child(child=child, dist=dist)


def get_basal(nodes, maxsize):
    """Identify `maxsize` most basal nodes from a list of sister nodes.

    Return 2 lists:
        - selected basal nodes,
        - excluded nodes (to detach)

    "Basal" means: divergence closer to the root.

    - `nodes`: a *list* of TreeNode instances;
    - `maxsize`: integer.

    Return ([],[]) if maxsize >= number of leaves.
    """

    nodedists = [(node, node.dist) for node in sorted(nodes,key=lambda n:n.dist)]

    # Not enough leaves
    if sum(len(n) for n,_ in nodedists) <= maxsize:  #minsize
        #return [], []
        return nodedists, []
        #return list(chain(n.get_leaves() for n in nodes))

    # Sort by distance from the original root:
    # An alternative sortkey could be
    #    sortkey = lambda nodedist: nodedist[0].get_closest_leaf(topology_only=False)[1]
    sortkey = lambda nodedist: nodedist[1]

    kept = []
    while len(kept) + len(nodedists) < maxsize:
        try:
            node, dist = nodedists.pop(0)  # Descend into the closest divergence.
        except IndexError:
            break

        nextnodes = node.children
        if nextnodes:
            nodedists.extend((nn, nn.dist + dist) for nn in nextnodes)
            nodedists.sort(key=sortkey)
        else:
            kept.append(node)

    descended_nodes = kept + [n for n,_ in nodedists]
    return descended_nodes[:maxsize], descended_nodes[maxsize:]


def thin(tree, node, leaves,
         get_data=lambda tree, data: [(ch, ch.dist) for ch in data.children],
         get_name=lambda tree, data: data[0].name
         ):
    """Keep only the intermediate nodes corresponding to the most recent common ancestors
    of the given leaves."""
    raise NotImplementedError

    #for leafpath in dclimb.iter_leaf_paths(tree, get_data, [(node,0)]):
    #    if

    newdata = []
    for child, chdist in get_data(tree, node, 0):
        unwanted_child = thin(tree, child, leaves, get_data, get_name)
    set_data(tree, node, data)
    return unwanted


def thin_ete3(root, keptleaves):
    """
    Thin down a tree to retain the minimal topology keeping the most recent
    ancestors of the keptleaves.

    Return the latest common ancestor of the found keptleaves if the root
    supports some keptleaves, else return None.

    This edits the structure (copying might be needed).
    """
    if root.is_leaf():
        return root if root in keptleaves else None
    for child in root.get_children():
        if thin_ete3(child, keptleaves) is None:
            root.remove_child(child)
    if not root.children:
        return None
    if len(root.children) == 1:
        newroot = root.children[0]
        newroot.dist += root.dist  # Transfer the branch length to the *child*
        root.delete(prevent_nondicotomic=False, preserve_branch_length=False)
        root = newroot

    return root


def thin_prottree(tree, root, rootdist, keptleaves):
    """
    Thin down a tree to retain the minimal topology keeping the most recent
    ancestors of the keptleaves.

    Return the latest common ancestor of the found keptleaves if the root
    supports some keptleaves, else return None.

    This edits the structure (copying might be needed).
    """
    if root not in tree.data:
        if root in keptleaves:
            return root, rootdist
        else:
            del tree.info[root]
            return None, None

    newnodedata = []
    for child,dist in tree.data[root]:
        thinned_child, thinned_dist = thin_prottree(tree, child, dist, keptleaves)
        if thinned_child is not None:
            newnodedata.append((thinned_child, thinned_dist))

    if not newnodedata:
        del tree.data[root]
        del tree.info[root]
        return None, None

    if len(newnodedata) == 1:
        del tree.data[root]
        del tree.info[root]
        root, newdist = newnodedata[0]
        rootdist += newdist  # Transfer the branch length to the *child*
    else:
        tree.data[root] = newnodedata

    return root, rootdist


def fuse_single_child_nodes_ete3(tree, copy=True):
    if copy:
        tree = tree.copy()
    count = 0
    for node in tree.iter_descendants(strategy='postorder'):
        if len(node.children) == 1:
            child, = node.children
            child.dist += node.dist
            # This does not preserve the order of the children,
            # if the parent had >=2 children.
            #node.delete(prevent_nondicotomic=False,
            #            preserve_branch_length=False)
            node_i = node.up.children.index(node)
            node.up.children[node_i] = child
            count += 1

    while len(tree.children) == 1:
        tree = tree.children[0].detach()
        count += 1

    logger.debug('Fused %d nodes', count)
    return tree


def collapse_clades(tree, get_items, set_items, root, clades, make_new_clade=None):
    """Note: this modifies tree **inplace**. Make a copy accordingly.

    make_new_clade: function to create the new node. Takes (clade, cladesize) as argument.
    By default, prefix the number of leaves to the clade name.
    Useful for PhylTree which does not support duplicated node names.

    Handles nested clades by collapsing the most basal.
    """
    if make_new_clade is None:
        def make_new_clade(clade, cladesize):
            if isinstance(clade, str):
                return '%d %s' % (cladesize, clade)
            elif hasattr(clade, 'name'):
                new = copy(clade)
                setattr(new, 'name', '%d %s' % (cladesize, clade.name))
                return new
    leaf_numbers = [0]*len(clades)
    # Iterate from root to leaves
    iter_tree = list(dfw_pairs_generalized(tree, get_items, queue=[((None,0), (root, 0))]))
    for _, (clade,_) in iter_tree:
    #iter_tree = list(dfw_descendants_generalized(tree, get_items, queue=[(tree.root, 0)], copy=True))
    #for _, items in iter_tree:
    #    for clade, _ in items:
        try:
            clade_i = clades.index(clade)
        except ValueError:
            continue
        leafdists = list(iter_distleaves(tree, get_items, root=clade))
        leaf_numbers[clade_i] = len(leafdists)
        maxdist = max(d for l,d in leafdists)
        #TODO: add the mindist information -> draw non ultrametric triangles
        set_items(tree, (clade, None), [(make_new_clade(clade, len(leafdists)), maxdist)])  # (clade, mindist)
                #TODO: some cleanup is needed in PhylTree (parents, dicLinks, ages)
    return leaf_numbers

