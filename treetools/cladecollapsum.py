#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Tool for phylogenetic comparative testing: Clades-collapse-summarize.

1 - Compute statistics in a phylogeny for each clade of a given rank
    (e.g, each family).

    Statistics are:
        - clade age
        - number of extant species
        - diversification rate
        - feature rate
    
    Return a table with values corresponding to each clade.

2 - Return the tree with collapsed clades.
"""

import sys
import os.path
import argparse
import ete3 # Will use PhyloTree, NCBITaxa
import numpy as np

from itertools import chain

RANKS = ['kingdom', 'subkingdom', 'phylum', 'subphylum', 'superclass', 'class', 'subclass', 'infraclass', 'superorder', 'order', 'suborder', 'infraorder', 'parvorder', 'superfamily', 'family', 'subfamily', 'genus', 'subgenus', 'species', 'subspecies']

def match_duplicate_taxid(taxids, node, taxid2name, ncbi):
    """Find the correct taxid corresponding to node name, by comparing the node
    ancestors with the taxid lineage."""
    lineages = ncbi.get_lineage_translator(taxids)
    named_lin = {t: set(taxid2name.get(lt) for lt in lin) for t,lin \
                                        in lineages.items()}
    ancestors = [n.name for n in node.get_ancestors()]
    
    descendants = [n.name for n in node.iter_leaves()]
    dname2taxids = ncbi.get_name_translator(descendants)
    _, dtaxids = dname2taxids.popitem()
    while len(dtaxids) > 1:
        _, dtaxids = dname2taxids.popitem()
    #dlineages = ncbi.get_lineage_translator(chain(*dtaxids.values()))
    dlineage0 = ncbi.get_lineage(dtaxids[0])

    candidates = sorted(taxids,
                    key=lambda t: dlineage0.index(t) if t in dlineage0 else -1)
    match = candidates[-1]

    # Find the correct taxid: the one whose lineage (ncbi) matches
    # first to the most recent ancestor (in tree).

    #match = None
    #while not match:
    #    try:
    #        anc = ancestors.pop()
    #    except IndexError:
    #        print(taxids, taxid2name, node.name, [n.name for n in node.get_ancestors()],
    #              named_lin, lineages, sep='\n', file=sys.stderr)
    #        raise
    #    for taxid, nl in named_lin.items():
    #        if anc in nl:
    #            if not match:
    #                match = taxid
    #            else:
    #                raise RuntimeError("Can't decide because previous "\
    #                                   "ancestor is in multiple lineages.")
    #            #break

    print('Choose %s' % match, file=sys.stderr)
    return match


def get_valid_tax_children(node, name2taxid):
    """return the list of node descendants such that all of them are in the
    taxonomy (name2taxid)."""
    nodelist = []
    for child in node.children:
        if child.name.replace('_', ' ') in name2taxid:
            nodelist.append(child)
        else:
            nodelist.extend(get_valid_tax_children(child, name2taxid))
    return nodelist


def def_group_feature_rate(stem_or_crown="crown"):
    """Return the computation function for either crown or stem group."""
    if stem_or_crown == "stem":
        iter_group = lambda node: node.traverse()
    elif stem_or_crown == "crown":
        iter_group = lambda node: node.iter_descendants()
    else:
        raise ValueError("'stem_or_crown' must be either 'stem' or 'crown'")

    def group_feature_rate(node, features):
        tot_dist = 0
        tot_ft = np.zeros(len(features))

        for descendant in iter_group(node):
            #if node.name == 'Carnivora': print('-', descendant.name)
            if descendant.dist > 0:
                tot_dist += descendant.dist
                try:
                    tot_ft += np.array([float(getattr(descendant, ft)) for ft in features])
                    #if node.name == 'Carnivora':
                    #    print(features)
                    #    print(np.array([float(getattr(descendant, ft)) for ft in features]))
                    #    print(tot_ft)
                    #    print('tot_dist:', tot_dist)
                except TypeError:
                    print([getattr(descendant, ft) for ft in features], file=sys.stderr)
                    raise
                except AttributeError:
                    print(node.name, 'root:', node.is_root(),
                          node.get_farthest_leaf(),
                          file=sys.stderr)
                    raise

        return tot_ft / tot_dist

    group_feature_rate.__name__ = stem_or_crown + '_group_feature_rate'
    return group_feature_rate


def main(inputtree, outbase, rank='family', div=True, features=None,
         stem_or_crown="crown", byage=None, bylist=None, bysize=None):
    """byage: collapse any node of age <= byage
       bylist: read list of nodes from file
       bysize: collapse oldest nodes with size < bysize"""
    group_feature_rate = def_group_feature_rate(stem_or_crown)

    tree = ete3.PhyloTree(inputtree, format=1, quoted_node_names=False)

    if bylist:
        rank = None
        with open(bylist) as s:
            nodelist = set(l.rstrip().split('\t')[0] for l in s)
        def is_leaf_fn(node):
            return node.name in nodelist
    elif byage:
        rank = None
        def is_leaf_fn(node):
            _, age = node.get_farthest_leaf()
            return age <= byage
    elif bysize:
        rank = None
        def is_leaf_fn(node):
            return len(node) <= bysize
    
    if rank or div:
        print("Loading taxonomy", file=sys.stderr)
        ncbi = ete3.NCBITaxa()

        name2taxid = ncbi.get_name_translator(
                            [node.name.replace('_', ' ') if node.is_leaf() \
                                else node.name for node in tree.traverse()])
        # Won't return anything for names not found

        if rank:
            taxid2rank = ncbi.get_rank(chain(*name2taxid.values()))
            taxid2name = ncbi.get_taxid_translator(chain(*name2taxid.values()))
            
            # Could also simply annotate the tree using ncbi.annotate_tree
            #ncbi.annotate_tree(tree, 'taxid')

            #def is_leaf_fn(node):
            #    return node.rank == rank

            def is_leaf_fn(node):
                """Stop at given rank"""
                taxids = name2taxid.get(node.name)
                if taxids is not None:
                    if len(taxids) > 1:
                        print('WARNING: non unique name %r: %s.' % (node.name, taxids),
                              end=' ', file=sys.stderr)
                        taxids = [match_duplicate_taxid(taxids, node, taxid2name, ncbi)]

                    # Non inclusive method (must be *exactly* rank, not above):
                    #noderank = taxid2rank[taxids[0]]
                    #return noderank == rank
                    # Must be included in the rank:
                    lineage = ncbi.get_lineage(taxids[0])
                    ranks = set(ncbi.get_rank(lineage).values())
                    #print('- %s:' % node.name, ' '.join(ranks), file=sys.stderr)
                    return rank in ranks
                else:
                    return False
    
    outsuffix = 'stem' if stem_or_crown == 'stem' else ''

    if byage:
        outsuffix += 'age%g' % byage
    elif bylist:
        outsuffix += 'list' + os.path.splitext(os.path.basename(bylist))[0]
    elif bysize:
        outsuffix += 'size%d' % bysize
    else:
        outsuffix += rank

    columns = [outsuffix, 'size', 'branches', 'age'] #'crown_age', 'stem_age']
    if div: columns.extend(('div_rate', 'ncbi_sp_sampling'))
    if features: columns.extend(features)

    with open(outbase + '-%s.tsv' % outsuffix, 'w') as outtsv, \
         open(outbase + '-%s.subtrees.nwk' % outsuffix, 'w') as outsub:

        outtsv.write('\t'.join(columns) + '\n')
        
        print("Iterating over found clades", file=sys.stderr)
        for node in tree.iter_leaves(is_leaf_fn):
            outsub.write(node.write(features, format=1, format_root_node=True) + '\n')
            
            # Collapse
            size = len(node)
            branches = len(node.get_descendants())
            _, age = node.get_farthest_leaf()
            if stem_or_crown == 'stem':
                age += node.dist
            values = [node.name, size, branches, age]
            if div:
                div_rate = float(size) / age if age else np.NaN
                try:
                    nodetaxids = name2taxid[node.name.replace('_', ' ')]
                    if len(nodetaxids) > 1:
                        nodetaxids = [match_duplicate_taxid(taxids, node,
                                                            taxid2name, ncbi)]

                except KeyError:
                    # This clade isn't in the taxonomy (example: Atlantogenata)
                    # take descendant nodes and join them
                    valid_tax_children = get_valid_tax_children(node, name2taxid)
                    vtc_names = [vtc.name.replace(' ', '_') for vtc in valid_tax_children]

                    print('WARNING: %r not found in NCBI Taxonomy. Merging the'
                          ' node children %s to get the descendant counts.' \
                          % (node.name, vtc_names), file=sys.stderr)
                    
                    nodetaxids = []
                    for vtc_n, vtc in zip(vtc_names, valid_tax_children):
                        vtc_taxids = name2taxid[vtc_n]
                        if len(vtc_taxids) == 1:
                            nodetaxids.append(vtc_taxids[0])
                        else:
                            nodetaxids.append(match_duplicate_taxid(vtc_taxids,
                                                      vtc, taxid2name, ncbi))

                try:
                    ncbi_sp = list(chain(*(ncbi.get_descendant_taxa(nt,
                                                    rank_limit='species') \
                                           for nt in nodetaxids)))
                                                    #collapse_subspecies=True))
                except:
                    print(node.name, nodetaxids, file=sys.stderr)
                    raise
                sp_sampling = float(size) / len(ncbi_sp)
                values.extend((div_rate, sp_sampling))

            if features:
                ft_rates = group_feature_rate(node, features)
                values += ft_rates.tolist()

            outtsv.write('\t'.join(str(v) for v in values) + '\n')

    tree.write(outfile=(outbase + '-%s.nwk' % outsuffix), format=1,
               is_leaf_fn=is_leaf_fn, format_root_node=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('inputtree')
    parser.add_argument('outbase',
                        help='basename for the 2 output files (.tsv, .nwk)')
    parser.add_argument('--nodiv', dest='div', action='store_false',
                        help='Do not compute diversification values')
    parser.add_argument('-S', '--stem', dest='stem_or_crown', default='crown',
                        action='store_const', const='stem',
                        help='include stem branch')
    parser.add_argument('-f', '--features', nargs='+',
                        help='compute the average rate of these features')
    method_parser = parser.add_mutually_exclusive_group()
    method_parser.add_argument('-r', '--rank', default='family', choices=RANKS,
                        help='taxonomic rank to collapse [%(default)s]')
    method_parser.add_argument('-a', '--byage', type=float,
                        help='collapse oldest clades with age <= byage')
    method_parser.add_argument('-l', '--bylist',
                        help='collapse by a list of node names (file)')
    method_parser.add_argument('-s', '--bysize', type=int,
                        help='collapse oldest clades with size <= bysize')
    
    
    args = parser.parse_args()
    main(**vars(args))
