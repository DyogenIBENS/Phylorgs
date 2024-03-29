#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Tools to handle "reconciled" phylogenetic gene trees

Reconciliation: mapping of a gene tree on a species tree.

Uses `ete3.Tree` objects.

This tools recognize different ways of annotating the species at each node:
    - TreeBest: the `S` tag in NHX comments.
    - custom (Ensembl gene trees): the start of the node name, e.g
                                   HomoPanENSGT0039000...
"""


from sys import stderr
from genomicustools.identify import ultimate_seq2sp
import logging
#logging.basicConfig(format='%(levelname)s:%(module)s l.%(lineno)d:%(funcName)s:%(message)s')
logger = logging.getLogger(__name__)
#ch = logging.StreamHandler()
#ch.setFormatter(logging.Formatter("%(levelname)s:%(message)s"))
#logger.addHandler(ch)


# For types
#import re, ete3


ENSEMBL_VERSION = 85


def make_ancgene2sp(ancestor, phyltree):
    return re.compile(r'('
                      + r'|'.join(re.escape(s) for s in
                                  list(phyltree.species[ancestor]) +
                                  sorted(phyltree.getTargetsAnc(ancestor),
                                         key=lambda a:len(a),
                                         reverse=True)).replace(' ', '.')
                      + r')([^a-z].*|)$')


def split_species_gene(nodename, ancgene2sp):
    match = ancgene2sp.match(nodename)
    try:
        taxon, genename = match.groups()
        taxon = taxon.replace('.', ' ')
    except AttributeError:
        taxon = None
        genename = nodename
    return taxon, genename


#def get_taxon(node: ete3.TreeNode, ancgene2sp : re._pattern_type,
#              ensembl_version : int = ENSEMBL_VERSION):
def get_taxon(node, ancgene2sp, ensembl_version=ENSEMBL_VERSION, seq_error='raise'):
    """from a gene name in my newick gene trees, find the taxon:
        either:
            - node is a leaf (e.g ENSMUSG00...)
            - node is internal (e.g Mus.musculusENSGT...)"""
    if node.is_leaf():
        try:
            taxon = ultimate_seq2sp(node.name, ensembl_version)
        except KeyError as err:
            if seq_error == 'raise':
                raise
            if seq_error == 'warn':
                logger.warning('ultimate_seq2sp: %s: %s', type(err).__name__, err)
                taxon = node.name
            elif seq_error == 'ignore':
                taxon = node.name
            else:
                raise ValueError('seq_error value must be "raise"|"warn"|"ignore"')
    else:
        try:
            taxon = ancgene2sp.match(node.name).group(1).replace('.', ' ')
        except AttributeError:
            raise ValueError("Cannot match species name in %r" % node.name)
    return taxon


def get_taxon_treebest(node, *args):
    """Get the taxon of a gene node, using a treebest output tree.
    
    *args are not used, they are here for compatibility with `get_taxon`"""
    try:
        return node.S.replace('.', ' ')
    except AttributeError as err:
        err.args += (node.name, node)
        raise


def eval_treebest_nhx_tag(node, tag='D'):
    value = getattr(node, tag, 0)
    try:
        value = int(value)
    except ValueError:
        if value in ('N', ''):
            value = False
        elif value == 'Y':
            value = True
        elif value.lower() in ('true', 'false', 'none'):
            value = eval(value.capitalize())
        else:
            try:
                value = float(value)  # handles D=0.0 as False
            except ValueError:
                logger.warning("Unexpected value %r for duplication tag 'D' at node %r", value, node.name)
    return value


def infer_gene_event(node, taxon, children_taxa):
    """Use taxon information to tell whether a gene tree node (Ete3 format) is:
    - a leaf,
    - a speciation,
    - or a duplication.
    - transfers may be annotated as duplications.
    
    Use TreeBest annotation if present.
    
    param: `children_taxa` must be a set (because the number of *uniq* elements
           is used).
    
    This is the original version from `genetree_drawer.py` (2018/11/19).
    """
    #if node.is_leaf():
    #import ete3
    #assert isinstance(node, ete3.TreeNode)
    #assert isinstance(taxon, str)
    #assert isinstance(children_taxa, set)

    if not children_taxa:
        return 'leaf'

    children_names = [ch.name for ch in node.children]
    treebest_isdup = eval_treebest_nhx_tag(node, 'D')
    treebest_istransfer = eval_treebest_nhx_tag(node, 'T')
    treebest_isdup = bool(treebest_isdup) or bool(treebest_istransfer)

    if (treebest_isdup) or taxon in children_taxa:
       #(len(children_taxa) == 1 and taxon in children_taxa): #\
                    #and len(node.children) > 1 \
       #or \
       #(len(node.children) == 1 and taxon in children_taxa):
        # Should rather check if taxon in children_taxa.
        if not treebest_isdup:
            logger.warning("The node %r -> %s "
                  "is marked as a *speciation* but looks like a duplication"
                  " after taxon information.",
                  node.name, children_names)
        elif (treebest_isdup) and taxon not in children_taxa:
            logger.warning("The node %r -> %s "
                    "is marked as a *duplication* but is ambiguous: "
                    "intermediate speciation nodes are missing.",
                  node.name, children_names)

        if len(children_taxa) > 1:  # or treebest_isdup == 'Y'
            logger.warning("The node %r -> %s "
                  "is a duplication + a speciation. Not "
                  "truly reconciled tree.",
                  node.name, children_names)
        elif len(node.children) == 1:
            logger.warning("The node %r -> %s "
                  " is a duplication/speciation with one "
                  "descendant.",
                  node.name, children_names)
        return 'dup'
    else:
        ### FIXME
        if len(children_taxa) == 1:
            msg = "The node %r -> %s"
            if len(node.children) > 1:
                msg += " is a duplication + a speciation. Not " \
                       "truly reconciled tree."
                #return 'dup'
            #else:
            #    msg += " is a duplication/speciation with one " \
            #           "descendant."
                #if not node.is_root(): event = 'dup'
                # TODO: check this in the duplication block above.
                logger.warning(msg, node.name, children_names)

        return 'spe'

# TODO:
# infer_gene_event_treebest
# infer_gene_event_taxa
# infer_gene_events that warns if mismatches between the 2.

def get_node_children_ete3(node, *args):
    return node.children

def get_node_children_dyogenprot(node, tree):
    """return children from a node in a LibsDyogen.myProteinTree format."""
    return [ch for ch,_ in tree.data.get(node, [])]

def get_node_name_ete3(node, *args):
    return node.name

def get_node_name_dyogenprot(node, tree):
    nodeinfo = tree.info[node]
    return nodeinfo.get('gene_name',
                        (nodeinfo['taxon_name'] +
                         nodeinfo.get('family_name', '').split('/')[0]))

def infer_gene_event_taxa(node, taxon, children_taxa,
                          get_node_children=get_node_children_ete3,
                          get_name=get_node_name_ete3, *args):

    """Use taxon information to tell whether a gene tree node is:
    - a leaf,
    - a speciation,
    - or a duplication.
    
    param: `children_taxa` must be a set (because the number of *uniq* elements
           is used).
    param: *args: extra arguments to be passed to `get_node_children`/`get_node_name`.
    """
    ### "ambiguous" or "dupspe" value should be returned in case of doubt.

    if not children_taxa:
        return 'leaf'

    children = get_node_children(node, *args)
    nodename = get_name(node, *args)
    
    if taxon in children_taxa:
       #(len(children_taxa) == 1 and taxon in children_taxa): #\
                    #and len(children) > 1 \
       #or \
       #(len(children) == 1 and taxon in children_taxa):
        # Should rather check if taxon in children_taxa.
        if len(children_taxa) > 1:
            logger.warning("The node %r -> %s "
                  "is a duplication + a speciation. Not "
                  "truly reconciled tree." %
                  (nodename, [get_name(ch, *args) for ch in children]))
        elif len(children) == 1:
            logger.warning("The node %r -> %s "
                  " is a duplication/speciation with one "
                  "descendant." %
                  (nodename, [get_name(ch, *args) for ch in children]))
        return 'dup'
    else:
        ### FIXME
        if len(children_taxa) == 1:
            msg = "The node %r -> %s"
            if len(children) > 1:
                msg += " is a duplication + a speciation. Not " \
                       "truly reconciled tree."
                #return 'dup'
            #else:
            #    msg += " is a duplication/speciation with one " \
            #           "descendant."
                #if not node.is_root(): event = 'dup'
                # TODO: check this in the duplication block above.
                logger.warning(msg, nodename, [get_name(ch, *args) for ch in children])

        return 'spe'


from collections import defaultdict
from dendro.bates import dfw_pairs_generalized, dfw_descendants_generalized, iter_leaves

def prottree_extract_genecounts(proteintrees, ancestor, phyltree,
                                speciesset=set(('Homo sapiens',)),
                                keeponly=None):
    """Walk each tree to find the ancestor taxon, then output the gene counts
    per descendant species, as well as the list of human genes.

    Param: `keeponly` should be None, "stem" or "crown".
    """
    
    # for myProteinTree class
    def get_children(tree, node):
        return [c for c,_ in tree.data.get(node, [])]

    ancestor_ancgenes   = []
    ancestor_genecounts = []  # In each species
    ancestor_spgenes    = []  # For genes from speciesset.
    ancestors = []  # Ancestors at which we got the node.

    clades_before_ancestor = set(phyltree.dicLinks[phyltree.root][ancestor])  # includes anc
    clades_after_ancestor = phyltree.allDescendants[ancestor]  # includes anc
    clades_outside_ancestor = (phyltree.outgroupSpecies[ancestor]
                               | ( phyltree.getTargetsAnc("/" + ancestor)
                                   - clades_before_ancestor )
                              )

    for tree in proteintrees:
        info = tree.info
        #if tree.root == 16401:
        #    import ipdb; ipdb.set_trace()

        if info[tree.root]['taxon_name'] in clades_before_ancestor:
            for parent, node in dfw_pairs_generalized(tree, get_children,
                                                      include_root=True):
                taxon_node = info[node]['taxon_name']
                if parent is not None:
                    taxon_parent = info[parent]['taxon_name']
                    isdup = (info[parent]['Duplication'] > 1)
                else:
                    # then the child node (i.e the root) should be kept if it
                    # is exactly equal to ancestor.
                    taxon_parent = taxon_node  # We know it's before (or equal).
                    isdup = False

                if taxon_parent in clades_before_ancestor and \
                        taxon_node in clades_after_ancestor:

                    if taxon_node != ancestor and taxon_parent != ancestor:
                        # The branch "jumps" over this ancestral species
                        # Process the most basal node. WHY? NO!
                        #node = parent
                        #taxon_node = taxon_parent
                        pass
                    elif taxon_node != ancestor and not isdup:
                        # The parent is an 'ancestor' speciation so this node should already 
                        # have been taken into account.
                        # Except if the parent is the root.
                        #assert info[parent]['family_name'] in ancestor_ancgenes,\
                        #    "At %d->%d: %s %s ->..." % (parent, node,
                        #                            taxon_parent,
                        #                            info[parent]['family_name'])
                        #
                        continue

                    if keeponly == 'crown' and info[node]['Duplication']>1:
                        # node is at the ancestor.
                        # We don't want to keep it if it is a duplication.
                        continue
                    
                    nodename = info[node]['family_name']
                    #assert nodename not in ancestor_ancgenes
                    ancestor_ancgenes.append(nodename)
                    ancestors.append((taxon_parent, taxon_node))
                    
                    spgenes = {sp: [] for sp in speciesset}
                    gene_counts = defaultdict(int)
                    for leaf in iter_leaves(tree, get_children, queue=[node]):
                        taxon_leaf = phyltree.officialName[info[leaf]['taxon_name']]
                        if taxon_leaf not in phyltree.listSpecies:
                            # Error because of tree.data.pop
                            #import ipdb; ipdb.set_trace()
                            errmsg = "%d '%s' is not a species! (tree %d, node %d)"\
                                          % (leaf, taxon_leaf, tree.root, node)
                            if taxon_leaf in phyltree.allNames:
                                if keeponly != 'stem':
                                    logger.error(errmsg)
                                continue
                            else:
                                raise RuntimeError(errmsg)
                        gene_counts[taxon_leaf] += 1
                        #if taxon_leaf == 'Homo sapiens':
                        #    spgenes.append(info[leaf]['gene_name'])
                        if taxon_leaf in speciesset:
                            spgenes[taxon_leaf].append(info[leaf]['gene_name'])

                    ancestor_genecounts.append(gene_counts)
                    #ancestor_spgenes.append(tuple(spgenes))  # tuple: Important for finding `()` in Series.
                    ancestor_spgenes.append({sp: tuple(genes) for sp,genes in
                                             spgenes.items()})

                    # Now do we wan't to score descendant ancestor nodes?
                    if keeponly == 'stem':
                        try:
                            tree.data.pop(node)
                        except KeyError:
                            # This node is not in data because it is a leaf:
                            assert node in tree.info
                elif taxon_node not in clades_before_ancestor:
                    # Avoid visiting outgroups and strict ingroups.
                    try:
                        # I don't know how, but this messes up the dfw iteration.
                        tree.data.pop(node)
                        logger.debug("pop data of %d (%s)", node, taxon_node)
                    except KeyError:
                        #assert 'gene_name' in info[node]
                        logger.debug("Ignore data of %d (%s)", node, taxon_node)
                        pass  # It's a leaf.

    assert len(ancestor_ancgenes) == len(ancestors)
    assert len(ancestor_ancgenes) == len(ancestor_spgenes)
    assert len(ancestor_ancgenes) == len(ancestor_genecounts)

    return ancestors, ancestor_ancgenes, ancestor_genecounts, ancestor_spgenes


def is_robust(tree, root, phyltree):
    info = tree.info
    def get_children(tree, node):
        return [c for c,_ in tree.data.get(node, [])]

    for parent, children in dfw_descendants_generalized(tree, get_children, queue=[root]):
                                                        #include_leaves=True):
        parent_taxon = info[parent]['taxon_name']
        children_taxa = [info[c]['taxon_name'] for c in children]
        expected_taxa = [anc for anc,_ in phyltree.items.get(parent_taxon, [])]
        if set(children_taxa) != set(expected_taxa):
            return False
    return True


def prottree_list_robusts(proteintrees, phyltree, ancestor=None, from_anc=True):
    if ancestor is None:
        ancestor = phyltree.root
    def get_children(tree, node):
        return [c for c,_ in tree.data.get(node, [])]

    clades_before_ancestor = set(phyltree.dicLinks[phyltree.root][ancestor])  # includes anc
    clades_after_ancestor = phyltree.allDescendants[ancestor]  # includes anc

    for tree in proteintrees:
        info = tree.info
        rootfamily = info[tree.root]['family_name']

        for parent, node in dfw_pairs_generalized(tree, get_children,
                                                  include_root=True):
            taxon = info[node]['taxon_name']
            if parent is not None:
                #logger.info('node %s%s', taxon, info[node]['family_name'])
                parent_taxon = info[parent]['taxon_name']
                isdup = (info[parent]['Duplication'] > 1)
            else:
                parent_taxon = None

            if parent_taxon not in clades_after_ancestor:
                if taxon in clades_after_ancestor:
                    if taxon != ancestor and from_anc:
                        continue
                    #if isdup: #TODO: handle this case
                    logger.info('Checking %s at %s', rootfamily, taxon)
                    if is_robust(tree, node, phyltree):
                        yield phyltree.fileName[taxon], info[node]['family_name']
                else:
                    continue

                #if parent in clades_after_ancestor
                    # Start checking for robustness.

