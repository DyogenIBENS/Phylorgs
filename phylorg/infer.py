#!/usr/bin/env python3


"""Phylogenetic inference tools based on sequences and trees."""


import numpy as np
from functools import reduce
from dendro.bates import rev_dfw_descendants
from seqtools.arrayal import presence_matrix

import logging
logger = logging.getLogger(__name__)


class Parsimony(object):
    #def parsimony_score(alint, tree, seqlabels, minlength=66, get_children=None,
    #                parts=None):
    """
    Computes the column-wise parsimony score based on the provided tree.
    
    The input alignment should be a matrix of integers (see seqtools.arrayal)

    parts: [tuples of indices] e.g. [(0,1,2)] represents leaves 0,1,2 as a clade.
    The outgroup must NOT be added into parts.
    """
    def __init__(self, alint, tree, seqlabels, minlength=66, get_children=None,
                 parts=None):

        self.alint = alint
        self.tree = tree
        self.seqlabels = seqlabels
        assert alint.shape[0] == len(seqlabels)
        self.minlength = minlength
        # Bio.Phylo by default.
        self.get_children = (lambda tree, node: node.clades) if get_children is None else get_children
        try:
            self.root = tree.clade
        except AttributeError:
            self.root = tree
        self.iter_tree = list(rev_dfw_descendants(tree, self.get_children, include_leaves=True,
                                        queue=[self.root]))

        # Holds the currently seen nodes, and their parsimony score and sequence.
        self.process_sequences = {}
        self.anc_states = {}  # These could be stored as the same object.

        # column-wise parsimony score
        self.score = np.zeros(alint.shape[1], dtype=int)
        
        self.parts = [set(p) for p in parts] if parts is not None else []
        self.node_to_part = {leaf_n: i for i,p in enumerate(self.parts)
                        for leaf_n in p}
        self.part_scores = [self.score.copy() for p in self.parts]
        self.part_branch_nbs = [0] * len(self.parts)

    def rootwards(self):
        # Index of encountered leaves/sequences
        leaf_nb = self.alint.shape[0] - 1

        process_sequences = self.process_sequences
        parts = self.parts
        merged_parts = [dict() for _ in parts]  # {node: set(seen_part_leaves)}
        # When merged_parts[i][node] == parts[i]: we found the MRCA of part.
        branch_nb = 0
        for parent, children in self.iter_tree:
            logger.debug('* %r -> %s', parent, children)
            if len(children) > 2:
                logger.warning("%d > 2 children at %r", len(children), parent)
            if not children:
                # Parent is a leaf. Obtain the sequence.
                assert parent.name == self.seqlabels[leaf_nb], \
                    "The alignment is not ordered as the tree. Seq %d: %s != leaf %s" \
                        % (leaf_nb, self.seqlabels[leaf_nb], parent.name)
                process_sequences[parent] = presence_matrix(self.alint[np.newaxis,leaf_nb,:], self.minlength)
                try:
                    p = self.node_to_part[parent] = self.node_to_part[leaf_nb]
                    merged_parts[p][parent] = set((leaf_nb,))
                    logger.info('Assigned leaf %d %r to part %d', leaf_nb, parent, p)
                except KeyError:  # When there is no part for this leaf.
                    pass

                leaf_nb -= 1
            else:
                # .pop(ch) ?
                try:
                    children_seqs = [process_sequences[ch] for ch in children]
                except KeyError as err:
                    #logger.debug('Processed sequences:\n%s',
                    #             '\n'.join('%r: %s' % (node, pseq.astype(np.int32))
                    #                       for node, pseq in process_sequences.items()))
                    logger.error('parent = %r; leaf_nb = %s', parent, leaf_nb)
                    raise
                children_inter = reduce(np.logical_and, children_seqs)
                children_union = reduce(np.logical_or, children_seqs)
                #print(children_inter, children_union, sep='\n')
                # Add one to each column where a substitution is needed
                empty_inter = ~children_inter.any(axis=0)

                # The new nucleotide set is the intersection if it's not empty,
                # otherwise the union
                process_sequences[parent] = children_inter
                process_sequences[parent][:, empty_inter] = children_union[:, empty_inter]
                #process_sequences[parent] = np.where(empty_inter,
                #                                     children_union,
                #                                     children_inter)

                if parts:
                    children_parts = list(set((self.node_to_part.get(ch) for ch in children)))
                    if len(children_parts) == 1:
                        # Still inside the same clade. This must be exclusive (monophyletic).
                        #assert all descendant leaves are in the same part.
                        p = children_parts[0]
                        try:
                            merged_parts[p][parent] = set.union(*(merged_parts[p][ch] for ch in children))
                            self.node_to_part[parent] = p
                            logger.info('Assigned node %r to part %d', parent, p)
                            # Just update the intra-clade score.
                            self.part_scores[p] += empty_inter
                            self.part_branch_nbs[p] += len(children)
                            # Still inside one clade, so skip updating the *global* score.
                            continue
                        except TypeError:
                            # p is None because no child was found in node_to_part
                            pass  # So we update the outgroup score...
                    else:
                        # We leave one or more clades.
                        # 1. assert that at least one is monophyletic.
                        part_oldest_node = []  # Monophyly check variable
                        # 2. TODO: Map substitutions to the descendent clades
                        for ch in children:
                            try:
                                p = self.node_to_part[ch]
                            except KeyError:
                                # this is the outgroup, skip. The global score should be updated though
                                continue
                            if merged_parts[p][ch] == parts[p]:
                                # The child is the exclusive MRCA (It contains exactly the part species).
                                # We update the score of the clade:
                                # NOTE: this can't be done on the first postorder traversal...
                                # Because we need to know the parent state to orient the substitution!

                                # Uncomment this to add the stem substitutions:
                                #part_scores[p] += empty_inter
                                #part_branch_nbs[p] += 1  # Add the stem branch
                                part_oldest_node.append(ch)

                        if not part_oldest_node:
                            raise ValueError("The given parts do not make monophyletic clades."
                                             + ';'.join('(%s)' % (
                                                     ','.join(self.seqlabels[i] for i in p))
                                                     for p in parts))
                        # Here, the current score will be updated (see just below)

                branch_nb += len(children)  # Except if children & merged_parts
                self.score += empty_inter
            #print(process_sequences[parent])

        # Number of branches of the **unrooted** tree:
        branch_nb -= 1

        if parts:
            return (self.score.astype(float), *self.part_scores), (branch_nb, *self.part_branch_nbs)

        self.branch_nb = branch_nb
        return self.score.astype(float) / branch_nb
    
    def leafwards(self):
        # Retraverse leafwards to get the ancestral states. Must have traversed rootwards first.
        self.score_leafward = np.zeros(self.alint.shape[1], dtype=int)
        anc_states = self.anc_states
        anc_states[self.root] = self.process_sequences.pop(self.root)
        for parent, children in reversed(self.iter_tree):
            parent_state = anc_states[parent]
            for ch in children:
                child_possibles = self.process_sequences.pop(ch)
                branch_inter = np.logical_and(parent_state, child_possibles)  # No substitution needed.
                branch_union = np.logical_or(parent_state, child_possibles)  # choice needed
                non_empty_inter = branch_inter.any(axis=0)
                self.score_leafward += (~non_empty_inter).sum()  # How many sites with empty intersection (i.e. change)

                # If non empty intersection, the child nucleotides are intersected with those of the parent
                # otherwise, unchanged.
                anc_states[ch] = child_possibles
                anc_states[ch][:, non_empty_inter] = branch_inter[:, non_empty_inter]
        #TODO: change self.part_scores
        #FIXME
        return self.score_leafward / self.branch_nb

    def __call__(self):
        self.rootwards()
        return self.leafwards()


# For backward compatibility:
def parsimony_score(*args, **kwargs):
    #FIXME: should return after leafwards: #return Parsimony(*args, **kwargs)()
    return Parsimony(*args, **kwargs).rootwards()
