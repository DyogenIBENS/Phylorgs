#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Measure indel frequencies and length distributions on each branch of the,
phylogeny, based on a given alignment.

Quick and crude parsimonious estimation.
"""

from functools import reduce
import numpy as np
from numpy import array
from Bio import AlignIO
from dendro.bates import rev_dfw_descendants, dfw_descendants_generalized
from dendro.parsers import parserchoice
from dendro.any import methodchoice

import logging
logger = logging.getLogger(__name__)


def map_gaps(align, codon=False):
    """Return the positions of gaps in each sequence (0-indexed, end included)"""
    step = 3 if codon else 1
    gap = '-' * step
    gaps = {}
    N = align.get_alignment_length()
    for record in align:
        rec_gaps = []
        gap_start = 0  # Also indicates if we are already inside a gap
        prev_is_gap = False
        for i in range(N // step):
            resid = record.seq[i*step:(i+1)*step]
            is_gap = (resid==gap)
            if is_gap != prev_is_gap:
                if is_gap:
                    # a gap opens
                    gap_start = i
                    gap_end = None
                else:
                    # a gap terminates
                    rec_gaps.append((gap_start, i-1))
                    gap_start = None
            prev_is_gap = is_gap
        if prev_is_gap:
            rec_gaps.append((gap_start, i))
        gaps[record.id] = rec_gaps
    return gaps


def combine_gaps(gaps1, gaps2):
    """List hypothetical parent gaps for children gaps `gaps1` and `gaps2`:
    - union of children gaps;
    - when some children gaps intersect, list all possible segments:
        e.g: (1,3) vs (1,2) -> (1,2), (3,3), (1,3)

    Gaps must be sorted.
    NOTE: gaps *include* the end coordinate.
    TODO: when a gap in 1 is next to a gap in 2, consider the fusion.
    """
    #starts1, ends1 = zip(*gaps1)
    #starts2, ends2 = zip(*gaps2)
    newgaps = []

    i = 0  # currently examined gap in gaps1
    j = 0  # currently examined gap in gaps2
    # Also spot intervening gaps (e.g. the segment between two gaps2 that is overlapped by a gap1)
    ends_under1 = []  # the last ended gaps2 that are overlapped by the current gap1
    ends_under2 = []  #

    while i < len(gaps1) and  j < len(gaps2):
        # At each iteration we process all the intersecting gaps
        s1, e1 = gaps1[i]
        s2, e2 = gaps2[j]
        logger.debug('# i,j = %d,%d ; (%d,%d), (%d, %d); ends_in1=%s, ends_in2=%s', i,j, s1,e1, s2,e2, ends_under1, ends_under2)
        if e1 < s2:
            newgaps.append((s1,e1))
            i += 1
            while ends_under1:
                newgaps.append((ends_under1.pop(), e1))
            continue
        elif e2 < s1:
            newgaps.append((s2, e2))
            j += 1
            while ends_under2:
                newgaps.append((ends_under2.pop(), e2))
            continue

        # Here they intersect.
        #newgaps.append((min(s1, s2), max(s1, s2)-1))
        # Add the preceeding segment
        if s1 < s2:
            newgaps.append((s1,s2-1))
            if e1 < e2:
                newgaps.append((s2, e1)) # intersecting segment
            # intervening segments:
            newgaps.extend((e_u1, s2-1) for e_u1 in ends_under1)
        elif s2 < s1:
            newgaps.append((s2,s1-1))
            if e2 < e1:
                newgaps.append((s1, e2)) # intersecting segment
            # intervening segments:
            newgaps.extend((e_u1, s2-1) for e_u1 in ends_under1)
        #else:
        #    newgaps.append((s1, min(e1, e2))) # intersecting segment

        if e1 < e2:
            newgaps.append((e1+1, e2))
            newgaps.append((s1, e1))  # Append the whole gap!
            # Store this end while gap2 is not closed.
            ends_under2.append(e1+1)
            i += 1
        elif e2 < e1:
            newgaps.append((e2+1, e1))
            newgaps.append((s2, e2))  # Append the whole gap!
            # Store this end while gap1 is not closed.
            ends_under1.append(e2+1)
            j += 1
        else:
            #if s1 < s2:
                #logger.debug('  ends_under1: %s; newgaps: %s', , ends_under1, newgaps)
                #while ends_under1:
                #    newgaps.append((ends_under1.pop(), e1))
            #elif s2 < s1:
                #logger.debug('  ends_under1: %s; newgaps: %s', , ends_under1, newgaps)
                #while ends_under2:
                #    newgaps.append((ends_under2.pop(), e2))
            newgaps.append((s1, e1))
            if s1 != s2:
                newgaps.append((s2, e2))
            i += 1
            j += 1
        logger.debug('  newgaps = %s', newgaps)

    logger.debug('# Final i,j = %d,%d ; gaps1=%s, gaps2=%s ; ends_in1=%s, ends_in2=%s', i,j, gaps1[i:], gaps2[j:], ends_under1, ends_under2)
    newgaps.extend(gaps1[i:] + gaps2[j:])
    newgaps.sort()  # Needed for recursive calls! #TODO: fill it in a sorted manner.
    return newgaps #, index1, index2 #(instead of the reindex function)


def reindex_states(states, coords, new_coords):
    """Apply state values in the new order.
    If new_coords contains *included* gaps that are overlapped by some gaps,
    their state is set to the overlapping gap state."""
    newstates = np.array([[1,0]]*len(new_coords))  # default state is not a gap
    coords_idx = {c:i for i,c in enumerate(coords)}
    for j,g in enumerate(new_coords):
        i = None
        try:
            i = coords_idx[g]
        except KeyError:
            for i, c in enumerate(coords):
                if (c[0] <= g[0] <= g[1] <= c[1]):
                    break
            else:
                continue
        newstates[j,:] = states[i,:]
    return newstates



def rootwards_gap_states(align, tree, get_children, get_label, root, codon=False):
    """Walk tree rootwards and compute each intermediate node state.
    
    State is a binary representation: 10='nogap' 01='gap' 11='unknown'
    """
    gap_coords = map_gaps(align, codon)  # gap positions. used to index states.

    #FIXME: use 2 for [1,0], 1 for [0,1] and 3 for [1,1]. any() becomes bitwise_inter>0.
    states = {sp: np.array([[0,1]]*len(gaps)) for sp, gaps in gap_coords.items()}
    #FIXME: what if new gaps occur (as fusion, extension, etc)

    # At each node, the parent gap index should contain all possible ancestral and child gaps

    score = 0

    for node, childnodes in rev_dfw_descendants(tree, get_children, queue=[root]):

        nodename = get_label(tree, node)
        children = [get_label(tree, chn) for chn in childnodes]

        combined_gaps = reduce(combine_gaps, [gap_coords[ch] for ch in children])

        #gaps[nodename] = list(sorted(all_gaps))

        #for pos in gaps[nodename]:
        ch_states = [reindex_states(states[ch], gap_coords[ch], combined_gaps)
                     for ch in children]
        state_union = reduce(np.logical_or, ch_states)
        state_inter = reduce(np.logical_and, ch_states)

        assert state_union.shape == (len(combined_gaps), 2)

        change = ~(state_inter.any(axis=1))
        assert change.shape == (len(combined_gaps),), 'change.shape=%s' % (change.shape,)
        score += change.sum()

        states[nodename] = np.where(change[:,None], state_union, state_inter)
        gap_coords[nodename] = combined_gaps
            
        assert states[nodename].shape == (len(combined_gaps), 2)
        #common_gaps = set.intersection(*(set(gaps[ch]) for ch in children))
        #specific_gaps = all_gaps - common_gaps
        ##intersecting_gaps
        #for i, ch in enumerate(children):
        #    if ch in uncertain:
        #        for g in uncertain[ch]:
        #            if all((g not in gaps[s])
        #                   for s in children[:i] + children[i+1:]):
        #                uncertain[ch].remove(g)
        #            elif all((g in gaps[s])
        #                     for s in children[:i] + children[i+1:]):
        #                common_gaps.add(g)
        #            else:
        #                specific_gaps.add(g)

        #gaps[nodename] = list(sorted(common_gaps))

        #uncertain[nodename] = list(sorted(specific_gaps))

    return gap_coords, states


def leafwards_branch_events(tree, get_children, get_label, root, gap_coords, states):

    rootname = get_label(tree, root)
    newstates = {rootname: states[rootname]}  # all indexed on gap_coords[rootname]
    rootcoords = gap_coords[rootname]

    branch_inser = {}
    branch_del = {}
    #branch_len
    for node, childnodes in dfw_descendants_generalized(tree, get_children, queue=[root],
            include_leaves=True):
        nodename = get_label(tree, node)
        if not childnodes:
            # This is a leaf. Restore the observed states
            newstates[nodename] = reindex_states(states[nodename], gap_coords[nodename], rootcoords)

        nodestate = newstates[nodename]
        assert nodestate.shape == (len(rootcoords), 2)

        for chn in childnodes:
            ch = get_label(tree, chn)
            branch_inser[(nodename, ch)] = []
            branch_del[(nodename, ch)] = []

            ch_state = reindex_states(states[ch], gap_coords[ch], rootcoords)
            assert ch_state.shape == nodestate.shape
            branch_inter = np.logical_and(nodestate, ch_state)
            branch_union = np.logical_or(nodestate, ch_state)
            assert branch_inter.shape == nodestate.shape
            change = ~(branch_inter.any(axis=1))
            #if nodename=='x' and ch=='a':
            #    logger.debug('nodestate=%s ch_state=%s change = %s', nodestate.tolist(), ch_state.tolist(), change.tolist())
            assert change.shape == (len(rootcoords),), 'change.shape=%s' % (change.shape,)

            newstates[ch] = np.where(change[:,None], branch_union, branch_inter)
            #TODO: at a leaf, states are constrained
            for i, gap in enumerate(rootcoords):
                if change[i]:
                    if nodestate[i,1]: # it was a gap
                        branch_inser[(nodename, ch)].append(gap)
                    else:
                        branch_del[(nodename, ch)].append(gap)


    return newstates, rootcoords, branch_inser, branch_del


def branch_len_distribs(branch_events):
    """branch_events are branch_inser/branch_del as outputted by leafwards_branch_events"""
    distribs = {}
    for branch, events in branch_events.items():
        distribs[branch] = np.array([e-s+1 for s,e in events])
    return distribs



