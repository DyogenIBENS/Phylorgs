#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Simulate gene trees with duplication/loss along a species tree."""


# 1. simple constant model:
# dup rate = loss rate; all trees have the same rate.

# 2. simple "site" model:
# dup rate = loss rate; gamma-distributed across gene trees.

# 3. unstable model:
# dup rate and loss rate are independent.

# 4. not-so-stable model:
# amplitude = dup rate + loss rate: follows a gamma distrib;
# (dup rate - loss rate) follows a normal distrib of std=amplitude/5
# NOTE: birth should be modeled.

# 5. family-varying & species-branch-varying dup/loss rates.


import sys
import os.path as op
import argparse as ap
import random
import numpy as np
#import ete3
from LibsDyogen import myPhylTree, myProteinTree
from dendro.bates import dfw_pairs_generalized
from dendro.any import myPhylTree as phyltree_methods
import logging
logger = logging.getLogger(__name__)


phyl_items = phyltree_methods.get_items

LOST = -1


def sim_notsostable(phyltree, n0=10000, gamma_shape=2., gamma_mean=10., 
                    root=None, overflowfactor=5000):
    """n0: starting number of gene trees.
    gamma_mean = gamma_shape * gamma_scale gives the mean number of duplication/loss per My;
    """
    if root is None:
        root = phyltree.root
    # Condition trees on their survival to the present?
    must_survive = False
    gamma_scale = gamma_mean/gamma_shape

    #for convenience, I gather all trees as sister branches of a "forest":
    #forest = myProteinTree.ProteinTree(data={0: [(i,None) for i in range(1, n0+1)]},
    #                                   root=0)
    forest = myProteinTree.ProteinTree(info={i: dict(taxon_name=root,
                                                     Duplication=0,
                                                     family_name='Fam%d' % i)
                                             for i in range(1, n0+1)})
    #parent_cache = {}  # child_id: (parent_id, item_index)  # Necessary to set a branch length after event.
    #myProteinTree.nextNodeID
    nextnode = n0+1
    nextloss = -1

    amplitudes = np.random.gamma(gamma_shape, gamma_scale, size=n0)
    netdivrates = np.random.normal(0, amplitudes/5, size=n0)
    #netdivrates = np.random.gamma(2, amplitudes/2, size=n0) - amplitudes
    duprates = (amplitudes + netdivrates)/2
    lossrates = (amplitudes - duprates)
    dupprobas = duprates/amplitudes  # proba that an event is a dup.
    #dupprobas = np.random.beta(size=n0)
    assert (duprates>=0).all() and (lossrates>=0).all()
    logger.debug('Amplitudes: mean=%g ±%g [%g, %g]', amplitudes.mean(),
                 amplitudes.std(), amplitudes.min(), amplitudes.max())
    logger.debug('net div rates: mean=%g ±%g [%g, %g]', netdivrates.mean(), netdivrates.std(), netdivrates.min(), netdivrates.max())

    # "heads" of the lineages: no node yet, but pointing back to the parent node id.
    ancestor_genes = {root: np.arange(1, n0+1)}  # anc: list(ids)
    ancestor_amplitudes = {root: amplitudes}  # sort?
    ancestor_dupprobas = {root: dupprobas}
    #ancestor_netdivrates = {root: netdivrates}
    #ancestor_duprates = {root: duprates}

    for (parent, _), (child, dist) in dfw_pairs_generalized(phyltree, phyl_items,
                                                   queue=[(None, (root, 0))],
                                                   include_root=False):
        anc_genes = ancestor_genes[parent].copy()  # COPY!
        n_genes = len(anc_genes)
        logger.info('Branch %s - %s (%g My); %d genes', parent, child, dist, n_genes)
        amplitudes = ancestor_amplitudes[parent].copy()
        #netdivrates = ancestor_netdivrates[parent].copy()
        #duprates = ancestor_duprates[parent].copy()
        dupprobas = ancestor_dupprobas[parent].copy()
        head_times = np.zeros(n_genes)  # /!\ lineage heads become nodes after an event.
        last_times_before_speciation = np.zeros(n_genes)

        # So the minimum waiting time (1st event), is a poisson-gamma mixture
        _rounds = 0
        while (head_times < dist).any():
            logger.debug('Events round %d: %d genes, net div rate %g ±%g, P(dup) %g ±%g',
                         _rounds, n_genes, amplitudes.mean(), amplitudes.std(),
                         dupprobas.mean(), dupprobas.std())
            logger.debug('Forest data = %s\nForest info = %s', forest.data, forest.info)
            _rounds += 1
            assert amplitudes.shape[0] == n_genes
            assert dupprobas.shape[0] == n_genes
            event_times = np.random.exponential(1./amplitudes, size=n_genes)
            head_times += event_times
            which_events = np.flatnonzero( head_times < dist )

            # Draw dup or loss
            isdup = np.random.random(len(which_events)) > dupprobas[which_events]
            which_dup = which_events[isdup]
            ndup = len(which_dup)
            nloss = len(which_events) - ndup

            notlost = sorted(set(range(n_genes)).difference(which_events[~isdup]))
            #notlost = np.flatnonzero(np.bincount(which_events[~isdup], minlength=n_genes) == 0)
            assert set(notlost) == set(np.flatnonzero(np.bincount(which_events[~isdup], minlength=n_genes) == 0))
            assert len(notlost) == n_genes - nloss
            
            # append duplicates
            for duploc in which_dup:
                parentg = anc_genes[duploc]
                # Append the duplication node.
                parentdata = forest.data.setdefault(parentg, [])
                parentdata.append((nextnode, event_times[duploc]))
                try:
                    family_name = forest.info[parentg]['family_name']
                except KeyError:
                    assert parent == root
                    family_name = 'Fam%d' % duploc
                forest.info[nextnode] = dict(taxon_name=child,
                                             Duplication=2,
                                             family_name=family_name)
                # Replace the heads
                anc_genes[duploc] = nextnode
                nextnode+=1

            logger.debug('parent of lost heads: %s', anc_genes[which_events[~isdup]])
            for lossloc in which_events[~isdup]:
                parentg = anc_genes[lossloc]
                parentdata = forest.data.setdefault(parentg, [])
                parentdata.append((nextloss, event_times[lossloc]))  # Append the event as a node.
                forest.info[nextloss] = {'Duplication': -1, 'taxon_name': child,
                                         'family_name': forest.info[parentg]['family_name']}
                nextloss -= 1

            transmitted = np.concatenate((notlost, which_dup))
            logger.debug('transmit:\n%s +\n%s', np.array(notlost), which_dup)
            anc_genes = anc_genes[transmitted]
            last_times_before_speciation = np.where(head_times<dist,
                                                    head_times,
                                                    last_times_before_speciation)[transmitted]
            amplitudes = amplitudes[transmitted]
            dupprobas  = dupprobas[transmitted]
            head_times = head_times[transmitted]

            n_genes += ndup - nloss
            assert n_genes == len(anc_genes)
            assert n_genes == len(head_times)
            if n_genes > n0*overflowfactor:
                logger.error('Overflow: n_genes > n0*%d. Return.', overflowfactor)
                break
                #return forest

        # After the duplications, we reach the speciation node
        ancestor_amplitudes[child] = amplitudes
        ancestor_dupprobas[child] = dupprobas
        for i, (parentg, last_t) in enumerate(zip(anc_genes, last_times_before_speciation)):
            parentdata = forest.data.setdefault(parentg, [])
            parentdata.append((nextnode, dist - last_t))
            family_name = forest.info[parentg]['family_name']
            forest.info[nextnode] = dict(taxon_name=child,
                                         Duplication=0,
                                         family_name=family_name)
            # Replace the heads
            anc_genes[i] = nextnode
            # Update the parent_cache
            nextnode+=1
        
        forest_nodes, all_nodes = set(forest.info), set(range(1, nextnode))
        assert forest_nodes >= all_nodes, all_nodes - forest_nodes

        if n_genes > n0*overflowfactor:
            break

        ancestor_genes[child] = anc_genes  # Needs to be copied later.

    return forest



def main():
    logging.basicConfig(format='%(levelname)s:%(funcName)s:%(message)s')
    logger.setLevel(logging.DEBUG)

    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('phyltreefile')
    parser.add_argument('-n', '--n0', type=int, default=100)
    parser.add_argument('-r', '--root', default='Simiiformes')
    parser.add_argument('--gamma-shape', '--gs', type=float, default=2)
    parser.add_argument('--gamma-mean', '--gm', type=float, default=0.1)
    parser.add_argument('--overflowfactor', '--of', type=int, default=10000)
    clargs = vars(parser.parse_args())

    phyltree = myPhylTree.PhylogeneticTree(clargs.pop('phyltreefile'))

    forest = sim_notsostable(phyltree, **clargs)
    for i in range(1, clargs['n0']+1):
        forest.printTree(sys.stdout, i)


if __name__ == '__main__':
    main()
