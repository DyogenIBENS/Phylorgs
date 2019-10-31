#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Various indices of tree imbalance.

See [Blum & François 2005,
On statistical tests of phylogenetic tree imbalance:
The Sackin and other indices revisited](http://membres-timc.imag.fr/Olivier.Francois/sackin.pdf)
"""

from collections import defaultdict
from queue import deque
from math import factorial
import numpy as np
from scipy.special import binom
from LibsDyogen import myTools
import logging
logger = logging.getLogger(__name__)


def sackin(tree):
    return sum(len(l.get_ancestors()) for l in tree.iter_leaves())

def expectation_sackin(tree):
    """For a dichotomic tree generated by the Yule model."""
    n = len(tree)  # Number of leaves
    return 2. * n * sum(1./j for j in range(2, n+1))

def normed_sackin(tree):
    # Beware. This is not a proper normalisation, as the standard deviation
    # of the Sackin index doesn't really depend on `n`
    n = len(tree)
    return (sackin(tree) - expectation_sackin(tree)) / float(n)


@myTools.memoize
def sackin_minmax(n):
    sackin_maximum = sum(range(1,n)) + n-1  # easy: caterpillar tree:

    # NOT true.
    ## For the minimum: decompose n into the largest packs of powers of 2.
    #n_binary = [int(x) for x in '{:b}'.format(n)]
    #n_binary.reverse()  # Largest powers to the right.
    #print('Binary repr = {:b}'.format(n))

    #sackin_minimum = 0
    #previous_depth = 0
    #max_p = len(n_binary)
    #for p in range(max_p-1, -1, -1):
    #    if n_binary[p]:
    #        depth = p
    #        remainder = any(n_binary[:p])
    #        if remainder:
    #            previous_depth += 1
    #        depth += previous_depth
    #        print('+= 2^{p} × ({p} + {pd}) [{d}]'.format(p=p,
    #                pd=previous_depth, d=depth))
    #        nleaves = 2**p
    #        sackin_minimum += nleaves*depth
    #        if not remainder:
    #            break

    sackin_minimum = 0
    nextnodes = deque([(1, n)])  # (depth, number of leaves)
    while nextnodes:
        depth, nleaves = nextnodes.popleft()
        if nleaves > 1:
            sackin_minimum += nleaves
            half_n = nleaves // 2
            nextnodes.extend(((depth+1, half_n),
                              (depth+1, half_n + nleaves%2)))

    return sackin_minimum, sackin_maximum


@myTools.memoize
def sackin_distribs(N):
    """N: number of leaves of the tree
    Uses the recursion: Sn = S(j) + S(n-j) + n for each internal node.
    """
    minmax = [sackin_minmax(1)]
    distribs = [defaultdict(int)]
    distribs[0][0] = 1  # 1 ways to make a tree of one node (sackin value of zero) 
    for n in range(2, N+1):
        m,M = sackin_minmax(n)
        logger.debug('# subtree of size %d; minmax=[%d, %d]', n, m, M)
        minmax.append((m,M))
        distrib = defaultdict(int)  #np.zeros(M - m + 1)
        distrib[M] = factorial(n) // 2  # Caterpillar tree. The single cherry can rotate.

        # For each left subtree of size k
        # if n is even: go up to n/2 (included)
        # if n is odd: go up to (n-1)/2 (included)
        for k in range(1, n//2 + 1):
            mk, Mk = minmax[k - 1]  # -1 because indexing of minmax starts at 0
            mn_k, Mn_k = minmax[n-k - 1]  # m(n-k), M(n-k)

            # Number of ways to put n leaves in 2 groups of sizes (k, n-k)
            ck = binom(n,k)  # Warning: the binom function may return a float approximation.
            if n%2 == 0 and k == n//2:
                # 2 splits are identical by rotation of the 2 subtrees.
                ck /= 2
            logger.debug('  - split k=%d: %g way(s)', k, ck)
            # For each possible value of the Sackin index of the tree of size n
            for v in range(m, M):
                # All possible ways to split the value between the 2 subtrees
                # Remember the recursive relation:
                pk = 0
                # Contraints:
                # 1)  mk         ≤ vk     ≤ Mk
                # 2)  mn_k       ≤ v-n-vk ≤ Mn_k
                #     v-n - Mn_k ≤     vk ≤ v-n - mn_k
                for vk in range(max(0, mk, v-n - Mn_k), min(v, Mk, v-n - mn_k) + 1):
                #for vk in range(max(1, mk), min(v, Mk) + 1):
                #for vk in range(1, v-n):
                    logger.debug('    value(n) = %d = n + value(k) + value(n-k) = %d + %d + %d',
                                 v, n, vk, v-n-vk)
                    logger.debug('    ways += %d × %d', distribs[k - 1][vk],
                                 distribs[n-k - 1][v-n - vk])
                    pk += distribs[k - 1][vk] * distribs[n-k - 1][v-n - vk]

                distrib[v] += ck*pk

        distribs.append(distrib)
    return distribs  # distrib


def test_sackin_minmax():
    assert sackin_minmax(1) == (0,0)
    assert sackin_minmax(2) == (2,2)
    assert sackin_minmax(3) == (5,5)
    assert sackin_minmax(4) == (8,9)
    assert sackin_minmax(5) == (12, 14)
    assert sackin_minmax(6) == (16, 20)
    #assert sackin_minmax(7)

# From Felsenstein, _Inferring Phylogenies_ (Sinauer and associates), Ch.3
n_rooted_labelled_topo = {1:1, 2:1, 3:3, 4:15, 5:105, 6:945, 7:10395,
        8:135135, 9:2027025, 10:34459425, 11:654729075, 12:13749310575,
        13:316234143225, 14:7905853580625, 15:213458046676875,
        16:6190283353629375, 17:191898783962510625, 18:6332659870762850625,
        19:221643095476699771875, 20:8200794532637891559375}


def test_sackin_distrib():

    # The sum of p should equal the number of possible labelled topologies
    d1, d2, d3, d4, d5, d6 = sackin_distrib(6)

    assert sum(d1.values()) == 1
    assert d1[0] == 1

    assert sum(d2.values()) == 1
    assert d2[2] == 1

    assert sum(d3.values()) == 3
    assert d3[5] == 3

    assert sum(d4.values()) == 3 + 4*3*2//2  # Balanced (4 choose 2)/2 + Caterpillar
    assert d4[8] == 3
    assert d4[9] == 12

    assert sum(d5.values()) == 105
    assert d5[14] == 60
    assert d5[12] == 30
    assert d5[13] == 15

    assert sum(d6.values()) == n_rooted_labelled_topo[6]  # 945
    assert d6[16] == 90+45  # 2 corresponding topologies
    assert d6[17] == 180
    assert d6[18] == 180
    assert d6[19] == 90
    assert d6[20] == factorial(6)/2  # 360


def n_cherries(tree):
    return len([n for n in tree.traverse() if all(ch.is_leaf() for ch in n.children)])


def colless(tree):
    s = 0
    for n in tree.traverse():
        children = n.children
        if len(children) < 2:
            continue
        elif len(children) > 2:
            raise ValueError("Tree should be dichotomic")
        ch0, ch1 = children
        s += abs(len(ch0) - len(ch1))


def normed_colless(tree):
    """Equals 1 for a caterpillar tree"""
    n = len(tree)
    return colless(tree) * 2. / ((n-1)*(n-2))


def mir(tree):
    """Mir et al 2013: Sum, for all pairs of leaves, of the depth of the
    lowest common ancestor."""
    pass