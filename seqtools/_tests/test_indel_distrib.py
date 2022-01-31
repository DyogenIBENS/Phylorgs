#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from seqtools.indel_distrib import *
#logging.basicConfig(format=logging.BASIC_FORMAT)
logger.setLevel(logging.DEBUG)
import pytest
from collections import Counter
from io import StringIO

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import ete3

from dendro.any import ete3 as ete3_methods

get_children = ete3_methods.get_children
get_label = ete3_methods.get_label
get_root = ete3_methods.get_root


def MSA(*seqs):
    return MultipleSeqAlignment([
        SeqRecord(Seq(seq), id=('%s' % chr(97+i)), name=('%s' % chr(97+i)))
        for i,seq in enumerate(seqs)])

class Test_map_indels_codons:
    def test_no_indel(self):
        al = MSA('ATGTTT')
        gaps = map_gaps(al, codon=True)
        assert len(gaps['a']) == 0
    def test_indel_middle(self):
        al = MSA('ATG------TTT')
        gaps = map_gaps(al, codon=True)
        assert len(gaps['a']) == 1
        assert gaps['a'][0] == (1,2)
    def test_indel_length1(self):
        al = MSA('ATG---TTT')
        gaps = map_gaps(al, codon=True)
        assert len(gaps['a']) == 1
        assert gaps['a'][0] == (1,1)
    def test_indel_beginning(self):
        al = MSA('---TTT')
        gaps = map_gaps(al, codon=True)
        assert len(gaps['a']) == 1
        assert gaps['a'][0] == (0,0)
    def test_indel_end(self):
        al = MSA('TTT---')
        gaps = map_gaps(al, codon=True)
        assert len(gaps['a']) == 1
        assert gaps['a'][0] == (1,1)

class Test_index_gaps:
    def test_no_overlap(self):
        g = combine_gaps([(1,2)], [(3,5)])
        assert set(g) == set(( (1,2), (3,5) ))
    def test_identical(self):
        g = combine_gaps([(1,2)], [(1,2)])
        assert set(g) == set(( (1,2), ))
        assert g == [(1,2)]
    def test_one_overlap_same_end(self):
        g = combine_gaps([(2,3)], [(1,3)])
        assert set(g) == set(( (1,1), (1,3), (2,3) ))
        assert len(g) == 3, "Wrong length g=%s" % g
    def test_one_overlap_same_start(self):
        g = combine_gaps([(1,3)], [(1,2)])
        assert set(g) == set(( (1,2), (3,3), (1,3) ))
        assert len(g) == 3, "Wrong length g=%s" % g
    def test_one_overlap_hanging(self):
        g = combine_gaps([(1,2)], [(2,3)])
        assert set(g) == set(( (1,2), (1,1), (2,2), (2,3), (3,3) ))
        assert len(g) == 5, "Wrong length g=%s" % g
    def test_one_overlap_contained(self):
        g = combine_gaps([(1,5)], [(2,3)])
        assert set(g) == set(( (1,1), (1,5), (2,3), (4,5) ))
        assert len(g) == 4, "Wrong length g=%s" % g
    def test_two_overlaps_same_extremities(self):
        g = combine_gaps([(1,5)], [(1,2), (4,5)])
        assert set(g) == set(( (1,2), (3,5), (1,5), (4,5), (1,3), (3,3) ))

    def test_two_overlaps_noduplicates(self):
        g = combine_gaps([(1,5)], [(1,2), (4,5)])
        print(len(g), len(set(g)))
        counter = Counter()
        for gap in g:
            counter[gap] += 1
        for k,v in counter.items():
            assert v == 1, "%s is found %d times" % (k, v)


class Test_reindex:
    def test_reindex_extras(self):
        states = np.array([[0,1]])
        newstates = reindex_states(states, [(4,4)], [(2,2), (4,4), (6,6)])
        # Only the original gap should have the [0,1] state
        assert (newstates == np.array([[1,0], [0,1], [1,0]])).all() 
    def test_reindex_included(self):
        states = np.array([[0,1]])
        newstates = reindex_states(states, [(3,4)], [(3,3), (4,4), (3,4)])
        # included gaps should inherit the state of the original overlapping gap.
        assert (newstates == np.array([[0,1], [0,1], [0,1]])).all() 


class Test_anc_gaps:
    def test_derived1(self):
        al = MSA('ATG---TTT',
                 'ATGCCCTTT',
                 'ATGCCCTTT')
        tree = ete3.Tree('((a,b)x,c)r;', format=1)

        gap_coords, states = rootwards_gap_states(al, tree, get_children, get_label, root=tree, codon=True)
        assert gap_coords['r'] == [(1,1)]

        # We should just have made the unions of states, if they don't intersect:
        assert ( states['x'] == np.array([[1,1]]) ).all()
        assert ( states['r'] == np.array([[1,0]]) ).all()

        newstates, rootcoords, branch_inser, branch_del, branch_evts = leafwards_branch_events(tree, get_children, get_label, tree, gap_coords, states)

        assert rootcoords == [(1,1)]
        assert (newstates['x'] == np.array([[1,0]])).all()

        # Also check that leaf states are identical to the observed states
        assert (newstates['a'] == np.array([[0,1]])).all()
        assert (newstates['b'] == np.array([[1,0]])).all()

    def test_derived2(self):
        al = MSA('ATG---TTT',
                 'ATG---TTT',
                 'ATGCCCTTT')
        tree = ete3.Tree('((a,b)x,c)r;', format=1)
        gap_coords, states = rootwards_gap_states(al, tree, get_children, get_label, root=tree, codon=True)

        assert gap_coords['x'] == [(1,1)]
        assert ( states['x'] == np.array([[0,1]]) ).all()

        newstates, rootcoords, branch_inser, branch_del, branch_evts = leafwards_branch_events(tree, get_children, get_label, tree, gap_coords, states)

        assert rootcoords == [(1,1)]
        assert (newstates['x'] == np.array([[0,1]])).all()
        #assert (newstates['r'] == np.array([[1,1]])).all()

    def test_1insertion(self):
        al = MSA('ATG---TTT',
                 'ATG---TTT',
                 'ATGCCCTTT',
                 'ATG---TTT',
                 'ATG---TTT')
        #tree = ete3.Tree('(((a,b)x,c)y,d)r;', format=1)
        tree = ete3.Tree('((((a,b)x,c)y,d)z,e)r;', format=1)

        gap_coords, states = rootwards_gap_states(al, tree, get_children, get_label, root=tree, codon=True)

        assert gap_coords['x'] == [(1,1)]
        assert gap_coords['y'] == [(1,1)]

        assert (states['x'] == np.array([[0,1]])).all()
        assert (states['y'] == np.array([[1,1]])).all()

        newstates, rootcoords, branch_inser, branch_del, branch_evts = leafwards_branch_events(tree, get_children, get_label, tree, gap_coords, states)

        assert rootcoords == [(1,1)]
        assert (newstates['x'] == np.array([[0,1]])).all()
        assert (newstates['y'] == np.array([[0,1]])).all()  # 'y' was a gap
        assert (newstates['z'] == np.array([[0,1]])).all()  # 'z' was a gap
        assert (newstates['r'] == np.array([[0,1]])).all()  # root was a gap (should be undertermined)

    def test_1extension(self):
        al = MSA('ATG------TTT',
                 'ATG---CCCTTT',
                 'ATG---CCCTTT')
        tree = ete3.Tree('((a,b)x,c)r;', format=1)

        gap_coords, states = rootwards_gap_states(al, tree, get_children, get_label, root=tree, codon=True)
        assert gap_coords['r'] == [(1,1), (1,2), (2,2)]

        assert ( states['x'] == np.array([[0,1], [1,1], [1,1]]) ).all()
        assert ( states['r'] == np.array([[0,1], [1,0], [1,0]]) ).all() # because it intersected states['x'] with states['c']

        newstates, rootcoords, branch_inser, branch_del, branch_evts = leafwards_branch_events(tree, get_children, get_label, tree, gap_coords, states)

        assert rootcoords == [(1,1), (1,2), (2,2)]
        assert (newstates['x'] == np.array([[0,1], [1,0], [1,0]])).all()
        assert (newstates['r'] == np.array([[0,1], [1,0], [1,0]])).all()

        # Check that leaf states are restricted to the observed states, 
        # including any subgap.
        assert (newstates['a'] == np.array([[0,1], [0,1], [0,1]])).all()
        assert (newstates['b'] == np.array([[0,1], [1,0], [1,0]])).all()

        #assert score == 1 # As implemented, there is an increment for each gap in
        # the index, not taking overlaps into account.

    def test_one_gap_distribs(self):
        al = MSA('ATG---TTT',
                 'ATGCCCTTT',
                 'ATGCCCTTT')
        tree = ete3.Tree('((a,b)x,c)r;', format=1)

        gap_coords, states = rootwards_gap_states(al, tree, get_children, get_label, root=tree, codon=True)
        newstates, rootcoords, branch_inser, branch_del, branch_evts = leafwards_branch_events(tree, get_children, get_label, tree, gap_coords, states)

        # As in test_derived1
        assert (newstates['a'] == np.array([[0,1]])).all()
        assert (newstates['b'] == np.array([[1,0]])).all()
        assert (newstates['x'] == np.array([[1,0]])).all()
        assert (~(np.logical_and(newstates['x'], newstates['a']).any(axis=1))).tolist() == [True]
        print(branch_del)
        assert branch_del[('x', 'a')] == [(1,1)]
        assert branch_inser[('x', 'a')] == []

        del_lengths = dict(branch_len_distribs(branch_del))
        assert del_lengths[('x', 'a')].tolist() == [1]
        assert del_lengths[('x', 'b')].shape == (0,)
        assert del_lengths[('r', 'x')].shape == (0,)
        assert del_lengths[('r', 'c')].shape == (0,)
        
        inser_lengths = dict(branch_len_distribs(branch_inser))
        for br,distrib in inser_lengths.items():
            assert distrib.shape == (0,), "distrib %s : %s" % (br, distrib)

