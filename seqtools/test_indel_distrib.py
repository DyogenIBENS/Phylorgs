#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from seqtools.indel_distrib import *
#logging.basicConfig(format=logging.BASIC_FORMAT)
logger.setLevel(logging.DEBUG)
import pytest
from collections import Counter

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
    def test_one_overlap_same_end(self):
        g = combine_gaps([(2,3)], [(1,3)])
        assert set(g) == set(( (1,1), (1,3), (2,3) ))
    def test_one_overlap_same_start(self):
        g = combine_gaps([(1,3)], [(1,2)])
        assert set(g) == set(( (1,2), (3,3), (1,3) ))
    def test_one_overlap_hanging(self):
        g = combine_gaps([(1,2)], [(2,3)])
        assert set(g) == set(( (1,2), (1,1), (2,2), (2,3), (3,3) ))
    def test_one_overlap_contained(self):
        g = combine_gaps([(1,5)], [(2,3)])
        assert set(g) == set(( (1,1), (1,5), (2,3), (4,5) ))
    def test_two_overlaps_same_extremities(self):
        g = combine_gaps([(1,5)], [(1,2), (4,5)])
        assert set(g) == set(( (1,2), (3,5), (1,5), (4,5), (1,3), (3,3) ))

class Test_anc_gaps:
    def test_anc_gaps_derived1(self):
        al = MSA('ATG---TTT',
                 'ATGCCCTTT',
                 'ATGCCCTTT')
        tree = ete3.Tree('((a,b)x,c)r;', format=1)

        gap_coords, states = rootwards_gap_states(al, tree, get_children, get_label, root=tree, codon=True)
        assert gap_coords['r'] == [(1,1)]

        assert ( states['x'] == np.array([[1,0]]) ).all()
        assert ( states['r'] == np.array([[1,0]]) ).all()

    def test_anc_gaps_derived2(self):
        al = MSA('ATG---TTT',
                 'ATG---TTT',
                 'ATGCCCTTT')
        tree = ete3.Tree('((a,b)x,c)r;', format=1)
        gap_coords, states = rootwards_gap_states(al, tree, get_children, get_label, root=tree, codon=True)

        #assert gap_coords['x'] == [(1,1)]
        assert ( states['x'] == np.array([[0,1]]) ).all()

#    def test_anc_1insertion(self):
#        al = MSA('ATG---TTT',
#                 'ATG---TTT',
#                 'ATGCCCTTT',
#                 'ATG---TTT',
#                 'ATG---TTT')
#        #tree = ete3.Tree('(((a,b)x,c)y,d)r;', format=1)
#        tree = ete3.Tree('((((a,b)x,c)y,d)z,e)r;', format=1)
#
#        anc_gaps, uncertain = get_anc_gaps(al, tree, get_children, get_label, root=tree, codon=True)
#
#        print(anc_gaps)
#        print(uncertain)
#        assert (1,1) in anc_gaps['x']
#        assert (1,1) in anc_gaps['y']

    def test_two_overlaps_noduplicates(self):
        g = combine_gaps([(1,5)], [(1,2), (4,5)])
        print(len(g), len(set(g)))
        counter = Counter()
        for gap in g:
            counter[gap] += 1
        for k,v in counter.items():
            assert v == 1, "%s is found %d times" % (k, v)
