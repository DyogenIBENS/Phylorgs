#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from seqtools.compo_freq import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

import pytest


class Test_get_seq_counts:
    def test_accept_unaligned(self):
        records = [SeqRecord(Seq('A')), SeqRecord(Seq('AA'))]
        length, *counts = get_seq_counts(records)
        assert len(length) == 2

    @pytest.mark.parametrize('input_type', [list, MultipleSeqAlignment])
    @pytest.mark.parametrize('is_aa', [False, True], ids=['DNA', 'protein'])
    def test_invalid_residues_valueerror(self, is_aa, input_type):
        records = input_type([SeqRecord(Seq('A$'))])
        with pytest.raises(ValueError):
            get_seq_counts(records, is_aa)

    @pytest.mark.parametrize('input_type', [list, MultipleSeqAlignment])
    @pytest.mark.parametrize('is_aa', [False, True], ids=['DNA', 'protein'])
    def test_valid_residues(self, is_aa, input_type):
        alphabet = AA if is_aa else NUCLEOTIDES
        records = input_type([SeqRecord(Seq(alphabet))])
        length, seq_resid, seq_gap, seq_N, *_ = get_seq_counts(records, is_aa)
        # Check dimensions
        assert seq_resid.shape == (len(alphabet), 1)
        assert seq_gap.shape == (len(records),)
        assert seq_N.shape == (len(records),)
        # Check counts
        assert np.all(seq_resid == 1)
        assert seq_gap[0] == 0
        assert seq_N[0] == 0

    @pytest.mark.parametrize('input_type', [list, MultipleSeqAlignment])
    @pytest.mark.parametrize('is_aa', [False, True], ids=['DNA', 'protein'])
    def test_gaps(self, is_aa, input_type):
        records = input_type([SeqRecord(Seq(GAPS))])
        length, seq_resid, seq_gap, seq_N, *_ = get_seq_counts(records, is_aa)
        assert np.all(seq_resid == 0)
        assert seq_gap[0] == len(GAPS)
        assert seq_N[0] == 0


class Test_get_seq_freqs:
    @pytest.mark.parametrize('extra_counts', [None, np.zeros(1), np.ones(1)], ids=['None_extra', 'zeros_extra', 'ones_extra'])
    @pytest.mark.parametrize('length', [5, [5]], ids=['scalar', 'list'])
    def test_one_sequence_accept_input(self, length, extra_counts):
        oneseq_counts = np.ones((4, 1)), np.zeros(1), np.zeros(1)
        oneseq_counts += (extra_counts,) * 2
        get_seq_freqs(length, *oneseq_counts)


class Test_get_al_compo_summary:
    @pytest.mark.parametrize('extra_counts', [None, np.zeros(1), np.ones(1)], ids=['None_extra', 'zeros_extra', 'ones_extra'])
    @pytest.mark.parametrize('length', [5, [5]], ids=['scalar', 'list'])
    def test_one_sequence_accept_input(self, length, extra_counts):
        oneseq_counts = np.ones((4, 1)), np.zeros(1), np.zeros(1)
        oneseq_counts += (extra_counts,) * 2
        get_al_compo_summary(length, oneseq_counts, get_seq_freqs(length, *oneseq_counts))


MSA = MultipleSeqAlignment([SeqRecord(Seq(GAPS*2 + NUCLEOTIDES)),
                            SeqRecord(Seq(GAPS + NUCLEOTIDES[:2] + NUCLEOTIDES))])

unaligned = [SeqRecord(Seq(NUCLEOTIDES)),
             SeqRecord(Seq(NUCLEOTIDES*2))]

# Valid inputs
valid_input_ids, valid_inputs = zip(*[('alignment_bysite', (MSA, False)),
                                      ('alignment_byseq',  (MSA, True)),
                                      ('list_aligned_bysite', (list(MSA), False)),
                                      ('list_aligned_byseq', (list(MSA), True)),
                                      ('list_unaligned_byseq', (unaligned, True))])

print(valid_input_ids)

class Test_make_al_compo:
    @pytest.mark.parametrize('records,byseq', valid_inputs, ids=valid_input_ids)
    def test_correct_output_shape(self, records, byseq):
        all_stats = []
        for is_aa in (False, True):
            stat_names, stats = make_al_compo(records, byseq, is_aa)
            if byseq:
                assert len(stat_names) == len(records)
            assert len(stats) == len(stat_names)
            all_stats.append(stats)

        stats_nucl, stats_aa = all_stats
        assert len(stats_nucl[0]) - len(NUCLEOTIDES) > len(stats_aa[0]) - len(AA), 'Expect more special features counted for nucleotides (GC, CpG, stops)'

