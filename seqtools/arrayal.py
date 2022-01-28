#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Module to load Multiple Sequence Alignments into numpy arrays.

The 'al2int' function encodes residues as integers.
"""


import re
import numpy as np
from itertools import product, combinations
from seqtools.IUPAC import gaps, unknown, nucleotides, ambiguous


STOPS = ['TAA', 'TAG', 'TGA']
CODONS = [''.join(codon) for codon in product(*[nucleotides]*3)
          if ''.join(codon) not in STOPS]
NACODON = '---'
#NCODONS = 
CODON2INT = {codon:i for i,codon in enumerate([NACODON] + CODONS + STOPS)}
#NUCL2INT  = {'-': 0, 'A': 1, 'C': 2, 'G': 3, 'T': 4}
NUCL2INT = {symbol: i for i, symbol in enumerate(gaps[0] + nucleotides)}


def make_unif_codon_dist():
    """Return a uniform codon distance matrix:
    each different nucleotide counts as one difference.
    Stop codons have a NaN distance.
    Invalid codons (containing N nucleotides) as well.
    """
    all_values = [NACODON] + CODONS
    size = len(all_values) + len(STOPS) + 1  # +1 because of invalid codons e.g 'NNN'
    uniform_codon_dist = 1 - np.diagflat(np.ones(size))
    for i,j in combinations(range(size - len(STOPS) - 1), 2):
        # sum of nucleotide differences
        cod0, cod1 = all_values[i], all_values[j]
        pair_dist = (np.array(list(cod0)) != np.array(list(cod1))).sum()
        uniform_codon_dist[i, j] = uniform_codon_dist[j, i] = pair_dist
    
    # Stops and invalid codons (e.g codons containing N's) get a NaN distance.
    uniform_codon_dist[:, -4:] = np.NaN
    uniform_codon_dist[-4:, :] = np.NaN

    return uniform_codon_dist

UNIF_CODON_DIST = make_unif_codon_dist()


#def read_codon_matrix(filename):
#    # ~/install/jprime-extra/arvecodon.txt


def al2list(align_col):
    """Convert an alignment object to a list of strings"""
    return [str(record.seq) for record in align_col]


def iter_strseq(align):
    """Iterate over an alignment object: yield sequences as strings"""
    for record in align:
        yield str(record.seq)


def al2array(align, nucl=False):
    """"""
    nseq = len(align)
    al_len = align.get_alignment_length()
    npos = al_len if nucl else al_len // 3
    step = 1      if nucl else 3
    #alarray = np.array([['']*npos]*len(align))

    return np.array([[seq[(j*step):((j+1)*step)] for j in range(npos)]
                                                 for seq in iter_strseq(align)])


def al2int(align, nucl=False, allow_N=False):
    """Converts a matrix of categorical values to integers"""
    alarray = al2array(align, nucl)
    converter_dict = NUCL2INT if nucl else CODON2INT
    return category2int(alarray, converter_dict, allow_N)


def category2int(array, converter_dict, allow_N=False):
    """Convert an array of categorical values using a dictionary"""
    if not allow_N:
        ufunc = converter_dict.__getitem__
    else:
        #ufunc = lambda residue: converter_dict.get(residue, np.NaN)
        Ncodon_regex = re.compile(r'[' + unknown + nucleotides + ''.join(ambiguous) + ']+$')
        Ncodon_int = max(converter_dict.values()) + 1
        
        def ufunc(residue):
            try:
                return converter_dict[residue]
            except KeyError:
                if Ncodon_regex.match(residue):
                    return Ncodon_int
                else:
                    raise

    return np.vectorize(ufunc)(array)


def count_matrix(vint, minlength=66):
    """column-wise matrix of counts of integers"""
    #assert np.issubdtype(vint.dtype, int  # 
    assert vint.dtype in (int, np.int64, np.int32, np.int16)
    return np.stack([np.bincount(vint[:,i], minlength=minlength) \
                        for i in range(vint.shape[-1])], axis=-1)

#count_matrix = np.vectorize(lambda array_1D: np.bincount(array_1D, minlength=minlength)
# Nope because this processes columns


def freq_matrix(vint, minlength=66):
    """Convert matrix of integers to frequencies (per column)"""
    return count_matrix(vint, minlength).astype(float) / np.alen(vint)


def presence_matrix(vint, minlength=66):
    """Convert a sequence of integers into a boolean matrix of presence of each value.
    Zero encodes gaps. >minlength encodes ambiguous/invalid residues."""

    assert len(vint.shape) == 2 and vint.shape[0] == 1, \
            "Can't input more than one sequence " + str(vint.shape)
    count_mat = count_matrix(vint, minlength)
    # Convert gaps (integer zero) and 'N' (integer minlength-1) to unknown: True everywhere.
    presence = count_mat[1:(minlength-1), :].astype(bool)
    # Broadcast row of gaps on all rows of `presence`.
    presence |= count_mat[0, :].astype(bool)
    # Broadcast row of 'N' on all rows of `presence`.
    presence |= count_mat[(minlength-1), :].astype(bool)

    return presence


    ### WTF
def presence_matrice(vint, minlength=66):
    assert len(vint.shape) == 2 and vint.shape[0] == 1, \
            "Can't input more than one sequence " + str(vint.shape)
    count_mat = count_matrix(vint, minlength)
    presence = count_mat.astype(bool)
    # Convert gaps (integer zero) and 'N' (integer minlength-1) to unknown: True everywhere.
    presence[0, :] = True
    return presence


def np_proportions(values):
    # Convert values to integers
    vunique, vcounts = np.unique(values, return_counts=True)
    return vunique, vcounts.astype(float) / len(values)


def np_entropy_subfunc(value_unique, value_prop, na=None):
    if na is None:
        na_pos = np.full(len(value_unique), True)
    else:
        na_pos = value_unique == na

    return - (value_prop * np.log2(value_prop))[~na_pos].sum()


def np_entropy(values, na=None):
    value_unique, value_prop = np_proportions(values)
    return np_entropy_subfunc(value_unique, value_prop, na)


def freqs2entropy(freqs):
    """Compute the entropy of columns from the frequency matrix.

    Return a vector of entropy of the same length as the number of columns.
    """
    freqs = freqs.copy()
    freqs[freqs == 0] = 1
    return - (freqs * np.log2(freqs)).sum(axis=0)


def get_position_stats(align, nucl=False, allow_N=False):
    """Compute gap proportion and entropy value for each position of the alignment"""
    al_len = align.get_alignment_length()
    if al_len % 3 and not nucl:
        logger.error("Not a codon alignment!")

    gap = '-'     if nucl else '---'
    npos = al_len if nucl else al_len // 3
    step = 1      if nucl else 3

    alint = al2int(align, nucl, allow_N)
    alint = alint[:, alint.sum(axis=0) > 0]

    #gap_prop = np.zeros(npos)
    #al_entropy = np.full(npos, np.NaN)

    #for i in range(npos):
    #    values = al2list(align[:, (step*i):(step*(i+1))])
    #    value_unique, value_prop = np_proportions(values)
    #    #gap_prop[i]   = value_prop[np.argmax(value_unique == gap)] or 0
    #    al_entropy[i] = np_entropy_subfunc(value_unique, value_prop, gap)

    #gap_prop = (alint == 0).sum(axis=0) / alint.shape[0]
    al_freqs = freq_matrix(alint, minlength=(6 if nucl else 66))
    gap_prop = al_freqs[0,:]
    # Convert zeros to ones so that x * log(x) = 0
    al_freqs = al_freqs[1:,:]
    # scipy.stats.entropy
    al_freqs[al_freqs == 0] = 1
    al_entropy = - (al_freqs * np.log2(al_freqs)).sum(axis=0)

    logger.info(alint.shape)
    return gap_prop, al_entropy, alint


def reorder_al(align, records=None, record_order=None):
    """
    param: `align`: Bio.Align object
    param: `record_order`: sequence of integers
    param: `records`: sequence of record names
    """
    # Mutually exclusive args
    if not (record_order is None) ^ (records is None):
        raise ValueError('`record_order` and `records` are mutually exclusive. Provide exactly one.')

    if records is not None:
        recnames = [rec.name for rec in align]
        record_order = [recnames.index(recordname) for recordname in records]

    aligntype = type(align)
    return aligntype([align[r] for r in record_order])


