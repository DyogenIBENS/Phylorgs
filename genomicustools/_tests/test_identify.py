#!/usr/bin/env python
# -*- coding: utf-8 -*-


import pytest
from genomicustools.identify import *
from glob import glob
from bz2 import BZ2File
from LibsDyogen import myProteinTree


def myopen(filename, *args, **kwargs):
    if filename.endswith('.bz2'):
        return BZ2File(filename, *args, **kwargs)
    else:
        return open(filename, *args, **kwargs)


def verify_gene_info(gene_info, ensembl_version):
    gene_info = op.expanduser(gene_info)
    splist_module = set(PROT2SP[ensembl_version].values())
    sp_from_filename = re.compile(gene_info.replace('%s', '([A-Za-z0-9.]+)'))
    gene_info_files = glob(gene_info.replace('%s', '*'))
    assert gene_info_files, "No files found with search path %s" % gene_info.replace('%s', '*')
    splist_file = set()
    for fn in gene_info_files:
        try:
            splist_file.add(sp_from_filename.match(fn).group(1).replace('.', ' '))
        except AttributeError as err:
            err.args += ('%s VS %r' % (sp_from_filename, fn),)
            raise

    if not splist_module == splist_file:
        raise FileNotFoundError('Differences in lists of species:\n' +
                             'module (%d): %s\n' % (len(splist_module),
                                                    splist_module - splist_file) +
                             'files  (%d): %s' % (len(splist_file),
                                                  splist_file - splist_module))
    return gene_info, splist_file


class Test_convert_prot2species:

    @pytest.mark.parametrize('ensembl_version', [85,93])
    def test_reject_wrong_strings_fails(self, ensembl_version):
        for wrong in ('xululul', '0000000', 'ENSXXXP', 'ENSG000'):
            with pytest.raises(KeyError):
                predicted_sp = convert_prot2species(wrong, ensembl_version)

    @pytest.mark.parametrize('ensembl_version', [85,93])
    def test_reject_wrong_strings_defaults(self, ensembl_version):
        for wrong in ('xululul', '0000000', 'ENSXXXP', 'ENSG000'):
            predicted_sp = convert_prot2species(wrong, ensembl_version, False)
            assert predicted_sp is False, "%r predicted %r" % (wrong, predicted_sp)

    def test_every_modernID(self, ensembl_version, gene_info, cprot):
        """
        test the 'convert_prot2species' function for every modernID using
        external gene_info tables
        """
        gene_info, splist_file = verify_gene_info(gene_info, ensembl_version)
        # Check every valid prot ID in the given files
        for sp in splist_file:
            filename = gene_info % sp.replace(' ', '.')
            print("Checking %s in %r" % (sp, op.basename(filename)), file=stderr)
            # Check that each species protein return the correct species.
            with myopen(filename) as IN:
                for line in IN:
                    prot = line.rstrip('\r\n').split('\t')[cprot]
                    try:
                        predicted_sp = convert_prot2species(prot, ensembl_version)
                        assert sp == predicted_sp, "%s: %r ≠ %r" % (prot, sp, predicted_sp)
                    except KeyError as err:
                        err.args = err.args[:-1] + \
                                   (err.args[-1] + ' '.join((sp, prot, "Not found")),)
                        raise


class Test_convert_gene2species:
    @pytest.mark.parametrize('ensembl_version', [85,93])
    def test_reject_wrong_strings(self, ensembl_version):
        for wrong in ('xululul', '0000000', 'ENSXXXG', 'ENSP000'):
            with pytest.raises(KeyError):
                predicted_sp = convert_gene2species(wrong, ensembl_version)

    def test_every_modernID(self, ensembl_version, gene_info, cgene):
        """test the above function for every modernID"""

        gene_info, splist_file = verify_gene_info(gene_info, ensembl_version)

        for sp in splist_file:
            filename = gene_info % sp.replace(' ', '.')
            print("Checking %s in %r" % (sp, op.basename(filename)), file=stderr)
            # Check that each species gene return the correct species.
            with open(filename) as IN:
                for line in IN:
                    gene = line.rstrip('\r\n').split('\t')[cgene]
                    try:
                        predicted_sp = convert_gene2species(gene, ensembl_version)
                    except KeyError as err:
                        err.args = err.args[:-1] + \
                                   (err.args[-1] + ' '.join((sp, gene, "Not found")),)
                        raise
                    assert sp == predicted_sp, "%s: %r ≠ %r" % (gene, sp, predicted_sp)


class Test_convert2species:
    @pytest.mark.parametrize('ensembl_version', [85,93])
    def test_reject_wrong_strings(self, ensembl_version):
        for wrong in ('xululul', '0000000', 'ENSXXXP', 'ENSG000'):
            predicted_sp = convert_prot2species(wrong, ensembl_version, False)
            assert predicted_sp is False, "%r predicted %r" % (wrong, predicted_sp)

    def test_every_modernID(self, forestfile, ensembl_version):
        """test the 'convert_prot2species' function for every modernID"""
        expected_species = set(GENE2SP[ensembl_version].values())
        expected_species_p = set(PROT2SP[ensembl_version].values())
        assert expected_species == expected_species_p

        for tree in myProteinTree.loadTree(op.expanduser(forestfile % ensembl_version)):
            for tip in (set(tree.info) - set(tree.data)):
                tipinfo = tree.info[tip]
                sp = tipinfo['taxon_name']
                gene = tipinfo['gene_name']
                prot = tipinfo['protein_name']

                assert sp in expected_species, 'Unexpected species %r' % sp

                try:
                    predicted_sp = convert_gene2species(gene, ensembl_version)
                except KeyError as err:
                    err.args = err.args[:-1] + \
                               (err.args[-1] + ' '.join((sp, gene, "Not found")),)
                    raise
                assert sp == predicted_sp, "%s: %r ≠ %r" % (gene, sp, predicted_sp)

                try:
                    predicted_sp = convert_prot2species(prot, ensembl_version, default=None)
                except KeyError as err:
                    err.args = err.args[:-1] + \
                               (err.args[-1] + ' '.join((sp, prot, "Not found")),)
                    raise
                assert sp == predicted_sp, "%s: %r ≠ %r" % (prot, sp, predicted_sp)

