#!/usr/bin/env python3

"""Configure tests for Pytest"""


# Add a command-line option
def pytest_addoption(parser):
    parser.addoption('--gene_info', help='Template for gene_info tabular files')
    parser.addoption('--cgene', default=1, help='Column for the gene name [1]')
    parser.addoption('--cprot', default=2, help='Column for the protein name [2]')
    parser.addoption('--ensembl_version', action='append', default=[])
    parser.addoption('--forestfile', help='example: ~/GENOMICUS%%d/tree.1.ensembl.bz2')

# Example config:
# --gene_info '~/ws7/DUPLI_data93/gene_info/%s_fromtree.tsv'
# --ensembl_version 93
# --forestfile '~/GENOMICUS%d/tree.1.ensembl.bz2'

# This will configure tests requiring a gene_info value
def pytest_generate_tests(metafunc):
    if 'gene_info' in metafunc.fixturenames:
        gene_info_param = metafunc.config.getoption('gene_info')
        metafunc.parametrize('gene_info', [gene_info_param] if gene_info_param else [])
        if gene_info_param:
            metafunc.parametrize('ensembl_version', [int(e) for e in metafunc.config.getoption('ensembl_version')])
            if 'cprot' in metafunc.fixturenames:
                metafunc.parametrize('cprot', [metafunc.config.getoption('cprot')])
            if 'cgene' in metafunc.fixturenames:
                metafunc.parametrize('cgene', [metafunc.config.getoption('cgene')])
    if 'forestfile' in metafunc.fixturenames:
        forestfile_param = metafunc.config.getoption('forestfile')
        metafunc.parametrize('forestfile', [forestfile_param] if forestfile_param else [])
        metafunc.parametrize('ensembl_version', [int(e) for e in metafunc.config.getoption('ensembl_version')])

