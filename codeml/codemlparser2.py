#!/usr/bin/env python3

from __future__ import print_function
from archiparse import ParseUnit, ParseBloc, floatlist, strlist, OrderedDict


mlc_parser = ParseBloc('mlc',
                 [ParseUnit('version',
                      r'^CODONML \(in paml version ([0-9a-z.-]+), ([A-Z][a-z]+ [0-9]+)\)',
                      keys=['number', 'date']),
                  ParseUnit('model',
                      r'^Model: (.*),\s*$'),
                  ParseUnit('Codon_freq_model',
                      r'^Codon frequency model: (.*)$'),
                  ParseUnit('nsls',
                      r'^ns =\s+([0-9]+)\s+ls =\s+([0-9]+)$',
                      keys=['ns', 'ls'], types=[int, int]),
                  ParseUnit('Codon usage in sequences',
                      r'^Codon usage in sequences$'),
                 ParseUnit('Codon position x base (3x4) table for each sequence.',
                     r'^Codon position x base \(3x4\) table for each sequence\.$'),
                 ParseUnit('Sums of codon usage counts',
                     r'^Sums of codon usage counts$'),
                 ParseUnit('Codon position x base (3x4) table, overall',
                     r'^Codon position x base \(3x4\) table, overall$'),
                 ParseUnit('Codon frequencies under model, for use in evolver',
                     r'^Codon frequencies under model, for use in evolver'),
                 ParseBloc('Nei & Gojobori',
                           units=[
                            ParseUnit('Description',
r'''^Nei & Gojobori 1986\. dN/dS \(dN, dS\)$'''),
#\(Pairwise deletion\)
#\(Note: This matrix is not used in later ML\. analysis\.
#Use runmode = -2 for ML pairwise comparison\.\)$'''),
                            ParseUnit('matrix',
                                '^(\w+)\ \ \ ?(.*)\n',
                                types=[str, str],
                                repeat=True)
                            ]),
                 ParseBloc('output', 
                            units=[
                        ParseUnit('numbered tree',
                                  r'^TREE # *\d+: +(.*;) +MP score: (-?\d+)$',
                                  keys=['tree', 'MP score'], types=[str, int]),
                        ParseUnit('lnL',
                            r'^lnL\(ntime: *(\d+)\s+np: *(\d+)\):\s+([0-9.+-]+)\s+([0-9.+-]+)$',
                        keys=['ntime', 'np', 'loglik', 'error'],
                        types=[int, int, float, float]),
                        ParseUnit('branches', r'((?: +[0-9]+\.\.[0-9]+)+) *$',
                                  types=[strlist]),
                        ParseUnit('branch lengths + parameters', r'((?: +[0-9]+\.[0-9]+)+)$',
                                  types=[floatlist]),
                        ParseUnit('tree length', r'^tree length =\s+([0-9.-]+)$',
                                  types=[float]),
                        ParseUnit('Detailed output identifying parameters', 
                                  r'^Detailed output identifying parameters'),
                        ParseUnit('kappa', r'^kappa \(ts/tv\) =\s+([0-9.-]+)$',
                                  types=[float]),
                        ParseUnit('omega',
                                  r'^w \(dN/dS\) for branches:\s+([ 0-9.-]+)$',
                                  types=[floatlist]),
                        ParseBloc('dNdS', #start='dN & dS for each branch',
                                  # after follows the table
                                  units=[ParseUnit('header',
                                          r'^\s*(branch +t +N +S +dN/dS +dN +dS +N\*dN +S\*dS)$',
                                          types=[strlist]),
                                      ParseUnit('rows',
                                          (r'^ *(\d+\.\.\d+)((?: +[0-9.-]+){8})$'),
                                          types=[str, floatlist],
                                          repeat=True,
                                          convert=OrderedDict)]),
                        ParseUnit('tree length for dN',
                                  r'^tree length for dN:\s+([0-9.-]+)$', 
                                  types=[float]),
                        ParseUnit('tree length for dS',
                                  r'^tree length for dS:\s+([0-9.-]+)$',
                                  types=[float]),
                        ParseUnit('dS tree', r'^dS tree:\n(.*;)$'),
                        ParseUnit('dN tree', r'^dN tree:\n(.*;)'),
                        ParseUnit('w ratios as labels for TreeView',
                                  r'^w ratios as labels for TreeView:\n(.*;)')]),
                        ParseUnit('Time used', r'^Time used:\s+(.*)$')])


def parse_mlc(mlcfile):
    with open(mlcfile) as mlc:
        return mlc_parser.parse(mlc.read())[0]