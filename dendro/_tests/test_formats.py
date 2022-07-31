#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from io import StringIO
import pytest
import ete3
from Bio import Phylo
from Bio.Phylo.Newick import Clade, Tree as NwkTree
from Bio.Nexus import Nexus
import Bio.Nexus.Trees as NexTrees
from dendro.formats import *


NHX_COMMENT = "&&NHX:varX=5:varY=6"
NEXUS_COMMENT = "&varX=5,varY=6"


@pytest.mark.parametrize('bracketed', [False, True])
def test_nhx_comment_parser(bracketed):
    comment = NHX_COMMENT
    if bracketed:
        comment = '[' + comment + ']'
    variables = nhx_comment_parser(comment)
    assert set(variables) == set(('varX', 'varY'))
    assert variables['varX'] == 5
    assert variables['varY'] == 6



class Test_apply_node_formats:
    def setup_method(self):
        a = Clade(name='a', comment=NEXUS_COMMENT)
        b = Clade(name='b', comment='&varX=9,varZ=uvw')
        self.tree = NwkTree(root=Clade(name='r', comment='&varX=0',
                                       confidence=50,
                                       clades=[a, b]))
        apply_node_formats(self.tree)

    def test_comments_from_beast_to_nhx(self):
        tree, (a, b) = self.tree, self.tree.get_terminals()

        # Check the comment conversion
        assert tree.root.comment == '&&NHX:varX=0'
        assert a.comment == NHX_COMMENT
        assert b.comment == '&&NHX:varX=9:varZ=uvw'

    def test_confidence_to_name_for_internal_nodes(self):
        tree, (a, b) = self.tree, self.tree.get_terminals()

        assert tree.root.name == '50'
        assert a.name == 'a'
        assert b.name == 'b'

    def test_readable_by_ete3(self):
        newick = self.tree.format('newick')
        r = ete3.Tree(newick, format=1)
        assert r.name == '50'
        assert len(r.children) == 2
        a, b = r.children
        assert a.name == 'a'
        assert b.name == 'b'
        # NHX comments correctly parsed:
        assert a.varX == '5' and a.varY == '6'
        assert b.varX == '9' and b.varZ == 'uvw'
        assert r.varX == '0'


class Test_apply_nexus_node_formats:
    def test_apply(self):
        a = NexTrees.Nodes.Node(NexTrees.NodeData(taxon='a', comment=NEXUS_COMMENT))
        b = NexTrees.Nodes.Node(NexTrees.NodeData(taxon='b', comment='&varX=9,varZ=uvw'))
        
        self.tree = NexTrees.Tree()
        root = self.tree.node(0)
        root.set_data(NexTrees.NodeData(taxon='r', comment='&varX=0', support=50))
        self.tree.add(a, 0)
        self.tree.add(b, 0)

        apply_nexus_node_formats(self.tree)


# internal names
NWK_BASE     = "((a:0.7[%s],b:0.5)x:0.2,c:1):0;"
NWK_BASE_NUM = "((1:0.7[%s],2:0.5)4:0.2,3:1):0;"  # requires Translate


NEXUS_TXT = """#NEXUS
Begin Trees;
Tree NAME_0 = %s
End;
""" % (NWK_BASE % NEXUS_COMMENT)

NHX_TXT = NWK_BASE % NHX_COMMENT


def test_read_newick2nexus():
    nx = read_newick2nexus(StringIO(NHX_TXT))
    assert len(nx.trees) == 1
    tree = nx.trees[0]
    assert len(tree.get_terminals()) == 3
    id_a = tree.search_taxon('a')
    comment_a = tree.node(id_a).data.comment
    assert comment_a.strip('[]') == NHX_COMMENT


def standardize_nexus(txt):
    """Uniformize capitalization, spacing and indentation"""
    special_words = {'begin': 'Begin', 'end': 'End', 'end;': 'End;',
                     'trees': 'Trees', 'tree': 'Tree', 'translate': 'Translate',
                     'taxa': 'Taxa', 'taxlabels': 'TaxLabels', 'dimensions': 'Dimensions'}
    return '\n'.join(' '.join(special_words.get(w.lower(), w) for w in line.strip().split())
                     for line in txt.splitlines() if line.strip()) + '\n'


class Test_write_nexus_trees:
    def test_simplest(self):
        nx = Nexus.Nexus(StringIO(NEXUS_TXT))
        apply_nexus_node_formats(nx.trees[0])
        out = StringIO()
        write_nexus_trees(nx, out)
        output = standardize_nexus(out.getvalue().rstrip())
        assert output == NEXUS_TXT

    def test_translated(self):
        nexus_txt = """#NEXUS
        Begin Trees;
        Translate 1 a,\n2 b,\n3 c,\n4 x;
        Tree NAME_0 = %s
        End;
        """ % (NWK_BASE_NUM % NEXUS_COMMENT)

        nx = Nexus.Nexus(StringIO(nexus_txt))
        apply_nexus_node_formats(nx.trees[0])

        out = StringIO()
        write_nexus_trees(nx, out, format_confidence='%d')
        output = standardize_nexus(out.getvalue().rstrip())
        
        # Do not compare the tree line, it's work in progress
        lines = output.splitlines()
        (treelineno, treeline), = [(i, line) for i, line in enumerate(lines) if line.startswith('Tree NAME_0')]

        obtained = '\n'.join(lines[:treelineno])
        expected = standardize_nexus(nexus_txt)[:len(obtained)]
        assert expected == obtained

        assert treeline == 'Tree NAME_0 = ' + NWK_BASE_NUM % NEXUS_COMMENT



class Test_nexus2nhx:
    def test_simple(self):
        out = StringIO()
        nexus2nhx(StringIO(NEXUS_TXT), out)
        output = out.getvalue().rstrip()
        assert output == NHX_TXT


def is_int_repr(txt):
    try:
        k = int(txt)
        return True
    except (ValueError, TypeError):
        return False


#class Test_write_nexus_trees_to_bayestraits:

class Test_nhx2bayestraits:
    def test_simple(self):
        out = StringIO()
        nhx2bayestraits(StringIO(NHX_TXT), out)
        output = standardize_nexus(out.getvalue().rstrip())
        assert 'Translate' in output
        treeline, = [line for line in output.splitlines()
                     if line.startswith('Tree ')]
        nx = Nexus.Nexus(output)
        
        all_labels = []
        all_supports = []
        for node_id, node in nx.trees[0].chain.items():
            all_labels.append(node.data.taxon)
            all_supports.append(node.data.support)
            print('node %s comment:' % all_labels[-1], node.data.comment)
            if all_labels[-1] == 'a' or nx.translate.get(all_labels[-1]) == 'a':
                comment_a = node.data.comment

        print(output)
        print('Labels:', all_labels)
        print('Supports', all_supports)
        assert comment_a.strip('[]') == NEXUS_COMMENT
        # Skip the check that labels are integers.
        #assert all(map(is_int_repr, all_labels)), "labels should be integers"


class Test_nexus:
    def test_unchanged(self):
        out = StringIO()
        nexus_rewrite(StringIO(NEXUS_TXT), out, float_fmt='%s')
        output = standardize_nexus(out.getvalue().rstrip())


### Also test some converter utilities

from dendro.converters import BioPhylo_to_BioNexusTree, BioNexusTrees_to_BioPhylo


def test_BioPhylo_to_BioNexus():
    tree = Phylo.read(StringIO(NHX_TXT), 'newick')
    ntree = BioPhylo_to_BioNexusTree(tree)


class Test_BioNexus_to_BioPhylo:
    def test_notranslate(self):
        nx = Nexus.Nexus(StringIO(NEXUS_TXT))
        trees = BioNexusTrees_to_BioPhylo(nx.trees)
    def test_translate(self):
        nx = Nexus.Nexus(StringIO(NEXUS_TXT))
        trees = BioNexusTrees_to_BioPhylo(nx.trees, nx.translate)

