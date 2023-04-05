#!/usr/bin/env python3
# -*- coding:utf8 -*-


"""Interactively display a tree using Ete3."""


from sys import stdin, stderr
import os.path as op
import ete3
import argparse


def fileinputtrees():
    newicks = ['']

    for line in stdin:

        line_chunks = line.rstrip().split(';')

        newicks[-1] += line_chunks.pop(0)
        for chunk in line_chunks:
            newicks[-1] += ';'
            newicks.append(chunk)  # Works even when chunk == ''

    if not newicks[-1]:
        newicks.pop()
    return newicks


def make_treename(newickfile):
    try:
        newickfile = newickfile.name
    except AttributeError:
        pass

    if op.isfile(newickfile):
        return op.splitext(op.basename(newickfile))[0]


PARSER_OPTIONS = {'q': ('quoted_node_names', True),
                  **{i: ('format', int(i)) for i in list('0123456789') + ['100']}}


def main(newickfiles, rootnode=None, ladderize=False,
         show_internal=False, show_dists=False, nodesize=1, text=False,
         compact=False, output=None, parser='ete3:1'):
    parser, _, parser_optiontext = parser.partition(':')
    parser = parser.strip().lower()
    if parser == 'ete3':
        # try to not depend on the dendro.converters module
        parser_options = parser_optiontext.split(',')
        parser_kwargs = {}
        for opt in parser_options:
            try:
                key, value = PARSER_OPTIONS[opt]
            except KeyError:
                raise ValueError('Invalid Ete3 parser option %r (not in q,0..9,100)' % opt)
            parser_kwargs[key] = value
        def parse_trees(newickfiles):
            if not newickfiles:
                newickfiles = fileinputtrees()
            for newick in newickfiles:
                try:
                    tree = ete3.Tree(newick, **parser_kwargs)
                except:
                    print("newick=%r" % newick, file=stderr)
                    raise
                yield tree, make_treename(newick)
    else:
        from dendro.parsers import chooseparser, eval_optiontext
        from dendro.converters import converterchoice
        treeparser = chooseparser(parser)
        treeconverter = converterchoice[parser]['ete3']
        parser_kwargs = eval_optiontext(parser_optiontext)
        def parse_trees(treefiles):
            if not treefiles:
                treefiles = [stdin]
            for treefile in treefiles:
                for tree in treeparser(treefile, **parser_kwargs):
                    yield treeconverter(tree), make_treename(treefile)

    for i, (tree, treename) in enumerate(parse_trees(newickfiles), start=1):
        print('Tree %d' % i)
        display_onetree(tree, treename, rootnode, ladderize, show_internal,
                        show_dists, nodesize, text, compact, output)


def display_onetree(tree, treename=None, rootnode=None, ladderize=False,
                    show_internal=False, show_dists=False, nodesize=1,
                    text=False, compact=False, output=None):

    if rootnode:
        tree = tree.search_nodes(name=rootnode)[0]

    if ladderize:
        tree.ladderize()

    if text:
        show_attr=['dist'] if show_dists else None
        print(tree.get_ascii(show_internal=show_internal, compact=compact,
                             attributes=show_attr))
    else:
        mynodestyle = ete3.NodeStyle(size=nodesize)

        ts = ete3.TreeStyle()
        if show_dists:
            ts.show_branch_length = True

        def mybasiclayout(node):
            node.set_style(mynodestyle)

        if show_internal:
            nameface = ete3.faces.AttrFace("name")
            def mylayout(node):
                mybasiclayout(node)
                if not node.is_leaf():
                    ete3.faces.add_face_to_node(nameface, node, column=0)
        else:
            mylayout = mybasiclayout

        if output:
            tree.render(output, tree_style=ts, layout=mylayout)
        else:
            tree.show(tree_style=ts, layout=mylayout, name=treename)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('newickfiles', nargs='*')
    parser.add_argument('-r', dest='rootnode', help='rootnode: name of an ' \
                        'internal node to take as root.')
    parser.add_argument('-l', dest='ladderize', action='store_true',
                        help='Ladderize tree before display')
    parser.add_argument('-i', dest='show_internal', action='store_true', 
                        help='Display internal nodes names')
    parser.add_argument('-d', dest='show_dists', action='store_true', 
                        help='Display branch lengths')
    parser.add_argument('-s', '--nodesize', type=int, default=1)
    parser.add_argument('-t', '--text', action='store_true', 
                        help='Print as text')
    parser.add_argument('-c', '--compact', action='store_true',
                        help='Compact output (only with --text)')
    parser.add_argument('-o', '--output', help='save to output file')
    parser.add_argument('-p', '--parser', default='ete3:1',
                        help='Implementation to use (ete3/prottree/phyltree + optional args) [%(default)s]')

    args = parser.parse_args()

    main(**vars(args))

