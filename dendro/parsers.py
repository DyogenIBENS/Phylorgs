#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Module centralizing multiple phylogenetic tree parsers"""

from sys import stdin
import os.path as op
import re
from collections import OrderedDict
import logging
logger = logging.getLogger(__name__)


def read_multinewick(lines, stripchars='\r\n'):
    """lines must be a file type or an iterable of lines"""
    newick = ''
    for line in lines:
        parts = line.strip(stripchars).split(';')
        newick += parts.pop(0)
        while parts:
            yield newick + ';'
            newick = parts.pop(0)
    if newick.strip(stripchars):
        yield newick + ';'


def nhx_comment_parser(text):
    variables = OrderedDict()
    if not text:
        return variables

    if text.startswith('&&NHX:'):
        text = text[6:]
    elif text.startswith('[&&NHX:') and text.endswith(']'):
        text = text[7:-1]
    else:
        logger.warning("NHX comment appears misshaped; lacking '&&NHX': %r", text)


    for keyval in text.split(':'):
        try:
            key, val = keyval.split('=', maxsplit=1)
        except ValueError as err:
            err.args += ('from "%s"' % keyval,)
            raise

        try:
            val = int(val)
        except ValueError:
            try:
                val = float(val)
            except ValueError:
                if val.lower() in ('none', 'false', 'true'):
                    val = eval(val.capitalize())
        variables[key] = val
    return variables


def beast_comment_parser(text):
    """A parser for Nexus tree comments formatted as '[&var1=value1,var2=value2]'

    More specific than Nexus, because it parses interval values
    '{lower_bound,higher_bound}' into a list.
    """
    variables = OrderedDict()
    #reg_interval = re.compile('^{-?\d+')
    if not text:
        return variables

    if not (text.startswith('[&') and text.endswith(']')) and not text.startswith('&'):
        logger.warning("comment appear misshaped; lacking '[& ]': %r", text)

    text = re.sub(r'^&|^\[&|\]$', '', text)
    while text:
        try:
            var, text = text.split('=', maxsplit=1)
        except ValueError as err:
            err.args += ('from "%s"' % text,)
            raise

        if text.startswith('{'):
            end = text.find('}')
            value = [float(x) for x in text[1:end].split(',')]
            text = text[(end+1):].lstrip(',')
        elif text.startswith('"{'):
            end = text.find('}"')
            value = [float(x) for x in text[2:end].split(',')]
            text = text[(end+2):].lstrip(',')
        else:
            try:
                value, text = text.split(',', 1)
            except ValueError:
                value = text
                text = ''
            try:
                value = float(value)
            except ValueError:
                pass
        variables[var] = value

    return variables


def read_newick2nexus(handle):
    from Bio import Phylo
    from Bio.Nexus import Nexus
    from dendro.converters import BioPhylo_to_BioNexusTree
    nx = Nexus.Nexus()
    for tree in Phylo.parse(handle, 'newick'):
        #nx.trees.append(Nexus.Tree(tree.format('newick')))  # Sorry, I'm relying on serialization into a newick string.
        nx.trees.append(BioPhylo_to_BioNexusTree(tree))
    return nx


def iter_as_phyltree(treefile, *args, **kwargs):
    from LibsDyogen import myPhylTree
    yield myPhylTree.PhylogeneticTree(treefile, *args, **kwargs)

def iter_as_prottree(treefile, *args, **kwargs):
    from LibsDyogen import myProteinTree
    for tree in myProteinTree.loadTree(treefile, *args, **kwargs):
        yield tree

def iter_as_ete3(treefile, *args, **kwargs):
    import ete3
    if treefile is not stdin and not op.exists(treefile):
        yield ete3.Tree(treefile, *args, **kwargs)
        return

    f = treefile if treefile is stdin else open(treefile)
    try:
        for newick in read_multinewick(f):
            try:
                yield ete3.Tree(newick, *args, **kwargs)
            except ete3.parser.newick.NewickError as err:
                err.args = (err.args[0] + 'ERROR with treefile %r ...' % treefile[:50],)
                raise
    finally:
        if f is not stdin:
            f.close()

def iter_as_biophylo(treefile, *args, **kwargs):
    from Bio import Phylo
    for tree in Phylo.parse(treefile, *args, **kwargs):
        yield tree

def iter_as_skbio(treefile, *args, **kwargs):
    import skbio
    raise NotImplementedError()



parserchoice = {'phyltree': iter_as_phyltree,
                'PhylTree': iter_as_phyltree,
                'Phyltree': iter_as_phyltree,
                'prottree': iter_as_prottree,
                'ProtTree': iter_as_prottree,
                'Prottree': iter_as_prottree,
                'ete3':     iter_as_ete3,
                'Ete3':     iter_as_ete3,
                'ete3_f1':  lambda treefile: iter_as_ete3(treefile, format=1),
                'Ete3_f1':  lambda treefile: iter_as_ete3(treefile, format=1),
                'skbio':    iter_as_skbio,
                'biophylo': iter_as_biophylo}
