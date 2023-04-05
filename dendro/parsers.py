#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Module centralizing multiple phylogenetic tree parsers"""

from sys import stdin
import os.path as op
import re
from collections import OrderedDict
from functools import partial
import logging
logger = logging.getLogger(__name__)


CONTROL = re.compile(r'[;"\']')

def read_multinewick(lines, stripchars='\r\n'):
    """lines must be a file type or an iterable of lines"""
    newick = ''
    i = 0
    for ln, line in enumerate(lines):
        line = line.strip(stripchars)
        logger.debug('stripped line: %r/%r/%r' % (line, line.strip(stripchars), line.strip('\r\n')))
        while line:
            match = CONTROL.search(line)
            if not match:
                logger.debug('full line %d: %r' % (ln, line))
                newick += line
                line = ''
            elif match.group() in ('"', "'"):
                try:
                    closing_pos = line[match.end():].index(match.group())
                    quoted_end = match.end() + closing_pos + 1
                except ValueError:
                    quoted_end = len(line)
                newick += line[:quoted_end]
                line = line[quoted_end:]
            else:
                logger.debug('break line %d: %r' % (ln,line))
                newick += line[:match.end()-1]
                yield newick + ';'
                i += 1
                newick = ''
                line = line[match.end():]
    if newick.strip('\r\n;'):
        logger.debug('Remainder of line %d: %r' % (ln, newick))
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


def read_nexus_trees(handle):
    """Return the tree strings"""
    for line in handle:
        if line.lstrip()[:4].lower() == 'tree':
            yield line.split(maxsplit=1)[1].rstrip()


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
                err.args = (err.args[0] + ' ERROR with treefile %r ...' % treefile[:50],)
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



PARSERS = {'phyltree': iter_as_phyltree,
           'prottree': iter_as_prottree,
           'proteintree': iter_as_prottree,
           'ete3':     iter_as_ete3,
           'ete3_f1':  partial(iter_as_ete3, format=1),
           'skbio':    iter_as_skbio,
           'biophylo': iter_as_biophylo}


def eval_optiontext(argtext):
    """Build list of parser arguments from a code string (comma separated)"""
    kwargs = {}
    for item in argtext.split(','):
        try:
            key, value = item.split('=', maxsplit=1)
        except ValueError:
            if item:
                kwargs[item.strip()] = True
            continue

        key = key.strip()
        value = value.strip()

        if value.lower() in ('none', 'true', 'false'):
            value = eval(value.capitalize())
        else:
            try:
                value = int(value)
            except ValueError:
                try:
                    value = float(value)
                except ValueError:
                    pass

        kwargs[key] = value
    return kwargs


PARSERSPEC_HELP = """
Parser specification code
=========================

# Parser name

Name of a python module for parsing tree files. It can be one of (case insensitive):

    Ete3, PhylTree, ProtTree, BioPhylo, skbio.

# Parser options

Arbitrary reader/writer options can be appended, after a colon, separated by commas:

* "key=value" form, such as "format=2"
* "flag" form, which sets "flag=True"

Ete3 short option codes:

* reader/writer options:
    f: format  q: quoted_node_names  0: format=0  1: format=1  ... 100: format=100
* writer-only options (ignored for parsing):
    r: format_root_node  d: dist_formatter  s: support_formatter  n: name_formatter

# Examples

    "ete3:1"   -> "ete3:format=1"
    "ete3:q"   -> "ete3:quoted_node_names=True
    "ete3:1,q" -> "ete3:format=1,quoted_node_names=True
"""

# Short option names for the ete3 parser specification:
LONG_KWARGS_ETE3 = dict(f='format', q='quoted_node_names', r='format_root_node',
                        d='dist_formatter', s='support_formatter', n='name_formatter',
                        **{i: ('format', int(i)) for i in list('0123456789')+['100']})


def expand_short_options(options: dict, long_codes: dict):
    for key, val in list(options.items()):
        try:
            new = long_codes[key]
        except KeyError:
            continue
        newval = options.pop(key)
        if isinstance(new, tuple):
            new, newval = new
        options[new] = newval


def chooseparser(parserspec, silent=False):
    parsername, _, argtext = parserspec.partition(':')
    parsername = parsername.strip().lower()

    if not argtext:
        return PARSERS[parsername]

    kwargs = eval_optiontext(argtext)

    if parsername == 'ete3':
        # Replace the short option names
        expand_short_options(kwargs, LONG_KWARGS_ETE3)
        # discard options that are only used in output:
        for out_option in ('format_root_node', 'dist_formatter', 'support_formatter', 'name_formatter'):
            try:
                del kwargs[out_option]
                if not silent:
                    logger.warning('Ete3 writer option %r ignored for parsing' % out_option)
            except KeyError:
                pass

    return partial(PARSERS[parsername], **kwargs)
