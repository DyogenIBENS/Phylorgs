#!/usr/bin/env python3

"""Parse results of codeml (.mlc files) to save dN, dS values in a table"""


from __future__ import print_function

import sys
import os.path
import re
from collections import OrderedDict, defaultdict
import numpy as np
#from numpy import array as a, concatenate as c
import ete3
import argparse
import logging

import LibsDyogen.myPhylTree as PhylTree # my custom python3 version

from genchron.MPL import *
from IOtools import Stream

mpllogger = logger
logger = logging.getLogger(__name__)


#from pamliped.codemlparser import mlc_parser

np.set_printoptions(formatter={"float_kind": lambda x: "%g" %x})

ENSEMBL_VERSION = 85
PHYLTREEFILE = "/users/ldog/glouvel/ws_alouis/GENOMICUS_SVN/data{0:d}/PhylTree.Ensembl.{0:d}.conf"

BEAST_MEASURES = set('%s%s' % (v,s) for v in ('height', 'length', 'rate')
                     for s in ('', '_median', '_95%_HPD', '_range')).union(('posterior',))
BEAST_DEFAULTS = {'%s%s' % (v,s): np.NaN for v in ('height', 'length', 'rate')
                  for s in ('', '_median')}
BEAST_DEFAULTS.update(posterior=np.NaN, dist=np.NaN)
#TODO? convert string "{0.0..42.42}" into low and up floats.
CODEML_MEASURES = set(('dN', 'dS', 't', 'N*dN', 'S*dS'))


def parse_beast_range(r):
    return [np.NaN, np.NaN] if r is None else [float(x) for x in r.strip('{}').split('..', 2)]


def reassign_beast_range(r, name):
    low, up = parse_beast_range(r)
    return {name+'_low': low, name+'_up': up}

BEAST_REASSIGNS = {'%s_%s' % (v,s): reassign_beast_range
                   for v in ('height', 'length', 'rate')
                   for s in ('95%_HPD', 'range')}


def ls(array):
    return str(array).replace('\n', '')

def printtree(tree, indent='', features=None, **kwargs):
    line = indent + tree.name
    if features:
        line += ': ' + ' '.join('%s=%s' % (ft,getattr(tree, ft))
                                for ft in features
                                if getattr(tree, ft, None))
    print(line, **kwargs)
    for ch in tree.children:
        printtree(ch, (indent+'  '), features)


def showtree(fulltree, ages):
    """Default `showtree` function: do nothing. This function can be
    redefined using `def_showtree`"""
    pass



def def_showtree(measures, show=None):  # ~~> genomicus.reconciled
#def def_showreconciledgenetree
    """Depending on the boolean argument 'show', return the function
    `showtree`:
    if show=None, this function does nothing;
    else, showtree display an interactive tree using ete3."""
    if show:
        if show == 'gui':
            def show_func(tree, ts=None):
                print("Showing tree.")
                tree.show(tree_style=ts, name=tree.name)
        elif show == 'notebook':
            def show_func(tree, ts=None):
                #print("Rendering tree.", file=sys.stderr)
                print(("Defining global FULLTREE and TS. Display with:\n"
                       "    FULLTREE.render('%%inline', tree_style=TS)"))
                #tree.render('%%inline', w=500, tree_style=ts) # 
                global FULLTREE, TS
                FULLTREE = tree
                TS = ts
        else:
            # Save to an image file.
            def show_func(tree, ts=None):
                print("Saving tree in %r." % show)
                tree.render(show, tree_style=ts)

        # define duplication node style:
        default_color = ete3.NodeStyle()['fgcolor']
        colors = {'spe': default_color,
                  'dup' :'red',
                  'leaf':'green',
                  'root':'black'}
            
        header = measures + ['age_'+m for m in measures] + ['type']
        if 'dist' in header: header.remove('dist')
        measure_faces = {feature: ete3.AttrFace(feature,
                                                text_prefix='%s: ' % feature,
                                                fsize=7,
                                                fgcolor='grey',
                                                formatter='%.4g')
                            for feature in header[:-1]}

        def mylayout(node):
            for feature in header[:-1]:
                if hasattr(node, feature):
                    node.add_face(measure_faces[feature], column=0,
                                  position='branch-bottom')
            node.set_style(ete3.NodeStyle(
                fgcolor=colors[getattr(node, 'type', 'root')],
                shape=('square' if getattr(node, 'cal', 0) else 'circle')))

        #leaf_info = dict(zip(header, [np.NaN]))

        def showtree(fulltree):
            # TODO: this should be a layout function. Then change ts.layout_fn.
            # define a tree style:
            ts = ete3.TreeStyle()
            ts.show_branch_length = True
            ts.layout_fn = mylayout

            show_func(fulltree, ts)
            return fulltree
    else:
        def showtree(fulltree):
            pass
    return showtree


def load_fulltree(resultfile, replace_nwk='.mlc', replace_by='.nwk'):
    """Return the ete3.Tree object corresponding to a given .mlc file.

    It simply finds the newick file by replacing the extension of the mlc file.
    Catch errors (file does not exist / wrong format)"""
    nwkfile = re.sub(replace_nwk, replace_by, resultfile)
    rootname = os.path.splitext(resultfile)[0]
    tree_alt_name = os.path.basename(re.sub(replace_nwk, '', resultfile))
    try:
        fulltree = ete3.Tree(nwkfile, format=1)
    except ete3.parser.newick.NewickError as e:
        if os.path.exists(nwkfile):
            logger.error("\nNewickError: Malformed newick tree structure in %r",
                          nwkfile)
        else:
            logger.error("\nNewickError: Unexisting tree file %r" % nwkfile)
        sys.exit(2)

    fulltree.add_features(basename=os.path.basename(rootname),
                          dirname=os.path.dirname(rootname),
                          tree_alt_name=tree_alt_name)
    return fulltree


def eat_lines_uptomatch(file_obj, regex, line=None, jump=0):
    if line is None:
        line = file_obj.readline()
    match = regex.match(line.rstrip())
    while not match:
        for _ in range(jump):
            file_obj.readline()

        line = file_obj.readline()
        if not line:
            raise LookupError("Reached end of file without matching '%s'" % regex.pattern)
        match = regex.match(line.rstrip())
    return line #, match


def yield_lines_whilematch(file_obj, regex, line=None, jump=0):
    if line is None:
        line = file_obj.readline()
    match = regex.match(line.rstrip())
    while match:
        yield line
        for _ in range(jump):
            file_obj.readline()

        line = file_obj.readline()
        match = regex.match(line.rstrip())


def update_tree_nodes(targettree, srctree, leaf1_2=None, update_features=['name'],
                      defaults=None, reassigns=None):
    """Appropriate for matching nodes only when the differences are polytomies
    or extra single-child nodes."""
    srcleaf_sets = OrderedDict()  # Conserve the tree post-ordering
    #clade2node = OrderedDict()
    if defaults is None:
        defaults = {}
    if reassigns is None:
        reassigns = {}

    for srcnode in srctree.traverse('postorder'):
        if srcnode.is_leaf():
            srcleaf_sets[srcnode] = set((srcnode.name,))
            #clade2node[(node.name,)] = node
        else:
            clade = set.union(*(srcleaf_sets[ch] for ch in srcnode.children))
            #cladekey = tuple(sorted(clade))
            #if cladekey not in clade2node:
            #    # Because in post-order traversal, overwrite can happen with 
            #    # single-child nodes, and we only want the MRCA.
            #    clade2node[cladekey] = node
            srcleaf_sets[srcnode] = clade

    leaf_sets = {}
    src_match = defaultdict(set)

    for node in targettree.traverse('postorder'):
        if node.is_leaf():
            clade = leaf_sets[node] = set((node.name,))
        else:
            clade = leaf_sets[node] = set.union(*(leaf_sets[ch]
                                                    for ch in node.children))
        
        # Iterate from the smallest clades to the largest. Do a `searchsorted`
        # Find a corresponding source node.
        for srcnode, srcclade in srcleaf_sets.items():
            # Find the first source clade containing or equal to the target clade.
            if srcclade >= clade:
                logger.info('Adding features from %r dist=%.2f -> %r', srcnode.name, srcnode.dist, node.name)
                newfeatures = dict((ft[0], getattr(srcnode, ft[1], defaults.get(ft[0])))
                                     if isinstance(ft, tuple)
                                     else (ft, getattr(srcnode, ft, defaults.get(ft)))
                                     for ft in update_features)
                logger.debug('srcnode features = %s; updated features = %s; target node features = %s', srcnode.features, update_features, node.features)
                for rawfeature, convert in reassigns.items():
                    try:
                        newfeatures.update(
                                **convert(newfeatures.pop(rawfeature,
                                                          getattr(srcnode, rawfeature)),
                                          rawfeature)
                                )
                    except AttributeError:
                        if not srcnode.is_root():
                            raise

                node.add_features(**newfeatures)
                if srcclade > clade:
                    logger.warning('Non exact clade match: %r != %r',
                                   srcnode.name, node.name)
                if srcnode in src_match:
                    logger.warning('Source %r already used for: %s (polytomy?)',
                                   srcnode.name, [n.name for n in src_match[srcnode]])
                src_match[srcnode].add(node)
                break
        else:
            # there is no match.
            logger.warning('Unmatched target node: %r', node.name)

    unmatched_leaves = leaf_sets[targettree] ^ srcleaf_sets[srctree]
    if unmatched_leaves:
        logger.warning('Unmatched leaves: %s', '; '.join(unmatched_leaves))

    unmatched_srcnodes = set(srcleaf_sets).difference(src_match)
    if unmatched_srcnodes:
        logger.warning('%d unmatched source nodes: %s', len(unmatched_srcnodes),
                        '; '.join('%r(%s)' % (n.name, id(n)) for n in unmatched_srcnodes))
                        #'; '.join(n.name for n in unmatched_srcnodes))



def branch2nb(mlc, fulltree):  # ~~> pamliped.codeml_parser?
    """Parse the codeml result file (.mlc) to return 2 trees:
    tree_nbs: the tree with node labelled as numbers.
    tree_ids: the tree with original node labels (including inner nodes).
    
    Also return the dictionaries to convert nb -> id and conversely.

    Arguments:
        - mlc     : an opened file (codeml result file)
        - fulltree: the ete3.Tree for the gene tree (with reinserted nodes)
    """
    regex_tree = re.compile(r'^(.*);$')
    regex_lnL  = re.compile(r'^lnL\(')
    regex_w    = re.compile(r'^w \(dN/dS\) for branches:')
    regex_freqtable = re.compile('^Codon position x base \(3x4\) table for each sequence\.$')
    regex_seqid = re.compile('^#(\d+): (.*)$')

    logger.debug('Starting mlc stream at byte %s.', mlc.tell())

    line = eat_lines_uptomatch(mlc, regex_freqtable)

    seqids = {}

    mlc.readline()
    line = mlc.readline()
    while True:
        m = regex_seqid.match(line.rstrip())
        if not m:
            break
        seqids[m.group(1)] = m.group(2)
        for _ in range(6):
            line = mlc.readline()
        assert line
        #if not line:
        #    break

    assert seqids

    # Get the line listing all branches
    while not regex_lnL.match(line):
        line = mlc.readline()
        assert line, "regex_lnL not matched before EOF"
        line = line.rstrip()
    branches_line = mlc.readline().rstrip()
    branches = branches_line.split()
    # get translation by looking at the newick tree lines.
    # Line containing the lengths *and* the additional parameter values.
    lengths_line = mlc.readline().rstrip()
    lengths = [float(l) for l in lengths_line.split()]

    # Get the tree labelled with numbers
    line = mlc.readline().rstrip()
    while not regex_tree.match(line):
        line = mlc.readline()
        assert line, "regex_tree (numbered) not matched before EOF"
        line = line.rstrip()

    tree_nbs = ete3.Tree(line)

    # Get the tree with original labels
    line = mlc.readline().rstrip()
    while not regex_tree.match(line):
        line = mlc.readline()
        assert line, "regex_tree (labeled) not matched before EOF"
        line = line.rstrip()
    tree_ids = ete3.Tree(line)

    # Get the omega (dN/dS) values
    checkpoint = mlc.tell()
    model = 1  # Indicator variable: are we parsing the `model = 1` output ?
    
    line = mlc.readline().rstrip()
    while not regex_w.match(line):
        line = mlc.readline()
        if not line:
            if model == 1:
                model = 5
                #"regex_w not matched before EOF"
                regex_w = re.compile(r'^Parameters in M5 \(gamma\):$')
                mlc.seek(checkpoint)
                continue
            elif model == 5:
                #raise LookupError('regex_w (model 5) not matched before EOF')
                model = 0
                regex_w = re.compile(r'^omega \(dN/dS\) =\s+([0-9.-]+)$')
                mlc.seek(checkpoint)
                continue
            else:
                raise LookupError('regex_w (model 1/model 5/model 0) not matched before EOF')
        line = line.rstrip()

    if model == 1:
        omegas = [float(w) for w in regex_w.sub('', line).split()]
    elif model == 0:
        omegas = [float(regex_w.match(line).group(1))] * len(branches)
    elif model == 5:
        line = mlc.readline()
        try:
            gamma_params = re.match(r'^ +a= +([0-9.-]+) +b= +([0-9.-]+)$', line).groups()
        except AttributeError as err:
            err.args += ('line = %r' % line,)
            raise
        regex_MLEw = re.compile(r'^MLEs of dN/dS \(w\) for site classes \(K=(\d+)\)')
        line = mlc.readline()
        while not regex_MLEw.match(line.rstrip()):
            line = mlc.readline()
            assert line, "regex_MLEw not matched before EOF"

        mlc.readline()  # known blank line
        # Parse Omegas for each site fraction
        sitefrac = [float(x) for x in
                 re.match('p:( +0\.[0-9]+)+$', mlc.readline()).group(1).split()]
        ws = [float(x) for x in
            re.match('w:( +[0-9]+\.[0-9]+)+$', mlc.readline()).group(1).split()]
        omegas = [np.average(ws, weights=sitefrac)] * len(branches)

    assert len(branches) == len(omegas) and len(branches) <= len(lengths), \
            (len(branches), len(omegas), len(lengths))
    branch_tw = list(zip(branches, lengths, omegas))

    id2nb = dict(zip(tree_ids.get_leaf_names(), tree_nbs.get_leaf_names()))
    nb2id = dict(zip(tree_nbs.get_leaf_names(), tree_ids.get_leaf_names()))

    assert set(seqids) == set(nb2id) and all((seqid == seqids[nb]) for nb,seqid in nb2id.items())

    # Remember nodes that were detached from fulltree, and all their descendants
    detached_subtrees = set()

    # get internal nodes nb-to-id conversion (fulltree: tree with internal node
    # annotations)
    for leafid in fulltree.get_leaves():
        try:
            leafid.add_feature('nb', id2nb[leafid.name])
        except KeyError as err:
            # This leaf is not a species node. skip.
            pass

    while branches:
        br = branches.pop()
        base, tip = br.split('..')
        debug_msg = "%-8s " % br
        # Update the tree_nbs
        base_nb_node = tree_nbs.search_nodes(name=tip)[0].up
        base_nb_node.name = base
        # Then find the matching node in fulltree
        try:
            tip_id = nb2id[tip]
            found_tips = fulltree.search_nodes(name=tip_id)
            if len(found_tips) > 1:
                logger.error('Fulltree has nodes with identical names. The correct conversion is not guaranteed.')
            try:
                base_node = found_tips[-1].up
            except IndexError as err:
                logger.error('Node %s:%r not found in fulltree', tip, tip_id)
                      #file=sys.stderr)
                # TODO: search in tree_ids, then detached_subtrees.add()
                detached_subtrees.add(tip)
                continue

                #while len(base_node.children) == 1:
                #    base_node = base_node.up
            except AttributeError as err:
                logger.warning('root (%r: %r) cannot have ancestors', tip, tip_id)
                #print(fulltree.search_nodes(name=tip_id), err)
                continue

            base_id = base_node.name
            if 'nb' not in base_node.features:
                # Add number in the fulltree:
                base_node.add_feature('nb', base)
                # Avoid duplicate IDs.
                duplicated = 0
                base_id_fmt = base_id + '%d'
                while base_id in id2nb:
                    logger.warning('base_id %r already seen. Increment.', base_id)
                    duplicated += 1
                    base_id = base_id_fmt % duplicated
                    base_node.name = base_id
                nb2id[base] = base_id
                id2nb[base_id] = base

            debug_msg += "%s -> %s  " % (base_id, tip_id)
            logger.debug(debug_msg + 'Ok')
        except KeyError as e:
            #if base in detached_subtrees:
            # I assume a progression from leaves to root, otherwise I will miss nodes
            logger.debug(debug_msg + 'Detached (%r not in nb2id -> %r)', tip, e)
            #else:
            # Not found now, put it back in the queue for later
            #    branches.insert(0, (base, tip))
            #    logger.debug('KeyError')

    logger.debug(tree_nbs.get_ascii())
    logger.debug('tree_ids:' + tree_ids.get_ascii())
    logger.debug('fulltree: ' + fulltree.get_ascii())
    #logger.debug(fulltree.get_ascii())
    #printtree(fulltree, features=['nb'])
    #showtree(fulltree)
    return id2nb, nb2id, tree_nbs, branch_tw, model

# Needs a unit test for the above (are the missing nodes properly re-inserted?)

def get_dNdS(mlc, skiptrees=False):  # ~~> pamliped.codeml_parser?
    """Parse table of dN/dS from codeml output file.
    
    mlc: filehandle

    WARNING: Numbers are rounded!!! (e.g, only 3 decimals for `t`, 4 for `dS`).
    """
    #reg_dNdS = re.compile(r'dN & dS for each branch')
    reg_dNdS = re.compile(r'\s+branch\s+t\s+N\s+S\s+dN/dS\s+dN\s+dS\s+N\*dN' \
                          r'\s+S\*dS')
    reg_dStree = re.compile(r'dS tree:$')
    reg_dNtree = re.compile(r'dN tree:$')
    
    line = mlc.readline().rstrip()
    while not reg_dNdS.match(line):
        line = mlc.readline()
        assert line, "reg_dNdS not matched before EOF"
        line = line.rstrip()

    assert mlc.readline().rstrip() == ''  # skip blank line
    dNdS = {}  # OrderedDict
    dNdS['colnames'] = line.split()[1:]
    line = mlc.readline().rstrip()
    while line != '':
        #print line
        fields = line.split()
        dNdS[fields[0]] = [float(x) for x in fields[1:]]
        line = mlc.readline().rstrip()

    if skiptrees:
        return dNdS, None, None

    line = mlc.readline().rstrip()
    while not reg_dStree.match(line):
        line = mlc.readline()
        assert line, "reg_dStree not matched before EOF"
        line = line.rstrip()
    dStreeline = mlc.readline().rstrip()
    
    line = mlc.readline().rstrip()
    while not reg_dNtree.match(line):
        line = mlc.readline()
        assert line, "reg_dNtree not matched before EOF"
        line = line.rstrip()
    dNtreeline = mlc.readline().rstrip()
    
    #TODO: np.array(dNdS)
    return dNdS, dStreeline, dNtreeline
    

def tree_nb_annotate(tree, id2nb, tree_nbs):  # ~~> pamliped.codeml_parser?
    """Add internal node names (numbers used by codeml) in the tree structure."""
    parent_nbs = {tuple(sorted(ch.name for ch in n.children)): n.name \
                for n in tree_nbs.traverse() if not n.is_leaf()}
    for node in tree.traverse('postorder'):
        if node.is_leaf():
            node.add_feature('nb', id2nb[node.name])
        else:
            node.add_feature('nb',
                    parent_nbs[tuple(sorted(ch.nb for ch in node.children))])


def dNdS_precise(dNdS, br_tw, id2nb, tree_nbs, dSnwk=None, dNnwk=None):  # ~~> codeml_parser?
    """Make the dNdS table from the codeml output file more precise:
    dS and dN values have more decimal in the tree string than in the table."""
    # First, update the dNdS table with the more accurate values br_lengths
    # (br_len is actually 't' in the table)
    colindex = {colname: i for i, colname in enumerate(dNdS['colnames'])}
    for br, br_len, br_w in br_tw:
        dNdS[br][colindex['t']] = br_len # 't'
        dNdS[br][colindex['dN/dS']] = br_w   # 'dN/dS'

    # Then update the dNdS table with the more accurate values in dStree dNtree
    col_dS, col_dN = colindex['dS'], colindex['dN']
    #col_S, col_N = colindex['S'], colindex['N']

    colname_and_tree = []
    if dSnwk is not None:
        dStr = ete3.Tree(dSnwk)
        tree_nb_annotate(dStr, id2nb, tree_nbs)
        colname_and_tree.append((col_dS, dStr))
    if dNnwk is not None:
        dNtr = ete3.Tree(dNnwk)
        tree_nb_annotate(dNtr, id2nb, tree_nbs)
        colname_and_tree.append((col_dN, dNtr))

    for col, tr in colname_and_tree:
        for node in tr.traverse():
            if not node.is_leaf():
                for child in node.children:
                    br = '%s..%s' % (node.nb, child.nb)
                    dNdS[br][col] = child.dist
    # N*dN and S*dS could also be more precise by calculating the product, 
    # but this amount of precision might not be relevant.

    if dSnwk is not None or dNnwk is not None:
        for branch, row in dNdS.items():  # Isn't it time for vectorized operations?
            if branch != 'colnames':
                #print(row, file=sys.stderr)
                row[colindex['S*dS']] = row[colindex['S']] * row[col_dS]
                row[colindex['N*dN']] = row[colindex['N']] * row[col_dN]
    return dNdS


#def set_data_fulltree(fulltree, datatable, raise_at_intermediates=True,
#                      proportional_var=('t', 'dS', 'dN', 'dist', 'N*dN', 'S*dS',
#                                        'length', 'length_median'):
def set_dNdS_fulltree(fulltree, id2nb, dNdS, raise_at_intermediates=True):
    """Add codeml dS on each branch of the complete tree (with missing species
    nodes).
    When a node has a single child: take the branch dS on which it is, 
    and divide proportionnally to each segment dist."""

    colindex = {colname: i for i, colname in enumerate(dNdS['colnames'])}
    branchnames = list(dNdS)
    branchnames.remove('colnames')

    # Set up the root values: N and S
    if 'S' in fulltree.features:
        fulltree.add_feature('Sp', fulltree.S)
    fulltree.add_features(S=dNdS[branchnames[0]][colindex['S']],
                          N=dNdS[branchnames[0]][colindex['N']])

    for node in fulltree.get_descendants('postorder'): # exclude root
        if len(node.children) != 1:
            parent = node.up
            intermediates = []
            #at root, up is None
            while parent.up and getattr(parent, 'reinserted', None):
                intermediates.append(parent)
                parent = parent.up
            #print("%s:%s -> %s:%s\n   %s" % (id2nb[parent.name],
            #                          parent.name,
            #                          id2nb[node.name],
            #                          node.name,
            #                          [ch.name for ch in parent.children]), file=sys.stderr)
            
            if parent.is_root():
                # Disable dating of nodes immediately below the root, as the
                # tree should be considered unrooted and those values aren't
                # reliable.
                logger.info("Node below root: %s. Discard measures.", node.name)
                totals = {valname: np.NaN for valname in colindex}
            else:
                try:
                    br = id2nb[parent.name] + '..' + id2nb[node.name]
                    totals = {valname: dNdS[br][i] for valname, i in colindex.items()}
                    #totals = {valname: datatable[node.name] for valname, i in colindex.items()}

                except KeyError as err:
                    logger.debug('id2nb = %s', id2nb)
                    err.args = (err.args[0] + ' (parent=%r node=%r)' % (parent.name, node.name),)
                    raise
    
            if intermediates: # this `if` is not necessary, maybe more efficient
                if raise_at_intermediates:
                    raise AssertionError('Node %r has a single child.' % parent.name)

                dist_tot = sum(n.dist for n in (intermediates + [node]))

                if dist_tot == 0:
                    logger.warning("WARNING: dS = 0 between %r and %r",
                                    parent.name, node.name)
                    dist_tot = 1. # should set to NaN to be sure.

                # Constant rates extrapolation.
                for inter_node in (intermediates + [node]):
                    dist_ratio = inter_node.dist / dist_tot
                    for valname, val in totals.items():
                        if valname in ('t', 'dN', 'dS', 'N*dN', 'S*dS'):
                            inter_node.add_feature(valname, val * dist_ratio)
                        else:
                            inter_node.add_feature(valname, val)
            else:
                for valname, val in totals.items():
                    if valname in ('N', 'S'):  # and not node.is_root()
                        continue
                    if valname in node.features:
                        # Do not override existing features (e.g: 'S' is the species)
                        logger.warning("Feature %r already in node, rename as '%s0'",
                                       valname, valname)
                        node.add_feature('%s0' % valname, getattr(node, valname))
                        # Better solution: prefix all new features with `mlc_`
                    node.add_feature(valname, val)


def rm_erroneous_ancestors(fulltree, phyltree):  # ~~> genomicustools
    """Some ancestors with only one child are present in the species phylogeny.
    They must be removed when they have an age of zero"""
    infinite_dist = 100
    for node in fulltree.iter_descendants():
        if node.dist > infinite_dist:
            logger.warning('DETACH node with dist>%d : %s, dist=%d',
                  infinite_dist, node.name, node.dist)
            # Detaching/deleting would create a mismatch between fulltree and
            # tree_nbs/tree_ids
            #node.add_feature('skip', True)
            parent = node.up
            node.detach()

        if not node.is_leaf():
            try:
                taxon = ANCGENE2SP.match(node.name).group(1).replace('.', ' ')
            except AttributeError:
                raise ValueError("Can not match species name in %r [children: %s]" % \
                                   (node.name, ' '.join(c.name for c in node.children)))
            if len(node.children) <= 1:
                if hasattr(node, 'reinserted'):
                    if phyltree.ages[taxon] > 0:
                        logger.warning('DELETE reinserted node with age>0: %s',
                                        node.name)
                        ### TODO: see if I should actually keep these nodes
                        ###       (most often known ancestors, so might be
                        ###        useful)
                else:
                    logger.warning('DELETE node with single child: %s',
                                    node.name)
                node.delete(prevent_nondicotomic=False,
                            preserve_branch_length=True)


def sum_average_dNdS(dNdS, nb2id, tree_nbs):
    """
    DEPRECATED: use `bound_average` with unweighted=True and original_leading_paths=False.

    Compute the average dS depth of each node (between this node and
    leaves).
    
    Probably the roughest method. Use rather `bound_average_dS`"""
    dS, dN = {}, {}
    for n in tree_nbs.traverse('postorder'):
        node = n.name
        if n.is_leaf():
            dS[nb2id[node]] = 0
            dN[nb2id[node]] = 0
            dS[node] = 0
            dN[node] = 0
        else:

            dS_branch = dNdS[node + '..' + c.name][5]
            dN_branch = dNdS[node + '..' + c.name][4]
            dS_children = [dS[nb2id[c.name]] + dS_branch for c in n.children]
            dN_children = [dN[nb2id[c.name]] + dN_branch for c in n.children]
            #dS_children = [dS[c.name]] for c in n.children]
            #dN_children = [dN[c.name]] for c in n.children]
            average_dS_ch = sum(dS_children) / len(dS_children)
            average_dN_ch = sum(dN_children) / len(dN_children)

            dS[nb2id[node]] = average_dS_ch
            dN[nb2id[node]] = average_dN_ch
            dS[node] = average_dS_ch
            dN[node] = average_dN_ch

    return dN, dS



def del_singletons(tree):
    for node in tree.iter_descendants('postorder'):
        #print(node.name, node.children)
        if len(node.children) == 1:
            child = node.children[0]
            #parent = node.up
            for feature in ['branch_dS', 'branch_dN', 'branch_t']:
                if hasattr(child, feature):
                    setattr(child, feature,
                            getattr(child, feature) + getattr(node, feature, np.nan))
            logger.info('DELETE %r before saving', node.name)
            node.delete(prevent_nondicotomic=False, preserve_branch_length=True)
    while len(tree.children) == 1:
        tree = tree.children[0].detach()

    return tree


def save_fulltree(fulltree, opened_outfile):
    # first delete nodes with single child, because R Ape cannot read it.
    #print("Save_fulltree")
    #print(fulltree.get_ascii())
    fulltree = del_singletons(fulltree)  # Move that outside the script (unix philosophy)

    features = (set.union(*(set(n.features) for n in fulltree.traverse()))
                - set(('name', 'dist', 'support', 'basename', 'dirname', 'ndataset', 'N*dN', 'S*dS')))
    # features = ['nb', 'cal', 'type', 'branch_dist', 't', 'dS', 'dN',
    #             'age_dist', 'age_t', 'age_dS', 'age_dN'],
    if fulltree.children:
        opened_outfile.write(fulltree.write(features=features,
                                            format=1,
                                            format_root_node=True) + '\n')

def save_subtrees(subtrees, opened_outfile):
    #for stree in subtrees:
    #    opened_outfile.write(stree.write(features=['type'], format=1) + '\n')
    #opened_outfile.write('\n'.join(subtrees) + '\n')
    for stree in subtrees:
        save_fulltree(stree, opened_outfile)


savefunctions = {'ages': save_ages,
                 'fulltree': save_fulltree,
                 'subtrees': save_subtrees} #save_subtrees}


def setup_fulltree(resultfile, phyltree, replace_nwk='.mlc', replace_by='.nwk',
                   measures=['dS']):
    """Put result values on the corresponding nodes of the target tree.

    Example measures at nodes made by Beast2+TreeAnnotator from the tree log:
        - posterior
        - height height_95%_HPD height_median height_range
        - length length_95%_HPD length_median length_range
        - rate rate_95%_HPD rate_median rate_range
    """
    orig_fulltree = load_fulltree(resultfile, replace_nwk, replace_by)
    rm_erroneous_ancestors(orig_fulltree, phyltree)
    if CODEML_MEASURES.intersection(measures):
        regex_dataset = re.compile(r'^Data set \d+$')
        with open(resultfile) as mlc:
            ndataset = 1
            while True:
                logger.debug('# Dataset %d from %s', ndataset, resultfile)
                fulltree = orig_fulltree.copy()
                fulltree.add_feature('ndataset', ndataset)
                logger.debug('fulltree %r: %d nodes; %d leaves; newick: %s',
                             fulltree.name, len(list(fulltree.traverse())), len(fulltree),
                             fulltree.write(outfile=None, format=1, format_root_node=True))
                id2nb, nb2id, tree_nbs, br_tw, model = branch2nb(mlc, fulltree)
                dNdS, dStreeline, dNtreeline = get_dNdS(mlc, skiptrees=(model!=1))

                dNdS = dNdS_precise(dNdS, br_tw, id2nb, tree_nbs, dStreeline, dNtreeline)
                set_dNdS_fulltree(fulltree, id2nb, dNdS)  #set_data_fulltree()
                yield fulltree
                try:
                    eat_lines_uptomatch(mlc, regex_dataset)
                    ndataset += 1
                except LookupError:
                    logger.debug('Parsed %d results from %s', ndataset, resultfile)
                    break
    elif BEAST_MEASURES.union(('beast:dist',)).intersection(measures):
        try:
            tree = ete3.Tree(resultfile, format=1)  # consensus tree from Beast + TreeAnnotator.
        except ete3.parser.newick.NewickError as e:
            if os.path.exists(resultfile):
                e.args = ("Malformed newick tree structure in %r" % resultfile,)
            else:
                e.args = ("Unexisting tree file %r" % resultfile,)
            raise

        fulltree = orig_fulltree.copy()
        update_tree_nodes(fulltree, tree,
                          update_features=['dist' if m=='beast:dist' else m
                                           for m in measures],
                          defaults=BEAST_DEFAULTS,
                          reassigns=BEAST_REASSIGNS)
        yield fulltree


#def process_tonewick(mlcfile, phyltree, replace_nwk='.mlc', measures=['dS'], 
#                     **kwargs):
#    fulltree = setup_fulltree(mlcfile, phyltree, replace_nwk, measures)
#
#    # replace the branch distances by the chosen measure
#    measure = measures[0]
#    if measure != 'dist':
#        for node in fulltree.traverse():
#            node.add_feature('old_dist', node.dist)
#            node.dist = getattr(node, measure, np.NaN)
#
#    return fulltree

def set_subtrees_distances(subtrees, measure):
    """Set the new branch lengths to the given measure (dS, dN or t)"""
    if measure != 'dist':
        for subtree in subtrees:
            for node in subtree.traverse():
                node.add_feature("branch_dist", node.dist)
                try:
                    newvalue = getattr(node, measure)
                except AttributeError as err:
                    if node.is_root():
                        newvalue = np.NaN
                    else:
                        raise

                node.dist = 0 if np.isnan(newvalue) else newvalue


#def process_toages(resultfile, phyltree, replace_nwk='.mlc', measures=['dS'],
def process(resultfile, ensembl_version, phyltree, replace_nwk='.mlc', replace_by='.nwk',
            measures=['dS'], todate="isdup", unweighted=False,
            original_leading_paths=False, correct_unequal_calibs='default',
            fix_conflict_ages=True, keeproot=False, dataset_fmt='{basename}'):

    # Convert argument to function
    todate_funcs = {'isdup': retrieve_isdup, 'd': retrieve_isdup,
            'isinternal': isinternal, 'isint': isinternal, 'i': isinternal,
            'taxon': def_is_any_taxon, 't': def_is_any_taxon,
            'true': true, 'false': false}

    todate = combine_boolean_funcs(todate, todate_funcs)
    logger.debug('todate function: %s', todate.__name__)

    def this_get_taxon(node):
        return get_taxon(node, ensembl_version)
    
    def get_eventtype(node, subtree):
        if node.is_leaf():
            subtree[node.name]['isdup'] = False
            return 'leaf'
        #return 'dup' if isdup(node, subtree) else 'spe'
        # Uses the side effect of `isdup_cache`
        return 'dup' if isdup_cache(node, subtree) else 'spe'
    
    def is_outgroup(node):
        return getattr(node, 'is_outgroup', 0)
    
    #def get_event_type(node, subtree):
    #    return 'leaf' if node.is_leaf() else 'dup' if isdup(node, subtree) else 'spe'
    for fulltree in setup_fulltree(resultfile, phyltree, replace_nwk, replace_by, measures):

        if BEAST_MEASURES.union(('beast:dist',)).intersection(measures):
            # Handle reassignments into splitted variables
            assigned_measures = []
            for m in measures:
                if m in BEAST_REASSIGNS:
                    assigned_measures.extend((m+'_low', m+'_up'))
                elif m=='beast:dist':
                    assigned_measures.append('dist')
                else:
                    assigned_measures.append(m)

            ages = tabulate_ages_from_tree(fulltree, todate,
                                           assigned_measures,
                                           keeproot=keeproot,
                                           node_info=[('taxon', this_get_taxon),
                                                      ('is_outgroup', is_outgroup)],
                                   node_feature_setter=[('type', get_eventtype)],
                                   dataset_fmt=dataset_fmt)
            subtrees = None
        else: #if CODEML_MEASURES.intersection(measures):
            node_info = [('taxon', this_get_taxon), ('is_outgroup', is_outgroup)]
            ages, subtrees = bound_average(fulltree, phyltree.ages, todate,
                                       measures,
                                       unweighted,
                                       original_leading_paths,
                                       correct_unequal_calibs,
                                       fix_conflict_ages=fix_conflict_ages,
                                       keeproot=keeproot,
                                       calib_selecter='taxon',
                                       node_info=node_info,
                                       node_feature_setter=[('type', get_eventtype)],
                                       dataset_fmt=dataset_fmt)

        showtree(fulltree)
        if not keeproot:
            ingroup_nodes = [n for n in fulltree.children
                                if int(getattr(n, 'is_outgroup', 0)) == 0]
            basename = getattr(fulltree, 'basename', None)
            dirname = getattr(fulltree, 'dirname', None)
            ndataset = getattr(fulltree, 'ndataset', None)
            try:
                fulltree, = ingroup_nodes
                if basename is not None:
                    fulltree.add_feature('basename', basename)
                if dirname is not None:
                    fulltree.add_feature('dirname', dirname)
                if ndataset is not None:
                    fulltree.add_feature('ndataset', ndataset)
            except ValueError:
                logger.warning("Ingroup node not found (%d candidates)",
                                len(ingroup_nodes))

        yield ages, fulltree, subtrees


def main(outfile, resultfiles, ensembl_version=ENSEMBL_VERSION,
         phyltreefile=PHYLTREEFILE, measures=['t', 'dN', 'dS', 'dist'],
         unweighted=False, original_leading_paths=False,
         correct_unequal_calibs='default', fix_conflict_ages=True, verbose=0,
         show=None, replace_nwk='.mlc', replace_by='.nwk', ignore_errors=False,
         saveas='ages', todate='isdup', keeproot=False, dataset_fmt='{basename}'):
    nb_results = len(resultfiles)
    
    if verbose > 1:
        loglevel = logging.DEBUG
    elif verbose > 0:
        loglevel = logging.INFO
    else:
        loglevel = logging.WARNING
    logger.setLevel(loglevel)
    mpllogger.setLevel(loglevel)

    logger.debug("outfile  %s\n"
                  "     resultfiles      %s %s\n"
                  "     ensembl       v.%d\n"
                  "     measures      %s\n"
                  "     verbose       %s\n"
                  "     show          %s\n"
                  "     unweighted    %s\n"
                  "     original_leading_paths %s\n"
                  "     correct_unequal_calibs %s\n"
                  "     fix_conflict_ages %s\n"
                  "     replace_nwk   %s\n"
                  "     replace_by    %s\n"
                  "     ignore_errors %s\n"
                  "     keeproot      %s\n"
                  "     dataset_fmt   %r\n",
                  outfile,
                  resultfiles[:5], '...' * (nb_results>5),
                  ensembl_version,
                  measures,
                  verbose,
                  show,
                  unweighted,
                  original_leading_paths,
                  correct_unequal_calibs,
                  fix_conflict_ages,
                  replace_nwk,
                  replace_by,
                  ignore_errors,
                  keeproot,
                  dataset_fmt)

    saveas_indices = {'ages': 0, 'fulltree': 1, 'subtrees': 2}
    saveas_i = saveas_indices[saveas]
    save_result = savefunctions[saveas]

    global showtree
    showtree = def_showtree(measures, show)

    phyltree = PhylTree.PhylogeneticTree(phyltreefile.format(ensembl_version))
    
    if 'codeml' in measures:
        measures.remove('codeml')
        measures.extend(sorted(CODEML_MEASURES))
    elif 'beast' in measures:
        measures.remove('beast')
        measures.extend(sorted(BEAST_MEASURES))

    ## Globbing suggestion
    #for m in measures:
    #    if '*' in m and m not in (CODEML_MEASURES | BEAST_MEASURES):
    #        expanded_m = [em for em in (BEAST_MEASURES | CODEML_MEASURES)
    #                           if re.search(m.replace('*', '.*'), em)]
    #       if not expanded_m:
    #           raise NotFoundError()
    #   else:
    #       .append(m)
    #measures = expanded_measures

    with Stream(outfile, 'w') as out:
        if saveas == 'ages':
            # There's a pb if 'dist' is in measures -> specify 'beast:dist' or codeml is implied
            logger.debug('measures: %s', measures)
            logger.debug('CODEML_MEASURES: %s', CODEML_MEASURES)
            logger.debug('BEAST_MEASURES: %s', BEAST_MEASURES)
            if set(measures) & CODEML_MEASURES.union(('dist',)):
                stats = ('branch', 'age', 'p_clock')
                if fix_conflict_ages:
                    stats += ('fixed_age',)
                measure_outputs = ['%s_%s' % (s,m)
                                    for s in stats for m in measures]
            elif set(measures) & BEAST_MEASURES.union(('beast:dist',)):
                measure_outputs = []
                for m in measures:
                    if m in BEAST_REASSIGNS:
                        measure_outputs.extend((m+'_low', m+'_up'))
                    else:
                        measure_outputs.append(m)

            header = ['name'] + measure_outputs + \
                     ['calibrated', 'parent', 'taxon', 'is_outgroup'] + \
                     ['type', 'root', 'subgenetree']
            out.write('\t'.join(header) + '\n')

        for i, resultfile in enumerate(resultfiles, start=1):
            percentage = float(i) / nb_results * 100
            print("\r%5d/%-5d (%3.2f%%) %s" % (i, nb_results, percentage, resultfile),
                  end=' ')
            try:
                ndataset = 1
                for result in process(resultfile, ensembl_version, phyltree,
                                 replace_nwk, replace_by, measures, todate,
                                 unweighted, original_leading_paths,
                                 correct_unequal_calibs, fix_conflict_ages,
                                 keeproot=keeproot, dataset_fmt=dataset_fmt):
                    if saveas=='fulltree':
                        set_subtrees_distances([result[1]], measures[0])
                    elif saveas=='subtrees':
                        set_subtrees_distances(result[2], measures[0])

                    save_result(result[saveas_i], out)
                    ndataset += 1
            except BaseException as err:
                print()
                if not isinstance(err, KeyboardInterrupt) and ignore_errors:
                    logger.error("Skip %r #%d: %r", resultfile, ndataset, err)
                else:
                    raise
    print()


def readfromfiles(filenames):  # ~~> CLItools
    lines = []
    for filename in filenames:
        with open(filename) as f:
            lines += [line.rstrip() for line in f if not line.startswith('#')]
    return lines


if __name__=='__main__':
    logging.basicConfig(format='%(levelname)s:%(funcName)-16s:%(message)s')
    parser = argparse.ArgumentParser(description=__doc__, 
                                formatter_class=argparse.RawTextHelpFormatter)
    
    ### Input options
    gi = parser.add_argument_group('INPUT PARAMETERS')

    gi.add_argument('outfile')
    gi.add_argument('resultfiles', nargs='+', help='.mlc file from codeml/annotated consensus tree from beast')
    gi.add_argument('--fromfile', action='store_true',
                    help='Take result files from a file (commented lines omitted).')
    gi.add_argument('-e', '--ensembl-version', type=int,
                    default=ENSEMBL_VERSION,
                    help='[%(default)s]')
    gi.add_argument('-p', '--phyltreefile', default=PHYLTREEFILE,
                    help='[%(default)s]')
    gi.add_argument('-r', '--replace-nwk', default='\.mlc$',
                    help='string to be replaced by REPLACE_BY to find the'\
                         ' tree file [%(default)s]')
    gi.add_argument('-R', '--replace-by', default='.nwk',
                    help='replacement string file [%(default)s]')
    
    ### Run options
    gr = parser.add_argument_group('RUN PARAMETERS')
    
    gr.add_argument('-d', '--todate', default="isdup", metavar='TESTS',
                    help=combine_boolean_funcs.__doc__ + """
    - tests may be any of: 'isdup', 'isinternal', 'taxon'.
       The 'taxon' test must be written with taxa as arguments:
           'taxon:mytaxon1,mytaxon2'.

Example to date all internal nodes except a calibrated speciation:
'isdup | isinternal !taxon:Simiiformes'

[%(default)r]""")
    gr.add_argument('-u', '--unweighted', action='store_true', 
                    help='average weighted by the size of clusters')
    gr.add_argument('-2', '--orig-leading-path', dest='original_leading_paths',
                    action='store_true')
    gr.add_argument('-c', '--correct-uneq-calib', dest='correct_unequal_calibs',
                    choices=['default', 'mine', 'pathd8', 'ignore'], default='default',
                    help='Formula to adjust path lengths when they end at ' \
                         'different calibrations [%(default)s]')
    gr.add_argument('-F', '--nofix-conflicts', dest='fix_conflict_ages',
                    action='store_false',
                    help='Do not force children to be younger than parent.')
    gr.add_argument('-k', '--keeproot', action='store_true',
                    help="Date the nodes immediately following the root " \
                         "(not recommended).")
    #gr.add_argument('--allow-uneq-ch-age', type=float, default=0,
    #                dest='allow_unequal_children_age',
    #                help="tolerance threshold to output the Mean Path Length "\
    #                     "value of a node, when the next calibrated children "
    #                     "have different ages"\
    #                     "(Use carefully) [%(default)s].")
    gr.add_argument("-i", "--ignore-errors", action="store_true", 
                    help="On error, print the error and continue the loop.")
    
    ### Output options
    go = parser.add_argument_group('OUTPUT PARAMETERS')
    
    go.add_argument('-m', '--measures', nargs='*',
                    default=['t', 'dN', 'dS', 'dist'],
                    choices=(['codeml', 'beast', 'dist', 'beast:dist'] +
                             list(CODEML_MEASURES) + list(BEAST_MEASURES)),
                    help='Which distance measure: dist (from the newick ' \
                         'tree) or dS,dN,t,S*dS,N*dN (from codeml), or ' \
                         'length,length_median (from beast)')
    go.add_argument('-f', '--dataset-fmt', default='{basename}',
                    help=('Format for the `subgenetree` column. Available '
                          'keys: basename,dirname,ndataset (relative to '
                          'resultfile) or tree_alt_name (after replace_nwk '
                          'removal). [%(default)s]'))
    go.add_argument('-t', '--tofulltree', dest='saveas', action='store_const',
                    const='fulltree', default='ages',
                    help='Do not compute the table, but save trees in one'\
                        ' newick file with the chosen measure as distance')
    go.add_argument('-s', '--tosubtrees', dest='saveas', action='store_const',
                    const='subtrees', default='ages',
                    help='Do not compute the table, but save subtrees in '\
                        'one newick file with the chosen measure as ' \
                        'distance. "subtree" means the gene subtree '\
                        'contained between two speciations.')
    
    ### Display options
    gd = parser.add_argument_group('Display parameters')

    gd.add_argument('--show',
            help=('"gui": start the interactive ete3 tree browser\n'
                  '"notebook": create global var FULLTREE and TS for rendering\n'
                  'other value: save to file'))
    gd.add_argument('-v', '--verbose', action='count', default=0,
                    help='print progression along tree')

    args = parser.parse_args()
    dictargs = vars(args)
    if dictargs.pop('fromfile'):
        dictargs['resultfiles'] = readfromfiles(dictargs['resultfiles'])

    main(**dictargs)
