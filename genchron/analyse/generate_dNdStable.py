#!/usr/bin/env python3

"""Parse results of codeml (.mlc files) to save dN, dS values in a table"""


from __future__ import print_function

import sys
import os.path
import re
from collections import OrderedDict, defaultdict
from queue import deque
import numpy as np
#from numpy import array as a, concatenate as c
import ete3
import argparse
import logging
import scipy.stats as st

import LibsDyogen.myPhylTree as PhylTree # my custom python3 version

from genomicustools.identify import convert_gene2species
from IOtools import Stream

logger = logging.getLogger(__name__)


#from pamliped.codemlparser import mlc_parser

np.set_printoptions(formatter={"float_kind": lambda x: "%g" %x})

ENSEMBL_VERSION = 85
PHYLTREEFILE = "/users/ldog/glouvel/ws_alouis/GENOMICUS_SVN/data{0:d}/PhylTree.Ensembl.{0:d}.conf"
ANCGENE2SP = re.compile(r'([A-Z][A-Za-z0-9_.-]+)ENS')
# ~~> genomicus.my_identify?

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
    #nwkfile = resultfile.replace(replace_nwk, '.nwk')
    rootname = re.sub(re.escape(replace_nwk) + '$', '', resultfile)
    nwkfile = rootname + replace_by
    try:
        fulltree = ete3.Tree(nwkfile, format=1)
    except ete3.parser.newick.NewickError as e:
        if os.path.exists(nwkfile):
            logger.error("\nNewickError: Malformed newick tree structure in %r",
                          nwkfile)
        else:
            logger.error("\nNewickError: Unexisting tree file %r" % nwkfile)
        sys.exit(2)

    fulltree.add_feature('treename', os.path.basename(rootname))
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
                newfeatures = dict((ft[0], getattr(srcnode, ft[1], defaults.get(ft[0])))
                                     if isinstance(ft, tuple)
                                     else (ft, getattr(srcnode, ft, defaults.get(ft)))
                                     for ft in update_features)
                # TODO: convert beast string ranges
                for rawfeature, convert in reassigns.items():
                    newfeatures.update(**convert(newfeatures.pop(rawfeature), rawfeature))

                logger.debug('Adding features %r -> %r', srcnode.name, node.name)
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
            else:
                raise ValueError('regex_w (model 5) not matched before EOF')
        line = line.rstrip()

    if model == 1:
        omegas = [float(w) for w in regex_w.sub('', line).split()]
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

    # branches follow the order of the newick string.
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
            try:
                base_node = found_tips[0].up
            except IndexError as err:
                logger.warning('Node %s:%r not found in fulltree', tip, tip_id)
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
            debug_msg += "%s -> %s  " % (base_id, tip_id)
            nb2id[base] = base_id
            id2nb[base_id] = base
            # Add number in the fulltree:
            base_node.add_feature('nb', base)
            logger.debug(debug_msg + 'Ok')
        except KeyError as e:
            #if base in detached_subtrees:
            # I assume a progression from leaves to root, otherwise I will miss nodes
            logger.debug(debug_msg + 'Detached')
            #else:
            # Not found now, put it back in the queue for later
            #    branches.insert(0, (base, tip))
            #    logger.debug('KeyError')

    logger.debug(tree_nbs.get_ascii())
    logger.debug(tree_ids.get_ascii())
    logger.debug(fulltree)
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
                raise ValueError("Can not match species name in %r" % \
                                   node.name)
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


# Functions used to determine whether a node should be calibrated in `bound_average`
#@memoize
def true(node, subtree=None):
    return True

def false(node, subtree=None):
    return False

def isdup(node, subtree):
    """Checks whether node is a duplication.

    Also see `genetree_drawer.infer_gene_event`.
    """
    taxon = subtree[node.name]['taxon']
    children_taxa = set((subtree[ch.name]['taxon'] for ch in 
                         node.children))
    #logger.debug('    %r isdup %s -> %s', node.name, taxon, children_taxa)
    return len( children_taxa & set((taxon,)) ) == 1


def isdup_cache(node, subtree):
    r = isdup(node, subtree)
    subtree[node.name]['isdup'] = r
    #logger.debug("Cached 'isdup': %s at %r", r, node.name)
    return r

def retrieve_isdup(node, subtree):
    return subtree[node.name]['isdup']


def isinternal(node, subtree=None):
    """True if node is neither a leaf nor the root"""
    return not node.is_leaf() and not node.is_root()

#def is_descendant_of(node, taxon=None, subtree=None):

# Example of a check for a given speciation node. Useless because cannot calibrate dups
#def is_murinae_spe(node, taxon, subtree):
#    return taxon=='Murinae' and not isdup(node, taxon, subtree)

def def_is_target_speciation(target_sp):
    """define the speciation node check function"""
    def is_target_spe(node, subtree):
        taxon = subtree[node.name]['taxon']
        return taxon==target_sp and not isdup(node, subtree)
    return is_target_spe

def def_is_any_taxon(*target_taxa):
    """define the speciation node check function"""
    def is_any_taxon(node, subtree):
        taxon = subtree[node.name]['taxon']
        return taxon in target_taxa
    is_any_taxon.__name__ += '_' + '|'.join(target_taxa)
    return is_any_taxon


## ~> function combiner module?
def negate(testfunc):
    def negated(*args,**kwargs):
        return not testfunc(*args, **kwargs)
    negated.__name__ = '!' + testfunc.__name__
    return negated

def all_of(tests):
    def all_tests(*args, **kwargs):
        return all(test(*args, **kwargs) for test in tests)
    all_tests.__name__ = '(' + '&'.join(t.__name__ for t in tests) + ')'
    return all_tests

def any_of(tests):
    def any_test(*args, **kwargs):
        return any(test(*args, **kwargs) for test in tests)
    any_test.__name__ = '(' + '|'.join(t.__name__ for t in tests) + ')'
    return any_test


def combine_boolean_funcs(key, funcs):
    """'Parse' the key (string) argument to combine several tests together.
    
    Syntax:
        - any whitespace means AND.
        - '|' means OR. AND takes precedence over OR.
        - tests can be negated by prefixing with '!'.
    """
    ## TODO: parse parentheses for precedence.
    alt_tests = [] # Split by 'OR' operators first.
    for alt_teststr in key.split('|'):
        tests = []
        for teststr in alt_teststr.split():
            applynegate = False
            if teststr[0] == '!':
                teststr = teststr[1:]
                applynegate = True
            if ':' in teststr:
                teststr, testargs = teststr.split(':')
                testargs = testargs.split(',')
                test = funcs[teststr](*testargs)
            else:
                test = funcs[teststr]

            if applynegate:
                test = negate(test)

            tests.append(test)

        alt_tests.append(all_of(tests))

    return any_of(alt_tests)


def get_taxon(node, ensembl_version=ENSEMBL_VERSION):
    if node.is_leaf():
        try:
            return convert_gene2species(node.name, ensembl_version)
        except KeyError as err:
            logger.warning(err)
    try:
        return ANCGENE2SP.match(node.name).group(1).replace('.', ' ')
    except AttributeError:
        raise ValueError("Can not match species name in %r" % node.name)


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


def rec_average(node, subtree, tmp_m='tmp_m', child_weight='cal_leaves', allow_unequal_children_age=0):
    """Perform the averaging of one node, given its descendants (unweighted).
    
    This operation must be repeated along the tree, from the leaves to the root.

    Update the subtree dictionary (which holds the temporary measures for nodes).

    :param: `child_weight`
    """
    # normkey could also be: age.

    # Store the age of the **next calibrated** node (propagate down).
    # Since it is a duplication, this should be the same for
    # children[0] and children[1]:
    ch_ages = [subtree[ch.name]['age'] for ch in node.children]

    if max(ch_ages) - min(ch_ages) > allow_unequal_children_age:

        logger.warning("At %r: unequal children's ages (next "
                        "calibrated ages): %s", node.name, ch_ages)

        subtree[node.name]['age'] = np.NaN
        #TODO: Rescale measure rates here.

    else:
        # This 'age' will *later* be replaced (taken from (dS) path lengths)
        subtree[node.name]['age'] = ch_ages[0]


    # Array operation here: increment tmp_m with current distances.
    # One child per row.
    children_ms = np.stack([subtree[ch.name]['br_m'] + subtree[ch.name][tmp_m]
                            for ch in node.children])
    logger.debug('%r: children_ms:\n%s', node.name, children_ms)
    # Weighting:
    if child_weight:
        children_ws = np.array([subtree[ch.name][child_weight] for ch in
                                node.children])
        logger.debug('%r: children_ws:\n%s', node.name, children_ws)

        subtree[node.name][child_weight] = children_ws.sum()
    else:
        children_ws = None

    # weighted mean of each measure (across children)
    subtree[node.name][tmp_m] = np.average(children_ms, axis=0, weights=children_ws)
    #return children_ms.sum(axis=0)


def rec_average_w(node, subtree, tmp_m='tmp_m', **kw):
    """Perform the averaging of one node, given its descendants.
    Weighted (implicitly here) by the number of descendant leaves.
    
    This operation must be repeated along the tree, from the leaves to the root

    UNUSED.
    """
    # Array operation here: increment tmp_m with current measure.
    # One child per row.
    children_ms = np.stack([subtree[ch.name]['br_m'] + subtree[ch.name][tmp_m]
                            for ch in node.children])
    # mean of each measure (across children)
    subtree[node.name][tmp_m] = children_ms.mean(axis=0)
    #return children_ms.mean(axis=0)

def rec_combine_weighted_paths(node, subtree, unweighted=True):
    # For the weighted algo.
    frac = 1. / len(node.children)
    
    subtree[node.name]['tmp_m'] = np.concatenate([
                            subtree[ch.name]['br_m'] + subtree[ch.name]['tmp_m']
                            for ch in node.children])
    # 'tmp_' ages/paths/leaves
    # cal_ages is the age of reference (calibrated), so one value for all measures.
    subtree[node.name]['cal_ages'] = np.concatenate([
                                            subtree[ch.name]['cal_ages']
                                            for ch in node.children])
    # cal_m
    subtree[node.name]['cal_paths'] = np.concatenate([
                                            subtree[ch.name]['cal_paths']
                                            for ch in node.children])  # 'cal_' + tmp_m
    subtree[node.name]['cal_leaves'] = np.concatenate([
                                            subtree[ch.name]['cal_leaves']
                                            for ch in node.children])
    # For the computation of confidence intervals (not good with internal calibrations)
    children_leaves = np.array([subtree[ch.name]['leaves'] for ch in node.children])[:,None]
    #if path_frac:
    #    children_leaves *= np.array([
    children_sds = np.array([subtree[ch.name]['sd'] for ch in node.children])
    children_brlen = np.array([subtree[ch.name]['br_m'] for ch in node.children])
    subtree[node.name]['leaves'] = children_leaves.sum()

    subtree[node.name]['w'] = np.concatenate([frac * subtree[ch.name]['w']
                                              for ch in node.children])
    #children_leaves *= 
    if unweighted:
        #logger.debug('%r update sd: children_sds=%s; children_brlen=%s',
        #             node.name, ls(children_sds), ls(children_brlen))
        subtree[node.name]['sd'] = np.sqrt(((children_sds**2).sum(axis=0) +
                                            (children_brlen).sum(axis=0)
                                           )) * frac
    else:
        #logger.debug('%r update sd: children_sds=%s; children_brlen=%s; children_leaves=%s',
        #             node.name, ls(children_sds), ls(children_brlen), ls(children_leaves))
        subtree[node.name]['sd'] = np.sqrt(
                                    (((children_sds*children_leaves)**2).sum()
                                     + (children_brlen*children_leaves**2).sum(axis=0)
                                    ) / sum(children_leaves)**2)

    # Shape checks
    expected_shape = (sum(subtree[ch.name]['tmp_m'].shape[0] for ch in node.children),
                      subtree[node.children[0].name]['tmp_m'].shape[1])
    assert subtree[node.name]['tmp_m'].shape == expected_shape, subtree[node.name]['tmp_m'].shape
    assert subtree[node.name]['cal_paths'].shape == expected_shape, subtree[node.name]['cal_paths'].shape

    assert subtree[node.name]['cal_ages'].shape == (expected_shape[0],), subtree[node.name]['cal_ages'].shape
    assert subtree[node.name]['sd'].shape == (expected_shape[1],)
    assert subtree[node.name]['w'].shape == (expected_shape[0],), subtree[node.name]['w'].shape


def old_timescale_paths(previous_node, previous_path, nodename, subtree,
                        tmp_m='tmp_m', weight_key=None, original_leading_paths=False):
    previous_age = subtree[previous_node]['age']

    if np.isnan(subtree[nodename]['age']):
        pass  #TODO

    m = subtree[nodename][tmp_m]  # descendant mean path length
    dist_to_calib = previous_age - subtree[nodename]['age']
    return previous_age - (1 - m/previous_path) * dist_to_calib


# Then root-to-leaf computations (preorder)
def timescale_paths(previous_node, previous_path, nodename, subtree,
                    unweighted=True, correct_unequal_calibs='default'):
    previous_age = subtree[previous_node]['age']
    cal_ages = subtree[nodename]['cal_ages'][:,None]  # transposing the vector.
    cal_paths = subtree[nodename]['cal_paths']  # 'cal_' + tmp_m
    paths = subtree[nodename]['tmp_m']
    
    logger.debug('        previous_path=%s cal_paths=%s', ls(previous_path), ls(cal_paths))
    
    assert (cal_ages.shape[1] == 1) and (paths.shape[1] == cal_paths.shape[1])
    assert cal_ages.shape[0] == paths.shape[0] == cal_paths.shape[0]
    
    cal_leaves = subtree[nodename]['cal_leaves']
    weights = cal_leaves if unweighted else subtree[nodename]['w']
    meanpath = np.average(paths, axis=0, weights=weights)

    logger.debug('        Expect broadcastable shapes: cal_ages %s; paths %s; '
                 'cal_paths %s; meanpath %s',
                 cal_ages.shape, paths.shape, cal_paths.shape, meanpath.shape)

    # Correction when the path end is at an age>0
    if correct_unequal_calibs=='default':
        ages = cal_ages + (previous_age-cal_ages) * (meanpath-cal_paths) / (previous_path - cal_paths)
        #if (previous_path - cal_paths < 0).any():
        #if (paths-cal_paths < 0).any():
    elif correct_unequal_calibs=='mine':
        durations = previous_age - cal_ages
        corrections = durations.max() / durations  # Or rescale with the min ?
        ages = (previous_age - cal_ages.min()) * (meanpath - cal_paths) * corrections
        ages /= (previous_path - cal_paths)
    elif correct_unequal_calibs=='pathd8':
        # Classical MPL. See PATHd8 (Britton et al. 2007, appendix 1)
        # But it's theoretically wrong (non homogeneous formula)
        ages = (cal_ages + (previous_age-cal_ages) * (meanpath-cal_paths)) / (previous_path - cal_paths)
    elif correct_unequal_calibs=='ignore':
        ages = previous_age * paths / previous_path
    elif correct_unequal_calibs=='rates':
        rate = previous_path/previous_age
        cal_rates = np.where(cal_ages>0, cal_paths/cal_ages, rate)

        # previous_paths must be rescaled as well.
        ages = previous_age * (paths - cal_paths + cal_paths*rate/cal_rates) / previous_path
    else:
        raise ValueError('`correct_unequal_calibs` should be default/mine/pathd8/ignore/rates.')

    logger.debug('    -> Inputs: mpl=%s weights=%s ages=%s',
                 meanpath, ls(weights), ls(ages))
    age = np.average(ages, axis=0, weights=weights)  # Potentially applying twice the weights.

    assert age.shape == (paths.shape[1],), "Shapes: age %s; paths %s" %(
            age.shape, paths.shape)

    return age


# What I should do:  # NOPE.
# class tree_store(object):  (squirrel hide)
#     """Store required variables for each node, for the dynamic prog algorith"""
#
#    def run_rootward(self, tree, root=None):
#       """Iterate from the current leaves to the next 'init_at' point"""
#       for node, children in rev_dfw_descendants(tree, get_children,
#                                                 queue=[root]):
#           for ch in children:
#           
#               if self.init_at(ch):
#                    self.init(ch)
#               
#               self.update_rootward(node, children)
#           if self.init_at(node):
#               break
#
#    def run_leafward(self, tree, root=None):
#       """Iterate from the current leaves to the next 'init_at' point"""
#       for node, children in rev_dfw_descendants(tree, get_children,
#                                                 queue=[root]):
#
#           self.update_leafward(node, children)
#
#
#    def run_leafward(self, tree, root=None):
#        """Iterate from the current node to the next 'init_at' leaves"""
#
#    def __init__(self, order=('r', 'l'), init_at, stop_at, initvars,
#                 update_rootward, update_leafward):
#       self.run_first = self.run_rootward if order[0] == 'r' else run_leafward
#       self.run_then = self.run_rootward if order[1] == 'r' else run_leafward
#       
#
#    def run_all(self, tree):
#        self.store = {}
#
#        running = True
#        for next_stop in self.run_first(tree):
#            for node in self.run_then(tree, next_stop):
#                self.emit(node)
#
# update_rootward(**kw)

##TODO:
def clock_test(node, subtree, measures=['dS'], unweighted=True): #weight_key=None):
    """Assuming total number of substitution > 30 for normal approx"""
    weight_key = 'cal_leaves' if unweighted else 'w'
    children_ws = [subtree[ch.name][weight_key] for ch in node.children]

    children_vars = np.array([subtree[ch.name]['sd']**2 + subtree[ch.name]['br_m'][0]
                             for ch in node.children])
    assert children_vars.shape == (len(node.children), len(measures)), children_vars.shape

    sd_d = np.sqrt(children_vars.sum(axis=0))

    children_ms = [subtree[ch.name]['tmp_m']+subtree[ch.name]['br_m'] \
                   for ch in node.children]
    children_ms = np.array([np.average(ch_m, axis=0, weights=ch_w) for ch_m,ch_w
                            in zip(children_ms, children_ws)])
    assert children_ms.shape == (len(node.children), len(measures)), children_ms.shape
    assert sd_d.shape == (len(measures),)

    logger.debug('    %r clock test with children_ms=%s; sd_d=%s', node.name,
                 ls(children_ms), ls(sd_d))

    # Since the process is assumed to be Poisson, variance is equal to mean.
    # standard deviation of the difference:
    #sds = np.sqrt(children_vars.sum(axis=0))

    if len(node.children)==2:
        d = children_ms[0] - children_ms[1]
        # Proba of being more extreme: 2-tailed.
        p = 2 * st.norm.sf(np.abs(d), 0, sd_d)
        # Note: it probably breaks the t-test assumption that group variances are equal...
    else:
        # Do ANOVA (We don't know the degrees of freedom).
        # Chi2 in the PATHd8 publication?
        #raise NotImplementedError('Multifurcation Clock test')

        children_subst = children_ms.copy()
        #weighting
        root = node.get_tree_root()
        for i,m in enumerate(measures):
            if m in ('dN', 'dS'):
                children_subst[:,i] *= getattr(root, m[1])
            elif m in ('t', 'dist'):
                try:
                    children_subst[:,i] *= (getattr(root, 'N') + getattr(root, 'S'))
                    if m=='dist':
                        # Because 'dist' is in unit/nucleotide
                        children_subst[:,i] *= 3
                except AttributeError:
                    pass
        p = np.array([st.chisquare(x) for x in children_subst.T])[:,1]

    return p



def bound_average(fulltree, calibration, todate=isdup,
                  measures=['dS'],  # should default to 'dist'
                  unweighted=False,
                  original_leading_paths=False,
                  correct_unequal_calibs='default',  # 'mine','pathd8','ignore'
                  fix_conflict_ages=True,
                  keeproot=False,
                  calib_selecter=None, node_info=None, node_feature_setter=None):

    """
    Normalize duplication node position between speciation nodes.

     Scaling:
     --------
       * d / (D' + d) with D' modified (Mean-Path-Length)
       * original_leading_paths: d / (D' + d) with D' as the original measure (supporting branch).
       - d: average (WPGMA/UPGMA-like) of branch length leading to next
            speciation, in units of 'measure'.
            'measure' must be an attribute present in fulltree ('dS', 'dist').
       - D': branch length from the previous speciation to this duplication.

    Arguments:
    ----------
      - calibration: dictionary of `calib_selecter` -> age
      - todate: function returning True when the node should be dated.
                Takes 2 params: (node, subtree).
                By default: check whether node is a duplication.
      - calib_selecter: the entries of `calibration`. Typically `taxon`.
                   This is a node field to search, or the attr from `node_attr_getter`.
      - node_info: paired list of functions retrieving additional info
                    (ex: ("taxon", get_taxon). Tuple of length 1 will just use `getattr`)
      - node_feature_setter: can be used to set new features, possibly based on other parts
                 of the tree (arguments: node, subtree)
    """
    # BEFORE 12/06/19
    #rec_rootward = rec_average
    #weight_key = 'cal_leaves' if unweighted else None  #/'leaves'
    #timescale_paths = old_timescale_paths

    reset_at_calib = False

    rec_rootward = rec_combine_weighted_paths
    calibration[None] = np.NaN  # When select_calib_id returns None

    node_info = [] if node_info is None \
                else [(getinfo[0], lambda node: getattr(node, attr, None))
                      if len(getinfo)==1 else getinfo
                      for getinfo in node_info]
    if node_feature_setter is None:
        node_feature_setter = []

    ## Future outputs
    ages = []
    # list of ete3 subtrees of fulltree (between consecutive calibrations)
    subtrees = []

    subtree = {} # temporary subtree while traversing from one calibration to
                 # another

    if calib_selecter is None:
        def select_calib_id(node):
            return node
    else:
        def select_calib_id(node):
            try:
                return getattr(node, calib_selecter)
            except AttributeError:
                return subtree[node.name][calib_selecter]

    n_measures = len(measures)
    # initial measure while passing a speciation for the first time
    measures_zeros = np.zeros((1, n_measures))
    is_next_cal = lambda node: (not node.is_root() and node.cal)

    for node in fulltree.traverse('postorder'):
        scname = node.name
        debug_msg = "* %s: " % scname
        ### DISCARD all dup just after the root speciation/duplication.
        if node.is_root() and not keeproot:
            logger.debug(debug_msg + "Root (discard)")
            continue

        subtree[scname] = {attr: getinfo(node)
                           for attr,getinfo in node_info}
        node.add_features(**{ft: setft(node, subtree) for ft,setft in node_feature_setter})

        #try:
        subtree[scname]['br_m'] = branch_measures = \
                np.array([getattr(node, m, np.NaN) for m in measures])
        #except AttributeError as err:
        #    if not node.is_root():
        #        err.args += ("Node: %r (is_root: %s)" % (node, node.is_root()),)
        #        raise
        if node.is_leaf():
            node.add_feature('cal' ,1)
            
            leaf_age = calibration.get(select_calib_id(node), 0)

            subtree[scname].update({
                               'tmp_m': measures_zeros,
                               #'total_m': measures_zeros,
                               'age': leaf_age,
                               'w': np.array([1]),   # Weights (path fraction)
                               'sd': np.zeros(n_measures),
                               'cal_ages': np.array([leaf_age]),
                               'cal_paths': measures_zeros,
                               'cal_leaves': np.array([1]), # number of leaves below the next calibrated nodes
                               'leaves': 1}) # number of descendant leaves
            #ages[scname] = 0
            ages.append([scname] + branch_measures.tolist() +
                        [leaf_age] * n_measures +
                        [np.NaN] * n_measures +
                        [1,
                         getattr(node.up, 'name', None)] +
                        [subtree[scname][attr] for attr,_ in node_info] +
                        [getattr(node, ft) for ft,_ in node_feature_setter] +
                        [fulltree.name, getattr(fulltree, 'treename', '')])
            logger.debug(debug_msg + "Leaf")
        else:
            #try:
            #    assert len(node.children) > 1
            #except AssertionError as err:
            #    if node.is_root():
            #        pass
            #    else:
            #        raise

            # Compute the temporary measure of the node.
            rec_rootward(node, subtree, unweighted=unweighted)
            #rec_rootward(node, subtree, 'total_m')
            #rec_average(node, subtree, measures, weight_key='ncal')

            subtree[scname]['p_clock'] = clock_test(node, subtree, measures, unweighted)

            if todate(node, subtree):
                # it is uncalibrated:
                node.add_feature('cal', 0)
                logger.debug(debug_msg + "Uncal.; m=%s; sd=%s; cal_leaves=%s",
                             ls(subtree[scname]['tmp_m']),
                             ls(subtree[scname]['sd']),
                             ls(subtree[scname]['cal_leaves']))
            else:
                # It is calibrated.
                node.add_feature('cal', 1)
                logger.debug(debug_msg + "Calibrated. sd=%s", ls(subtree[scname]['sd']))

                # store the age of this taxon
                node_age = subtree[scname]['age'] = calibration[select_calib_id(node)]
                discarded_nodes = {scname: subtree[scname]}

                ages.append([scname] + branch_measures.tolist() +
                            [node_age] * n_measures +
                            list(subtree[scname]['p_clock'].flat) +
                            [1, # calibrated
                             getattr(node.up, 'name', None)] +
                            [subtree[scname][attr] for attr,_ in node_info] +
                            [getattr(node, ft) for ft,_ in node_feature_setter] +
                            [fulltree.name, getattr(fulltree, 'treename', '')])

                node.add_features(**{'age_'+m: node_age for m in measures})

                # Climb up tree until next speciation and assign an age to
                # each duplication met on the way

                # Path weights.
                scaling_weights = subtree[scname]['cal_leaves'] if unweighted else subtree[scname]['w']
                next_tmp_m = np.average(subtree[scname]['tmp_m'], axis=0, weights=scaling_weights)
                #scaling_m = subtree[scname]['tmp_m'] * dist_to_calib/mean_dist_to_calib
                if original_leading_paths:
                    # Defined at each node
                    scaling_m = None
                else:
                    scaling_m = next_tmp_m

                logger.debug("    climb up to next calibrations: "
                             "scaling_m=%s m=%s sd=%s w=%s leaves=%s "
                             "cal_leaves=%s",
                             ls(scaling_m),
                             ls(subtree[scname]['tmp_m']),
                             ls(subtree[scname]['sd']),
                             None if unweighted else ls(subtree[scname]['w']),
                             subtree[scname]['leaves'],
                             ls(subtree[scname]['cal_leaves']))
                ### This is where version _2 is different
                # walk descendants until speciation.
                # need to get the age of next speciation and compute the
                # time between the two speciation.

                # Dynamic list during the climbing:
                # store (node, measures for this node's branch)
                # TODO: deque
                nextnodes = deque([(ch,subtree[ch.name]['br_m'])
                                   for ch in node.children])
                while nextnodes:
                    nextnode, next_path_m = nextnodes.popleft()

                    try:
                        nextnode_m = subtree[nextnode.name]['tmp_m']
                    except KeyError as err:
                        err.args += ("Error: Node exists twice in the "
                          "tree. You may need to rerun `prune2family.py`",)
                        raise
                    logger.debug("    - %s: measure(to calib)=%s",
                                 nextnode.name, ls(nextnode_m))
                    # or: any(nextnode_m > measures_zeros) ?
                    #if nextnode_m is measures_zeros:

                    if not nextnode.cal:
                        if original_leading_paths:
                            logger.debug("          measure(from calib)=%s",
                                         ls(next_path_m))
                            scaling_weights = (subtree[nextnode.name]['cal_leaves']
                                               if unweighted else
                                               subtree[nextnode.name]['w'])
                            scaling_m = np.average(next_path_m + nextnode_m,
                                                   axis=0,
                                                   weights=scaling_weights)

                        logger.debug('        Previous shapes: previous_paths %s; scaling_weights %s',
                                     scaling_m.shape,
                                     scaling_weights.shape)

                        if (scaling_m == 0).any():
                            logger.warning("Scaling measure = %s (cannot divide) at %r",
                                           ls(scaling_m), nextnode.name)

                        age = timescale_paths(scname, scaling_m,
                                              nextnode.name, subtree,
                                              unweighted,
                                              correct_unequal_calibs)
                        parent_age = discarded_nodes[nextnode.up.name]['age']
                        if fix_conflict_ages:
                            child_age = subtree[nextnode.name]['cal_ages'].max()
                            too_recent = (child_age > age)
                            too_old = (age > parent_age)
                            if isinstance(parent_age, (int, float)):
                                age[too_old] = parent_age
                            else:
                                age[too_old] = parent_age[too_old]
                            age[too_recent] = child_age

                            for i,m in enumerate(measures):
                                if too_recent[i]:
                                    node.add_feature('fixed_age_%s' % m, -1)
                                if too_old[i]:
                                    node.add_feature('fixed_age_%s' % m, 1)

                        subtree[nextnode.name]['age'] = age
                        logger.debug("    -> Shape of age %s", age.shape)

                        nextnode_measures = list(subtree[nextnode.name]['br_m'].flat)
                        ages.append([nextnode.name] + nextnode_measures +
                                    list(age.flat) +
                                    list(subtree[nextnode.name]['p_clock'].flat) +
                                    [0, #dated
                                     nextnode.up.name] +
                                    [subtree[nextnode.name][attr]
                                        for attr,_ in node_info] +
                                    [getattr(nextnode, ft)
                                        for ft,_ in node_feature_setter] +
                                    [fulltree.name,
                                     getattr(fulltree, 'treename', '')])
                        for i, m in enumerate(measures):
                            nextnode.add_feature('age_'+m, age[i])
                        # When climbing up to next spe, need to increment
                        # the lengths 'next_path_m' from the previous spe
                        nextnodes.extend(
                                (nch, subtree[nch.name]['br_m'] + next_path_m)
                                for nch in nextnode.children)
                    else:
                        logger.debug('        calibrated.')
                    discarded_nodes[nextnode.name] = subtree.pop(nextnode.name)

                # copy the node so that it becomes a rooted tree
                nodecopy = node.copy()  # Could be avoided?
                cal_nodes = nodecopy.get_leaves(is_leaf_fn=is_next_cal)
                # No need to save subtrees that only contain calibrated nodes.
                logger.debug('=> next cal: %s',
                             ' '.join(n.name for n in cal_nodes))
                if set(cal_nodes) != set(nodecopy.children):
                    logger.debug('=> subtree %s', nodecopy.name)
                    for clf in cal_nodes:
                        logger.debug('=> clf children: %s',
                                ' '.join(c.name for c in clf.children))
                        for clfch in clf.get_children():  # Copy!!
                            logger.debug('=> detach %s', clfch.name)
                            clfch.detach()
                
                    subtrees.append(nodecopy)

                next_tmp_m = np.array([next_tmp_m])
                logger.debug('    After calibrated node: next_tmp_m=%s',
                             ls(next_tmp_m))
                subtree[scname].update(
                            cal_ages=np.array([node_age]),
                            tmp_m=next_tmp_m,
                            cal_paths=next_tmp_m,
                            cal_leaves=np.array([subtree[scname]['cal_leaves'].sum()]),
                            w=np.array([1]))
                # then reset measure (dS, dist) to zero
                if reset_at_calib:
                    subtree[scname].update(tmp_m=measures_zeros, # and no need for total_m
                                           cal_paths=measures_zeros,
                                           leaves=1,
                                           cal_leaves=np.array([1]))
    #logger.debug(subtree)
    return ages, subtrees


def tabulate_ages_from_tree(fulltree, todate=true,
                            measures=['height'],
                            keeproot=True,
                            node_info=None,
                            node_feature_setter=None):

    node_info = [] if node_info is None \
                else [(getinfo[0], lambda node: getattr(node, attr, None))
                      if len(getinfo)==1 else getinfo
                      for getinfo in node_info]
    if node_feature_setter is None:
        node_feature_setter = []

    # Future outputs
    ages = []
    
    # Some cached info
    subtree = {}

    for node in fulltree.traverse('postorder'):
        debug_msg = "* %s: " % node.name
        if node.is_root() and not keeproot:
            logger.debug(debug_msg + "Root (discard)")
            continue

        subtree[node.name] = {attr: getinfo(node) for attr,getinfo in node_info}

        new_node_features = {ft: setft(node, subtree) for ft,setft in node_feature_setter}
        subtree[node.name].update(new_node_features)
        node.add_features(**new_node_features)
        branch_measures = [getattr(node, m, np.NaN) for m in measures]

        logger.debug('subtree[%r] = %s', node.name, subtree[node.name])
        cal = 0 if todate(node, subtree) else 1

        node.add_feature('cal', cal)

        ages.append([node.name] + branch_measures +
                    [cal,
                     getattr(node.up, 'name', None)] +
                    [subtree[node.name][attr] for attr,_ in node_info] +
                    [getattr(node, ft) for ft,_ in node_feature_setter] +
                    [fulltree.name, getattr(fulltree, 'treename', '')])
    return ages


def save_ages(ages, opened_outfile):
    for row in ages:
        opened_outfile.write("\t".join(str(elem) for elem in row)
        #                                if not np.isnan(elem) else '')
                             + "\n")

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
                - set(('name', 'dist', 'support', 'treename', 'N*dN', 'S*dS')))
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
    fulltree = load_fulltree(resultfile, replace_nwk, replace_by)
    rm_erroneous_ancestors(fulltree, phyltree)
    if CODEML_MEASURES.intersection(measures):
        with open(resultfile) as mlc:
            id2nb, nb2id, tree_nbs, br_tw, model = branch2nb(mlc, fulltree)
            dNdS, dStreeline, dNtreeline = get_dNdS(mlc, skiptrees=(model!=1))

        dNdS = dNdS_precise(dNdS, br_tw, id2nb, tree_nbs, dStreeline, dNtreeline)
        set_dNdS_fulltree(fulltree, id2nb, dNdS)  #set_data_fulltree()
    elif BEAST_MEASURES.union(('beast:dist',)).intersection(measures):
        try:
            tree = ete3.Tree(resultfile)  # consensus tree from Beast + TreeAnnotator.
        except ete3.parser.newick.NewickError as e:
            if os.path.exists(resultfile):
                e.args = ("Malformed newick tree structure in %r" % resultfile,)
            else:
                e.args = ("Unexisting tree file %r" % resultfile,)
            raise

        update_tree_nodes(fulltree, tree,
                          update_features=['dist' if m=='beast:dist' else m
                                           for m in measures],
                          defaults=BEAST_DEFAULTS,
                          reassigns=BEAST_REASSIGNS)
    return fulltree


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
            fix_conflict_ages=True, keeproot=False):

    fulltree = setup_fulltree(resultfile, phyltree, replace_nwk, replace_by, measures)

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
                               node_feature_setter=[('type', get_eventtype)])
        subtrees = None
    else: #if CODEML_MEASURES.intersection(measures):
        ages, subtrees = bound_average(fulltree, phyltree.ages, todate,
                                   measures,
                                   unweighted,
                                   original_leading_paths,
                                   correct_unequal_calibs,
                                   fix_conflict_ages=fix_conflict_ages,
                                   keeproot=keeproot,
                                   calib_selecter='taxon',
                                   node_info=[('taxon', this_get_taxon),
                                              ('is_outgroup', is_outgroup)],
                                   node_feature_setter=[('type', get_eventtype)])

    showtree(fulltree)
    if not keeproot:
        ingroup_nodes = [n for n in fulltree.children
                            if int(getattr(n, 'is_outgroup', 0)) == 0]
        treename = getattr(fulltree, 'treename', None)
        try:
            fulltree, = ingroup_nodes
            if treename:
                fulltree.add_feature('treename', treename)
        except ValueError:
            logger.warning("Ingroup node not found (%d candidates)",
                            len(ingroup_nodes))

    return ages, fulltree, subtrees


def main(outfile, resultfiles, ensembl_version=ENSEMBL_VERSION,
         phyltreefile=PHYLTREEFILE, measures=['t', 'dN', 'dS', 'dist'],
         unweighted=False, original_leading_paths=False,
         correct_unequal_calibs='default', fix_conflict_ages=True, verbose=False,
         show=None, replace_nwk='.mlc', replace_by='.nwk', ignore_errors=False,
         saveas='ages', todate='isdup', keeproot=False):
    nb_results = len(resultfiles)
    
    loglevel = logging.DEBUG if verbose else logging.WARNING
    logger.setLevel(loglevel)

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
                  "     keeproot      %s\n",
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
                  keeproot)

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
            if set(measures) & CODEML_MEASURES:
                measure_outputs = ['%s_%s' % (s,m)
                                    for s in ('branch', 'age', 'p_clock')
                                      for m in measures]
            elif set(measures) & BEAST_MEASURES.union(('beast:dist',)):
                measure_outputs = []
                for m in measures:
                    if m in BEAST_REASSIGNS:
                        measure_outputs.extend((m+'_low', m+'up'))
                    else:
                        measure_outputs.append(m)

            header = ['name'] + measure_outputs + \
                     ['calibrated', 'parent', 'taxon', 'is_outgroup', 'type',
                      'root', 'subgenetree']
            out.write('\t'.join(header) + '\n')

        for i, resultfile in enumerate(resultfiles, start=1):
            percentage = float(i) / nb_results * 100
            print("\r%5d/%-5d (%3.2f%%) %s" % (i, nb_results, percentage, resultfile),
                  end=' ')
            try:
                result = process(resultfile, ensembl_version, phyltree,
                                 replace_nwk, replace_by, measures, todate,
                                 unweighted, original_leading_paths,
                                 correct_unequal_calibs, fix_conflict_ages,
                                 keeproot=keeproot)
                if saveas=='fulltree':
                    set_subtrees_distances([result[1]], measures[0])
                elif saveas=='subtrees':
                    set_subtrees_distances(result[2], measures[0])

                save_result(result[saveas_i], out)
            except BaseException as err:
                print()
                if not isinstance(err, KeyboardInterrupt) and ignore_errors:
                    logger.error("Skip %r: %r", resultfile, err)
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
    logging.basicConfig(format=logging.BASIC_FORMAT)
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
                    help='Do not constrain estimates within calibration boundaries.')
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
    gd.add_argument('-v', '--verbose', action='store_true', 
                    help='print progression along tree')

    args = parser.parse_args()
    dictargs = vars(args)
    if dictargs.pop('fromfile'):
        dictargs['resultfiles'] = readfromfiles(dictargs['resultfiles'])

    main(**dictargs)
