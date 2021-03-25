#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from sys import stdout, stderr
import re
import os.path as op
import argparse as ap
import csv
from queue import deque
import numpy as np
import scipy.stats as st

from dendro.parsers import parserchoice
from dendro.converters import converterchoice
from genomicustools.identify import convert_gene2species, GENE2SP

import logging
logger = logging.getLogger(__name__)


ENSEMBL_VERSION = max(GENE2SP.keys())
ANCGENE2SP = re.compile(r'([A-Z][A-Za-z0-9 _.-]+?)(ENS|$)')  #WARNING, introducing space here
#TODO: a test suite for each time new labelling format are introduced, to check backward compat.

def ls(array):
    return str(array).replace('\n', '')

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
                  calib_selecter=None, node_info=None, node_feature_setter=None,
                  dataset_fmt='{basename}'):

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
                   This is a node field to search, or the attr from `node_info`.
      - node_info: paired list of functions retrieving additional info
                    (ex: ("taxon", get_taxon). Tuple of length 1 will just use `getattr`)
                    (a first time visiting the node, i.e. pre-processing)
      - node_feature_setter: can be used to set new features, possibly based on other parts
                 of the tree (arguments: node, subtree)
                 (post-processing)
    """
    # BEFORE 12/06/19
    #rec_rootward = rec_average
    #weight_key = 'cal_leaves' if unweighted else None  #/'leaves'
    #timescale_paths = old_timescale_paths

    reset_at_calib = False

    rec_rootward = rec_combine_weighted_paths
    calibration[None] = np.NaN  # When select_calib_id returns None

    node_info = [] if node_info is None \
                else [(getinfo[0], lambda node: getattr(node, getinfo[0], None))
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
            node.add_feature('cal', 1)
            
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
                               'leaves': 1}                 # number of descendant leaves
                               )
            #ages[scname] = 0
            ages.append([scname] + branch_measures.tolist() +
                        [leaf_age] * n_measures +
                        [np.NaN] * n_measures +
                        ([""] * n_measures if fix_conflict_ages else []) +
                        [1,
                         getattr(node.up, 'name', None)] +
                        [subtree[scname][attr] for attr,_ in node_info] +
                        [getattr(node, ft) for ft,_ in node_feature_setter] +
                        [fulltree.name, dataset_fmt.format(**vars(fulltree))])
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
                            ([""] * n_measures if fix_conflict_ages else []) +
                            [1, # calibrated
                             getattr(node.up, 'name', None)] +
                            [subtree[scname][attr] for attr,_ in node_info] +
                            [getattr(node, ft) for ft,_ in node_feature_setter] +
                            [fulltree.name, dataset_fmt.format(**vars(fulltree))])

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
                                    nextnode.add_feature('fixed_age_%s' % m, -1)
                                if too_old[i]:
                                    nextnode.add_feature('fixed_age_%s' % m, 1)
                            subtree[nextnode.name]['fixed_age'] = (-1)*too_recent + 1*too_old

                        subtree[nextnode.name]['age'] = age
                        logger.debug("    -> Shape of age %s", age.shape)

                        nextnode_measures = list(subtree[nextnode.name]['br_m'].flat)
                        ages.append([nextnode.name] + nextnode_measures +
                                    list(age.flat) +
                                    list(subtree[nextnode.name]['p_clock'].flat) +
                                    list(subtree[nextnode.name].get('fixed_age', [])) +
                                    [0, #dated
                                     nextnode.up.name] +
                                    [subtree[nextnode.name][attr]
                                        for attr,_ in node_info] +
                                    [getattr(nextnode, ft)
                                        for ft,_ in node_feature_setter] +
                                    [fulltree.name,
                                     dataset_fmt.format(**vars(fulltree))])
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
                        logger.debug('=> Calibrated leaf %r: children: %s',
                                clf.name,
                                ' '.join('%r' % c.name for c in clf.children))
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
                            node_feature_setter=None, dataset_fmt='{basename}'):

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
                    [fulltree.name, dataset_fmt.format(**vars(fulltree))])
    return ages


def get_taxon(node, ensembl_version=ENSEMBL_VERSION):
    if node.is_leaf():
        try:
            return convert_gene2species(node.name, ensembl_version)
        except KeyError as err:
            logger.warning(','.join(err.args))
    try:
        return ANCGENE2SP.match(node.name).group(1).replace('.', ' ')
    except AttributeError:
        raise ValueError("Can not match species name in %r" % node.name)


def save_ages(ages, opened_outfile):
    for row in ages:
        opened_outfile.write("\t".join(str(elem) for elem in row)
        #                                if not np.isnan(elem) else '')
                             + "\n")


def save_mpl_tree(tree, anc_ages, opened_outfile):
    #TODO: understand why the age_* features do not appear.
    for node in tree.iter_descendants():
        node.dist = anc_ages[node.up.name] - anc_ages[node.name]
    opened_outfile.write(tree.write(format=1, format_root_node=True) + '\n')


def main():
    logging.basicConfig(format='%(levelname)s:%(funcName)-16s %(message)s')
    parser = ap.ArgumentParser(description=__doc__)
    parser.add_argument('treefile')
    parser.add_argument('calibfile', #nargs='?',
                        help='column file')
    #parser.add_argument('-p', '--phyltreefile')
    parser.add_argument('-p', '--parser', default='ete3_f1')
    parser.add_argument('-D', '--delimiter', default='\t',
                        help='Column separator in the calibfile [%(default)r]')
    parser.add_argument('-a', '--ages-column', default=1,
                        help='Column containing the ages [%(default)r]')
    parser.add_argument('-H', '--header-row', action='store_true',
                        help='Whether to discard one header row.')
    parser.add_argument('-m', '--measures', default='dist',
                        help='comma separated list [%(default)s]')
    parser.add_argument('-e', '--ensembl-version', type=int, default=ENSEMBL_VERSION)
    parser.add_argument('-t', '--totree', action='store_true',
                        help='Do not compute the table, but save trees in one'\
                            ' newick file with MPL lengths as new distances.')

    gr = parser.add_argument_group('RUN PARAMETERS')
    
    gr.add_argument('-d', '--todate', default="isinternal", metavar='TESTS',
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
    parser.add_argument('-f', '--dataset-fmt', default='{basename}',
                        help=('Format for the `subgenetree` column. Available '
                          'keys: basename,dirname (of the treefile). [%(default)s]'))
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='print progression along the tree')

    args = parser.parse_args()
    if args.verbose:
        logger.setLevel(logging.DEBUG)

    measures = args.measures.split(',')
    ensembl_version = args.ensembl_version

    #NOTE: copied from genchron.analyse.generate_dNdStable.process()
    # Convert argument to function
    todate_funcs = {'isdup': retrieve_isdup, 'd': retrieve_isdup,
            'isinternal': isinternal, 'isint': isinternal, 'i': isinternal,
            'taxon': def_is_any_taxon, 't': def_is_any_taxon,
            'true': true, 'false': false}

    todate = combine_boolean_funcs(args.todate, todate_funcs)
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
    
    #NOTE: here starts genchron.MPL custom code.
    calib = {}
    with open(args.calibfile) as cf:
        csvreader = csv.reader(cf, delimiter=args.delimiter)
        if args.header_row:
            calib_header = next(csvreader)
        else:
            calib_header = list(range(50))
        try:
            col = int(args.ages_column)
        except ValueError:
            col = calib_header.index(args.args_column)
        namecol = 1 if (col==0) else 0
        
        for row in csvreader:
            calib[row[namecol]] = float(row[col])

    stats = ('branch', 'age', 'p_clock')
    if args.fix_conflict_ages:
        stats += ('fixed_age',)
    header = ['name'] + ['%s_%s' %(s,m) for s in stats for m in measures] + \
             ['calibrated', 'parent', 'taxon', 'is_outgroup'] + \
             ['type', 'root', 'subgenetree']
    if not args.totree:
        print('\t'.join(header))

    parse = parserchoice[args.parser]
    convert = converterchoice[args.parser]['ete3']
    for tree_obj in parse(args.treefile):
        tree = convert(tree_obj)
        tree.add_features(dirname=op.dirname(args.treefile),
                          basename=op.basename(args.treefile))

        node_info = [('taxon', this_get_taxon), ('is_outgroup', is_outgroup)]
        ages, subtrees = bound_average(tree, calib, todate,
                                       measures,
                                       args.unweighted,
                                       args.original_leading_paths,
                                       args.correct_unequal_calibs,
                                       fix_conflict_ages=args.fix_conflict_ages,
                                       keeproot=args.keeproot,
                                       calib_selecter='taxon',
                                       node_info=node_info,
                                       node_feature_setter=[('type', get_eventtype)],
                                       dataset_fmt=args.dataset_fmt)

        if args.totree:
            age_col = header.index('age_'+measures[0])
            # Better have unique node names.
            anc_ages = {row[0]: row[age_col] for row in ages}
            save_mpl_tree(tree, anc_ages, stdout)
        else:
            save_ages(ages, stdout)


if __name__ == '__main__':
    main()

