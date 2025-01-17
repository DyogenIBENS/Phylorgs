#!/usr/bin/python3
# coding: utf-8

# flake8-in-file-ignores: noqa: F403,E261,E266

import sys
import os.path as op
from itertools import chain
import textwrap
from datetime import datetime as dt

import warnings
from copy import deepcopy
import numpy as np
from numpy.lib.arraysetops import setdiff1d
import pandas as pd
import matplotlib as mpl
#mpl.use('TkAgg', warn=False)
#try:
#    get_ipython().magic('matplotlib inline')
#except UnknownBackend:
import matplotlib.pyplot as plt
import seaborn as sb

from objectools import Args, as_args, generic_update, generic_remove_items
from seqtools.compo_freq import weighted_std
from dendro.bates import dfw_pairs_generalized, dfw_pairs
from datasci.graphs import scatter_density, \
                           intersection_plot, \
                           heatmap_cov, \
                           plot_loadings, \
                           plottree, \
                           plottree_set_xlim, \
                           stackedbar, \
                           cathist, \
                           fade_color_hex, \
                           darken_hex
from datasci.stats import r_squared, multicol_test, multi_vartest,\
                          rescale_groups, iqr, iqr90, iqr95, ci95, mad, trimstd, trimmean, \
                          mean_absdevmed, VIF, partial_r_squared_fromreduced
from datasci.routines import *
from datasci.dataframe_recipees import *

from dendro.any import myPhylTree as phyltree_methods, ete3 as ete3_methods

from LibsDyogen import myPhylTree

from scipy import stats
from sklearn.decomposition import PCA, FactorAnalysis
import statsmodels.api as sm
import statsmodels.stats.multitest as smm
import statsmodels.graphics as smg  # smg.gofplots.qqplot

from IPython.display import display_html, display_markdown
from IPython.utils.io import capture_output
from datasci.savior import HtmlReport, CellReport, DoubleReport, \
                      css_dark_style, generate_slideshow, slideshow_generator,\
                      reroute_loggers
import logging

try:
    # BUG in Jupyter Notebook: I must catch the following error:
    from joblib.externals.loky.process_executor import TerminatedWorkerError
except ImportError:
    pass

_LOGFMT = "%(levelname)-7s:l.%(lineno)3s:%(funcName)-20s:%(message)s"
_LOGF = logging.Formatter(_LOGFMT)

try:
    from UItools import colorlog
    _CLOGFMT = "$LVL%(levelname)-7s:${RESET}l.%(lineno)2s:$RESET%(message)s"
    _COLORLOGF = colorlog.ColoredFormatter(_CLOGFMT)
    _CLOGFMT2 = "$LVL%(levelname)-7s:${RESET}%(module)s:$BLACK%(funcName)-20s$RESET:%(message)s"
    _COLORLOGF2 = colorlog.ColoredFormatter(_CLOGFMT2)
except ImportError:
    _COLORLOGF = _COLORLOGF2 = _LOGF

# Notebook setup
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
# Loggers needed:
logger_stats = logging.getLogger('datasci.stats')  # == datasci.stats.logger
logger_graphs = logging.getLogger('datasci.graphs') # == datasci.graphs.logger
logger_routines = logging.getLogger('datasci.routines')
logger_dfr = logging.getLogger('datasci.dataframe_recipees')
loggers = (logger, logger_stats, logger_graphs, logger_routines, logger_dfr)

sh = logging.StreamHandler(sys.stdout)
sh.setFormatter(_COLORLOGF)
for hdlr in logger.handlers:
    hdlr.close()
logger.handlers = []
logger.addHandler(sh)

sh2 = logging.StreamHandler(sys.stdout)
sh2.setFormatter(_COLORLOGF2)
for logr in loggers[1:]:  # From imported modules
    for hdlr in logr.handlers:
        hdlr.close()
    logr.handlers = []
    logr.addHandler(sh2)
    #logging.basicConfig(handlers=[sh])
    #logging.basicConfig()

myslideshow_js = op.expanduser('~/mydvpt/atavistic-doc-tools/myslideshow.js')
myslideshow_css = myslideshow_js.replace('.js', '.css')


mpl.style.use("softer")  # See "softdark" as well.
pd.set_option("display.max_columns", 50)
pd.set_option("display.width", 115)
pd.set_option("display.max_colwidth", 50)
pd.set_option("display.show_dimensions", True)  # even on non truncated dataframes

# wide_screen_style
mpl.rcParams['figure.figsize'] = (14, 10) # width, height

# Calculation: fit 2 figures/page. margins: side=2.5cm top/bottom=3cm (+header)
# --> imposes a ratio height/width < 0.6
# available width = 7.28, height<4.43
thesisfigsize = (8.74, 5.25)  # visual fontsize = 10
thesisfigsize = (9.71, 5.85)  # visual fontsize = 9
#thesisfigsize = (10.92, 6.59)  # visual fontsize = 8

# Useful variables: the name of statistics of interest:
## Piping of the metrics:
## measures: '{m}'

## chronos ages: 'age.{m}'     with m in ('PL1', 'PL100', 'R1', ..., 'MPL')

## age_measures -> 'age_{m}    m in ('dS', 't', 'dN', 'dist') (represents substitutions/site)

## Original branch lengths in substitutions
## dist_measures/branch_measures: 'branch_{m}'
## same but with the parent branch: += '_parent'

## rate_measures -> 'rate_{m}' computed by dividing the branch_m by the real length (My)
## rate_std_measures -> 'rate_{m}_std' std dev of the rate in each tree.

## consecutive_null_{m}
## freq_null_{m}
## null_{m}_after
## null_{m}_before

## median_measures

MEASURES = ['dS', 'dN', 't', 'dist']
DIST_MEASURES = ['branch_%s' % m for m in MEASURES]
RATE_MEASURES = ['%s_rate' % m for m in MEASURES]
RATE_STD_MEASURES = [r + '_std' for r in RATE_MEASURES]


# ~~> numbertools? timetools? converters?
def time2seconds(time_str):
    """Convert "time used" into seconds."""
    factors = [1, 60, 3600, 3600*24]
    s = 0
    for factor, n_units in zip(factors, reversed(time_str.split(':'))):
        s += factor * int(n_units)
    return s


def load_stats_al(alfile):
    """Load stats table computed from alignments"""
    return pd.read_csv(alfile, index_col=0, sep='\t')

load_any_stats = load_stats_al


def load_stats_tree(treefile):
    """Load stats table computed from trees"""
    ts = pd.read_csv(treefile, index_col=0, sep='\t',
                     dtype={'leaves_robust':      bool,
                            'single_child_nodes': bool,
                            'nodes_robust':       bool,
                            'only_treebest_spe':  int,
                            'aberrant_dists':     int})
                            #'root2tip_var': float
    ts['really_robust'] = ts.leaves_robust & ts.nodes_robust & ~ts.only_treebest_spe
    try:
        ts['unlike_clock'] = ts.root2tip_sd / ts.root2tip_mean
    except AttributeError:
        logger.warning('"root2tip_mean" or "root2tip_sd" not in treestats columns.')
    return ts

def load_stats_codeml(codemlfile):
    """Load stats table computed from Codeml output"""
    cs = pd.read_csv(codemlfile, index_col=0, sep='\t')
    cs['seconds'] = cs['time used'].apply(time2seconds)
    return cs


def load_stats_chronosruns(runsfile):
    """Load stats table computed from APE::chronos calculations"""
    rs = pd.read_csv(runsfile, sep='\t', index_col=0)
    rs.set_axis([c.replace('.', '_') for c in rs.columns], axis=1, inplace=True)
    for c in rs.columns:
        if c.endswith('_message'):
            rs[c] = rs[c].astype('category')
    return rs

def load_stats_chronoslogs(logsfile):
    ls = pd.read_csv(logsfile, sep='\t', index_col=0,
                     true_values='T', false_values='F',
                     dtype={'action': 'category'})
    return ls

stat_loaders = {'al':     load_stats_al,
                'tree':   load_stats_tree,
                'codeml': load_stats_codeml,
                'cleaning': load_any_stats,
                'chronos.runs': load_stats_chronosruns,
                'chronos.logs': load_stats_chronoslogs}

stat_loaders['treeI'] = stat_loaders['tree']
stat_loaders['codemlfsahmmc'] = stat_loaders['codemlfsa'] = stat_loaders['codemlI'] = stat_loaders['codeml']
stat_loaders['alfsahmmc'] = stat_loaders['alfsa'] = stat_loaders['alI'] = stat_loaders['al']
stat_loaders['cleaningfsa'] = stat_loaders['cleaning']


def load_subtree_stats(template, stattypes=('al', 'tree', 'codeml')):
    """
    param: template to the csv/tsv files containing al/tree/codeml stats.

    Example template: 'subtreesRawTreeBestO2_{stattype}stats_Simiiformes.tsv'
    """

    # ## Load tree/alignment statistics
    output_tables = []
    for stattype in stattypes:
        filename = template.format(stattype=stattype)
        logger.info('Load %s', filename)

        try:
            loader = stat_loaders[stattype]
        except KeyError:
            loader = load_any_stats
            logger.warning('No specific parameters for loading %r',
                    stattype)
        output_tables.append(loader(filename))

    return tuple(output_tables)


#def check_load_subtree_stats(aS, ts, cs, clS=None, figsize=(10,3)):
def check_subtree_stats(subtrees_stats, figsize=(10,3)):

    for stype, ss in subtrees_stats:
        logger.info("%s shape: %s. Has dup: %s", stype, ss.shape, ss.index.has_duplicates)
        only_stype = ss.index.difference(set.union(*(set(oss.index)
                                            for ostype, oss in subtrees_stats
                                            if ostype != stype))
                                        )
        l_only = len(only_stype)
        if l_only:
            logger.warning("%d only in %s stats: %s", l_only, stype,
                           list(only_stype)[:min(5, l_only)])

    common_subtrees = set.intersection(*(set(ss.index) for _, ss in subtrees_stats))
    logger.info("%d common subtrees", len(common_subtrees))
    intersect_kwargs = {stype: ss.index for stype, ss in subtrees_stats}

    ax = intersection_plot(**intersect_kwargs)
    ax.figure.set_size_inches(figsize)
    return ax


# ## Function to merge additional subgenetree information into `ages`
def merge_criterion_in_ages(criterion_serie, ages=None, ages_file=None,
                            criterion_name=None):
    """Merge a column into the *node ages* table: the common field is the *subgenetree*."""

    assert (ages is not None or ages_file) and not (ages_file and ages is not None), "At least `ages` (dataframe) or `ages_file` (filename) must be given."
    if ages is None:
        ages = pd.read_csv(ages_file, sep='\t')

    logger.info("Input shape: %s", ages.shape)
    criterion = criterion_serie.name if not criterion_name else criterion_name
    ages_subgenetrees = ages.subgenetree.unique()

    if len(criterion_serie) < ages_subgenetrees.size \
            or set(criterion_serie.index) < set(ages_subgenetrees):
        logger.error("Not all genetrees have a criterion value: %d < %d" % (len(criterion_serie), ages_subgenetrees.size))

    #criterion_df = pd.DataFrame({criterion: criterion_serie})
    criterion_df = criterion_serie.to_frame(criterion)
    assert (criterion_df.index == criterion_serie.index).all()

    ages_c = pd.merge(ages, criterion_df, how='left', left_on='subgenetree', right_index=True)
    # If all subgenetrees have a criterion, equivalent to:
    #ages[criterion] = criterion[ages.subgenetree]
    return ages_c



def load_prepare_ages(ages_file, ts=None, measures=['dist', 'dS', 'dN', 't'], compute_extra_ns=None):
    """Load ages dataframe, join with parent information, and compute the
    'robust' info.

    `compute_extra_ns`: a function that takes the ages.groupby('subgenetree')
                object and output a dataframe to be joined with `ns` (dataset
                specific statistics, one row per subgenetree).
                if it produces a column named "treelink", this will be used to join `ts`.
                (otherwise it joins on its index, i.e. the subgenetree values).
    """

    beast_renames = {'height': 'age_beast',
                     'height_median': 'age_beastmedian',
                     'length': 'branch_beast',
                     'length_median': 'branch_beastmedian',
                     'beast:dist': 'dist'}

    ages = pd.read_csv(ages_file, sep='\t', index_col=0)\
            .rename(columns=lambda n: n.replace('.', '_'))\
            .rename(columns=beast_renames)\
            .rename_axis('name')  # If coming from date_dup.R
    if measures is None:
        measures = [c[4:] for c in ages.columns if c.startswith('age_')]

    # Type checks:
    dtypes = {} if measures is None else {'%s_%s' % (s, m): float
                                          for s in ('age', 'branch')
                                          for m in measures}
    dtypes.update(calibrated=bool, parent=str, taxon=str, is_outgroup=bool,
                  type=str, root=str, subgenetree=str)
    try:
        ages = ages.astype(dtypes, copy=False, errors='raise')
    except KeyError as err:
        err.args += ('%r in `dtypes` not in columns.' % set(dtypes).difference(ages.columns),)
        raise

    # Convert beast string ranges to floats, in 2 columns:
    rangevars = ['height_95%_HPD', 'height_range', 'length_95%_HPD', 'length_range', 'rate_95%_HPD', 'rate_range']
    if not ages.columns.intersection(rangevars).empty:
        logger.debug('ages[rangevars].dtypes = %s', ages[rangevars].dtypes)
        for var in rangevars:
            parsedrange = ages[var].replace("None", "NaN").astype(str)\
                          .str.strip('{}').str.split(r'\.\.', n=2, expand=True)\
                          .astype(float)
            logger.debug('%s: parsedrange.shape = %s; parsedrange.columns = %s',
                         var, parsedrange.shape, parsedrange.columns.tolist())
            ages[var+'_low'] = parsedrange[0]
            ages[var+'_up'] = parsedrange[1]
        ages.drop(columns=rangevars, inplace=True)

    logger.info("Shape ages: %s; has dup node names: %s; has dup (name,subgenetree): %s",
                ages.shape, ages.index.has_duplicates, ages.set_index('subgenetree', append=True).index.has_duplicates)
    #FIXME: upstream code should always index on (name,subgenetree), thus not be affected by duplicated names.
    n_nodes = ages.shape[0]
    ages_int = ages[ages.type != 'leaf']
    logger.info("Shape ages internal nodes: %s; has dup node names: %s; has dup (node,subgenetree): %s",
                ages_int.shape,
                ages_int.index.has_duplicates,
                ages_int.set_index('subgenetree', append=True).index.has_duplicates)
    n_nodes_int = (ages.type != 'leaf').sum()

    #TODO:
    # Really important here:
    # 1. index (=node.name) should be UNIQUE.
    # 2. the "parent" column should intersect with the index!
    #
    # Therefore the design: a unique index + uniq parent
    # based on (subgenetree, name) and (subgenetree, parent).

    # Fetch parent node info  ## FIXME: handle missing column 'is_outgroup'
    ages_p = pd.merge(ages.reset_index(),
                      ages.loc[ages.type != 'leaf',
                            ['taxon', 'type', 'calibrated', 'is_outgroup',
                             'subgenetree']
                            + ['%s_%s' % (s,m) for s in ('age', 'branch')
                               for m in measures]],
                      how="left", left_on=["parent", "subgenetree"],
                      right_on=["name", "subgenetree"],
                      suffixes=('', '_parent'), indicator=True,
                      validate='many_to_one')\
                .set_index('name')
    logger.info("Shape ages with parent info: %s", ages_p.shape)
    n_nodes_p = ages_p.shape[0]

    if n_nodes_p < n_nodes:
        logger.warning('%d nodes were lost when fetching parent information.',
                       n_nodes - n_nodes_p)

    orphans = ages_p._merge=='left_only'
    ages_orphans = ages_p[orphans]
    # if --keeproot in generate_dNdStable: orphans are the root nodes;
    # else they are the direct children of the roots that could be dated with MPL,
    #      or if tabulate_ages_from_tree was used, any child of the root.

    logger.info("Orphans: %d", orphans.sum())
    #display_html(ages_orphans.head())

    #expected_orphans = ages.is_outgroup_parent.isna()
    recognized_outgroups = (ages_orphans.is_outgroup | ages_orphans.is_outgroup_parent)
                            #| ages_orphans.is_outgroup_parent == 1)
    orphan_taxa = ages_orphans[~recognized_outgroups].taxon.unique()

    #print("All orphans are expected (leaf, or child of the root).")
    logger.info('Computed ages of orphan taxa:\n%s\n%s',
                ages_orphans.loc[~recognized_outgroups,
                                 ['age_%s' %m for m in measures]].head(),
                ages_orphans.loc[~recognized_outgroups,
                                 ['age_%s' %m for m in measures]].isna().sum()\
                            .rename('Sum of NaNs').to_frame().T)
    logger.warning('CHECK those orphan taxa (should all be unrecognized outgroups sequences (could be duplicates from the ingroup taxa): %s.',
                   ', '.join(orphan_taxa))

    e_nochild = set(ages.index) - set(ages.parent)  # expected
    parent_nodata = set(ages[ages.type!='leaf'].parent) - set(ages_p.index)
    n_nochild = len(e_nochild)
    logger.info("Expected nodes without children (leaves): %d", n_nochild)
    logger.info("Observed nodes not found as parents: %d",
                len(set(ages.index) - set(ages_p[ages_p._merge=='both'].index)))
    logger.info("Parent nodes without data: %d", len(parent_nodata))
    #assert len(nochild) == n_nochild, \
    #    "Found %d unexpected nodes without child:\n%s" % \
    #        (len(nochild) - n_nochild,
    #         ages.loc[nochild - e_nochild])

    # Drop outgroups! (they create duplicated index values and might confuse stuff
    # HOWEVER DON'T drop the ingroup root.
    # (doesn't have an is_outgroup_parent info)

    child_of_root = ages_p.root == ages_p.parent

    # Ad hoc formatting: to find the ingroup root, assume that the subgenetree name is the ingroup name...
    ingroup_root = ages_p.index.to_series() == ages_p.subgenetree
    if ingroup_root.sum() == 0:
        logger.warning('ZERO ingroup root found (based on matching the subgenetree name)')
    #recognized_outgroups = (ages_p.is_outgroup == 1 | ages_orphans.is_outgroup_parent == 1)
    #to_remove = (~child_of_root & orphans
    #             | (ages_p.is_outgroup == 1)
    #             | ages_p.is_outgroup_parent == 1)
    #to_remove = (orphans
    #             | (ages_p.is_outgroup == 1)
    #             | ages_p.is_outgroup_parent == 1)
    to_remove = (~ingroup_root & orphans
                 | (ages_p.is_outgroup == 1)
                 | ages_p.is_outgroup_parent == 1)
    if (to_remove == child_of_root).all():
        logger.error('ALL root children will be removed')

    ages_p = ages_p.loc[~to_remove].copy(deep=False)
    #ages_p = ages_p.loc[~to_remove].set_index('subgenetree', append=True, drop=False)
    if ages_p.index.has_duplicates:
        debug_columns = ['parent', 'subgenetree', 'taxon', 'taxon_parent']
        debug_columns += [s+measures[0] for s in ('branch_', 'age_')]
        logger.error("Failed to remove index duplicates in 'ages_p':\n%s\n...",
                     ages_p.loc[ages_p.index.duplicated(), debug_columns].head(10))

    ages_spe2spe = ages_p[(ages_p.type.isin(('spe', 'leaf')))
                          & (ages_p.type_parent == 'spe')]

    logger.info("\nShape ages speciation to speciation branches (no dup): %s",
                ages_spe2spe.shape)

    #ages_robust = ages_treestats[ages_treestats.really_robust & \
    #                             (ages_treestats.aberrant_dists == 0)]\
    #                        .drop(['really_robust', 'aberrant_dists'], axis=1)
    return add_robust_info(ages_p, ts, measures, compute_extra_ns)

def orig_extra_ns(sgg):
    return sgg.subgenetree.first().rename('ensembl_version', inplace=False)\
             .str.extract(r'^[A-Z][a-zA-Z_]+ENSGT00(\d\d)0.*').astype(int)

def add_robust_info(ages_p, ts=None, measures=['dist', 'dS', 'dN', 't'], compute_extra_ns=None):
    """
    Compute Number of duplications/speciation per tree,
    and additional tree specific statistics.
    """
    branch_measures = ['branch_'+m for m in measures]  # global dist_measures
    #NOTE: already computed by `subtrees_stats` in current version.
    for m in measures:
        logger.debug('measure: %s' %m)
        ages_p['consecutive_null_%s' %m] = \
                ~ages_p[['branch_%s_parent' %m, 'branch_%s' %m]].any(axis=1,
                                                                  skipna=False)
        #NOTE: there's a pitfall with this method as the `branch_dist` of the ingroup node is not documented by `generate_dNdS.py`!
        #-> with skipna=False, the NaN are True.

    sgg = subgenetree_groups = ages_p.groupby('subgenetree', sort=False)

    logger.info('Aggregating `ns` ("new stats" specific to this dataset. ages_p.shape = %s)...', ages_p.shape)
    # This is a slow operation. Needs optimization.
    
    # tree statistics making use of each node.
    ns = pd.concat((sgg.type.agg(lambda v: sum(v == "dup")),
                    sgg.type.agg(lambda v: sum(v == "spe"))),
                   axis=1, keys=['Ndup', 'Nspe'])

    logger.debug('ns original shape = %s', ns.shape)
    #print(ns.Ndup.describe())
    #print(ns.Nspe.describe())

    def freq_of_null(v):
        return (v==0).mean()

    # Already computed by codeml.subtrees_stats.get_codeml_stats
    agg_funcs = {'consecutive_null_%s' %m: 'mean' for m in measures}

    #
    agg_funcs.update({'branch_%s' %m: freq_of_null for m in measures})
    ns = ns.join(sgg.agg(agg_funcs))\
           .rename(columns={'branch_%s' %m: 'freq_null_%s' %m for m in measures})

    # Now, count null branches **around** dated nodes (before/after)

    sgg_before = ages_p[ages_p.calibrated==0].groupby('subgenetree', sort=False)
    sgg_after = ages_p[ages_p.calibrated_parent==0].groupby('subgenetree', sort=False)

    ns = ns.join(sgg_before[branch_measures].agg(freq_of_null), sort=False)\
           .rename(columns={'branch_%s' %m: 'null_%s_before' %m
                            for m in measures})\
           .join(sgg_after[branch_measures].agg(freq_of_null), sort=False)\
           .rename(columns={'branch_%s' %m: 'null_%s_after' %m
                            for m in measures})
           #.join(sgg_after[['branch_dN', 'branch_dS']]\
           #      .apply(lambda r: (r==0).all(axis=1))
    logger.debug('ns joined shape = %s', ns.shape)

    if compute_extra_ns is not None:
        nscolumns = ns.columns.tolist()
        ns = ns.join(compute_extra_ns(sgg), sort=False)
        logger.debug('ns.columns = %s', ns.columns.tolist())
        if set(nscolumns) == set(ns.columns):
            logger.error('No new column added by `compute_extra_ns`. Head is:\n%s', compute_extra_ns(sgg).head())

    # merge tree stats to select robust trees
    merged_ns = ['Ndup', 'Nspe']
    linktree = None  # merge on index by default
    if 'linktree' in ns:
        linktree = 'linktree'
        merged_ns.append(linktree)

    if ts is None:
        nsts = ns[merged_ns]
    else:
        nsts = ns[merged_ns].join(ts[['really_robust', 'aberrant_dists',
                                      'rebuilt_topo']],
                                  on=linktree, sort=False)  # Is this useful?
        if np.all(nsts.really_robust.isna()):
            msg = 'Merging ns with ts failed: zero common tree?'
            linked_ts = ts[linktree] if linktree else ts.index.to_series()
            pasted_heads = '\n'.join(line1 + '\t' + line2 for line1, line2 in
                                     zip(str(ns.index.to_series().head()).split('\n'),
                                         str(linked_ts.head()).split('\n')))
            logger.error('%s Check:\nns.index\tts[%s]\n%s', msg, linktree, pasted_heads)
            raise ValueError(msg)
    ages_treestats = pd.merge(ages_p.drop('_merge', axis=1),
                              nsts,
                              how='left',
                              left_on='subgenetree',
                              right_index=True,
                              indicator=True, validate='many_to_one')
    logger.info("Ndup %s; Nspe %s; ages_treestats %s",
                (~ns.Ndup.isna()).sum(), (~ns.Nspe.isna()).sum(),
                ages_treestats.shape)
    logger.info('merge types:\n%s',
                ages_treestats._merge.value_counts(dropna=False))

    return ages_treestats, ns
    #return ages_p, ns


# for averaging by taking into account branch length: with Omega.
# NOTE: columns will be reordered following `var`.
def old_group_average(g, var, weight_var="median_brlen"):
    g = g[~g[weight_var].isna()]
    #g = g.dropna(subset=[weight_var])
    #if not g.shape[0]:
    #    return pd.Series([np.NaN]*len(var))
    return pd.Series(np.average(g[var], axis=0, weights=g[weight_var]))
    # TODO: append the count of removed rows (NA in weight_var)

def raw_group_average(g, var, weight_var):
    gv = g[var].values
    gw = g[weight_var].values
    keep = ~np.isnan(gw)
    return np.average(gv[keep], axis=0, weights=gw[keep])

def npraw_group_average(gvalues, var_idx, weight_idx):
    #values = np.ma.array(g[var].values, mask=g[var].isna().values)
    gv = gvalues[:,var_idx]
    gw = gvalues[:,weight_idx]
    keep = ~np.isnan(gw)
    return np.average(gv[keep], axis=0, weights=gw[keep])

def npraw_group_average_w0(gvalues):
    """Vectorized weighted average from numpy input (2D) to numpy output (1D).

    From an array `gvalues` of dim (N, M) (lines, columns):
    - take the row weights in column 0;
    - ignore rows with NaN weight;
    - compute the weighted average of the other columns.
    """
    #values = np.ma.array(g[var].values, mask=g[var].isna().values)
    gv = gvalues[:,1:]
    gw = gvalues[:,0]
    keep = ~np.isnan(gw)  # What about isnan(gv) (when branch==0 and branchtime==0)
    try:
        return np.average(gv[keep], axis=0, weights=gw[keep])
    except ZeroDivisionError:
        return np.NaN

def group_average(g, var, weight_var):
    #values = np.ma.array(g[var].values, mask=g[var].isna().values)
    gv = g[var].values
    gw = g[weight_var].values
    keep = ~np.isnan(gw)
    return pd.Series(np.average(gv[keep], axis=0, weights=gw[keep]))


def group_weighted_std(g, var, weight_var="median_brlen"):
    return pd.Series(weighted_std(g[var], axis=0, weights=g[weight_var]))

# for calculating rates: with dS, dN, t
def tree_dist_2_rate(g, dist_var, norm_var="median_brlen"):
    # in pandas, sum is done per columns.
    return g[dist_var].sum() / g[norm_var].sum()


def add_control_dates_lengths(ages, phyltree, control_ages_CI=None,
                              measures=['dist', 'dS', 'dN', 't'],
                              control_condition='really_robust & aberrant_dists == 0',
                              control_names=['timetree']  # 'dosreis'
                              ):
    # 'consecutive_zeros_dS==0'
    # Merge control dates
    ages_forcontrol = ages if not control_condition else ages.query(control_condition)
    logger.info("%d nodes from control trees", ages_forcontrol.shape[0])

    age_measures = ['age_'+m for m in measures]
    median_age_measures = ['median_age_'+m for m in measures]
    control_ages = ages_forcontrol[ages_forcontrol.type.isin(("spe", "leaf"))]\
                                   .groupby("taxon", sort=False)[age_measures]\
                                   .median()\
                                   .rename(columns=dict(zip(age_measures,
                                                            median_age_measures)))
    logger.debug('median_taxon_ages.columns = [%s]',
                 ' '.join(control_ages.columns))
    # FIXME: WHY does `median_taxon_ages` only have the 1st age_measure?
    #TODO: change for each measure as in `join_extra_ages`

    if control_ages_CI is not None:
        toconcat = (control_ages,)
        if 'timetree_age' not in control_ages_CI.columns:
            timetree_age = control_ages.index.to_series()\
                                             .apply(phyltree.ages.get)\
                                             .rename('timetree_age')  # age_timetree
            toconcat += (timetree_age,)

        toconcat += (control_ages_CI,)
        control_ages = pd.concat(toconcat, axis=1, sort=False)
    logger.debug('Control ages (columns): ' + ' '.join(control_ages.columns))

    # Add potentially missing ages (0 for extant species).
    for ctl in control_names:
        na_ages = control_ages['%s_age' % ctl].isna()
        if na_ages.any():
            for taxon in control_ages.index.values[na_ages]:
                # PhylTree ages are those from TimeTree
                if taxon in phyltree.listSpecies or ctl=='timetree':
                    control_ages.loc[taxon, '%s_age' % ctl] = phyltree.ages.get(taxon, np.NaN)

    #print(control_ages.sort_values('timetree_age', ascending=False))

    ages_controled = pd.merge(ages, control_ages,
                              left_on="taxon", right_index=True,
                              validate="many_to_one")

    # Merge control branch lengths
    invalid_taxon_parent = ages_controled.taxon_parent.isna()
    invalid_measures = ages_controled[age_measures].isna().all(1)
    invalid_nodes = invalid_taxon_parent & ~invalid_measures
    # calibrated | branch_dS.isna()
    # Should be nodes whose parent node is the root.
    debug_columns = ['parent', 'subgenetree', 'taxon', 'taxon_parent'] \
                    + median_age_measures
    if 'dist' in measures:
        debug_columns += ['branch_dist', 'age_dist']
    elif 'dS' in measures:
        debug_columns += ['branch_dS', 'age_dS']
    elif 'beast' in measures:
        debug_columns += ['branch_beast', 'age_beast']
    if invalid_nodes.any():
        #assert (ages_controled[invalid_taxon_parent].parent == \
        #        ages_controled[invalid_taxon_parent].root).all()
        logger.warning("%d invalid 'taxon_parent':head:\n%s\n"
                     "The following taxa have no parent taxa information, "
                     "please check:\n%s\n**DROPPING** this data!",
                       invalid_nodes.sum(),
                       ages_controled[invalid_nodes][
                           debug_columns
                           ].head(),
                       ', '.join(ages_controled[invalid_nodes].taxon.unique()))
        ages_controled.dropna(subset=['taxon_parent'], inplace=True)

    # Do not include length from/to duplication.
    ages_controled_spe2spe = ages_controled.dropna(subset=['taxon_parent']).query('type != "dup" & type_parent != "dup"')
    same_sp = ages_controled_spe2spe.taxon == ages_controled_spe2spe.taxon_parent
    if same_sp.any():
        logger.error("Failed to filter out duplications:\n%s",
                     ages_controled_spe2spe[same_sp].head(15))

    for ctl in control_names:
        ages_controled[ctl+'_brlen'] = \
            ages_controled_spe2spe.taxon_parent\
                                  .apply(control_ages[ctl+'_age'].get) \
            - ages_controled_spe2spe[ctl+'_age']

    #TODO: NaN for data not in ages_forcontrol?
    for m,am,mm in zip(measures, age_measures, median_age_measures):
        ages_controled['median_brlen_'+m] = \
            ages_controled_spe2spe.taxon_parent.apply(control_ages[mm].get) \
            - ages_controled_spe2spe[mm]
    #control_ages.reindex(ages_controled.taxon_parent)\
    #        .set_axis(ages_controled.index, inplace=False)
    # would be more 'Pandas-like'.

    # An approximation of the real branch length
    for m,am in zip(measures, age_measures):
        ages_controled['approx_%s' %m] = (ages_controled[am+'_parent']
                                            - ages_controled[am])
        # Also include an approx that can't be zero even with 0 substitutions:
        #almost 0:
        #min_val = ages_controled.loc[ages_controled['approx_%s' % m]>0, 'approx_%s' % m].min()
        # Can be as low as 8e-19
        #ages_controled['approx_%s_nonzero' % m] = ages_controled['approx_%s' % m].replace(0, 1e-19)

    # Resulting branch lengths
    branch_info = ["taxon_parent", "taxon"]
    #control_brlen = ages_controled.loc[
    #                    ~ages_controled.duplicated(branch_info),
    control_brlen = ages_controled.query('type != "dup" & type_parent != "dup"')\
                    .groupby(branch_info, sort=False)\
                    [['median_%s_%s' % (typ, m) for typ in ('age', 'brlen')
                        for m in measures]
                     + ["%s_%s" %(ctl, typ) for typ in ('age', 'brlen')
                         for ctl in control_names]]\
                    .first()#agg(lambda s: s[0])
                #.sort_values('timetree_age')
                    #].sort_values("taxon_parent", ascending=False)
                    #.reset_index(branch_info, drop=True)
                    #.set_axis(pd.MultiIndex.from_arrays(
                    #                control_brlen[branch_info].values.T,
                    #                names=branch_info),
                    #          inplace=False)\
                    #.drop(branch_info, axis=1)
    ###FIXME: Still some branches that are not in PhylTree (jump over taxa)
    #display_html(control_brlen)
    return ages_controled, control_ages, control_brlen


def check_control_dates_lengths(control_brlen, phyltree, root,
                                measures=MEASURES, out=None):
    """Check NA values and if the branches fit the phylogenetic tree.

    Return: (unexpected_branches, lost_branches)"""
    get_phylchildren = lambda phyltree, ancdist: phyltree.items.get(ancdist[0], [])

    expected_branches, expected_dists = zip(*(((p[0],ch[0]),ch[1]) for p,ch in
                             dfw_pairs_generalized(phyltree,
                                                   get_phylchildren,
                                                   queue=[(None, (root,0))])))

    #logger.debug(expected_branches)
    #logger.debug(expected_dists)

    na_brlen = control_brlen.isna().any(axis=1)
    if na_brlen.any():
        print("MISSING branch lengths:\n" + str(control_brlen[na_brlen]), file=out)

    median_measures = ['median_brlen_%s' % m for m in measures]
    median_brlen_sum = control_brlen[median_measures].sum()
    print("Sum of median branch lengths =", median_brlen_sum, "My", file=out)
    if 'timetree_brlen' in control_brlen:
        timetree_brlen_sum = control_brlen.timetree_brlen.sum()
        print("Sum of timetree branch lengths =", timetree_brlen_sum, "My", file=out)
        real_timetree_brlen_sum = sum(expected_dists)
        print("Real sum of TimeTree branch lengths (in phyltree) =",
              real_timetree_brlen_sum, "My", file=out)
    else:
        print("No branch lengths from TimeTree", file=out)

    unexpected_branches = set(control_brlen.index) - set(expected_branches)
    if unexpected_branches:
        logger.error("Extraneous branches not seen in phyltree:\n%s",
                     unexpected_branches)
    lost_branches = set(expected_branches) - set(control_brlen.index)
    if lost_branches:
        logger.error("Forgotten branches in phyltree:\n%s",
                     lost_branches)

    median_treelen_phyltree = control_brlen.reindex(list(expected_branches))[median_measures].sum()
    if 'timetree_brlen' in control_brlen:
        timetree_treelen_phyltree = control_brlen.reindex(list(expected_branches)).timetree_brlen.sum()
        print("Sum of median branch lengths for branches found in phyltree =\n",
              str(median_treelen_phyltree).replace('\n', '\t\n'), file=out)
        print("Sum of timetree branch lengths for branches found in phyltree =",
              timetree_treelen_phyltree, file=out)
    return unexpected_branches, lost_branches


def compute_dating_errors(ages_controled, control='median', measures=['dS'],
                          rescale=None):
    """:param: `control` in median/timetree"""
    if measures is None:
        age_vars = ['age']
        brlen_vars = ['brlen']
    else:
        age_vars = ['age_'+m for m in measures]
        brlen_vars = ['brlen_'+m for m in measures]

    control += '_'
    prefix = '' if control=='median_' else control

    ### Inplace
    for age_var, brlen_var in zip(age_vars, brlen_vars):
        #ages_controled = ages_controled.assign(**{

        if control == 'median_':
            ctl_age, ctl_brlen = control+age_var, control+brlen_var
        else:
            # control in ('timetree_', 'dosreis_')
            ctl_age, ctl_brlen = control+'age', control+'brlen'

        logger.debug('Control=%r age_var=%r => subtract %r', control, age_var, ctl_age)

        ages_controled[prefix+"signed_dev_"+age_var] = (ages_controled[age_var]
                                                   - ages_controled[ctl_age])
        ages_controled[prefix+"signed_dev_"+brlen_var] = (ages_controled[age_var+'_parent']
                                                     - ages_controled[age_var]
                                                     - ages_controled[ctl_brlen])
        if rescale is None:
            ages_controled[prefix+"abs_dev_"+age_var] = ages_controled[prefix+"signed_dev_"+age_var].abs()
            ages_controled[prefix+"abs_dev_"+brlen_var] = ages_controled[prefix+"signed_dev_"+brlen_var].abs()
        elif rescale=='sym':
            # 1. Symetrical method
            raise NotImplementedError
        elif rescale=='asym':
            # 2. Asymetrical method. Let c be the real age, c0 and c1 the calibrations before and after, and x the estimation:
            # if x - c > 0 (the estimation is older) => normalize by c0 - c
            # if x - c < 0 (the estimation is younger) => normalize by c1 - c
            #ages_controled["sym_dev_" + age_var] = (ages_controled["signed_dev_" + age_var]
            #                                        /( - ages_controled[control+age_var]))
            raise NotImplementedError
        else:
            raise ValueError('`rescale` not in (None, "sym", "asym")')

    # Compute the mean only for nodes that were not calibrated.
    sgg = ages_controled.groupby("subgenetree", sort=False)
    #dev_measures = ['abs_age_dev', 'signed_age_dev', 'abs_brlen_dev', 'signed_brlen_dev']
    mean_errors = pd.concat((
                    sgg[[prefix+dev+age_var for age_var in age_vars
                            for dev in ('abs_dev_', 'signed_dev_')]].sum().div(
                        sgg['calibrated'].agg(lambda v: (1-v).sum()), axis=0),
                    sgg[[prefix+dev+brlen_var for brlen_var in brlen_vars
                            for dev in ('abs_dev_', 'signed_dev_')]].sum().div(
                        sgg[['calibrated', 'calibrated_parent']]\
                           .apply(lambda df: ((1-df.calibrated) | (1-df.calibrated_parent)).sum()), axis=0)
                    ),
                    axis=1)
    #TODO: median_errors

    return mean_errors


def join_extra_ages(new_ages_file, ages_data,
                    control_condition='really_robust & aberrant_dists == 0',
                    suffixes=('_extra', ''),
                    control_names=['timetree']):
                    #subset=None):
    """Re-use outputs of `analyse_age_errors`/`load_prepare_ages`/
    `add_control_dates_lengths`
    """
    extra_ages = pd.read_csv(new_ages_file, sep='\t', index_col=0)\
                   .rename(columns=lambda n: n.replace('.', '_'))\
                   .rename_axis(index='name')\
                   .reset_index()\
                   .assign(subgenetree=lambda df: df.subgenetree.fillna(method='ffill'))
    age_measures = [colname for colname in extra_ages.columns
                        if colname.startswith('age_')]

    logger.info('Joining new age measures: ' + ' '.join(age_measures))

    # 1. Needs to concat the parent data!! (for compute_mean_errors of 'brlen')
    extra_ages = pd.merge(extra_ages,
                          extra_ages.loc[extra_ages.type!='leaf',
                                         ['name', 'subgenetree'] + age_measures],
                      how="left",
                      left_on=["parent", "subgenetree"],
                      right_on=["name", "subgenetree"],
                      suffixes=('', '_parent'),
                      indicator=True, validate='many_to_one')

    # 2. Join to the previous data.
    ages = pd.merge(
               extra_ages.drop(columns=['calibrated', 'type', 'parent', 'taxon']),
               ages_data.ages_controled_withnonrobust.reset_index(),
               'inner',  # to remove what has been removed (outgroups)
               on=['name', 'subgenetree'],
               suffixes=suffixes,
               sort=False)\
           .set_index('name')
    logger.debug('New ages shape: %s', ages.shape)
    logger.debug('ages.columns = %s', ' '.join(ages.columns))

    if set(ages_data.ages_controled.columns).intersection(age_measures):
        logger.warning('New ages already in the data: %s (joined with suffix %r)',
                       age_measures, suffixes[0])
        ages.rename(columns={am+'_parent'+suffixes[0]: am+suffixes[0]+'_parent'
                             for am in age_measures
                             if am in ages_data.ages_controled.columns},
                    inplace=True)
        #logger.debug('extra_ages.columns = %s', ' '.join(extra_ages.columns))
        age_measures = [am + suffixes[0]
                        if am in ages_data.ages_controled.columns
                        else am
                        for am in age_measures]
    measures=[am[4:] for am in age_measures]

    # 3. Rerun a small part of `add_control_dates_lengths`
    ages_forcontrol = ages.query(control_condition)
    logger.info("%d nodes from robust trees", ages_forcontrol.shape[0])

    # 3.1. Recompute the median ages.
    median_taxon_ages = ages_forcontrol[ages_forcontrol.type\
                                                       .isin(("spe", "leaf"))]\
                                   .groupby("taxon", sort=False)[age_measures].median()\
                                   .rename(columns={am: 'median_'+am
                                                    for am in age_measures})
    control_ages = pd.concat((median_taxon_ages, ages_data.control_ages),
                             axis=1, sort=False)

    #print(control_ages.sort_values('timetree_age', ascending=False))

    ages_controled = pd.merge(ages,
                              control_ages[['median_'+am for am in age_measures]],
                              left_on="taxon", right_index=True,
                              validate="many_to_one")

    # 3.2. Check duplication nodes (optional?)
    ages_controled_spe2spe = ages_controled.dropna(subset=['taxon_parent']).query('type != "dup" & type_parent != "dup"')
    same_sp = ages_controled_spe2spe.taxon == ages_controled_spe2spe.taxon_parent
    if same_sp.any():
        logger.error("Failed to filter out duplications:\n%s",
                     ages_controled_spe2spe[same_sp].head(15))

    # 3.3. Median branch lengths from the given measures.
    ages_controled = ages_controled.assign(**{
        ('median_brlen_%s' % m):
            (ages_controled_spe2spe.taxon_parent\
                        .apply(control_ages['median_%s' % am].get)
             - ages_controled_spe2spe['median_%s' % am])
        for m,am in zip(measures, age_measures)})

    branch_info = ["taxon_parent", "taxon"]
    #control_brlen = ages_controled.loc[
    #                    ~ages_controled.duplicated(branch_info),
    control_brlen = ages_data.control_brlen.join(
                        ages_controled.query('type != "dup" & type_parent != "dup"')\
                        .groupby(branch_info, sort=False)\
                        [['median_%s' % am for am in age_measures]]\
                        .first())
    #check_control_dates_lengths(control_brlen, phyltree, root, measures)

    ages_controled_cond = ages_controled.query(control_condition).copy(deep=False)
    # copy to silence the `SettingWithCopyWarning`. Alternatively, set `_is_copy=None`
    mean_errors = pd.concat(
                    [compute_dating_errors(ages_controled_cond, ctl, measures)
                     for ctl in ['median']+control_names],
                    axis=1, sort=False, copy=False)

    return age_analysis_data(ages_controled_cond,
                             ages_controled,
                             None,  # ns  #TODO: add a dict of functions to aggregate columns.
                             control_ages,
                             control_brlen,
                             mean_errors)


def lineage_evolutionary_rates(anc, ages_file,
                               phyltree,
                               control_ages_CI,
                               stats_tmpl='subtreesGoodQualO2_%sstats-%s.tsv',
                               control_condition='really_robust & aberrant_dists == 0',
                               measures=MEASURES,
                               control_names=['timetree'],  # 'dosreis'
                               saveas=None):
    # Also see `analyse_age_errors`
    ts = load_stats_tree(stats_tmpl % ('tree', anc))
    aS = load_stats_al(stats_tmpl % ('al', anc))

    ages, ns = load_prepare_ages(ages_file, ts, measures)
    ages = pd.merge(ages, aS.join(ns, how='outer'),
                    left_on='subgenetree', right_index=True, sort=False,
                    validate='many_to_one')

    ages_controled_withnonrobust, control_ages, control_brlen =\
        add_control_dates_lengths(ages, phyltree, control_ages_CI, measures,
                                  control_condition, control_names)

    ages_controled = ages_controled_withnonrobust.query(control_condition).copy(deep=False)
    unexpected_branches, lost_branches = check_control_dates_lengths(
                                                control_brlen, phyltree, anc,
                                                measures)
    control_brlen.drop(unexpected_branches, inplace=True)
    age_analysis = age_analysis_data(ages_controled,
                            ages_controled_withnonrobust,
                            ns, control_ages, control_brlen, mean_errors=None)

    dist_measures = ['branch_%s' % m for m in measures]
    return (age_analysis,) + lineage_evolutionary_rates_fromdata(age_analysis, dist_measures, control_names)


def lineage_evolutionary_rates_fromdata(age_analysis, dist_measures,
                                        control_names=['timetree']):
    def wmean(g, wkey='ingroup_glob_len'):
        gdata = g.drop(columns=wkey)
        return pd.Series(np.average(gdata, axis=0, weights=g[wkey]), index=gdata.columns)
        #return np.average(gdata, axis=0, weights=g[wkey])

    #def func(g, wkey='ingroup_glob_len'):
    #    gdata = g.drop(columns=wkey)
    #    return pd.DataFrame(
    #                np.average(gdata, axis=0, weights=g[wkey]),
    #                     index=gdata.columns)

    lineage_groups = age_analysis.ages_controled.groupby(['taxon_parent', 'taxon'], sort=False)
    lineage_brlen = lineage_groups[dist_measures]\
            .agg(['median', 'mean', 'std'])\
            .join(lineage_groups[dist_measures + ['ingroup_glob_len']]\
                        .apply(wmean)\
                        .set_axis(pd.MultiIndex.from_product([dist_measures, ['wmean']]),
                            axis=1, inplace=False),
                  sort=False)\
            .sort_index(1)
    ctl = control_names[0]

    print('lineage_brlen =\n', lineage_brlen.head().iloc[:,:4])
    print('control_brlen =\n', age_analysis.control_brlen.head().iloc[:,:4])
    lineage_rates = lineage_brlen.div(age_analysis.control_brlen[ctl+'_brlen'],
                                      axis='index')

    branch_na_rates = lineage_rates.isna().any(axis=1)
    if branch_na_rates.any():
        logger.info('DROP NA rates rows:\n%s', lineage_rates.index[branch_na_rates].values)
        lineage_rates.dropna(inplace=True)

    return lineage_rates, lineage_brlen


def display_lineage_evolutionary_rates(lineage_rates, lineage_brlen,
                                       age_analysis, phyltree, anc,
                                       measures=MEASURES,
                                       control='timetree',
                                       figsize=None,
                                       extra_text='',
                                       cmap_minmax=None):

    ages_controled = age_analysis.ages_controled
    control_brlen = age_analysis.control_brlen
    dist_measures = ['branch_%s' % m for m in measures]
    rate_measures = ['%s_rate' % m for m in measures]

    ordered_branches = ['%s--%s' % br for br in dfw_pairs(phyltree, queue=[(None, anc)], closest_first=True)]

    ordered_branches_bylen = ['%s--%s' % v
                              for v in control_brlen.sort_values(control+'_brlen').index]

    # Table of summary values
    styled_rates = lineage_rates\
            .style.background_gradient(cmap='PRGn')#,
                                       #subset=[(d, stat)
                                       #    for d in dist_measures
                                       #    for stat in ('median', 'mean', 'std')])
    display_html(styled_rates)
    outputs = [styled_rates]

    m = measures[0]
    br_m = 'branch_'+m
    rates = lineage_rates[(br_m, 'median')]
    print('Min of median %s rate: %s (%s)' %(m, rates.min(), rates.idxmin()))
    print('Max of median %s rate: %s (%s)'%(m, rates.max(), rates.idxmax()))

    # Violin plot of rates
    ages_data = ages_controled\
                    .assign(
                        branchtaxa=(ages_controled.taxon_parent + '--' + ages_controled.taxon)
                    )\
                    .join(ages_controled[dist_measures].div(
                                ages_controled.timetree_brlen, axis=0),
                          how='outer', rsuffix='_rate')

    ax = sb.violinplot('branchtaxa', br_m+'_rate',
                       data=ages_data,
                       width=1, order=ordered_branches_bylen)

    values = ages_data.groupby('branchtaxa', sort=False)[br_m+'_rate']
    #rate_q = values.quantile([0.01, 0.99]).values
    #rate_range = rate_q[1] - rate_q[0]
    #ylim = (rate_q[0] - 0.01*rate_range, rate_q[0] + 0.01*rate_range)
    ndev = 2
    get_ylow = lambda v: max(v.min(), v.mean() - ndev*v.std())
    get_yup = lambda v: min(v.max(), v.mean() + ndev*v.std())

    ylim = (values.apply(get_ylow).min(),
            values.apply(get_yup).max())

    logger.debug('Setting new ylim: %s', ylim)
    ax.set_ylim(ylim)
    ax.set_ylabel(m+' rate')
    plt.setp(ax.get_xticklabels(), rotation=45, va='top', ha='right')
    ax.set_title('Variation of synonymous substitutions/site/My' +
                 (('. '+extra_text) if extra_text else ''))
    outputs.append(ax.figure)

    # Tree with branch length in median number of substitutions
    fig, ax = plt.subplots(figsize=figsize)
    plottree(phyltree,
             get_items=(lambda t,nd: [(ch,
                                       lineage_brlen.loc[(nd[0], ch),
                                                         (br_m, 'median')])
                                     for ch,d in t.items.get(nd[0], [])]),
             get_label=(lambda t,n: n),
             root=anc, ax=ax)
    ax.set_xlabel('Median %s branch length' % m)
    ax.set_title(extra_text)
    outputs.append(fig)

    # Tree with colored branches
    fig, (cax,ax) = plt.subplots(1, 2, figsize=figsize,
                                 gridspec_kw={'width_ratios': [1,20]})
    cax.set_position(cax.get_position().anchored('SE').shrunk(1, 0.5).shrunk_to_aspect(20))

    lines, coords, subaxes = plottree(phyltree,
                 lambda t,nd: [(ch, control_brlen.loc[(nd[0], ch), control+'_brlen'])
                                 for ch,d in t.items.get(nd[0], [])],
                 phyltree_methods.get_label,
                 root=anc,
                 edge_colors=lineage_rates.reset_index('taxon_parent', drop=True)[
                             (br_m, 'median')],
                 edge_cmap='afmhot', add_edge_axes=None, style='squared',
                 ax=ax)
    #cax = ax.inset_axes((0, 0.5, 0.05, 0.5))
    cbar = fig.colorbar(lines, cax=cax)
    cbar.set_label('Median '+br_m)
    cax.yaxis.set_ticks_position('left')
    cax.yaxis.set_label_position('left')
    ax.set_title('Median evolutionary %s rate across %s descendants' % (measures[0], anc)
                 + (('. '+extra_text) if extra_text else ''))
    ax.set_xlabel('Age (My)') #+ control
    outputs.append(fig)
    #fig.tight_layout()
    return outputs


def compute_branchrate_std(ages_controled, dist_measures,
                           branchtime='median_brlen_dS', taxon_age=None,
                           mean_condition=None,
                           std_condition=None,
                           weighted=False, poisson=True, debug=False):
    """
    Example filter_condition for approximated dS (nonlocal):
    '(calibrated==1) & (calibrated_parent==1)'

    Summary:
    - create cs_rates, containing mean_rate for each subgenetree
    - repeat the tree mean for each branch (construct aligned_csrates: not saved)
    - this allows computing the rate deviation, for each branch,
      then the std dev for each subgenetree
    """

    groupby_cols = ["subgenetree", "taxon_parent", "taxon",
                    branchtime] + dist_measures
    if taxon_age is not None:
        groupby_cols.append(taxon_age)  ## "median_age"/"median_age_dS"
    #ages_controled["omega"] = ages_controled.branch_dN / ages_controled.branch_dS

    ages_controled_full = ages_controled
    if mean_condition:
        ages_controled = ages_controled.query(mean_condition)
        if not ages_controled.shape[0]:
            logger.warning('Selected 0 rows with `mean_condition`.')
    elif not ages_controled.shape[0]:
        logger.warning('0 rows in data')
    sgg = subgenetree_groups = ages_controled[groupby_cols].groupby('subgenetree', sort=False)

    ### Average (substitution) rates over the tree:
    #     sum of all branch values / sum of branch lengths

    if debug:
        print('%d rows queried for mean (ages_controled)' % ages_controled.shape[0])
        check_nan(sgg[branchtime].sum(), 'sgg[branchtime].sum()')

    if poisson:
        # Sum aggregation + division broadcasted on columns
        # This is the maximum likelihood estimate for a constant Poisson rate.
        cs_rates = sgg[dist_measures].sum().div(sgg[branchtime].sum(), axis=0)
        # Problem: directly proportional to the tree length. Because branch times same across genes
    elif not weighted:
        cs_rates = ages_controled[dist_measures]\
                    .div(ages_controled[branchtime], axis=0)\
                    .join(ages_controled[['subgenetree']], sort=False)\
                    .groupby('subgenetree', sort=False)\
                    .mean()
    else:
        # Will be slow.
        cs_rates = sgg.apply(lambda g: g[dist_measures].mul(g[branchtime], axis=0))\
                        .div(sgg[branchtime].sum()**2, axis=0)
    #cs_rates["omega"] = (sgg.branch_dN / sgg.branch_dS).apply()

    if debug:
        check_nan(cs_rates, 'rate means')
        check_inf(cs_rates, 'rate means')

    # This weird operation fills the subgenetree-mean-rate into each node name.
    aligned_csrates = ages_controled[['subgenetree']]\
                            .join(cs_rates, on='subgenetree', sort=False)\
                            .drop(columns='subgenetree')

    if debug:
        check_nan(cs_rates, 'rate means (aligned to ages_controled)')
        check_inf(cs_rates, 'rate means (aligned to ages_controled)')

        diff_subgenetrees = set(ages_controled.subgenetree).difference(cs_rates.index)
        if diff_subgenetrees:
            logger.warning('%d unmatched subgenetrees in ages_controled: %s',
                           len(diff_subgenetrees), ' '.join(diff_subgenetrees))
            #The NaN should correspond to genetrees with zero node matching mean_condition?
    if debug == 'aligned_csrates':
        return aligned_csrates

    #check_nan(cs_rates, 'rate means (aligned to ages_controled)')
    #check_inf(cs_rates, 'rate means (aligned to ages_controled)')

    ### Weighted standard deviation of substitution rates among branches
    # the squared deviations from the mean rate:
    rate_dev = (ages_controled[dist_measures]\
                    .div(ages_controled[branchtime], axis=0)\
                    .sub(aligned_csrates, axis=0)**2)\
                .join(ages_controled[['subgenetree', branchtime]], sort=False)
    if debug:
        age_measures = [dm.replace('branch_', 'age_') for dm in dist_measures]
        debug_cols = (groupby_cols + age_measures +
                      [v+'_parent' for v in dist_measures + age_measures])

        #logger.debug('rate_dev.shape = %s; columns = %s', rate_dev.shape, rate_dev.columns)
        rate_dev_na_rows = check_nan(rate_dev, 'rate_dev')
        rate_dev_inf_rows = check_inf(rate_dev, 'rate_dev')
        if rate_dev_na_rows.any() or rate_dev_inf_rows.any():
            # Possible divisions 0/0
            for dist_m in dist_measures:
                div_0by0 = (ages_controled[[dist_m, branchtime]]==0).all(axis=1)
                print('%s/%s -> %d divisions 0/0' % (dist_m, branchtime,
                        div_0by0.sum()))
                print('data 0/0:')
                display_html(ages_controled.loc[div_0by0, debug_cols].head(20))

                if (rate_dev_na_rows & ~div_0by0).any():
                    print('%d Other NaN' % (rate_dev_na_rows & ~div_0by0).sum())
                # Possible divisions by 0
                div_by0 = (ages_controled[branchtime]==0) & ~div_0by0
                if div_by0.any():
                    print('%d divisions by 0 (by %s)' % (div_by0.sum(), branchtime))
                    print('data /0:')
                    display_html(ages_controled.loc[div_by0, debug_cols].head(20))
                    if (rate_dev_inf_rows & ~div_by0).any():
                        print('%d Other Inf' % (rate_dev_inf_rows & ~div_by0).sum())

    if debug == 'rate_dev':
        return rate_dev

    cs_rates.rename(columns=lambda d: d.replace('branch_', '')+'_rate', inplace=True)
    rate_measures = cs_rates.columns.tolist()  #[(m.replace('branch_', '') + '_rate') for m in dist_measures]

    #tmp = pd.merge(ages_controled[["subgenetree", branchtime] + dist_measures],
    #               cs_rates, left_on="subgenetree", right_index=True, sort=False)

    ### subtract branch rate with mean rate, then square.
    #rate_dev_dict = {d: (tmp[d] / tmp[branchtime] - tmp[r])**2
    #                 for d,r in zip(dist_measures, rate_measures)}

    #rate_dev = pd.DataFrame(rate_dev_dict)\
    #                .join(tmp[['subgenetree', branchtime]], sort=False)

    # Caching the column indices to apply raw numpy computations (faster).
    #dist_measures_idx = rate_dev.columns.get_indexer(dist_measures)
    #branchtime_idx = rate_dev.columns.tolist().index(branchtime)
    rsgg = rate_dev.groupby("subgenetree", sort=False)
    #rsgg_g0name = list(rsgg.groups.keys())[0]
    #rsgg_g0 = rsgg.get_group(rsgg_g0name)
    #logger.debug('group 0 "%s" shape = %s; columns = %s; name = %s',
    #             rsgg_g0name, rsgg_g0.shape, rsgg_g0.columns, getattr(rsgg_g0, 'name', None))
    #logger.debug('group 0 mean rate =\n%s;\nvalues = %s',
    #             cs_rates.loc[rsgg_g0name], rsgg_g0.values)
    #logger.debug('group 0 group_average = %s', npraw_group_average_w0(rsgg_g0.values))

    if weighted:
        #NOTE: The applied function omits NaN values in weights, but propagates NaNs in rates.
        #FIXME: with Pandas v1, fails if measures=['dS', 'dN', 't', 'dist'] (more than one measure).
        #       -> TypeError: Series.name must be hashable type
        cs_stds = rsgg[[branchtime] + dist_measures]\
                        .apply(lambda g: pd.Series(np.sqrt(
                                                    npraw_group_average_w0(g.values)),
                                                  name=g.name)
                        )
    else:
        cs_stds = np.sqrt(rsgg[dist_measures].mean())
    cs_stds.set_axis([r+'_std' for r in rate_measures], axis=1, inplace=True)

    #cs_stds = pd.DataFrame(
    #            np.stack(
    #                rsgg.apply(
    #                (lambda x:
    #                        sqrt(
    #                            npraw_group_average_w0(x.values)
    #                            )
    #                        )
    #                )
    #            ),
    #            columns=[r+'_std' for r in rate_measures])

    cs_rates = pd.concat((cs_rates, cs_stds), axis=1, sort=False, verify_integrity=True)
    # Checks
    inf_cols = np.isinf(cs_rates).any(axis=0)
    inf_rows = np.isinf(cs_rates).any(axis=1)  #check_inf(cs_rates)

    if inf_rows.any():
        logger.warning('%d Inf values in columns %s. DROPPING rows!',
                       inf_rows.sum(),
                       cs_rates.columns[inf_cols].tolist())
        cs_rates = cs_rates[~inf_rows].copy(deep=False)
    return cs_rates


def triplet_aggfunc(func, func_args, func_kwargs, parentgrouped, row):
    """To use within `apply` only! (because of .name attribute)"""
    try:
        children_rows = parentgrouped.get_group(row.name)
    except KeyError:
        #logger.debug('parent %r not found', row.name)
        return pd.Series(np.NaN, index=row.index)
    return pd.concat((row, children_rows), sort=False).agg(func, *func_args, **func_kwargs)


def raw_triplet_std(row, parentgrouped):
    try:
        children_rows = parentgrouped.get_group(row.name)
    except KeyError:
        #logger.debug('parent %r not found', row.name)
        return np.full(row.shape[0], np.NaN)
    return np.nanstd(np.vstack((row, children_rows)), axis=0, ddof=0)


def compute_correlated_rate(ages_controled, dist_measures,
                            branchtime='median_brlen_dS', taxon_age=None,
                            mean_condition=None,
                            std_condition=None):
    rates = ages_controled[dist_measures].div(ages_controled[branchtime], axis=0)
    sister_rates = rates.groupby(ages_controled.parent, sort=False)
    #compute_triplet_std = partial(triplet_aggfunc, 'std', (), {'ddof': 0}, sister_rates)
    #compute_triplet_std = lambda row: pd.Series(raw_triplet_std(row, sister_rates))
    #logger.debug('rates: type=%s; shape=%s; columns=%s', type(rates), rates.shape, rates.columns)
    #logger.debug('rates.head() =\n%s', rates.head(20))

    #TODO: subset ages_controled where 'type!="leaf"'
    triplet_rate_corrstds = rates.apply(raw_triplet_std, axis=1, raw=False,
                                    result_type='expand', args=(sister_rates,))

    groupby_cols = ["subgenetree", "taxon_parent", "taxon",
                    branchtime] + dist_measures
    if taxon_age is not None:
        groupby_cols.append(taxon_age)  ## "median_age"/"median_age_dS"

    if mean_condition:
        triplet_rate_corrstds = triplet_rate_corrstds.loc[ages_controled.query(mean_condition).index]
        if not triplet_rate_corrstds.shape[0]:
            logger.warning('Selected 0 rows with `mean_condition`.')
    elif not triplet_rate_corrstds.shape[0]:
        logger.warning('0 rows in data')

    #return triplet_rate_corrstds
    return triplet_rate_corrstds.groupby(ages_controled['subgenetree'], sort=False).mean()


def subset_on_criterion_tails(criterion_serie, ages=None, ages_file=None,
                              outbase=None, criterion_name=None, nquantiles=4,
                              thresholds=None, save=False):
    """From the input data, output two files:
    - one with the lower quantile of criterion values,
    - one with the upper quantile.

    Otherwise thresholds can be used as a tuple (lower_tail_max, upper_tail_min)

    output files will be named as `outbase + ("-lowQ%d" % nquantiles) + ".tsv"`
    """
    ages_c = merge_criterion_in_ages(criterion_serie, ages, ages_file, criterion_name)

    print(ages.columns)

    if not outbase and ages_file:
        outbase, _ = op.splitext(ages_file)
    elif not outbase and not ages_file:
        outbase = "dataset"

    if not thresholds:
        q = 1 / nquantiles
        low_lim = criterion_serie.quantile(q)
        high_lim = criterion_serie.quantile(1. - q)
        outlow = "%s-%slowQ%d.tsv" % (outbase, criterion_name, nquantiles)
        outhigh = "%s-%shighQ%d.tsv" % (outbase, criterion_name, nquantiles)
    else:
        low_lim, high_lim = thresholds
        outlow = outbase + "-%slow%1.3f.tsv" % (criterion_name, low_lim)
        outhigh = outbase + "-%shigh%1.3f.tsv" % (criterion_name, high_lim)

    print("Output %s values outside %s (from %d quantiles)" % (criterion_name,
                                    [low_lim, high_lim],
                                    (nquantiles if not thresholds else None)))
    if save:
        ages_c[ages_c[criterion_name] <= low_lim].to_csv(outlow, sep='\t')
        ages_c[ages_c[criterion_name] >= high_lim].to_csv(outhigh, sep='\t')
    print("Output files:", outlow, outhigh, sep='\n')


def get_tails_on_criterion(df, criterion_name, nquantiles=4):
    # WARNING: unused?
    q = 1. / nquantiles
    low_lim, high_lim = df[criterion_name].quantile([q, 1. - q])
    df_low = df[df[criterion_name] <= low_lim].copy()
    df_high = df[df[criterion_name] >= high_lim].copy()
    print(low_lim, high_lim)
    return df_low, df_high


def annot_quantiles_on_criterion(ages, criterion_serie, criterion_name=None,
                                 nquantiles=4, transform=None):
    assert criterion_name not in ages.columns

    # Exclude the points 0. and 1.
    q = np.linspace(1./nquantiles, 1, nquantiles-1, endpoint=False)

    if transform is not None:
        criterion_serie = transform(criterion_serie)
        criterion_name = transform.__name__ + '_' + criterion_name

    quantiles = criterion_serie.quantile(q)
    print("Quantiles (n=%d) values for %s:" % (nquantiles, criterion_name),
          quantiles, sep='\n')

    ages_c = merge_criterion_in_ages(criterion_serie, ages, criterion_name=criterion_name)
    ages_c["Q_" + criterion_name] = ages_c[criterion_name].apply(
                                lambda v: np.searchsorted(quantiles.values, v))
    return ages_c


def _violin_spe_ages_vs_criterion_quantiles(annot_df, criterion_name, isin=None,
                                           split=True, order=None, cut=0,
                                           points=False,
                                           **kwargs):
    Q_col = "Q_" + criterion_name
    if isin is None:
        # Look at extreme quantiles only
        isin = (annot_df[Q_col].min(), annot_df[Q_col].max())
        logger.info(isin)

    ax = sb.violinplot(x="taxon", y="age_dS", hue=Q_col,
                       data=annot_df[(annot_df.type == "spe")
                                     & annot_df[Q_col].isin(isin)],
                       order=order,
                       split=split,
                       cut=cut,
                       **kwargs)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, va='top', ha='right')
    if points:
        #violins = ax.get_children()
        orig_legend_args = ax.get_legend_handles_labels()
        orig_legend_title = ax.legend_.get_title().get_text()

        sb.pointplot(x="taxon", y="age_dS", hue=Q_col,
                     data=annot_df[(annot_df.type == "spe") & annot_df[Q_col].isin(isin)],
                     order=order,
                     estimator=np.nanmedian,
                     #ci='sd',
                     join=False,
                     dodge=(1./(len(isin)+1)),
                     palette=['#b8b8b8'],
                     alpha=0.8,
                     ax=ax,
                     legend=False)
        # Reset the legend.
        #print(ax.get_children())
        #plt.setp(set(ax.get_children()).difference(violins), label='')
        ax.legend_ = None
        ax.legend(*orig_legend_args, title=orig_legend_title)
    return ax


def violin_spe_ages_vs_criterion(ages, criterion_serie, criterion_name=None,
                                 nquantiles=10, split=False, order=None, **kwargs):
    criterion_name = criterion_name or criterion_serie.name
    annot_ages = annot_quantiles_on_criterion(ages, criterion_serie,
                                              criterion_name, nquantiles)
    isin = None if split else list(range(nquantiles))
    return _violin_spe_ages_vs_criterion_quantiles(
                                annot_ages[annot_ages.taxon != 'Simiiformes'],
                                criterion_name, isin, split, order, **kwargs)


class age_analysis_data(object):
    """Just a namespace to hold several objects of one analysis."""
    def __init__(self,
                 ages_controled,
                 ages_controled_withnonrobust,
                 ns,
                 control_ages,
                 control_brlen,
                 mean_errors):
        for k,v in locals().items():
            if k != 'self':
                setattr(self, k, v)


def analyse_age_errors(ages_file, root, phyltree, control_ages_CI, ts=None,
                       measures=['dS', 'dist', 'dN', 't'],
                       control_names=['timetree'], compute_extra_ns=None, out=None):
    #, aS, cs, clS=None
    ages, ns = load_prepare_ages(ages_file, ts, measures, compute_extra_ns)
    available_robust_info = set(ages.columns.intersection(('really_robust', 'aberrant_dists')))
    if len(available_robust_info) < 2:
        logger.warning('Available robustness info (%s) is lacking %s',
                       ','.join(available_robust_info),
                       ','.join(set(('really_robust', 'aberrant_dists')) - available_robust_info))

    control_condition_list = []
    if 'really_robust' in available_robust_info:
        control_condition_list.append('really_robust')
    if 'aberrant_dists' in available_robust_info:
        control_condition_list.append('(aberrant_dists == 0)')
    control_condition = ' & '.join(control_condition_list)

    ages_controled_withnonrobust, control_ages, control_brlen =\
        add_control_dates_lengths(ages, phyltree, control_ages_CI, measures,
                                  control_condition, control_names=control_names)
    if control_condition:
        ages_controled = ages_controled_withnonrobust.query(control_condition).copy(deep=False)
    else:
        ages_controled = ages_controled_withnonrobust
    unexpected_branches, lost_branches = check_control_dates_lengths(
                                                control_brlen, phyltree, root,
                                                measures, out=out)
    control_brlen.drop(unexpected_branches, inplace=True)
    # Rates
    #cs_rates = compute_branchrate_std(ages_controled, dist_measures)
    #cs_rates_approx
    #cs_rates_onlyAnc
    #cs_rates_withoutAnc

    mean_errors = pd.concat(
                [compute_dating_errors(ages_controled_withnonrobust, ctl, measures)
                 for ctl in ['median']+control_names],
                axis=1, sort=False, copy=False)

    return age_analysis_data(ages_controled, ages_controled_withnonrobust, ns,
                             control_ages, control_brlen, mean_errors)


def compare_params(x='taxon', y='age_dS', split=True, order=None, hue_order=None, width=1,
                   controls=('timetree', 'dosreis'), **datasets):
    """Splitted violinplot of taxon~age_dS with hue='data'.
    Use **datasets to name your datasets: ex um1=df1, um2=df2.

    From Jupyter Notebook `Compare_dating_methods.ipynb`.
    """
    stacked = pd.concat([df.assign(data=paramname) for paramname, df in datasets.items()], sort=False,
                         ignore_index=True)
    #display_html(stacked.head())
    # Show one value for x and data.
    hues = sorted(datasets.keys()) if hue_order is None else hue_order

    def kurtosis(a): return stats.kurtosis(a, axis=None, fisher=True, nan_policy='omit')

    scales = [iqr95, iqr90, iqr, mad, mean_absdevmed, 'std', trimstd]
    centers = ['median', 'median', 'median', 'median', 'median', 'mean', trimmean]

    subgroup_dispersion = stacked.groupby([x, 'data'])[y].agg(scales + ['skew', kurtosis, 'count']).unstack('data')
    subgroup_disptest = stacked.groupby(x)[[y, 'data']]\
        .apply(lambda g: pd.Series(stats.levene(*(g.loc[g.data==h, y].dropna() for h in hues))))\
        .set_axis(pd.MultiIndex.from_product([['BrownForsythe'], ['stat', 'pvalue']]), axis=1, inplace=False)
    subgroup_disptest_corr = pd.DataFrame(smm.multipletests(subgroup_disptest[('BrownForsythe', 'pvalue')],
                                                            alpha=0.01, method='fdr_bh')[:2],
                                          index=pd.MultiIndex.from_product([['correction'],
                                                                            ['rejectH0', 'pval_corr']]),
                                          columns=subgroup_disptest.index).T
    subgroup_disp = pd.concat((subgroup_dispersion, subgroup_disptest, subgroup_disptest_corr), axis=1, sort=False)

    #scalename = getattr(scale, '__name__', scale)
    display_html(
        subgroup_disp\
        .style.highlight_min(axis=1, subset=[('mad', d) for d in datasets])\
        .highlight_min(axis=1, subset=[('mean_absdevmed', d) for d in datasets], color='gold')\
        .highlight_min(axis=1, subset=[('iqr', d) for d in datasets], color='lemonchiffon')\
        .highlight_min(axis=1, subset=[('std', d) for d in datasets], color='khaki')\
        .highlight_min(axis=1, subset=[('trimstd', d) for d in datasets], color='moccasin')\
        .highlight_max(axis=1, subset=[('skew', d) for d in datasets], color='coral')\
        .highlight_max(axis=1, subset=[('kurtosis', d) for d in datasets], color='sandybrown'))#\
        #.apply(lambda values: ['color: red' if v is True else '' for v in values],
        #       subset=pd.IndexSlice[:, pd.IndexSlice[('correction', 'rejectH0')]]))

    #assert center in ('median', 'mean')
    ltests = []
    for center,scale in zip(centers, scales):
        ltests.append(multi_vartest(y, 'data', by=x, data=stacked, ref_hue=hues, center=center, scale=scale))
        print('Multi var test (%s, %s):' %(getattr(center, '__name__', center), getattr(scale, '__name__', scale)),
              ltests[-1], '\n-----')

    ax = sb.violinplot(x, y, hue='data', data=stacked, order=order, hue_order=hue_order,
                       split=split, width=width, cut=0)
    plt.setp(ax.get_xticklabels(), rotation=45, va='top', ha='right')
    violin_handles, violin_labels = ax.get_legend_handles_labels()

    # medians for each hue group.
    sb.pointplot(x, y, hue='data', data=stacked, order=order, hue_order=hue_order, dodge=width/(len(datasets)+1),
                 join=False, palette=['#4c4c4c'], estimator=np.nanmedian,
                 ci='sd', ax=ax)
    current_handles, current_labels = ax.get_legend_handles_labels()
    point_handles = [h for h in current_handles if h not in violin_handles]

    # Reference values
    ylim = ax.get_ylim()
    #ax.plot('timetree_age', 'r_', data=control_ages.loc[ordered_simii_anc],
#     ax.plot([phyltree.ages[anc] for anc in ordered_simii_anc], 'r_', label='TimeTree age',
#             markersize=24, linewidth=24)
#     ax.plot('timetree_CI_inf', 'y_', data=timetree_ages_CI.loc[ordered_simii_anc],
#             markersize=18)
#     ax.plot('timetree_CI_sup', 'y_', data=timetree_ages_CI.loc[ordered_simii_anc],
#             markersize=18)
    X = np.arange(len(stacked[x].unique()) if order is None else len(order))
    xlines = [x for xanc in X for x in [xanc-0.5, xanc+0.5, None]]

    for ctl in controls:
        ylines = [y for anc in ordered_simii_anc for y in [calibrations[ctl+'_age'].loc[anc]]*2 + [None]]
        color = 'y' if ctl == 'timetree' else 'b'
        lines, = ax.plot(xlines, ylines, color+'--', alpha=0.9, label=ctl+' age')
        #print(type(lines), lines)

        ax.bar(X, (calibrations[ctl+'_CI95_sup'] - calibrations[ctl+'_CI95_inf']).loc[order].values,
               bottom=calibrations[ctl+'_CI95_inf'].loc[order],
               color='y', alpha=0.5, zorder=(-1 if ctl=='timetree' else -2), width=0.95,
               label=ctl+' 95% interval')

    handles, labels = ax.get_legend_handles_labels()
    handles, labels = zip(*sorted([(h,l) for h,l in zip(handles, labels) if h not in point_handles],
                                  key=lambda hl: hl[0] not in violin_handles))
    #print(handles)
    ax.legend(handles, labels)
    return stacked, subgroup_disp, ltests


# Copied from Compare_dating_methods20200318.py [2020/05/25 11:09]
def plottree_add_hists(taxon_col, x, positions, data=None, nbins=100,
                       order=None, ax=None, zorder=None, labels=None,
                       invert=False, scale=0.95, alpha=0.7, recolor_samepos=True,
                       **kwargs):
    """
    x: variable whose histograms are drawn (column names from data).
    taxon_col: groups to make histograms;
    positions: histogram bottoms from each group;
    order: ordered list of taxa for positioning.

    invert: Whether the first histogram is plotted *under* the branch
    (useful for inverted yaxis.)"""
    # positions is a bad name. It refers to the bottom of the histograms.
    broadcasted = [x]
    if data is not None:
        broadcasted.append(data)
    if labels is not None:
        broadcasted.append(labels)
    input_lists = [isinstance(arg, (list, tuple)) for arg in broadcasted]
    if any(input_lists):
        if not all(input_lists):
            raise ValueError('(x, data, labels) should all be lists/tuples or '
                             'all single elements.')
        if data is None:
            data = [None]*len(x)
        if labels is None:
            labels = [None]*len(x)
    else:
        x = [x]
        data = [data]
        labels = [labels]

    if len(x) > 2:
        raise NotImplementedError("Can't handle more than 2 histograms.")

    if ax is None:
        ax = plt.gca()
    xbound = ax.get_xbound()

    assert positions is not None, "positions argument is mandatory (list/array)"
    positions = np.asarray(positions)
    if order is not None: order = np.asarray(order)

    # Set a different color to overlapping positions
    if recolor_samepos:
        step = 0.5
        # Group together positions closer than 0.5
        pos_bins = np.arange(round(min(positions)) - step/2.,
                             round(max(positions)) + step/2.,
                             step)
        #pos_hist, _ = np.histogram(positions, pos_bins)
        pos_group = np.digitize(positions, pos_bins, right=True)

        distinct_positions = [np.unique(pos_group, return_index=True)[1]]
        next_position_index = setdiff1d(np.arange(len(positions), dtype=int),
                                        distinct_positions[-1])
        # At each iteration: exclude overlapping positions into a new set of positions.
        while next_position_index.size:
            # Extract a set of positions, non-overlapping between them.
            # I.e. Get the first time each pos_group is seen:
            #uniq_pos, uniq_pos_index = np.unique(pos_group, return_index=True)
            #distinct_positions.append(positions[uniq_pos_index])
            filled_pos, filled_pos_index = np.unique(pos_group[next_position_index],
                                                     return_index=True)
            distinct_positions.append(next_position_index[filled_pos_index])
            # Delete them
            next_position_index = setdiff1d(next_position_index, distinct_positions[-1])
            #pos_hist -= 1
            #pos_group[uniq_pos_index] = np.NaN  # NaN is always different from itself, so is always unique
    else:
        distinct_positions = [np.arange(len(positions))]  # FIXME, should be integer

    #for i, (vtaxon, xi) in enumerate(zip(var_taxon, var_xs), start=1):
    bars = []
    for i, (xi, dat, lab) in enumerate(zip(x, data, labels), start=(int(invert))):
#         print('* Plot data %r: (shape: %s)\n'
#                   '  positions: %s\n'
#                   '  order: %s\n'
#                   '  range: %s\n'
#                   '  bins: %s' % (x[i-1], getattr(xi, 'shape',None),
#                                   positions, order, xbound, bins))
        for distinct_pos in distinct_positions:
            logger.debug('# i=%d. %d positions; distinct_pos=%s dtype %s; order=%s',
                    i, len(positions), distinct_pos, distinct_pos.dtype, order)
            logger.debug('Plotting taxa dupli:\n%s',
                 dat.groupby(taxon_col)[xi].count()[
                     slice(None) if order is None else np.asarray(order)[distinct_pos]])  #FIXME: distinct_pos invalid indexer when recolor_samepos=False
            _, bin_edges, b = cathist(
                    #var_taxon, xi,
                    taxon_col, xi, data=dat,
                    bins=nbins,
                    positions=positions[distinct_pos],
                    scale=(scale * (-1)**i),
                    order=None if order is None else order[distinct_pos],
                    range=xbound,
                   ax=ax, zorder=zorder, label=lab, alpha=alpha, **kwargs)
            bars.append(b)
    return bars
#        print('%d bars, %d bin_edges: %s' % (len(bars), len(bin_edges), bin_edges))


def display_errors(params, control='timetree', age_col='age_dS',
                   plotsize=(12, 5), save=False, plot=True):
    """From notebook `Speciations.ipynb`"""
    global ordered_simii_anc, phyltree, age_analysis
    
    mean_abs_errors            = np.full((len(ordered_simii_anc), len(params)), np.NaN)
    scaled_mean_abs_errors     = np.full((len(ordered_simii_anc), len(params)), np.NaN)
    asymscaled_mean_abs_errors = np.full((len(ordered_simii_anc), len(params)), np.NaN)
    
    std_abs_errors            = np.full((len(ordered_simii_anc), len(params)), np.NaN)
    scaled_std_abs_errors     = np.full((len(ordered_simii_anc), len(params)), np.NaN)
    asymscaled_std_abs_errors = np.full((len(ordered_simii_anc), len(params)), np.NaN)
    
    for i,anc in enumerate(ordered_simii_anc):
        
        #FIXME
        real_age = calibrations.loc[anc, "%s_age" % control]  #phyltree.ages[anc]
        
        # Asymmetric error measure
        next_calib = max(phyltree.ages[ch] for ch,_ in phyltree.items[anc])
        parent_calib = phyltree.ages[phyltree.parent[anc].name]

        # Symmetric error measure
        calib_scale = parent_calib - next_calib

        #fig, axes = plt.subplots(len(params), figsize=(plotsize, 5*len(params)), sharex=True, sharey=True)
        if plot:
            fig = plt.figure(figsize=plotsize)
            ax = None
    
        nax = sum(1 for p in params if p in age_analyses[anc])
        for j, param in enumerate(params):
            try:
                data = age_analyses[anc][param].ages_controled.query('calibrated==0')[age_col]
            except KeyError as err:
                logger.warning("! age_analysis[%r][%r] => %r" % (anc, param, err))
                continue
            
            #measure = 'dS'
            #dev_data = age_analyses[anc][param].mean_errors[('' if control=='median' else control) + 
            #                                                '_abs_dev_age_' + measure]
            #dev_data.mean() / calib_scale
            
            # Symmetric
            errs = (real_age - data).abs()
            err = errs.mean()
            scaled_err = err / calib_scale
            print('# %-25s %-8s:\tscaled (sym)  mean_abs_error = %.4f' % (anc, param, scaled_err))

            mean_abs_errors[i, j] = err
            scaled_mean_abs_errors[i, j] = scaled_err
            std_abs_errors[i, j] = errs.std()
            scaled_std_abs_errors[i, j] = std_abs_errors[i, j] / calib_scale
            # Asymmetric

            #dev_data = age_analyses[anc][param].mean_errors[('' if control=='median' else control) + 
            #                                                '_signed_dev_age_' + measure]
            #scaled_errs = np.where(dev_data<0, dev_data/(next_calib-real_age), dev_data/(parent_calib - real_age))
            
            errs = (data - real_age)
            scaled_errs = np.where(errs<0, errs/(next_calib-real_age), errs/(parent_calib - real_age))
            scaled_err = np.nanmean(scaled_errs)
            asymscaled_std_abs_errors[i, j] = np.nanstd(scaled_errs, ddof=1)
            
            print('# %-25s %-8s:\tscaled (asym) mean_abs_error = %.4f\t[%7.4f ← %7.4f → %7.4f]' % (anc, param,
                                    scaled_err, next_calib, real_age, parent_calib))

            asymscaled_mean_abs_errors[i, j] = scaled_err

            # Plot distrib
            if not plot:
                continue
            ax = fig.add_subplot(nax, 1, j+1, sharex=ax)
            data.hist(bins=50, ax=ax)
            
            y = ax.get_ylim()
            ax.plot([real_age]*2, y, label=('%s age' % control))
            ax.plot([data.median()]*2, y, label='median')
            ax.plot([data.mean()]*2, y, label='mean')
            ax.legend()
            ax.set_title(param)
            
        if plot:
            fig.suptitle(anc)
            fig.set_size_inches(plotsize[0], nax*plotsize[1])
            fig.show()
            if save:
                fig.savefig('../fig/ages_%s.pdf' % anc, bbox_inches='tight')
                
    return (pd.DataFrame(scaled_mean_abs_errors, index=ordered_simii_anc, columns=params),
            pd.DataFrame(scaled_std_abs_errors, index=ordered_simii_anc, columns=params),
            pd.DataFrame(asymscaled_mean_abs_errors, index=ordered_simii_anc, columns=params),
            pd.DataFrame(asymscaled_std_abs_errors, index=ordered_simii_anc, columns=params))


def check_nan(df, description=None, out=None):
    if description is None: description = ''
    if isinstance(df, pd.Series):
        df = df.to_frame()
    print('%s -> Any NA: %s' % (description, df.isna().sum(axis=0).any()), file=out)
    print('%s -> *All* NA: %s' % (description, df.columns[df.isna().all()].tolist()), file=out)
    print('Shape %s: %s' % (description, df.shape), file=out)
    na_rows = df.isna().any(axis=1)
    if na_rows.any():
        print(('%d NA rows' % na_rows.sum()) +
              (' (%s)' % description if description else ''), file=out)
        na_counts = df.isna().sum().sort_values(ascending=False)
        print('Columns with NA%s:\n%s'
              % ((' (%s)' % description if description else ''),
                 na_counts[na_counts>0]), file=out)
        #self.a_t.dropna(inplace=True)
    return na_rows

def check_inf(df, description=None, out=None):
    if description is None: description = ''
    if isinstance(df, pd.Series):
        df = df.to_frame()
    dfn = df.select_dtypes(np.number)  # Numeric columns only.
    inf_rows = np.isinf(dfn).any(axis=1)
    if inf_rows.any():
        print('%s -> %d Inf rows' % (description, inf_rows.sum()), file=out)
        print('Data with Inf: columns %s' % dfn.columns[np.isinf(dfn).any()], file=out)
    return inf_rows

def check_nonnumeric(df, description=None, out=None):
    non_numbool = df.select_dtypes(exclude=['number', 'bool'])
    if non_numbool.shape[1]:
        print('Non numeric/bool dtypes (please validate):', non_numbool.dtypes, file=out)
        print(non_numbool.head(), file=out)

def check_dupindex(df, description=None, out=None):
    if df.index.has_duplicates or df.columns.has_duplicates:
        n_dup_rows = df.index.duplicated().sum()
        dup_cols = df.columns[df.columns.duplicated()]
        logger.warning('df has duplicates in index: %d, columns: %d (%s)',
                       n_dup_rows, len(dup_cols), dup_cols.tolist())

def check_constants(df, description=None, out=None):
    #dfnum = df.select_dtypes('number')
    #dfbool = df.select_dtypes('bool').astype(np.uint8)
    df = df.select_dtypes(['number', 'bool'])
    with warnings.catch_warnings(record=True) as ws:
        constant_vars = df.columns[(df.nunique() == 0)].tolist()
    #for w in ws:
    #    if not (issubclass(w.category, FutureWarning)
    #            and ("Method .ptp is deprecated and will be removed"
    #                 in str(w.message))):
    #        # Then resend the warning
    #        warnings.warn(w.category(w.message))
    if constant_vars:
        logger.warning('Constant features: %s', ' '.join(constant_vars))
    return constant_vars


def sanity_check(df, description=None, out=None):
    check_dupindex(df, description, out)
    check_nonnumeric(df, description, out)
    check_inf(df, description, out)
    check_nan(df, description, out)
    #check_constants(df, description, out)

def standardize_dataframe(df, ref_df=None):
    if ref_df is None: ref_df = df
    dfstd = ref_df[df.columns].std(ddof=1)
    assert isinstance(dfstd, pd.Series), type(dfstd)
    assert (set(dfstd.index) == set(df.columns)), dfstd.index
    assert not dfstd.index.has_duplicates
    if (dfstd == 0).any():
        logger.warning('std = 0 at: %s', dfstd.index[dfstd==0].tolist())

    return df.div(dfstd, axis=1)

def zscore_dataframe(df, ref_zscore=None):
    if ref_zscore is None:
        ref_zscore = df
    #a = np.divide(
    #    np.subtract(df.values, np.nanmean(df.values, axis=0), axis=1),
    #    np.nanstd(df.values, axis=0),
    #    axis=1)
    dfnum = df.select_dtypes(np.number)  #exclude=
    if set(df.columns.difference(dfnum.columns)):
        #logger.error(
        raise ValueError("Non-numeric columns can't be zscored (=>dropped): %s"
                         % (df.columns.difference(dfnum.columns)))
        #df = dfnum

    dfmean = ref_zscore[df.columns].mean()
    assert isinstance(dfmean, pd.Series), type(dfmean)
    assert (set(dfmean.index) == set(df.columns)), df.columns.difference(dfmean.index)
    assert (not dfmean.index.has_duplicates), dfmean.index[dfmean.index.duplicated()]
    dfcentered = df.sub(dfmean, axis=1)
    assert isinstance(dfcentered, pd.DataFrame), type(dfcentered)
    assert (dfcentered.shape == df.shape), dfcentered.shape
    assert (set(dfcentered.index) == set(df.index)), dfcentered.index
    assert (set(dfcentered.columns) == set(df.columns)), dfcentered.columns
    assert not dfcentered.index.has_duplicates
    assert not dfcentered.columns.has_duplicates
    return standardize_dataframe(dfcentered, ref_zscore)


def lassoselect_and_refit(a, y, features, atol=1e-2, method='elastic_net',
                          alpha=0.01, L1_wt=1, cov_type='HC1', model=sm.OLS,
                          maxiter=50, out=None, **psummary_kw):
    exog = sm.add_constant(a[features])
    if 'const' not in exog:
        logger.warning('No constant added to `exog`: some features are already constant')
    regul_params = dict(method=method, alpha=alpha, start_params=None)
    if model.__name__ == 'OLS':
        regul_params.update(L1_wt=L1_wt)
        fit_params = dict(cov_type=cov_type)
    elif model.__name__ == 'Logit':
        fit_params = dict(method='bfgs', warn_convergence=True, maxiter=maxiter)
        # BFGS does not inverse the Hessian, so it will not fail if the matrix is singular.
        regul_params['maxiter'] = maxiter
        if np.isscalar(alpha):
            # The doc of Logit suggests that scalar alpha will not penalize anything.
            regul_params['alpha'] = np.full(len(features)+1, alpha)
    # NOTE:
    # valid methods for sm.Logit are 'l1' and 'l1_cvxopt_cp'

    modelled = model(a[y], exog)
    fitlasso = modelled.fit_regularized(**regul_params)
    if method=='sqrt_lasso':
        alpha = 1.1 * np.sqrt(exog.shape[0]) * stats.norm.ppf(1 - 0.05 / (2*exog.shape[1]))
    print('Regularized fit: %s alpha=%g L1_wt=%g' % (method, alpha, L1_wt), file=out)

    NaN_coefs = fitlasso.params.isna()
    if NaN_coefs.all():
        raise RuntimeError("ALL coefs estimated as NaN.")
    elif NaN_coefs.any():
        logger.warning("Some coefs estimated as NaN.")

    almost_zero = np.isclose(fitlasso.params.drop('const'), 0, atol=atol)
    print('Non zeros params: %d/%d (from %d input features).\n'
          'Almost zeros params (<%g): %d' % (
            (fitlasso.params.drop('const') != 0).sum(),
             fitlasso.params.drop('const').shape[0],
             len(features),
             atol,
             almost_zero.sum()), file=out)

    sorted_coefs = fitlasso.params.drop('const').abs().sort_values(ascending=False)
    sorted_features = sorted_coefs.index.tolist()
    multicol_test_cumul = pd.DataFrame(columns=['cumul_colinearity',
                                                'cumul_colinearity_noconst'],
                               index=['const'] + sorted_features,
                               dtype=float)
    multicol_test_cumul.loc['const', 'cumul_colinearity'] = 1
    for i, ft in enumerate(sorted_features, start=1):
        multicol_test_cumul.loc[ft] = [multicol_test(exog[['const'] + sorted_features[:i]]),
                                       multicol_test(exog[sorted_features[:i]])]
        if multicol_test_cumul.loc[ft].isna().any():
            logger.error("NaN multi-colinearity condition numbers at %s", ft)

    pslopes = sm_pretty_summary(fitlasso, multicol_test_cumul, out=out)
    try:
        out.html(pslopes)  #with out=HtmlReport()
    except AttributeError:
        display_html(pslopes)

    #scatter_density('ingroup_glob_len', y, data=a, alpha=0.5);
    #sb.violinplot('triplet_zeros_dS', 'abs_age_dev', data=a);

    print('\n#### %s refit of Lasso-selected variables (%g)' % (model.__name__, atol), file=out)
    #TODO: delete but make a refit.
    selected_features = sorted_coefs[sorted_coefs > atol].index.tolist()
    print('Removed variables:', ', '.join(sorted_coefs[sorted_coefs <= atol].index), file=out)
    fit = model(a[y], exog[['const'] + selected_features]).fit(**fit_params)

    bars_to_show = ['coef', 'Lasso coef', 'Simple regression coef']
    bars_to_show += psummary_kw.pop('bars', [])
    bars_mid_to_show = ['VIFs'] + psummary_kw.pop('bars_mid', [])
    info_to_join = psummary_kw.pop('join', None)
    # We **must** apply a multiple test correction, even if this is a multiple regression,
    # because we have first performed a variable selection step (post-selection inference).
    # Benjamini-Hochberg is not stringent enough (the number of selected features is too low)
    # We can instead do a Bonferroni by considering the total number of *initial* features.
    adj_pvalues = fit.pvalues.loc[selected_features] * len(features)
    simple_regressions =  [model(a[y], sm.add_constant(a[selected_features])[['const']]).fit(**fit_params)]
    simple_regressions += [model(a[y], sm.add_constant(a[ft])).fit(**fit_params)
                           for ft in selected_features]
    simple_reg_adj_pvals = [reg.pvalues[ft] * len(features) for ft, reg in
                            zip(selected_features, simple_regressions[1:])]
    partial_r2s = leave1out_eval(['const']+selected_features,
                                 partial(partial_r_squared_fromreduced, fit))
    vifs = [VIF(fit, ft) for i,ft in enumerate(['const']+selected_features)]
    r2_attr = 'prsquared' if model.__name__ == 'Logit' else 'rsquared'
    param_info = pd.concat((
                    adj_pvalues.rename('adj. p-value'),
                    pslopes.data.coef.rename('Lasso coef'),
                    pd.DataFrame([(reg.params.loc[ft],
                                   getattr(reg, r2_attr),
                                   reg.pvalues[ft],
                                   *reg.conf_int().loc[ft])
                                  for ft, reg in zip(['const'] + selected_features,
                                                     simple_regressions)],
                                 index=['const'] + selected_features,
                                 columns=['Simple regression coef',
                                          'Simple regression R2',
                                          'Simple regression p-value',
                                          'Simple regression [0.025',
                                          'Simple regression 0.975]']),
                    pd.Series(vifs, index=['const']+selected_features, name='VIFs'),
                    pd.Series(simple_reg_adj_pvals, index=selected_features, name='Simple regression adj. p-value')
                    ) + (() if info_to_join is None else (info_to_join,)),
                    axis=1, sort=False)

    preslopes = sm_pretty_summary(fit, param_info, bars=bars_to_show,
                                  bars_mid=bars_mid_to_show, out=out, **psummary_kw)
    try:
        out.html(preslopes)  #with out=HtmlReport()
    except AttributeError:
        display_html(preslopes)

    assert set(selected_features) == set(preslopes.data.drop('const').index)
    return fitlasso, fit, pslopes, preslopes


class fullRegression(object):
    init_vars = ['same_alls',
                 'responses',
                 'features',
                 'ref_suggested_transform',
                 'impose_transform',    # global _must_transform
                 'to_decorr',           # global _must_decorr
                 'must_drop_features',  # global _must_drop_features
                 'protected_features',  # global _protected_features
                 'must_drop_data',      # global _must_drop_data
                 'regul_kw',
                 'out', 'logger']  #, 'widget'

    init_defaults = {'ref_suggested_transform': dict,
                     'impose_transform': dict,
                     'must_drop_features': set,
                     'to_decorr': Args,
                     'protected_features': set,
                     'must_drop_data': dict,
                     'regul_kw': dict}

    def __init__(self, same_alls, responses, features,
                 ref_suggested_transform=None, impose_transform=None,
                 must_drop_features=None, to_decorr=None,
                 protected_features=None, must_drop_data=None, regul_kw=None,
                 out=None, logger=None, widget=None):
        for k,v in locals().items():
            if k != 'self':
                if k in self.init_defaults and v is None:
                    v = self.init_defaults[k]()  # initialize to the proper type.
                setattr(self, k, v)
        self.displayed = []  # Figures and styled dfs
        self.set_output()

    def set_output(self):
        try:
            self.show = self.out.show
        except AttributeError:
            logger.debug("Output %s does not have a .show method", self.out)
            self.show = plt.show
        try:
            self.display = self.out.display
        except AttributeError:
            logger.debug("Output %s does not have a .display method", self.out)
            self.display = display
        try:
            self.display_html = self.out.html
        except AttributeError:
            logger.debug("Output %s does not have a .html method", self.out)
            self.display_html = display_html
        #print('Setting logger: current __name__ = %r; ' % __name__,
        #      'current self.__module__ = %r' % self.__module__)
        if self.logger is None:
            # For some reason, here __name__ == 'builtins'.
            self.logger = logging.getLogger(self.__module__)

    @classmethod
    def from_other(cls, other_regression):
        """Instantiate by copying all attributes from another instance."""
        self = cls(**{ivar: getattr(other_regression, ivar)
                               for ivar in cls.init_vars})
        #vars() doesn't contain properties nor private attributes
        for k,v in vars(other_regression).items():
            if k not in cls.init_vars:
                try:
                    setattr(self, k, v)
                except AttributeError as err:
                    err.args += (k,)
                    raise
        self.set_output()
        return self


    @property  #@dependency
    def suggested_transform(self):
        try:
            return self._suggested_transform
        except AttributeError:
            responses, features = self.responses, self.features
            ref_suggested_transform = self.ref_suggested_transform
            impose_transform = self.impose_transform
            alls = self.alls

            suggested_transform = test_transforms(alls,
                                    [ft for ft in responses+features
                                        if ft not in impose_transform],
                                    out=self.out, widget=self.widget)
            suggested_transform.update(**dict(
                                        ( ft, t(alls[ft]) )
                                        if (t.__name__ == 'make_best_logtransform' and ft in alls)
                                        else (ft,t)
                                        for ft,t in impose_transform.items())
                                      )
            # List the features that were skipped because constant, or non numeric/bool
            self.features_constant = [ft for ft in responses+features if ft not in suggested_transform]


            if ref_suggested_transform:
                onlyref = []
                onlynew = []
                diff_funcs = []
                diff_args = []
                for k in set(ref_suggested_transform).union(suggested_transform):
                    if k not in ref_suggested_transform:
                        onlynew.append(k)
                    elif k not in suggested_transform:
                        onlyref.append(k)
                    else:
                        reffunc = ref_suggested_transform[k].__name__
                        newfunc = suggested_transform[k].__name__
                        if reffunc != newfunc:
                            reffunc0, *refargs = reffunc.split('(')
                            newfunc0, *newargs = newfunc.split('(')
                            if reffunc0 == newfunc0:
                                newargs = newargs[0] if newargs else ''
                                refargs = refargs[0] if refargs else ''
                                diff_args.append((k, '%s VS %s' % (newargs.rstrip(')'),
                                    refargs.rstrip(')'))))
                            else:
                                diff_funcs.append((k, '%s VS %s : *updating*.' % (reffunc, newfunc)))
                                # *Updating*
                                suggested_transform[k] = ref_suggested_transform[k]
                if onlynew:
                    print('Transforms only in new:', ', '.join(onlynew), file=self.out)
                if onlyref:
                    print('Transforms only in ref:', ', '.join(onlyref), file=self.out)
                if diff_funcs:
                    print('Different transforms:\n', '\n'.join('%35s\t%s' % t for t in diff_funcs), file=self.out)
                if diff_args:
                    print('Different transform args:\n', '\n'.join('%35s\t%s' % t for t in diff_args), file=self.out)

            for ft in list(suggested_transform.keys()):
                delete_hardcoded = []
                if ft not in alls.columns:
                    delete_hardcoded.append(ft)
                    suggested_transform.pop(ft)
            if delete_hardcoded:
                self.logger.warning('%d hardcoded features not available, delete: %s', len(delete_hardcoded), ', '.join(delete_hardcoded))

            self.features_delete_hardcoded = delete_hardcoded
            self._suggested_transform = suggested_transform
            return suggested_transform


    def prepare_design_matrix(self):
        """This method should be overriden for additional design matrix preprocessing."""
        self.alls = self.same_alls.copy()

    #def fix_decorr(self):
    #    """This method should be overriden for case-by-case customisation of the instance:
    #    from types import MethodType
    #    fullreg.fix_decorr = MethodType(mynewfixdecorr)
    #    """
    #    pass

    def do(self):
        self.do_undecorred_fit()
        self.do_decorr()
        # You can fix the decorrelation inbetween these steps
        self.do_decorred_fit()
        reslopes2 = self.do_bestfit()
        self.do_extreme_points()
        return reslopes2

    def do_transform(self, ref_rescale=None, rescale=True):
        logger = self.logger

        features = self.features
        ft_idx = pd.Index(features)
        assert (not ft_idx.has_duplicates), ft_idx[ft_idx.duplicated()]

        responses = self.responses
        y = responses[0]
        print('Variable Y :', y, file=self.out)

        self.prepare_design_matrix()
        alls = self.alls
        print('%d observation × %d columns; %d input features.' % (alls.shape + (len(set(features).intersection(alls.columns)),)),
              file=self.out)
        if set(features).difference(alls.columns):
            raise KeyError('Not in data: %s' % set(features).difference(alls.columns))

        self.na_amount = na_amount = alls.isna().sum(axis=0).sort_values(ascending=False)
        print('Amount of NA:', na_amount.head(10), sep='\n', file=self.out)
        print('Amount of Inf:', np.isinf(alls.select_dtypes(np.number)).sum(axis=0).sort_values(ascending=False).head(10), sep='\n', file=self.out)

        too_many_na = na_amount.index[(na_amount >= 0.9*alls.shape[0])].tolist()
        many_na = na_amount.index[(na_amount >= 0.5*alls.shape[0])].tolist()
        if too_many_na:
            logger.warning('**DROP** these columns with >90%% of NA: %s',
                           ' '.join(too_many_na))
            alls.drop(columns=too_many_na, inplace=True)
            self.features = features = [ft for ft in features if ft not in too_many_na]
            many_na = [ft for ft in many_na if ft not in too_many_na]
        if many_na:
            logger.warning('Columns with >50%% of NA (kept as is): %s',
                           ' '.join(many_na))
        self.features_too_many_na = too_many_na

    #def do_transforms(self):
        # Calling the cached property: depends on self.alls and self.features.
        suggested_transform = self.suggested_transform

        # Remove constant features
        self.features = features = [ft for ft in features if ft in suggested_transform]
        self.features_posttransform_constant = check_constants(alls[features], 'After suggest transforms', self.out)
        if self.features_posttransform_constant:
            logger.warning('Dropping constants %s', self.features_posttransform_constant)
            for ft in self.features_posttransform_constant:
                features.remove(ft)
        
        # All binary variables should **NOT** be z-scored!
        self.bin_features = bin_features = [ft for ft in responses+features
                            if (suggested_transform[ft].__name__.startswith('binarize') or
                                np.issubdtype(alls[ft].dtype, bool))]
        zscored_features = [ft for ft in responses+features if ft not in bin_features]

        if alls.isin((-np.Inf, np.Inf)).any(axis=None):
            logger.warning('Replace Inf,-Inf by NaN')
        alls_transformed = alls.replace([-np.Inf, np.Inf], np.NaN).transform(suggested_transform)

        #print('transformed scalar types:', alls_transformed.iloc[0].map(type).tolist(), file=self.out)
        check_nonnumeric(alls_transformed, 'after transform', self.out)

        self.a_t = a_t = alls_transformed

        self.na_rows_t = check_nan(a_t, 'after transform', out=self.out)
        self.inf_rows_t = check_inf(a_t, 'after transform', out=self.out)
        constant_vars_t = check_constants(a_t, 'after transform', out=self.out)

    #def do_norm(self):
        logger.debug('Will Z-score: from %r a_t %s %s', type(a_t), a_t.shape,
                     [ft for ft in responses+features if ft not in bin_features])
            #{ft: zscore for ft in responses+features if ft not in bin_features})
        logger.debug('Will join binary features: %s', bin_features)
        bin_responses = set(responses).intersection(bin_features)
        if bin_responses:
            # Special case: it is not desirable to standardize binary response variable,
            # for example the Logit regression requires all values to be 0 or 1.
            logger.info('Binary responses will NOT be standardized (divided by stddev): %s',
                        ','.join(bin_responses))

        # This deep shit of .transform() is creating a RecursionError.
        # Also tried: zscore(df) but specifically fails on a_t.
        # Also tried .apply(zscore, raw=True, result_type='broadcast')
        # Back to good ol' numpy:
        # Ok, the problem was duplicates in index. Well no.
        check_dupindex(a_t, 'after transform', self.out)

        if ref_rescale is None:
            ref_rescale = a_t

        if rescale:
            if bin_responses:
                for ft in bin_responses:
                    bin_features.remove(ft)
            a_n = zscore_dataframe(a_t[zscored_features], ref_rescale[zscored_features])\
                 .join(standardize_dataframe(a_t[bin_features], ref_rescale[bin_features]))\
                 .join(a_t[list(bin_responses)].copy())

        else:
            a_n = a_t[zscored_features+bin_features].copy()

        #a_n = pd.concat()
        self.na_rows_n = na_rows_n = check_nan(a_n, 'after zscore', out=self.out)

        print('Check a_n.head():', file=self.out)
        self.display_html(a_n.head())
        if na_rows_n.any():
            print('Drop %d NA rows (after zscoring)' % na_rows_n.sum(), file=self.out)
            a_n.dropna(inplace=True)

        a_nn = a_n.select_dtypes(np.number)  # Numeric columns only.
        self.inf_rows_n = inf_rows_n = np.isinf(a_nn).any(axis=1) # = check_inf(a_n, '', self.out)
        if inf_rows_n.any():
            print('Drop %d Inf rows' % inf_rows_n.sum(), file=self.out)
            print('Data with Inf: columns %s:\nhead:\n' % (
                        a_nn.columns[np.isinf(a_nn).any()],), file=self.out)
            self.display(alls[~na_rows_n][inf_rows_n].head(10))
            a_n = a_n[~inf_rows_n].copy(deep=False)
        
        constant_vars = check_constants(a_n, '', self.out)
        if constant_vars:
            logger.warning('**DROPPING** constant features')
            a_n.drop(columns=constant_vars, inplace=True)
            self.features = features = [ft for ft in features if ft not in constant_vars]
        self.a_n = a_n

        print('\nResponse transforms and normalisations:\n'
                + '\n'.join('%s -> %s : Mean=%.4f Std=%.4f'
                            %(r, suggested_transform[r].__name__,
                              a_t[r].mean(), a_t[r].std())
                          for r in responses), file=self.out)

    #def do_fitall(self):
    def do_undecorred_fit(self):
        self.do_transform()
        a_n = self.a_n
        y = self.responses[0]
        features = self.features
        print('\n### Fit of all features (lasso)', file=self.out)

        #fitlasso, fit, pslopes, displayed =
        *_, pslopes0, preslopes0 = lassoselect_and_refit(a_n, y, features,
                                                    out=self.out, **self.regul_kw)
        self.displayed.extend((pslopes0, preslopes0))

        #sb.violinplot('null_dS_before', 'abs_age_dev', data=a_n, cut=0);
        #scatter_density('r2t_dS_mean', 'abs_brlen_dev', data=a_n, alpha=0.5)
        #scatter_density('ingroup_glob_len', 'abs_age_dev', data=a_n, alpha=0.5)
        #scatter_density('dS_rate_std', 'abs_age_dev', data=a_n, alpha=0.5)
        #scatter_density(alls.null_dist_before, alls.null_dS_before)

    #def do_undecorred_fa(self):
        print('\n### Dimension reduction of features\n#### PCA', file=self.out)

        ft_pca, ft_pca_outputs = detailed_pca(a_n, features, abs_cov=True,
                                              make_corr=True, heat_dendro=True,
                                              dendro_pad=0.2, out=self.out)
        self.ft_pca = ft_pca
        self.displayed.extend(ft_pca_outputs)

        print('\n#### Factor analysis', file=self.out)

        # To take into account continuous and categorical variables

        FA, fa_outputs = detailed_pca(a_n, features, FA=True, out=self.out)
        self.FA = FA
        self.displayed.extend(fa_outputs)

        # Outlier study (can be separated by a line in PC1-PC2 space)
        # (example from Catarrhini exact m1w04_fsa)
        #scatter_density(transformed0[:,0], transformed0[:,1], alpha=0.4)
        #ax = plt.gca()
        #ax.plot([-3, 3], [0.25, -0.5])

        #xsep = [-3, 3]; ysep = [0.25, -0.5]
        #eq_sep = lambda x: -0.125 - 0.75/6*x
        ##x = np.array([-2, -1, 0, 1, 2])
        ##plt.plot(x, eq_sep(x), 'r.')

        #transformed0_out = transformed0[:,1] > eq_sep(transformed0[:,0])
        #plt.plot(transformed0[transformed0_out, 0], transformed0[transformed0_out, 1], 'r.')

    def do_decorr(self, ref_rescale=None, rescale=True):
        print('\n### Feature decorrelation: 1. drop features; 2. decorr.', file=self.out)
        suggested_transform = self.suggested_transform
        must_drop_features = self.must_drop_features
        a_t = self.a_t
        a_n = self.a_n
        na_rows_n = self.na_rows_n
        inf_rows_n = self.inf_rows_n
        y = self.responses[0]
        ### **TODO**!!! Check that decorrelation is done on the un-zscored data!!
        
        #FIXME: there should be no "renorm" step here, but just a general zscore/standardize later.
        decorr_names = {renorm_decorrelatelogs: 'renorm_decorrelatelogs',
                        decorrelatelogs: 'decorrelatelogs',
                        renorm_decorrelate: 'renorm_decorrelate',
                        decorrelate: 'decorrelate',
                        renorm_unregress: 'renorm_unregress',
                        renorm_unregress0: 'renorm_unregress0'}

        # First drop features that were excluded manually (after seeing the first PCA).
        try:
            a_n_inde = a_n.drop(must_drop_features, axis=1)
        except KeyError:
            logger.warning("Not in `self.a_n`:" + ' '.join(must_drop_features.difference(a_n.columns)))
            a_n_inde = a_n.drop(must_drop_features, axis=1, errors='ignore')
        logger.info('Manual drop: %s', ' '.join(a_n.columns.intersection(must_drop_features)))
        print('%d Independent columns (%d rows)' % a_n_inde.shape[::-1], file=self.out)

        # Check for unfound variables that were selected for decorr.
        self.missing_todecorr = [pair for pair in
                     chain(*(decorr_args.values() for decorr_args in self.to_decorr.values()))
                     if pair[0] not in a_t or pair[1] not in a_t]
        self.dropped_todecorr = [pair for pair in
                     chain(*(decorr_args.values() for decorr_args in self.to_decorr.values()))
                     if pair[0] in must_drop_features]
        if self.missing_todecorr:
            logger.warning('Unexpected missing pairs to renorm(log)decorr: %s',
                           ' '.join('%s~%s' % pair for pair in self.missing_todecorr))
        if self.dropped_todecorr:
            logger.warning('Dropped before renorm(log)decorr: %s',
                           ' '.join('%s~%s' % pair for pair in self.dropped_todecorr))
        # Filtered
        valid_item = lambda item: (item[1][0] in a_t
                                            and item[1][1] in a_t
                                            and item[1][0] not in must_drop_features)

        to_decorr = {decorr_func: Args.fromitems(
                                       *filter(valid_item, decorr_args.items()))
                    for decorr_func, decorr_args in self.to_decorr.items()}

        # Check if the proposed decorrelation is in line with the transforms.
        # (subtract for logs, divide for raw/sqrt)
        to_divide_items = list(to_decorr.get(renorm_decorrelate, {}).items()) + \
                                list(to_decorr.get(decorrelate, {}).items())
        # reversed() because popping from list.
        for func in (renorm_decorrelatelogs, decorrelatelogs):
            for k, (var, corrvar) in reversed(list(to_decorr.get(func, {}).items())):
                msg = []
                var_transform = suggested_transform[var].__name__
                corrvar_transform = suggested_transform[corrvar].__name__
                if 'log' not in var_transform and 'binarize' not in var_transform:
                    msg.append('suggested_transform[%r] = %s' % (var, var_transform))
                if 'log' not in corrvar_transform and 'binarize' not in corrvar_transform:
                    msg.append('/ suggested_transform[%r] = %s' % (corrvar, corrvar_transform))
                if msg:
                    logger.warning(' '.join(msg) + ': automatic switch from "decorrelatelogs" (subtract) to "decorr" (divide).')
                    dest_func = decorrelate if decorr_names[func] == 'decorrelate' else renorm_decorrelate
                    try:
                        to_decorr[dest_func].additem(k, to_decorr[func].pop(k))
                    except KeyError:
                        to_decorr[dest_func] = Args()
                        to_decorr[dest_func].additem(k, to_decorr[func].pop(k))

        # Conversely, if dividing is on the log, ask confirmation.
        # But NO auto-switch, because sometimes dividing by the log-scaled var makes sense.
        for k, (var, corrvar) in to_divide_items:
            msg = []
            var_transform = suggested_transform[var].__name__
            corrvar_transform = suggested_transform[corrvar].__name__
            if 'log' in var_transform or 'binarize' in var_transform:
                msg.append('suggested_transform[%r] = %s' % (var, var_transform))
            if 'log' in corrvar_transform or 'binarize' in corrvar_transform:
                msg.append('/ suggested_transform[%r] = %s' % (corrvar, corrvar_transform))
            if msg:
                logger.warning(' '.join(msg) + ': are you sure to "decorr" (divide logs/binary)?')
                #to_decorr[renorm_decorrelatelogs].additem(k, to_decorr[renorm_decorrelate].pop(k))
                #if isinstance(k, int):
                #    to_renormdecorr[0].remove((var, corrvar))
                #    to_renormlogdecorr[0].append((var, corrvar))
                #else:
                #    to_renormdecorr[1].pop(k)
                #    to_renormlogdecorr[1][k] = (var, corrvar)

        if a_n_inde.shape[0] != a_t.shape[0]:
            msg = ('a_n_inde and a_t have incompatible shapes: %s %s '
                   '(NA rows: %d/%d, Inf rows: %d/%d)') % (
                         a_n_inde.shape, a_t.shape,
                         na_rows_n.sum(), na_rows_n.shape[0],
                         inf_rows_n.sum(), inf_rows_n.shape[0])
            if a_t.shape[0] - a_n_inde.shape[0] != (na_rows_n | inf_rows_n).sum():
                logger.error(msg)
            else:
                logger.warning(msg)

        #a_n_inde = renorm_decorrelatelogs(a_n_inde, a_t[~na_rows_n],
        #                                 *to_renormlogdecorr[0],
        #                                 **to_renormlogdecorr[1])
        #a_n_inde = renorm_decorrelate(a_n_inde, a_t[~na_rows_n],
        #                              *to_renormdecorr[0],
        #                              **to_renormdecorr[1])

        #self.missing_zeros_todecorr = [pair for pair in self.to_subtract#[0]
        #                          if pair[0] not in a_t or pair[1] not in a_t]
        #if self.missing_zeros_todecorr:
        #    logger.warning('Unexpected missing pairs of zero measures to decorr: %s',
        #                   self.missing_zeros_todecorr)
        #zeros_to_decorrelate = [pair for pair in self.to_subtract
        #                        if pair[0] in a_t and pair[1] in a_t]

        #a_n_inde = decorrelatelogs(a_n_inde, a_t[~na_rows_n], *zeros_to_decorrelate)

        self.unregressed_coefs = {}  # newfeature: (intercept, slope)
        for decorr_func, decorr_args in to_decorr.items():
            #FIXME: UserWarning: Boolean Series key will be reindexed to match DataFrame index.
            a_n_inde = decorr_args(decorr_func, a_n_inde, a_t[~na_rows_n & ~inf_rows_n])
            # decorring from a_t or alls might yield very different results.
            if decorr_func in (check_unregress, renorm_unregress):
                for k, (var, corrvar) in decorr_args.items():
                    fit = sm.OLS(a_t[~na_rows_n & ~inf_rows_n][var],
                        sm.add_constant(a_t[~na_rows_n & ~inf_rows_n][corrvar])).fit()
                    if isinstance(k, int): k = 'R'+var
                    self.unregressed_coefs[k] = (fit.params['const'], fit.params[corrvar])
                    if not np.allclose(a_n_inde[k],
                            a_t.loc[~na_rows_n & ~inf_rows_n, var] - (
                                fit.params['const']
                                + fit.params[corrvar]*a_t.loc[~na_rows_n & ~inf_rows_n, corrvar])):
                        logger.error("INCONSISTENCY in unregress residuals: %s~%s", var, corrvar)
            #elif isinstance(decorr_func, partial) and decorr_func.args[0].__name__ == 'refdecorr':

        # FIXME: special hardcoded decorr
        if set(('ingroup_mean_CpG', 'ingroup_mean_GC')).intersection(self.alls.columns):
            CpG_odds = (self.alls.ingroup_mean_CpG / (self.alls.ingroup_mean_GC**2))[~na_rows_n & ~inf_rows_n]
            CpG_odd_transform = make_best_logtransform(CpG_odds)
            a_n_inde['CpG_odds'] = CpG_odd_transform(CpG_odds)

        print('%d independent columns (%d rows)' % a_n_inde.shape[::-1], file=self.out)

        #TODO: fix decorrs
        #self.fix_decorr()

        to_decorr_iter_args = [(decorr_func, *item)
                               for decorr_func, decorr_args in to_decorr.items()
                               for item in decorr_args.items()]
        if set(('ingroup_mean_CpG', 'ingroup_mean_GC')).intersection(self.alls.columns):
            to_decorr_iter_args.append(('CpG/GC^2', 'CpG_odds', ('ingroup_mean_CpG', 'ingroup_mean_GC')))
        if self.widget:
            to_decorr_iter_args = self.widget(to_decorr_iter_args)
        for decorr_func, *decorr_item in to_decorr_iter_args:
            fig = display_decorrelate(decorr_item, self.alls, a_t, a_n_inde, color_var=y)
            fig.suptitle('%s: %s ~ %s' %(decorr_names.get(decorr_func, decorr_func), *decorr_item[1]))
            self.show(fig); plt.close(fig)

        # Description of the transformation for the user.
        #self.decorr = {'R'+var1: ('-', var2) for var1, var2 
        #                in to_renormlogdecorr[0] + zeros_to_decorrelate}
        #self.decorr.update(
        #        CpG_odds=('ingroup_mean_CpG', '/', 'ingroup_mean_GC^2'),
        #        **{'R'+var1: ('/', var2) for var1, var2 in to_renormdecorr[0]},
        #        **{key: (v1, '/', v2) for key, (v1,v2)
        #              in chain(to_renormlogdecorr[1].items(),
        #                       to_renormdecorr[1].items())})
        decorr_symbol = {renorm_decorrelatelogs: '-',
                         decorrelatelogs: '-',
                         renorm_decorrelate: '/',
                         decorrelate: '/',
                         renorm_unregress: '~',
                         check_unregress: '~'}
        self.decorred = dict(('R'+var1, (decorr_symbol.get(decorr_func, str(decorr_func)), var2))
                                if isinstance(k, int) else
                                (k, (var1, decorr_symbol.get(decorr_func, str(decorr_func)), var2))
                             for decorr_func, decorr_args in to_decorr.items()
                             for k, (var1, var2) in decorr_args.items())
        if set(('ingroup_mean_CpG', 'ingroup_mean_GC')).intersection(self.alls.columns):
            self.decorred.update(
                CpG_odds=('ingroup_mean_CpG', '/', 'ingroup_mean_GC^2'))
        # Discriminate binary features from continuous ones, to choose "standardize" instead of "zscore".
        decorred_continuous = []
        decorred_binary = []
        for ft, decorred_item in self.decorred.items():
            corrvar = decorred_item[-1]
            var = ft[1:] if len(decorred_item)==2 else decorred_item[0]
            if 'binary' in self.suggested_transform[var].__name__ and 'binary' in self.suggested_transform[corrvar].__name__:
                decorred_binary.append(ft)
            else:
                decorred_continuous.append(ft)

        self.decorr_source = {key: (key[1:] if len(tup)==2 else tup[0])
                              for key, tup in self.decorred.items()}

        if self.decorred:
            self.a_decorred = a_n_inde[list(self.decorred)].copy()  # Used in self.predict()
            self.decorred_stats = a_n_inde[list(self.decorred)].agg(['mean', 'std'])
        else:
            self.a_decorred, self.decorred_stats = {}, {}

        if rescale:
            self.a_n_inde = a_n_inde.assign(
                    **dict(zscore_dataframe(a_n_inde[decorred_continuous], ref_rescale)),
                    **dict(standardize_dataframe(a_n_inde[decorred_binary], ref_rescale))
                    )
        else:
            self.a_n_inde = a_n_inde

    # Could be: @property(inde_features)
    def do_decorred_fa(self):
        a_n_inde = self.a_n_inde  # Could make this a property, triggering do_decorr.
        inde_features = [ft for ft in self.features if ft in a_n_inde]
        #print('inde_features', len(inde_features), file=self.out)
        inde_features += [colname for colname in a_n_inde.columns.difference(self.a_n.columns)]
        print('inde_features', len(inde_features), file=self.out)

        new_inde_features = set(inde_features) - set(self.features)
        self.inde_features = inde_features

        na_rows_decorr = check_nan(a_n_inde, 'after decorr', self.out)
        check_inf(a_n_inde, 'after decorr', self.out)
        #if na_rows_decorr.any():
        #    print('%d NA rows after decorr!' % na_rows_decorr.sum(), file=self.out)

        #ft_pca_inde, ft_pca_inde_outputs = detailed_pca(a_n_inde, inde_features, out=self.out)
        self.ft_pca_inde = ft_pca_inde = PCA(n_components=15)
        # The above line is necessary to get all attributes (e.g. covariance)
        ft_pca_inde.fit_transform(a_n_inde[inde_features]) # -> transformed data

        heatmap_cov(np.abs(ft_pca_inde.get_covariance()), inde_features,
                    make_corr=True).suptitle('Inde features absolute correlation (PCA)')
        self.displayed.append(plt.gcf())
        self.show(); plt.close()

        #FA_inde, FA_inde_outputs = detailed_pca(a_n_inde, inde_features, FA=True, out=self.out)
        self.FA_inde = FA_inde = FactorAnalysis(n_components=15)
        self.FA_inde_components = FA_inde_components = FA_inde.fit_transform(a_n_inde[inde_features])

        heatmap_cov(np.abs(FA_inde.get_covariance()), inde_features, make_corr=True)
        self.displayed.append(plt.gcf())
        self.displayed[-1].suptitle('Inde features absolute correlation (FA)')
        self.show(); plt.close()
        fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
        scatter_density(FA_inde_components[:,1], FA_inde_components[:,0], alpha=0.4,
                        ax=ax1)
        ax1.set_ylabel('PC1')
        ax1.set_xlabel('PC2')
        scatter_density(FA_inde_components[:,2], FA_inde_components[:,0], alpha=0.4,
                        ax=ax2)
        ax2.set_xlabel('PC3')
        self.displayed.append(fig)
        self.show(); plt.close()


    def do_decorred_fit(self):
        logger = self.logger
        self.do_decorred_fa()
        y = self.responses[0]

        print('\n### Fit of less colinear features (after internal decorrelation)', file=self.out)

        self.fitlasso, self.fit, pslopes, preslopes = lassoselect_and_refit(
                                self.a_n_inde, y, self.inde_features, atol=1e-2,
                                #method='elastic_net', L1_wt=1,
                                cov_type='HC1', out=self.out,
                                **self.regul_kw) #**psummary_kw
        self.slopes = pslopes.data
        self.reslopes = preslopes.data
        self.displayed.extend((pslopes, preslopes))

    #def do_randomforest(self):
        print('\n#### Random Forest Regression', file=self.out)
        try:
            RFcrossval_r2 = randomforest_regression(self.a_n_inde[self.inde_features],
                                                    self.a_n_inde[y],
                                                    out=self.out)
        except TerminatedWorkerError:
            logger.error("Can't compute RFcrossval_r2 from random forest due to parallel error (TerminatedWorkerError)")


    @property
    def suggest_multicolin(self):
        try:
            return self._suggest_multicolin
        except AttributeError:
            print('\n#### Finding the most colinear features', file=self.out)
            #multicol_test(a_n[features]), multicol_test(a_n_inde[inde_features])
            protected_features = [ft for ft in self.protected_features
                                    if ft in self.inde_features]
            #assert np.isclose(self.refitlasso.condition_number,
            #                  multicol_test(self.a_n_inde[self.inde_features]))
            self._suggest_multicolin, outputs = loop_leave1out(self.inde_features,
                                          lambda x: multicol_test(self.a_n_inde[x]),
                                          criterion='min',
                                          stop_criterion='<20',
                                          #start=self.refitlasso.condition_number,
                                          protected=protected_features,
                                          out=self.out, widget=self.widget)
            self.displayed.extend(outputs)
            return self._suggest_multicolin

    @property
    def suggest_multicolin2(self):
        try:
            return self._suggest_multicolin2
        except AttributeError:
            print('\n#### Finding the most colinear features (after bad-data removed)', file=self.out)
            #multicol_test(a_n[features]), multicol_test(a_n_inde[inde_features])
            protected_features2 = self.protected_features.intersection(self.inde_features2)
            #assert np.isclose(self.refitlasso.condition_number,
            #                  multicol_test(self.a_n_inde[self.inde_features]))
            self._suggest_multicolin2, outputs = loop_leave1out(self.inde_features2,
                                          lambda x: multicol_test(self.a_n_inde2[x]),
                                          criterion='min',
                                          stop_criterion='<20',
                                          protected=protected_features2,
                                          out=self.out, widget=self.widget)
            self.displayed.extend(outputs)
            return self._suggest_multicolin2


    def do_drop_bad(self, ref_rescale=None, rescale=True):
        """
        Drop colinear features (auto selection) and bad properties:

        - colinear features in columns;
        - bad properties in rows.
        """
        logger = self.logger
        a_n_inde = self.a_n_inde
        inde_features = self.inde_features
        y = self.responses[0]
        assert not a_n_inde[inde_features + [y]].isna().any(axis=None)

        print(self.suggest_multicolin, file=self.out)

        print('\n#### Refit without most colinear features', file=self.out)

        print('Drop trees with bad properties (before removing colinear features)', file=self.out)
        bad_props = {}
        missing_bad_props = set()
        bad_props_ttests = {}
        self.bad_props_counts = {}
        for badp, badval in self.must_drop_data.items():
            if badp not in a_n_inde.columns:
                missing_bad_props.add(badp)
                continue
            print('\n--------\n#', badp, file=self.out)
            vcounts = a_n_inde[badp].value_counts(dropna=False)
            if vcounts.shape[0] > 2:
                logger.warning('%r not binary.', badp)
            print('\n'.join(str(vcounts).split('\n')[:-1]), file=self.out)
            #is_badval = a_n_inde[badp].isin(badval)
            is_badval = self.a_t[badp].reindex(a_n_inde.index).isin(badval)
            if (is_badval.sum()) > 0.05 * vcounts.sum():
                logger.warning('Discarding %s.isin(%s) trees will remove >5%% of the subtrees',
                               badp, badval)
            if (~is_badval).sum() == 0 or is_badval.sum() == 0:
                logger.error('0 data available in one group [badval = %s]', badval)
            tt_badp = stats.ttest_ind(
                        a_n_inde.loc[is_badval, y],
                        a_n_inde.loc[~is_badval, y],
                        equal_var=False)
            print('T-test: Y(badval) VS Y(goodval): t=%g, p=%g' % tt_badp, file=self.out)
            bad_props[badp] = badval
            bad_props_ttests[badp] = tt_badp
            self.bad_props_counts[badp] = is_badval.sum()

        if missing_bad_props:
            logger.warning('Hard-coded Bad properties not in the columns: %s',
                           ' '.join(missing_bad_props))

        self.bad_data_rows = bad_data_rows = self.a_t.reindex(a_n_inde.index).isin(bad_props).any(axis=1)
        self.bad_props_counts['ANY'] = bad_data_rows.sum()
        a_n_inde2 = a_n_inde.loc[~bad_data_rows]
        tt_bad = stats.ttest_ind(
                    a_n_inde.loc[bad_data_rows, y],
                    a_n_inde2[y],
                    equal_var=False)
        print(("Removed VS kept data: t=%g; p=%g (two-tailed Welch's t-test; "
               +"sample sizes: %d, %d)") % (*tt_bad,
                  bad_data_rows.sum(), (~bad_data_rows).sum()), file=self.out)
        bad_props_ttests = pd.DataFrame(bad_props_ttests, index=['T', 'P']).T
        na_ttests = bad_props_ttests['P'].isna()

        try:
            bad_props_ttests['P_over'] = bad_props_ttests.apply(lambda row:
                                        (1 - row['P']/2 if row['T']<=0 else row['P']/2),
                                        axis=1)
            if bad_props_ttests.shape[0]:
                bad_props_ttests['Pcorr'] = smm.multipletests(
                                                bad_props_ttests['P_over'],
                                                alpha=0.05,
                                                method='fdr_bh')[1] # Benjamini/Hochberg
                # Even with *some* NaN p-values, the added column should not be all NaN:
                if bad_props_ttests.Pcorr.isna().all():
                    logger.error('Pcorr are all NaNs. There were %d input NaN pvalues (%s)',
                                 na_ttests.sum(), bad_props_ttests.index[na_ttests])
            else:
                bad_props_ttests['Pcorr'] = np.NaN
            # Add last row being the pooled features T-test.
            self.bad_props_ttests = pd.concat((bad_props_ttests,
                                        pd.DataFrame([tt_bad], columns=['T', 'P'], index=['ANY']))
                                        )
            self.bad_props_ttests['removed'] = pd.Series(self.bad_props_counts)
            self.display_html(self.bad_props_ttests.style.applymap(
                             lambda v: ('background: khaki' if v<0.01 else
                                        'background: lemonchiffon; alpha: 0.5' if v<=0.05 else ''),
                                     subset=['Pcorr']))
        except ValueError as err:
            # Cannot set a frame with no defined index and a value that cannot be converted to a Series
            logger.exception('Possibly no features removed.')

        print("New shape after removing bad data:", a_n_inde2.shape, file=self.out)
        print("New Mean,SD of response: %g; %g" % (
                a_n_inde2[y].mean(), a_n_inde2[y].std()), file=self.out)
        good_features = [ft for ft in inde_features if ft not in set(bad_props)]

        var_ranges = a_n_inde2[good_features].max(axis=0) - a_n_inde2[good_features].min(axis=0)

        #check_constants(a_n_inde2[inde_features2])
        self.inde_features2, self.features_postnorm_constant = [], []
        for ft in good_features:
            if var_ranges[ft]>0:
                self.inde_features2.append(ft)
            else:
                self.features_postnorm_constant.append(ft)
        print('%d/%d newly introduced constant variables:\n' % (len(self.features_postnorm_constant), len(good_features)),
              var_ranges.sort_values(), file=self.out)

        # Removing bad data probably shifted Y. We need to re-center and standardize.
        self.continuous_inde_features2 = []
        self.bin_inde_features2 = []
        self.bin_responses2 = []
        for ft in self.inde_features2:
            if self.suggested_transform[self.decorr_source.get(ft, ft)].__name__.startswith('binarize'):
                self.bin_inde_features2.append(ft)
            else:
                self.continuous_inde_features2.append(ft)
        for ft in self.responses:
            if self.suggested_transform[self.decorr_source.get(ft, ft)].__name__.startswith('binarize'):
                self.bin_responses2.append(ft)
            else:
                self.continuous_inde_features2.append(ft)

        if rescale:
            self.a_n_inde2 = a_n_inde2[self.bin_responses2].assign(
                **dict(
                    zscore_dataframe(a_n_inde2[self.continuous_inde_features2], ref_rescale)  #.apply(zscore, raw=True, result_type='broadcast')
                    ),
                **dict(
                    standardize_dataframe(a_n_inde2[self.bin_inde_features2], ref_rescale)
                    )
                )
            #TODO: drop some constant features.
        else:
            self.a_n_inde2 = a_n_inde2

        print(self.suggest_multicolin2, file=self.out)

    def do_bestfit(self):
        self.do_drop_bad()  # Generates a_n_inde2, inde_features2 and suggest_multicolin2
        a_n_inde2 = self.a_n_inde2
        self.inde_features2 = inde_features2 = [ft for ft in self.inde_features2
                                                if ft not in self.suggest_multicolin2]
        y = self.responses[0]

        #print('\n##### OLS fit (2)', file=self.out)
        #self.fit2 = fit2 = ols.fit()
        #fit2.summary()
        ##display_html(sm_pretty_slopes(fit2), file=self.out)
        #display_html(fit2.summary(), file=self.out)

        print('\n##### LASSO fit (2)', file=self.out)
        inde_features2_source = [self.decorr_source.get(ft, ft) for ft in inde_features2]
        param_info = pd.concat((
                        pd.Series([self.suggested_transform[sft].__name__
                                   for sft in inde_features2_source],
                                  index=inde_features2,
                                  name='transform'),
                        self.a_t.reindex(a_n_inde2.index)[inde_features2_source]\
                                .agg(['mean', 'std'])\
                                .rename(lambda c: 'transformed '+c).T\
                                .set_axis(inde_features2, axis=0),
                        pd.Series({ft: ' '.join(self.decorred.get(ft, ()))
                                   for ft in inde_features2}, name='decorr')
                        ),
                        axis=1, sort=False)

        print('Re-check newly introduced constant features:', file=self.out)
        print(a_n_inde2[inde_features2].agg(np.ptp).sort_values().head(), file=self.out)

        logger.debug('regul_kw = %s', self.regul_kw)
        fitlasso2, refitlasso2, slopes2_styled, reslopes2_styled = \
                lassoselect_and_refit(a_n_inde2, y, inde_features2, join=param_info, out=self.out,
                                      **self.regul_kw)
        self.displayed.extend((slopes2_styled, reslopes2_styled))
        self.fitlasso2 = fitlasso2
        self.refitlasso2 = refitlasso2
        self.slopes2_styled = slopes2_styled
        self.reslopes2_styled = reslopes2_styled

        self.slopes2 = slopes2_styled.data
        reslopes2 = self.reslopes2 = self.reslopes2_styled.data
        self.selected_features2 = reslopes2.drop('const').index.tolist()
        self.F_pval2 = getattr(refitlasso2, 'f_pvalue', None)  # no such attribute for LogitResults
        self.lL2 = refitlasso2.llf
        #print('fstat = %g\nP(F > fstat) = %g\nlog-likelihood = %g' %(
        #        refitlasso2.fvalue, self.F_pval2, self.lL2), file=self.out)

        # Plot Y~X with each coef (the 15 largest).
        if self.suggested_transform[y].__name__.startswith('binarize') or \
                np.issubdtype(a_n_inde2[y].dtype, bool):
            def responseplot(x, y, ax=None):
                ax = sb.violinplot(a_n_inde2, x=x, y=y, inner='stick', orient='h', ax=ax, color='.7')
                ax.autoscale(False, 'y')  # Stay focused on the data limits, do not adjust Y axis because of regression lines.
        else:
            def responseplot(x, y, ax=None):
                scatter_density(a_n_inde2[x], a_n_inde2[y], alpha=0.4, ax=ax)
                ax = ax if ax is not None else plt.gca()
                ymin, ymax = ax.get_ylim()
                ax.autoscale(False, 'y')
                ax.set_ylim(ymin, ymax)

        fit_figs = []
        iter_coeffs = range(1, min(16, len(self.reslopes2.index)))
        if self.widget is not None:
            iter_coeffs = self.widget(iter_coeffs)
        for coef_i in iter_coeffs:
            fig, ax = plt.subplots(num=777, clear=True) # If plt.show go crazy, this reuses the previous fig
            responseplot(self.slopes2.index[coef_i], y, ax=ax)
            #ax.relim(); ax.autoscale_view()
            x = np.array(ax.get_xlim())
            a, b = reslopes2.loc[['const', self.slopes2.index[coef_i]], 'Simple regression coef']
            ax.plot(x, 0 + b*x, '--', label='Simple Regression line (a=%g b=%g)' % (0,b))
            a, b = reslopes2.loc[['const', self.slopes2.index[coef_i]], 'Lasso coef']
            ax.plot(x, a + b*x, '--', label='Lasso regression line (a=%g b=%g)' % (a,b))
            ax.legend()
            ax.set_title('Response ~ %s (%s largest Lasso coef).'
                         % (self.slopes2.index[coef_i],
                            '1st' if coef_i==1 else
                            '2nd' if coef_i==2 else
                            '3rd' if coef_i==3 else '%dth' % coef_i))
            ax.set_ylabel(y)
            ax.set_xlabel(self.slopes2.index[coef_i])
            fit_figs.append(fig)
            self.show(fig); plt.close()
            #ax = sb.regplot(a_n_inde2[self.slopes2.index[coef_i]], a_n_inde2[y])
            #self.show(); plt.close()
        self.displayed.append(fit_figs)

        # Residual plot
        fig, axes = plt.subplots(2)
        ax = axes[0]
        scatter_density(a_n_inde2[y] - refitlasso2.fittedvalues,
                        refitlasso2.fittedvalues,
                        alpha=0.4, ax=ax)
        ax.set_title('Residuals plot')
        ax.set_ylabel('Residual error')
        ax.set_xlabel('Predicted response')

        try:
            residuals = refitlasso2.resid  # OLS fit
        except AttributeError:
            residuals = refitlasso2.resid_response  # Logit fit
        smg.gofplots.qqplot(residuals, line='r', ax=axes[1])
        self.displayed.append(fig)
        self.show(fig); plt.close()

        if hasattr(refitlasso2, 'cov_HC0'):
            fig = heatmap_cov(refitlasso2.cov_HC0, ['const']+self.selected_features2,
                        cmap='seismic', make_corr=True)
            fig.suptitle('cov_HC0')
            n = refitlasso2.cov_HC0.shape[0]
            w, h = fig.get_size_inches()
            fig.set_size_inches(min(w, w/20. * n), min(h, h/20 * n))
            self.displayed.append(fig)
            self.show(fig); plt.close()
        else:
            logger.warning("No attribute 'cov_HC0' in `refitlasso2`.")
        if hasattr(refitlasso2, 'cov_HC1'):
            fig = heatmap_cov(refitlasso2.cov_HC1, ['const']+self.selected_features2,
                        cmap='seismic', make_corr=True)
            n = refitlasso2.cov_HC1.shape[0]
            w, h = fig.get_size_inches()
            fig.set_size_inches(min(w, w/20. * n), min(h, h/20. * n))
            fig.suptitle('cov_HC1')
            self.displayed.append(fig)
            self.show(fig); plt.close()  # Not needed because this is the last figure before return.
        else:
            logger.warning("No attribute 'cov_HC1' in `refitlasso2`.")

        # See this page for all possible accessible attributes/methods:
        # https://www.statsmodels.org/stable/generated/statsmodels.regression.linear_model.OLSResults.html#statsmodels.regression.linear_model.OLSResults

        return reslopes2


    def reshow(self):
        for out in self.outputs:
            if isinstance(out, (pd.Series, pd.DataFrame, pd.io.formats.style.Styler)):
                self.display_html(out)
            elif isinstance(out, mpl.figure.Figure):
                self.display(out)
            else:
                print(out, file=self.out)


    def do_extreme_points(self):
        print('\n### Investigate the data points with highest %s' % self.responses[0],
              file=self.out)
        self.display_html(self.alls.sort_values(self.responses[0], ascending=False).head(50))
        print('\n### Investigate the data points with lowest %s' % self.responses[0],
              file=self.out)
        self.display_html(self.alls.sort_values(self.responses[0], ascending=True).head(50))


    def plot_coefs(self, fig_columns=['coef'], renames=None, signif_levels=[0.05, 0.01]):
        # fig_columns could also include: ['Lasso coef', 'Simple regression coef']
        if renames is None: renames = {}
        # Adjust significance levels if not given
        if signif_levels is None:
            raise NotImplementedError('Automatic choice of significance levels')
            signif_levels = [0.05]
            signif_pvals = pvalues[pvalues>0.05]
            signif_levels.append(10**round(np.log10(median(signif_pvals))))
        signif_levels = - np.asarray(signif_levels)  # negate pvalue cutoffs to sort by increasing significance
        signif_levels.sort()

        ax_titles = {'coef': 'Multiple regression coefficient',
                     'Lasso coef': 'Lasso regression',
                     'Simple regression coef': 'Simple regression'}
        ncols = len(fig_columns)
        coefs = self.reslopes2  # Should run .do() or .do_bestfit() first. DataFrame expected.

        coef_err = coefs[['[0.025', '0.975]']].subtract(coefs.coef, axis=0).abs().values.T
        nan_err = np.full(coef_err.shape, np.NaN)
        pvalues = coefs['adj. p-value'].values
        # Specify the colors explicitly using the dict form, otherwise only the first color of the list is used.
        barcolors = {colname: '#d65f5f' for colname in fig_columns}
        barcolors['coef'] = ['#d65f5fff' if p<=-signif_levels[0] else '#d65f5f88' for p in pvalues]

        axes = matplotlib_stylebar(coefs.rename(renames), fig_columns, barcolors,
                                   err=np.array([coef_err]+[nan_err]*(ncols-1)),
                                   float_fmt='% .4f')
        #ax0_xmax = axes[0].texts[0].get_window_extent().x1  # In figure pixels? Requires a renderer.
        ax0_xmax = coefs['0.975]'].max()
        for j, (ft, pval) in enumerate(coefs['adj. p-value'].items()):
            starring = '*' * np.searchsorted(signif_levels, - pval, side='right')
            txt = axes[0].annotate('%-4s' % starring, (ax0_xmax, j), va='center', ha='left',
                                   textcoords='offset points', xytext=(45,2))  # Ad hoc xytext, must go left of the coef texts.

        axes[0].set_xlabel(None)
        #axes[0].set_ylabel(ylabel)
        axes[0].set_title(ax_titles[fig_columns[0]])
        axes[0].set_xticks([-0.5, 0, 0.5])
        axes[0].tick_params(bottom=True, labelbottom=True)
        axes[0].grid(axis='x', color='.2', alpha=0.5)
        axes[0].set_xticklabels(['-0.5', '0', '0.5'])
        axes[0].annotate('\n'.join('%-4s: adj. p-value ≤ %-4g' % (('*' * i), -lvl)
                                   for i,lvl in enumerate(signif_levels, start=1)),
                        xy=(axes[0].get_xlim()[1]/2., axes[0].get_ylim()[0]), #xycoords='axes fraction',
                        xytext=(0,-18), textcoords='offset points',
                        va='top', ha='left', fontsize='xx-small')

        if len(axes)>1:
            axes[1].set_xlabel(None)
            axes[1].set_title(ax_titles.get(fig_columns[1], fig_columns[1]))
        if len(axes)>2:
            axes[2].set_xlabel(None)
            axes[2].set_title(ax_titles.get(fig_columns[2], fig_columns[2]))
        fig = axes[0].figure
        fig.tight_layout()
        fig.set_size_inches(3+3*len(fig_columns), 0.28 * (coefs.shape[0]+2))
        return fig, axes
        #for ext in ('pdf', 'png'):
        #    plt.savefig('../fig/coefs_fsahmmc+mpl.%s' % ext, bbox_inches='tight')

    def predict_extreme_values(self, percent=10):
        print('\n### Predicted extreme values', file=self.out)
        y = self.responses[0]
        N = self.a_n_inde2.shape[0]
        #print(N)
        sorted_fit = self.refitlasso2.fittedvalues.sort_values().rename('fitted_y')
        n = N * percent // 100
        fitted_lowest = sorted_fit.head(n).index  # Lowest error
        fitted_highest = sorted_fit.tail(n).index
        top_features = self.slopes2.index[1:11].tolist()
        for decorred_ft in set(top_features).intersection(self.decorr_source):
            print('Decorred %r: %s' % (decorred_ft, ' '.join(self.decorred[decorred_ft])),
                  file=self.out)

        top_features = [self.decorr_source.get(ft, ft) for ft in top_features]
        try:
            response_transform = self.suggested_transform[y].__name__
        except KeyError:
            response_transform = 'notransform'

        retro_trans = make_retro_trans(response_transform)

        #tr_mean, tr_std = self.reslopes2.loc[y,
        #                            ['transformed mean', 'transformed std']]
        #FIXME: don't remember if it's after dropping bad properties.
        tr_mean, tr_std = (self.a_t.loc[self.a_n_inde2.index, y].mean(),
                           self.a_t.loc[self.a_n_inde2.index, y].std())

        self.fitted_lowest_max = sorted_fit[fitted_lowest[-1]]
        self.fitted_highest_min = sorted_fit[fitted_highest[0]]
        self.fitted_lowest_max_tr = self.fitted_lowest_max*tr_std + tr_mean
        self.fitted_highest_min_tr = self.fitted_highest_min*tr_std + tr_mean
        self.fitted_lowest_max_orig = retro_trans(self.fitted_lowest_max*tr_std + tr_mean)
        self.fitted_highest_min_orig = retro_trans(self.fitted_highest_min*tr_std + tr_mean)
        print(('\n#### %g%% lowest *fitted* %s: <= %g\n'
                 '     %g%% Highest >= %g\n'
                 '#### In transformed scale: <= %g / >= %g\n'
                 '#### In original scale:     %g / %g [undoing %r]\n') % (
              percent, y, self.fitted_lowest_max,
              percent, self.fitted_highest_min,
              self.fitted_lowest_max_tr,
              self.fitted_highest_min_tr,
              self.fitted_lowest_max_orig,
              self.fitted_highest_min_orig,
              response_transform),
              file=self.out)

        print('\n#### Original data: (mean ± std)', file=self.out)
        summary_lowest = self.alls.loc[fitted_lowest,
                                    [y] + top_features
                                  ].agg(['median', 'mean', 'std', ci95])
        summary_highest = self.alls.loc[fitted_highest,
                                    [y] + top_features
                                  ].agg(['median', 'mean', 'std', ci95])
        self.predicted_extremes = pd.concat((summary_lowest, summary_highest),
                               keys=['Predicted %s%% lowest' % percent,
                                     'Predicted %s%% highest' % percent],
                               verify_integrity=True, sort=False, copy=False)
        # Format only the {mean} ± {std}
        self.predicted_extremes_txt = (
                self.predicted_extremes.xs('mean', level=1).applymap('{:.3g}'.format) +
                self.predicted_extremes.xs('std', level=1).applymap(' ± {:.3g}'.format))
        self.display(self.predicted_extremes_txt)

        iter_features = top_features if self.widget is None else self.widget(top_features)
        for ft in iter_features:
            colored = pd.Series('#7f7f7f', index=self.a_n_inde2.index, name='color')
            colored[fitted_highest] = '#991009'
            colored[fitted_lowest] = '#202099'
            fig, ax = plt.subplots()
            ax.scatter(ft, y, color='color', data=self.alls.join(colored, how='right'),
                       alpha=0.4, linewidths=0, edgecolors='none')
            ax.set_ylabel(y)
            ax.set_xlabel(ft)
            self.show(fig)
            negative_corr = (self.predicted_extremes[ft].loc[('Predicted %s%% lowest' % percent,
                                                             'mean')]
                            > self.predicted_extremes[ft].loc[('Predicted %s%% highest' % percent,
                                                               'mean')])
            sorted_ft_propbad = self.alls.join(sorted_fit)\
                                   .sort_values(ft, ascending=negative_corr)\
                                   .assign(frac_invalid=lambda df: (df.fitted_y>=self.fitted_highest_min).expanding().mean())
            # First index that falls under 50%
            ft_index_halfbad = np.searchsorted(sorted_ft_propbad.frac_invalid < 0.5, True)-1
            ft_threshold_halfbad = sorted_ft_propbad[ft].iloc[ft_index_halfbad]
            self.display_html(sorted_ft_propbad[[ft, 'fitted_y', 'frac_invalid']].head(20))
            print('Limit of %s such that 50%% of tested trees have a predicted'
                  ' %s >= %g: %g (%d trees, negative relation=%s)'
                  % (ft, y, self.fitted_highest_min,
                     ft_threshold_halfbad, ft_index_halfbad, negative_corr),
                  file=self.out)
            plt.close()


    def predict(self, out=None, logger=None, widget=None, **new_attrs):
        """Return predicted Y in transformed scale (unnormalised).
        Also attach the prediction object (of same class) to self."""
        selected_features = self.reslopes2.index[1:].tolist()
        ### Apply all the transforms + decorrelations
        self.prediction = prediction = self.__class__.from_other(self)
        prediction.out = out
        prediction.logger = logger
        prediction.widget = widget
        prediction.set_output()

        training_observations = self.a_n_inde.index
        training_observations2 = self.a_n_inde2.index
        for attrname, value in new_attrs.items():
            setattr(prediction, attrname, value)
        prediction.do_transform(rescale=False) #ref_rescale=self.alls.loc[training_observations])
        prediction.do_decorr(rescale=False) #ref_rescale=self.a_decorred)

        ref_rescale = self.a_t.assign(**dict(self.a_decorred)).loc[training_observations]
        prediction.do_drop_bad(rescale=False)  # a_t_inde
        a_n_inde2 = prediction.a_n_inde2
        if not training_observations2.difference(a_n_inde2.index).empty:
            logger.warning('Lost %d observations in the prediction.', len(training_observations2.difference(a_n_inde2.index)))
        if a_n_inde2.index.difference(training_observations2).empty:
            logger.error('Testing set is empty.')
        prediction.a_n_inde2 = a_n_inde2.assign(
                                **dict(zscore_dataframe(a_n_inde2[self.continuous_inde_features2],
                                                        a_n_inde2.loc[training_observations2])),
                                **dict(standardize_dataframe(a_n_inde2[self.bin_inde_features2],
                                                        a_n_inde2.loc[training_observations2]))
                                )

        #pred_reslopes2 = prediction.do_bestfit()
        common_features = set(self.inde_features2).intersection(prediction.inde_features2)
        logger.info('%d common features.', len(common_features))
        discarded_features = set(self.inde_features2).difference(prediction.inde_features2)
        if discarded_features:
            logger.warning('Prediction discarded %d features: %s' % (
                    len(discarded_features), '\n'.join(sorted(discarded_features))))
        if not a_n_inde2.columns.difference(common_features).empty:
            logger.warning('Extraneous features in prediction (NOT rescaled): %s',
                           ' '.join(a_n_inde2.columns.difference(common_features)))

        prediction_features = self.refitlasso2.params.drop('const').index.tolist()
        stds = a_n_inde2.std(ddof=1)
        not_scaled = (stds - 1).abs() > 0.0001
        if not_scaled.any():
            if not_scaled[prediction_features].any():
                logger.error('Not rescaled but used for prediction: %s', 
                    ' '.join(not_scaled[prediction_features][not_scaled[prediction_features]].index))
            logger.warning('Other not rescaled: %s',
                ' '.join(not_scaled[not_scaled].drop(prediction_features, errors='ignore').index))

        logger.info('training data shape: %s', self.refitlasso2.model.exog.shape)
        logger.info('training data used %d features for refitting.', len(prediction_features))

        predicted = self.refitlasso2.predict(
                        sm.add_constant(
                            prediction.a_n_inde2[prediction_features]))\
                    .rename('predicted')
        prediction.predicted_raw = predicted
        # Transform back to the original scale.
        y = self.responses[0]
        err_transform = self.suggested_transform[y]
        err_t = err_transform(prediction.alls[y])[training_observations2]
        self.err_t_center, self.err_t_scale = err_t.mean(), err_t.std(ddof=1)
        print('%r training data: mean = %.5f; std = %.5f' % (y, self.err_t_center,
               self.err_t_scale),
              file=out)
        prediction.predicted = predicted*self.err_t_scale + self.err_t_center
        return prediction.predicted
        # Check that predicted == fitted for the training set.

    def latex_feature_summary_tables(self, longnames=None):
        if longnames is None:
            longnames = {}

        pre_removal = []  # features removed prior to the first fit.
        pre_removal += ['hardcoded'] * len(self.features_delete_hardcoded)
        pre_removal += ['too many NA'] * len(self.features_too_many_na)
        # FIXME: constant features removed in suggested_transform (during test_transforms) are missing.
        pre_removal += ['constant'] * len(self.features_constant)
        pre_removal += ['post-transform constant'] * len(self.features_posttransform_constant)
        pre_removal += ['post-transform normalize constant'] * len(self.features_postnorm_constant)
        pre_features = self.features_delete_hardcoded \
                     + self.features_too_many_na \
                     + self.features_constant \
                     + self.features_posttransform_constant \
                     + self.features_postnorm_constant  # after transform & z-score

        feature_decorred_name = {src: new for new, src in self.decorr_source.items()}
        feature_removal = []
        for ft in self.features:
            if ft in self.must_drop_features:
                feature_removal.append('PCA drop')
            elif ft in self.dropped_todecorr:
                feature_removal.append('decorr denominator')
            else:
                ft = feature_decorred_name.get(ft, ft)  # use the name of the post-decorr var
                if ft in self.bad_props_counts:
                    feature_removal.append('bad prop')
                elif ft in self.suggest_multicolin2:
                    feature_removal.append('collinearity condition')
                elif ft not in self.inde_features2:
                    feature_removal.append('constant after bad prop')
                elif ft not in self.selected_features2:
                    # the feature was not retained by the last LASSO regression.
                    feature_removal.append('Lasso drop')
                else:
                    feature_removal.append('')

        self.feature_track = pd.DataFrame({'Description': [longnames.get(ft, '') for ft in pre_features+self.features],
                                           'Removal': pre_removal + feature_removal},
                                           index=pre_features+self.features)
        # see param_info[['transform', 'decorr']], based on inde_features2 though
        transformed = self.features_posttransform_constant + self.features_postnorm_constant + self.features
        feature_transform = pd.Series([self.suggested_transform[ft].__name__
                                       for ft in transformed],
                                      index=transformed, name='transform')
        feature_decorr = pd.DataFrame.from_records(
                                   [[src, ' '.join(self.decorred[new]), new]
                                    for new,src in self.decorr_source.items()],
                                   index='feature', columns=['feature', 'decorrelation', 'decorred_name'])
        self.feature_transform_decorr = pd.concat((feature_transform, feature_decorr), axis=1).fillna('')

        print('\n### LaTeX table of fitted variables', file=self.out)
        self.display("\n```latex\n" \
                     + self.reslopes2.style\
                        .format_index(escape='latex', axis=0)\
                        .format(precision=4)\
                        .applymap(lambda v: ('cellcolor:{Khaki}' if v<0.01 else
                                        'cellcolor:{LemonChiffon!50}' if v<=0.05 else ''),
                                         subset=['P>|z|', 'Simple regression p-value'])\
                        .to_latex() \
                    + "\n```\n")

        print('\n### LaTeX table of all variables', file=self.out)
        pandas_max_rows = pd.get_option('display.max_rows')
        pd.set_option('display.max_rows', None)
        self.display_html(self.feature_track)
        pd.set_option('display.max_rows', pandas_max_rows)
        self.display("\n```latex\n" \
            + self.feature_track.style.format_index(escape='latex', axis=0).to_latex(environment='longtable') + "\n```\n")

        # Escape special characters for LaTeX:
        latex_escape = str.maketrans({'%': '\\%', '\\': '\\\\', '_': '\\_', '~': '$\\sim$'})
        print('### LaTeX table of variable transforms and decorrelations', file=self.out)
        self.display("\n```latex\n" \
            + self.feature_transform_decorr.style\
                .format_index(escape='latex', axis=0)\
                .format(lambda txt: txt.translate(latex_escape))\
                .to_latex() + "\n```\n")


def refdecorr(a, b, v, cv):
    logger.info('refdecorr: coefficients: a=%s b=%s', a, b)
    return v - (a+b*cv)

def make_ref_todecorr(trained_to_decorr, unregressed_coefs):
    to_decorr = {func: Args.fromcontent(args.content)
                 for func,args in trained_to_decorr.items()}
    
    for i, func_unregress in enumerate((check_unregress, renorm_unregress)):
        unregress_args = to_decorr.pop(func_unregress)
        #unregress_args.pop('sitelnL')
        #except KeyError as err:
        #    logger.error('Key (%d) %r in to_decorr: %s\nOther keys: %s',
        #                 i, func_unregress, func_unregress in to_decorr,
        #                 ' '.join('%r' % k for k in to_decorr.keys()))
        #    raise
        for k, (var, corrvar) in unregress_args.items():
            if isinstance(k, int): k = 'R'+var
            try:
                a, b = unregressed_coefs[k]  # Intercept, slope
            except KeyError:
                logger.warning('Feature %r not in unregressed_coefs.', k)
                continue
            logger.info('make_ref_todecorr: %s ~ %s -> coefficients: %s + x*%s', var, corrvar, a, b)
            ref_unregress_var = partial(check_decorrelator, partial(refdecorr, a, b))
            to_decorr[ref_unregress_var] = Args.fromitems((k, (var, corrvar)))
    return to_decorr


# full_dating_regression init parameters:
# tailored for the Simiiformes *robust* trees dataset!
_must_transform = dict(
        ingroup_nucl_entropy_median=binarize,
        #ingroup_nucl_parsimony_median=binarize, #constant
        ingroup_codon_entropy_median=binarize,
        ingroup_codon_parsimony_median=binarize,
        ingroup_nucl_entropy_mean=sqrt,
        ingroup_nucl_entropy_std=sqrt,  # Should be the same because of decorrelation.
        ingroup_nucl_parsimony_mean=make_best_logtransform, # log (generates Inf for few data)
        ingroup_nucl_parsimony_std=make_best_logtransform,  # log
        rebuilt_topo=binarize,
        consecutive_zeros=binarize, # - triplet_zeros
        sister_zeros=binarize,      # - triplet_zeros
        triplet_zeros=binarize,
        bootstrap_min=notransform,
        prop_splitseq=binarize,
        convergence_warning=binarize,  # notransform
        codemlMP=notransform,
        consecutive_zeros_dS=binarize,
        sister_zeros_dS=binarize,
        triplet_zeros_dS=binarize,
        consecutive_zeros_dN=binarize,
        sister_zeros_dN=binarize,
        triplet_zeros_dN=binarize,
        r2t_dN_mean=make_best_logtransform,  # This one should be called on the data first.
        gb_Nblocks=notransform,
        hmmc_propseqs=notransform,
        freq_null_dS=binarize,
        null_dist_before=binarize,
        null_dS_before=binarize,
        null_dN_before=binarize,
        null_dist_after=binarize,
        null_dS_after=binarize,
        null_dN_after=binarize,
        #null_dN_before=sqrt
        **{'dN_rate'+('_'+setting if setting else ''): make_best_logtransform
           for setting in ('', 'global', 'local', 'nonlocal', 'global_approx',
                           'local_approx', 'global_beastS')#self.rate_settings
          }
       )

# Example parametrisation (seemed reasonable for most regressions so far)
_must_drop_features = {"ls", "seconds",  # ~ ingroup_glob_len
                       "ingroup_std_gaps", # ~ ingroup_std_len
                       "dS_treelen",       # ~dS_rate
                       "dN_treelen",
                       "treelen",
                       "ingroup_nucl_entropy_mean", # ~ ingroup_nucl_parsimony_mean
                       "ingroup_nucl_entropy_std",  # 
                       "ingroup_nucl_parsimony_median",  # 
                       "ingroup_codon_entropy_std",   # ~ ingroup_nucl_parsimony_std
                       "Ringroup_codon_entropy_std",   # ~ ingroup_nucl_parsimony_std
                       "ingroup_codon_entropy_mean",
                       "ingroup_codon_entropy_median",  # ~ ingroup_nucl_entropy_median
                       "ingroup_codon_parsimony_mean",
                       "ingroup_codon_parsimony_median",
                       "ingroup_codon_parsimony_std",
                       "Ringroup_codon_parsimony_std",
                       "r2t_t_mean", "r2t_dS_mean", "r2t_dN_mean",
                       "r2t_t_std",  "r2t_dS_std",  "r2t_dN_std",
                       "bootstrap_mean",  # ~ bootstrap_min
                       "brOmega_skew",
                       "ingroup_std_N",
                       "ingroup_std_CpG",  # ~ ingroup_std_GC
                       "NnonsynSites",  # ~ ls/3 - NsynSites
                       # decorrelated:
                       "ingroup_mean_CpG",
                       # "ingroup_codon_parsimony_std",
                       # "NnonsynSites", "Nsynsites", "brOmega_std",
                       # "ingroup_mean_CpG", "ingroup_std_N", "lnL",
                       # "dS_rate_std", "t_rate_std", "dN_rate_std", "dist_rate_std"
                       # BeastS
                       "treeL_12_mean", #~mean
                       "birthRateY_mean",
                       "ucldMean_12_mean",
                       "ucldStdev_12_mean",
                       "ucldMean_3_mean",
                       "TreeHeight_mean",
                       "rate_12_mean_mean",
                       "gammaShape_mean",
                       "gammaShape_12_mean",
                       "gammaShape_3_mean",
                       "rateAG_mean",
                       "rate_mean",
                       "kappa_12_mean",
                       "kappa_3_mean",
                       "beastclock_rate",  # Instead, use the 'median' estimates
                       "beastclock_rate_std",
                       "beast_rate",
                       "beast_rate_std",
                       "ucldMean_12_mean",
                       "ucldStdev_12_mean",
                       "ucldMean_3_mean",
                       "ucldStdev_3_mean",
                       "ucldMean_12_med", # Redundant with beastclock rate
                       "ucldStdev_12_med",
                       "ucldMean_3_med",
                       "ucldStdev_3_med",
                     }


def renorm_logsqdecorr(y, x):
    return zscore(y - 2*x)

# if data is log, equivalent to log(ey/(ex)**2)

_must_decorr = {
    renorm_decorrelatelogs: Args(
        ),
    decorrelatelogs: Args(
        ('brOmega_std', 'brOmega_mean'),
         #('ingroup_std_N', 'ingroup_mean_N')]
          # Normalise the rate deviations by the rate mean.
        ('treeL_3_stddev', 'treeL_3_med'),
        ('treeL_12_stddev', 'treeL_12_med'),
        *(('r2t_%s_std' % m, 'r2t_%s_mean' % m) for m in MEASURES),
        RsynSites=('NsynSites',    'ls'),
        RnonsynSites=('NnonsynSites', 'ls')  # But this one should be removed.
    ),
    #decorrelate: Args(
    renorm_decorrelate: Args(
    #    sitelnL=('lnL', 'ingroup_glob_len')
    ),
    decorrelate: Args(
    #    sitelnL=('lnL', 'ingroup_glob_len')
    ),
    #decorrelatelogs: Args(
        #[('%s_zeros%s' %(how, ('' if m=='dist' else '_'+m)),
        #  'triplet_zeros'+('' if m=='dist' else '_'+m))
        # for m in MEASURES
        # for how in ('sister', 'consecutive')],
        #[],  # No need anymore since this is now handled in treestats.
        #{}
    #),
    renorm_unregress: Args(),
    check_unregress: Args(
        ('ingroup_nucl_parsimony_std', 'ingroup_nucl_parsimony_mean'),
        ('ingroup_nucl_entropy_std', 'ingroup_nucl_entropy_mean'),
        ('ingroup_codon_parsimony_std', 'ingroup_codon_parsimony_mean'),
        ('ingroup_codon_entropy_std', 'ingroup_codon_entropy_mean'),
        *(('%s_rate_std%s' %(m, ('_'+setting if setting else '')),
            '%s_rate%s' %(m, ('_'+setting if setting else '')))
          for setting in ('', 'global', 'local', 'nonlocal', 'global_approx',
                           'local_approx', 'global_beastS')#self.rate_settings
          for m in MEASURES+['beastS', 'beast', 'beastmedian', 'beastclock', 'beastclockmedian']),
        sitelnL=('lnL', 'ingroup_glob_len')
        )#,
    #renorm_logsqdecorr = Args(
    #     CpG_odds=('ingroup_mean_CpG', 'ingroup_mean_GC')
    #    )
    }
    #renorm_bestlog_and_divide:
    #renorm_bestlog_and_unregress:


# variable name, variable values. Dropped by .isin()
# Must be the decorr variable name.
_must_drop_data = dict(prop_splitseq=(1,),
                       #really_robust=(0,),
                       **{'null_%s_%s' % (m,where): (1,)
                          for m in MEASURES
                          for where in ('before', 'after')},
                       **{'%s_zeros%s' %(what,('' if m=='dist' else '_'+m)): (1,)
                          for m in MEASURES
                          for what in ('triplet', 'Rsister', 'sister',
                                       'Rconsecutive', 'consecutive')})

_protected_features = {'RdS_rate_std', 'dS_rate_std',
                       'RdS_rate_std_local', 'dS_rate_std_local',
                       'RdS_rate_std_nonlocal', 'dS_rate_std_nonlocal',
                       'RbeastS_rate_std', 'beastS_rate_std',
                       'Rbeastmedian_rate_std', 'beastmedian_rate_std',
                       'ingroup_glob_len'}


# anc = 'Catarrhini'
# param = 'um1.new'
# nsCata = age_analyses[anc][param].ns
class full_dating_regression(fullRegression):

    init_vars = ['data',  # Not in fullRegression
                 'dataset_params',  # Not in fullRegression
                 'measures']  # Not in fullRegression. measures of branch lengths and ages.
    init_vars += fullRegression.init_vars
    #            ['same_alls',
    #             'responses',
    #             'features',
    #             'ref_suggested_transform',
    #             'impose_transform',    # global _must_transform
    #             'to_decorr',           # global _must_decorr
    #             'must_drop_features',  # global _must_drop_features
    #             'protected_features',  # global _protected_features
    #             'must_drop_data',      # global _must_drop_data
    #             'out', 'logger', 'widget']

    #init_defaults = {'ref_suggested_transform': dict,
    #                 'impose_transform': dict,
    #                 'must_drop_features': set,
    #                 'to_decorr': Args,
    #                 'protected_features': set,
    #                 'must_drop_data': dict}

    def __init__(self, data, same_alls, dataset_params, responses, features,
                 measures=MEASURES, *args, **kwargs):
        super().__init__(same_alls, responses, features, *args, **kwargs)
        self.data = data
        self.dataset_params = dataset_params
        self.measures = measures

    def do_rates(self, unnamed_rate_setting=None, **named_rate_settings):
        # Compute cs_rates
        dist_measures = ['branch_'+m for m in self.measures]
        if unnamed_rate_setting is not None:
            named_rate_settings[''] = unnamed_rate_setting
        self.rate_settings = named_rate_settings

        toconcat = []
        for setting, rate_args in named_rate_settings.items():
            print('### Compute rates with setting %r and measures %s.'
                    % (setting, ','.join(self.measures)), file=self.out)
            kwargs = {'branchtime': 'median_brlen_%s' % self.measures[0],
                      'taxon_age': 'median_age_%s' % self.measures[0],
                      'dist_measures': dist_measures,  # defaults kwargs
                      **rate_args}
        #display_html('<h3>Compute rates</h3>', raw=True)
            cs_rates = compute_branchrate_std(self.data.ages_controled_withnonrobust,  #TODO: ages_controled_withnonrobust
                                              **kwargs)
            if setting:
                cs_rates.rename(columns=lambda c: c+'_'+setting, inplace=True)
            toconcat.append(cs_rates)

        self.cs_rates = pd.concat(toconcat, axis=1,
                                  sort=False, verify_integrity=True)
        #self.features.extend(self.cs_rates.columns)


    def prepare_design_matrix(self):
        """Run self.do_rates() first."""
        cs_rates = self.cs_rates

        data = self.data
        same_alls = self.same_alls
        dataset_params = self.dataset_params
        responses = self.responses

        print('\n# Merge features', file=self.out)

        data_to_concat = [data.mean_errors, cs_rates, same_alls]
        if dataset_params is not None:
            data_to_concat.append(data.ns[dataset_params])

        self.alls = alls = pd.concat(data_to_concat,
                                     axis=1, join='inner', sort=False,
                                     verify_integrity=True)

        if alls.shape[0] == 0:
            msg = ('NULL inner join from:\nmean_errors\n------\n'
                         + str(mean_errors.iloc[:5, :5]) + '\n'
                         '\ncs_rates\n------\n'
                         + str(cs_rates.iloc[:5, :5]) + '\n'
                         '\nsame_alls\n------\n'
                         + str(same_alls.iloc[:5, :5]) + '\n'
                         )
            if dataset_params is not None:
                msg += ('\ndata.ns[dataset_params]\n------\n'
                        + str(data.ns[dataset_params].iloc[:5, :5]) + '\n')
            raise RuntimeError(msg)


def predict_nonrobust(reg_approx, same_alls, hr, renames, lang_fr=False, nfeatures=3):
    """Essentially compare the histogram of response in the training VS testing
    dataset. Then check feature distributions relative to the 90% decile of
    fitted response.
    """
    hr.mkd('## Values at extreme points of the fit')
    reg_approx.out = hr
    reg_approx.logger = logger
    reg_approx.widget = slideshow_generator(hr)
    reg_approx.set_output()
    model = reg_approx.model  # the regression model, usually StatsModel's OLS.

    reg_approx.predict_extreme_values()

    hr.mkd('## Prediction on new data')
    training_observations = reg_approx.a_n_inde.index
    training_observations2 = reg_approx.a_n_inde2.index
    print('### Check content of `training_observations(2)`', file=hr)
    print(('%d trees in training_observations\n'
           '%d trees in training_observations2\n'
           '%d in intersection\n'
           '%d only in 1\n'
           '%d only in 2\n') % (
               len(training_observations),
               len(training_observations2),
               len(training_observations.intersection(training_observations2)),
               len(training_observations.difference(training_observations2)),
               len(training_observations2.difference(training_observations))),
           file=hr)

    #to_decorr = {func: Args.fromcontent(args.content)
    #             for func,args in reg_approx.to_decorr.items()}
    #to_decorr[check_unregress].args.remove(('ingroup_nucl_parsimony_std', 'ingroup_nucl_parsimony_mean'))
    #to_decorr[check_unregress].pop('sitelnL')
    #ref_unregress_nuclparsstd = partial(check_decorrelator,
    #                                    lambda v,cv: v - (cv*0.508705-0.652472))
    #ref_unregress_lnL = partial(check_decorrelator,
    #                            lambda v,cv: v - (cv*(-0.557353)-0.962986))
    #to_decorr[ref_unregress_nuclparsstd] = Args(('ingroup_nucl_parsimony_std', 'ingroup_nucl_parsimony_mean'))
    #to_decorr[ref_unregress_lnL] = Args(sitelnL=('lnL', 'ingroup_glob_len'))
    to_decorr = make_ref_todecorr(reg_approx.to_decorr, reg_approx.unregressed_coefs)
    must_drop_data = {k: v for k,v in reg_approx.must_drop_data.items()}
    must_drop_data.pop('really_robust', None)
    predicted = reg_approx.predict(same_alls=same_alls,
                                   to_decorr=to_decorr,
                                   must_drop_data=must_drop_data,
                                   out=hr, logger=logger,
                                   widget=slideshow_generator(hr))
    hr.html(predicted.to_frame())  # In transformed scale (unnormalised)
    try:
        scatter_density(reg_approx.refitlasso2.fittedvalues[training_observations2],\
                        predicted[training_observations2], alpha=0.4, s=3)
    except np.linalg.LinAlgError:
        # Singular matrix: GOOD news, the prediction on training data matches the fit!
        plt.scatter(reg_approx.refitlasso2.fittedvalues[training_observations2],
                    predicted[training_observations2], alpha=0.4, s=3)
    ax = plt.gca()
    ax.set_ylabel('Predicted %s on training data' % reg_approx.responses[0])
    ax.set_xlabel('Fitted %s (trained)' % reg_approx.responses[0])
    hr.show()
    plt.close()

    prediction_data = reg_approx.prediction.a_n_inde2
    #prediction_data['data'] = 'testing'
    #prediction_data.loc[training_observations2, 'data'] = 'training'
    #err_threshold = 0.5  # In transformed scale.
    err_threshold = reg_approx.fitted_highest_min_tr
    # Or alternatively, the median of the total/testing set?

    low_err_trained = predicted[training_observations2] <= err_threshold
    low_err_tested = predicted.drop(training_observations2) <= err_threshold
    n_low_err_trained = low_err_trained.sum()
    n_low_err_tested = low_err_tested.sum()
    response_transform = reg_approx.suggested_transform[reg_approx.responses[0]]
    response_retro_trans = make_retro_trans(response_transform.__name__)
    retro_err_thr = response_retro_trans(err_threshold)

    cycled = mpl.rcParams['axes.prop_cycle'].by_key()['color']
    colorpairs = [func(c, frac)
                  for c in cycled
                  for func, frac in [(fade_color_hex, 0.7), (darken_hex, 0.3)]]

    for ft, i in zip(reg_approx.refitlasso2.params.drop('const').index, generate_slideshow(hr)):
        fig = plt.figure()
        gridspec = fig.add_gridspec(2, 2)

        gs00 = gridspec[0,0] if ft in reg_approx.decorr_source else gridspec[:,0]
        ax00 = fig.add_subplot(gs00)
        bicmap = mpl.colors.ListedColormap(['#393b79', '#e6550d'])
        ax00.scatter(ft, ft+'_pred',
                    data=reg_approx.a_n_inde2[[ft]].join(
                        prediction_data[[ft]],
                        rsuffix='_pred', sort=False).loc[training_observations2],
                    alpha=0.5, s=3, c=low_err_trained.loc[training_observations2].astype(float), cmap=bicmap)
        ax00.set_xlabel('Original %s' % ft)
        ax00.set_ylabel('Prediction %s' % ft)
        if ft in reg_approx.decorr_source:
            ax10 = fig.add_subplot(gridspec[1,0])
            srcft = reg_approx.decorr_source[ft]
            ax10.scatter(srcft, srcft+'_pred',
                        data=reg_approx.a_t[[srcft]].join(
                            reg_approx.prediction.a_t[[srcft]],
                            rsuffix='_pred', sort=False).loc[training_observations2],
                        alpha=0.5, s=3, c=low_err_trained.loc[training_observations2].astype(float), cmap=bicmap)
            ax10.set_xlabel('Original %s' % srcft)
            ax10.set_ylabel('Prediction %s' % srcft)

        trained_data = prediction_data[ft].loc[training_observations2]
        tested_data = prediction_data[ft].drop(training_observations2)
        n_trained = np.isfinite(trained_data.values).sum()
        n_tested = np.isfinite(tested_data.values).sum()
        ax1 = fig.add_subplot(gridspec[:,1])
        ax1.set_prop_cycle('color', colorpairs)
        heights, bins, _ = ax1.hist([trained_data[low_err_trained],
                  trained_data[~low_err_trained],
                  tested_data[low_err_tested],
                  tested_data[~low_err_tested]],
                 bins=50, histtype='barstacked',
                 label=['training (tr err ≤ %g) (n=%4d)' % (err_threshold, n_low_err_trained),
                        'training (tr err > %g) (n=%4d)' % (err_threshold, n_trained - n_low_err_trained),
                        'testing  (tr err ≤ %g) (n=%4d)' % (err_threshold, n_low_err_tested),
                        'testing  (tr err > %g) (n=%4d)' % (err_threshold, n_tested - n_low_err_tested)])
        ax1.step(bins, heights[1].tolist() + [heights[1][-1]],
             bins, heights[3].tolist() + [heights[3][-1]],
             where='post', color='k', alpha=0.8)
        ax1.legend()
        hr.show(fig)
        plt.close()

    # Check again the Ringroup_nucl_parsimony_std: why regression so different?
    fig, axes = plt.subplots(ncols=3)
    xvar, yvar = 'ingroup_nucl_parsimony_mean', 'ingroup_nucl_parsimony_std'
    prediction = reg_approx.prediction
    axes[0].set_ylabel(yvar)
    for data, ax, title in zip((reg_approx.a_t, prediction.a_t), axes,
                               ('Original fit', 'Prediction')):
        ax.scatter(xvar, yvar, data=data.drop(training_observations),
                    label='testing', alpha=0.4, s=3)
        ax.scatter(xvar, yvar, data=data.loc[training_observations],
                    label='training', alpha=0.4, s=3)
        x = np.array(ax.get_xlim())
        regdata = data.loc[training_observations, [yvar, xvar]]
        if regdata.isna().any(axis=None):
            logger.warning('%d NaN rows in reg0t data[[%r, %r]]',
                           regdata.isna().any(axis=1).sum(), yvar, xvar)
            regdata = regdata.dropna()
        reg0t = model(regdata[yvar],
                      sm.add_constant(regdata[xvar])).fit()
        a, b = reg0t.params.loc['const'], reg0t.params.loc[xvar]
        ax.plot(x, a + b*x, '--', label='training: %.5f + %.5f * x' % (a, b))
        regdata = data[[yvar, xvar]]
        if regdata.isna().any(axis=None):
            logger.warning('%d NaN rows in reg0 data[[%r, %r]]',
                           regdata.isna().any(axis=1).sum(), yvar, xvar)
            regdata = regdata.dropna()
        if (~np.isfinite(regdata)).any(axis=None):
            logger.warning('%d Inf rows in reg0 data[[%r, %r]]',
                           (~np.isfinite(regdata)).any(axis=1).sum(), yvar, xvar)
            regdata = regdata[np.isfinite(regdata).all(axis=1)]
        reg0 = model(regdata[yvar], sm.add_constant(regdata[xvar])).fit()
        a, b = reg0.params.loc['const'], reg0.params.loc[xvar]
        ax.plot(x, a + b*x, '--', label='complete: %.5f + %.5f * x' % (a, b))
        ax.legend()
        ax.set_title(title)
    ft = 'R'+yvar
    a, b = reg_approx.unregressed_coefs['R'+yvar]
    ax = axes[2]
    data = prediction.a_t.assign(**{ft: (lambda df: df[yvar] - a - b*df[xvar])})
    ax.scatter(xvar, ft, data=data.drop(training_observations),
                label='testing', alpha=0.4, s=3)
    ax.scatter(xvar, ft, data=data.loc[training_observations],
                label='training', alpha=0.4, s=3)
    ax.set_xlabel(xvar)
    ax.set_ylabel(ft)
    ax.set_title('Manual decorr of prediction data')

    fig.suptitle('Relation between source variables (before decorr, transformed)')
    hr.show(fig)
    plt.close()

    print('%s ~ %s: a=%s; b=%s' % (yvar, xvar, a, b), file=hr)
    fig, axes = plt.subplots(ncols=3, figsize=(14,7))
    axes[0].scatter(reg_approx.a_t.loc[training_observations, xvar],
                    prediction.a_t.loc[training_observations, xvar],
                    marker='.', alpha=0.5)
    axes[0].set_title('Original %s' % xvar)
    axes[1].scatter(reg_approx.a_t.loc[training_observations, yvar],
                    prediction.a_t.loc[training_observations, yvar],
                    marker='.', alpha=0.5)
    axes[1].set_title('Original %s' % yvar)
    axes[2].scatter((reg_approx.a_t.loc[training_observations, yvar]
                     - (a+b*reg_approx.a_t.loc[training_observations, xvar])),
                    (prediction.a_t.loc[training_observations, yvar]
                     - (a+b*prediction.a_t.loc[training_observations, xvar])),
                    marker='.', alpha=0.5)
    axes[2].set_title('Residuals of %s' % yvar)
    fig.suptitle('Compare original fit VS prediction features, step-by-step.')
    hr.show(fig)
    plt.close()

    fig, axes = plt.subplots(ncols=3, figsize=(14,7))
    axes[0].scatter(reg_approx.a_n_inde2[ft],
                    (prediction.a_t.loc[training_observations2, yvar]
                   - (a+b*prediction.a_t.loc[training_observations2, xvar])),
                  marker='.', alpha=0.5)
    axes[0].set_ylabel('Manual decorr %s' % yvar)
    axes[0].set_xlabel('Original fit %s' % ft)
    axes[0].set_title('prediction manually decorred\n%s' % yvar)
    axes[1].scatter(reg_approx.a_n_inde2[ft],
                  (reg_approx.a_t.loc[training_observations2, yvar]
                   - (a+b*reg_approx.a_t.loc[training_observations2, xvar])),
                  marker='.', alpha=0.5)
    axes[1].set_xlabel('Original %s' % ft)
    axes[1].set_ylabel('Manual decorr %s' % yvar)
    axes[1].set_title('Original fit manually decorred\n%s' % yvar)
    axes[2].scatter(reg_approx.prediction.a_n_inde2.loc[training_observations2, ft],
                  (reg_approx.prediction.a_t.loc[training_observations2, yvar]
                   - (a+b*reg_approx.prediction.a_t.loc[training_observations2, xvar])),
                  marker='.', alpha=0.5)
    axes[2].set_xlabel('Prediction %s' % ft)
    axes[2].set_ylabel('Manual decorr %s' % yvar)
    axes[2].set_title('prediction VS prediction,\nmanually decorred.')
    fig.suptitle('Compare original fit at different steps.')
    hr.show(fig)
    plt.close()

    predicted_data = pd.DataFrame({'predicted': predicted})
    predicted_data['data'] = 'testing'
    predicted_data.loc[training_observations2, 'data'] = 'training'

    fig, (ax0, ax1) = plt.subplots(ncols=2)
    _, bins, _ = ax0.hist([predicted_data.loc[predicted_data.data=='training', 'predicted'],
              predicted_data.loc[predicted_data.data=='testing', 'predicted']],
             bins=50, histtype='barstacked', label=['training', 'testing'])
    ax0.legend()
    ax0.set_xlabel('Predicted error (in %s scale)' % response_transform.__name__)
    ax0.set_ylabel('Number of trees')

    predicted_data['predicted_retro'] = response_retro_trans(predicted_data['predicted'])
    heights, bins, _ = ax1.hist(
             [predicted_data.loc[predicted_data.data=='training', 'predicted_retro'],
              predicted_data.loc[predicted_data.data=='testing', 'predicted_retro']],
             bins=response_retro_trans(bins), histtype='barstacked',
             label=['training', 'testing'],
             edgecolor='none')
    ax1.step(bins, heights[0].tolist() + [heights[0][-1]], 
             bins, heights[1].tolist() + [heights[1][-1]],
             where='post', color='k', alpha=0.8)
    ax1.legend()
    ax1.set_xscale('log')
    ax1.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax1.xaxis.set_minor_formatter(mpl.ticker.ScalarFormatter())
    xt_labels = [xtl.get_text() for xtl in ax1.get_xticklabels(minor=True)]
    for i, xt in enumerate(ax1.get_xticks(minor=True)):
        if xt in (0.5,2,3,5,20,50):
            xt_labels[i] = '%g' % xt
    ax1.set_xticklabels(xt_labels, minor=True)
    ax1.set_xlabel('Predicted error')
    hr.show(fig)
    plt.close()
    return predict_nonrobust_paperfig(reg_approx, renames, lang_fr, hr, nfeatures)

def predict_nonrobust_paperfig(reg_approx, renames, lang_fr=False, hr=None, nfeatures=3):
    """Figure to save for the paper/thesis."""
    training_observations2 = reg_approx.a_n_inde2.index
    prediction_data = reg_approx.prediction.a_n_inde2
    predicted = reg_approx.prediction.predicted
    err_threshold = reg_approx.fitted_highest_min_tr

    low_err_trained = predicted[training_observations2] <= err_threshold
    low_err_tested = predicted.drop(training_observations2) <= err_threshold
    n_low_err_trained = low_err_trained.sum()
    n_low_err_tested = low_err_tested.sum()

    response_transform = reg_approx.suggested_transform[reg_approx.responses[0]]
    response_retro_trans = make_retro_trans(response_transform.__name__)
    predicted_retro = response_retro_trans(predicted)  #TODO: Move to reg.predict()
    #predicted_data = pd.DataFrame({'predicted': reg_approx.prediction.predicted})
    #predicted_data['data'] = 'testing'
    #predicted_data.loc[training_observations2, 'data'] = 'training'
    #predicted_data['predicted_retro'] = response_retro_trans(predicted_data['predicted'])
    retro_err_thr = response_retro_trans(err_threshold)
    cycled = mpl.rcParams['axes.prop_cycle'].by_key()['color']
    colorpairs = [func(c, frac)
                  for c in cycled
                  for func, frac in [(fade_color_hex, 0.7), (darken_hex, 0.3)]]

    fig = plt.figure(figsize=thesisfigsize)  # contrained_layout=True
    gridspec = fig.add_gridspec(nfeatures+1,2, hspace=0.5, height_ratios=[0.15]+[0.85*1/nfeatures]*nfeatures)  #width_ratios=[2,1]
    ax = fig.add_subplot(gridspec[:,0])  # Entire first column
    retro_err_bins = response_retro_trans(np.histogram(predicted, 50)[1])
    heights, bins, _ = ax.hist(
            [predicted_retro.loc[training_observations2],
             predicted_retro.drop(training_observations2)],
            bins=retro_err_bins,
            histtype='barstacked',
            label=(['entraînement', 'test'] if lang_fr else ['training', 'testing']),
            linewidth=0, edgecolor='none', alpha=0.7) #, rwidth=1.1)
    bin_lim = np.searchsorted(retro_err_bins, retro_err_thr)
    print('len(retro_err_bins)=%d; len(heights)=%d' % (len(retro_err_bins), len(heights)), file=hr)
    ax.fill_between([retro_err_thr] + retro_err_bins[bin_lim:-1].tolist(),
                    0, heights[0][bin_lim-1:],
                    step='post', color=darken_hex(cycled[0])) #, rwidth=1.1)
    ax.fill_between([retro_err_thr] + retro_err_bins[bin_lim:-1].tolist(),
                    heights[0][bin_lim-1:], heights[1][bin_lim-1:],
                    step='post', color=darken_hex(cycled[1])) #, rwidth=1.1)
    ax.step(bins, heights[0].tolist() + [heights[0][-1]],
            bins, heights[1].tolist() + [heights[1][-1]],
            color='k', alpha=0.8, where='post')
    ax.axvline(retro_err_thr, linestyle='--', color='k', alpha=0.8,
               label='err=%.2f\n' % retro_err_thr + ('(10e décile du fit)' if lang_fr else '(10th decile of fit)'))
    #ax.annotate('err = %.2f' % retro_err_thr, (retro_err_thr, )
    print('err threshold = %g (%s-transformed) -> original scale: %g'
          %(err_threshold, response_retro_trans.__name__, retro_err_thr), file=hr)
    lx, ly = retro_err_thr, max(heights[1])/2.
    #_, xright = ax.get_xlim()
    #_, ytop = ax.get_ylim()
    #ax.legend(loc='upper left', bbox_to_anchor=(lx, ly, xright-lx, ytop-ly), bbox_transform=ax.transData)
    ax.set_xscale('log')
    ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    xt_labels = [xtl.get_text() for xtl in ax.get_xticklabels(minor=True)]
    for i, xt in enumerate(ax.get_xticks(minor=True)):
        if xt in (0.5, 2,3,4,5, 20,50):
            xt_labels[i] = '%g' % xt
    ax.set_xticklabels(xt_labels, minor=True)
    ax.set_xlabel('Erreur prédite (Ma/nœud)' if lang_fr else 'Predicted error (My/node)')
    ax.set_ylabel("Nombre d'arbres" if lang_fr else 'Number of trees')
    for i, ft in enumerate(reg_approx.refitlasso2.params.drop('const').index[:nfeatures]):
        axi = fig.add_subplot(gridspec[1+i, 1])
        trained_data = prediction_data[ft].loc[training_observations2]
        tested_data  = prediction_data[ft].drop(training_observations2)
        #TODO: convert to original units. At least not z-scored.
        n_trained = np.isfinite(trained_data.values).sum()
        n_tested  = np.isfinite(tested_data.values).sum()
        n_low_err_trained = low_err_trained.sum()
        n_low_err_tested  = low_err_tested.sum()
        axi.set_prop_cycle('color', colorpairs)
        heights, bins, _ = axi.hist(
                 [trained_data[low_err_trained],
                  trained_data[~low_err_trained],
                  tested_data[low_err_tested],
                  tested_data[~low_err_tested]],
                 bins=50, histtype='barstacked',
                 label=[('entraînement' if lang_fr else 'training')+ ' (err ≤ %.2f) (n=%4d)' % (retro_err_thr, n_low_err_trained),
                        ('entraînement' if lang_fr else 'training')+ ' (err > %.2f) (n=%4d)' % (retro_err_thr, n_trained - n_low_err_trained),
                        ('test   ' if lang_fr else 'testing')+ '  (err ≤ %.2f) (n=%4d)' % (retro_err_thr, n_low_err_tested),
                        ('test   ' if lang_fr else 'testing')+ '  (err > %.2f) (n=%4d)' % (retro_err_thr, n_tested - n_low_err_tested)])
        axi.step(bins, heights[1].tolist() + [heights[1][-1]], 
             bins, heights[3].tolist() + [heights[3][-1]],
             where='post', color='k', alpha=0.8)
        axi.set_xlabel(renames[ft])
        xticks = axi.get_xticks()
        sft = reg_approx.decorr_source.get(ft, ft)
        sft_center = reg_approx.a_t.loc[training_observations2, sft].mean()
        sft_scale = reg_approx.a_t.loc[training_observations2, sft].std(ddof=1)
        print('### %s (from %s)\ntransformed center, scale = %g, %g' % (ft, sft, sft_center, sft_scale), file=hr)
        print('%s transform: %s' % (sft, reg_approx.reslopes2.loc[ft,'transform']), file=hr)
        if sft == ft:  # Useful info only if the feature was not decorrelated.
            #axi.set_xticklabels(['%g' % (xt*sft_scale + sft_center) for xt in xticks])
            print('Ticks in transformed units (prior to center-rescale): ',
                  ' '.join('%g' % (xt*sft_scale + sft_center) for xt in xticks), file=hr)
        if ft in reg_approx.suggested_transform:
            ft_retro_trans = make_retro_trans(reg_approx.suggested_transform[ft].__name__)
            print('Ticks in original units (prior to transform): ',
                  ' '.join('%g' % ft_retro_trans(xt*sft_scale + sft_center) for xt in xticks), file=hr)
        #axi.legend()
    axl = fig.add_subplot(gridspec[0,1])
    axl.legend(*axi.get_legend_handles_labels(), loc='upper right')
    axl.axis('off')
    return fig


# Convert original code names to names displayed in the paper (see fig.1).
dataset_to_pipeline = {'fsahmmc': 'fsa+cleaned,branchMPL',   # (6)
                       'fsahmmcBeastS': 'fsa+cleaned,Beast', # (8)
                       'fsa': 'fsa,branchMPL', # (5)
                       'beastS': 'fsa,Beast'}  # (7)


def fix_decorr(reg):
    var, corrvar = 'dS_rate_std_local', 'dS_rate_local'
    logger.info('Fixing decorr: %s ~ %s (renorm_bestlog_and_divide)', var, corrvar)
    reg.a_n_inde = check_decorrelator(renorm_bestlog_and_divide, reg.a_n_inde,
                                      reg.alls[~reg.na_rows_n & ~reg.inf_rows_n],
                                              (var, corrvar))
    fig = display_decorrelate((0, (var, corrvar)), reg.alls, reg.a_t, reg.a_n_inde,
                              reg.responses[0])
    fig.suptitle('Fixed decorrelation %s ~ %s' % (var, corrvar))
    reg.decorred['RdS_rate_std_local'] = ('divided(bestlogs)', 'dS_rate_local')
    reg.show(fig); plt.close()


def fix_decorr2(reg):
    var, corrvar = 'dS_rate_std_local', 'dS_rate_local'
    logger.info('Fixing decorr: %s ~ %s (renorm_bestlog_and_unregress)', var, corrvar)
    reg.a_n_inde = check_decorrelator(renorm_bestlog_and_unregress, reg.a_n_inde,
                                      reg.alls[~reg.na_rows_n & ~reg.inf_rows_n],
                                              (var, corrvar))
    fig = display_decorrelate((0, (var, corrvar)), reg.alls, reg.a_t, reg.a_n_inde,
                              reg.responses[0])
    fig.suptitle('Fixed decorrelation %s ~ %s' % (var, corrvar))
    reg.decorred['RdS_rate_std_local'] = ('unregressed(bestlogs)', 'dS_rate_local')
    reg.show(fig); plt.close()


# How many layers of parametrisation can I stack?
class regressAncs(object):
    """
    2020/02/21. Output regressions for each anc in an HTML report.
    
    Based on the fsahmmc (fsa+cleaned,MPL) regressions for all ancestors,
    with rate variables local and nonlocal.
    
    Purpose: local & nonlocal are correlated; we started to add "local divided by nonlocal",
    + this class allow any rate reparametrisation
    + control over the protected variables.
    + automatically perform the post-regression joint analysis (display R², plot coeffs).

    """
    # I make it a class, so that internal variables may be retrieved, even if it crashes.

    # "Patching" some of full_dating_regression functionalities
    # or adding some functionality for a set of regressions.
    param_keys = ['rate_settings',
                  'ref_suggested_transform_anc', # = 'Catarrhini'
                  'rate_reprocess', # = dict(newname=<function(cs_rates)>)
                  'fix_decorr']     # = <function(full_dating_reg)>
    param_defaults = {'rate_settings': dict,  #(lambda: {'local': local_rate_args_weighted, 'nonlocal':})
                  'ref_suggested_transform_anc': (lambda: 'Catarrhini'),
                  'rate_reprocess': dict, # = dict(newname=<function(cs_rates)>)
                  'fix_decorr': (lambda: (lambda reg: None))} #sorry.
    # Attributes added by `do()` (not created during __init__)
    todo = ['ref_suggested_transform', 'reports', 'all_reports', 'all_adjR2', 'spec_coefs',
            'union_top5_coefs', 'union_top3_coefs', 'sorted_union_top5_coefs',
            'spec_simple_coefs', 'ref_anc', 'ref_anc_i', 'other_ancs']

    def __getitem__(self, anc):
        return self.all_regressions[anc]

    def __init__(self, description, ancestors, age_analyses, same_alls,
                 out_template='outputs/Regression_{anc}_{desc}.html',
                 **full_dating_kwargs):
        self.description = description
        self.ancestors = ancestors
        self.age_analyses = age_analyses
        self.parametrisation = {k: self.param_defaults.get(k, lambda: None)()
                                for k in self.param_keys}  # It's actually a reparametrisation, or param update.
        # {'rate_settings': } #ref_suggested_transform_anc = 'Catarrhini'
        self.all_regressions = {}
        self.all_coeffs = {}
        self.out_template = out_template
        self.logger_states = {}  # Hold the current handlers
        self.dataset = 'fsahmmc'  # fsa+cleaned

        #default_full_dating_args = Args(
        default_full_dating_kwargs = dict(
                                    #age_analyses[anc][self.dataset],
                                    same_alls=same_alls, #.copy()
                                    dataset_params=list(dataset_params_dS),
                                    responses=['abs_dev_age_dS', 'abs_dev_brlen_dS'],
                                    features=list(features),
                                    measures=['dS'],
                                    impose_transform=deepcopy(_must_transform),
                                    must_drop_features=_must_drop_features.union(
                                            ('gb_percent', 'gb_Nblocks')),
                                    to_decorr=dict((k, as_args(args)) for k,args in _must_decorr.items()),
                                    protected_features=set(_protected_features),
                                    must_drop_data=deepcopy(_must_drop_data),
                                    regul_kw=dict(method='elastic_net', alpha=0.001, L1_wt=1),
                                    logger=logger)
        self.full_dating_kwargs = {**default_full_dating_kwargs,
                                   **full_dating_kwargs}
        self.added_to_params = {'to_decorr': {}}
        self.removed_from_params = {'to_decorr': {}}

    @classmethod
    def from_other(cls, other):
        new = cls(other.description, other.ancestors, other.out_template,
                  **other.full_dating_kwargs)
        for attr in ('parametrisation', 'all_regressions', 'all_coeffs',
                'logger_states', 'dataset', 'ref_suggested_transform',
                'added_to_params', 'removed_from_params'):
            setattr(new, attr, getattr(other, attr))
        for attr in cls.todo:
            try:
                setattr(new, attr, getattr(other, attr))
            except AttributeError:
                pass
        return new

    def parametrise(self, **kwargs):
        for k in self.param_keys:
            try:
                self.parametrisation[k] = kwargs.pop(k)
            except KeyError:
                pass
        deprecated_reparams = ('add_to_drop_features', 'add_to_decorr', 'add_to_protect', 'deprotect', 'add_to_drop_data', 'add_to_must_transform')
        if set(kwargs).intersection(deprecated_reparams):
            warnings.warn(DeprecationWarning("Deprecated parametrisation keys: "
                                             + ' '.join(deprecated_reparams)))
        if set(kwargs).difference(deprecated_reparams):
            logger.warning('Unknown parametrisation keys: %s' % set(kwargs).difference(deprecated_reparams))
        self.parametrisation.update(**kwargs)

        #self.add_to_params(features=list(self.parametrisation['rate_reprocess']))

    def add_to_param(self, **kwargs):
        """Update/union any valid full_dating_regression parameter."""
        for k,v in kwargs.items():
            try:
                current = self.full_dating_kwargs[k]
            except KeyError:
                raise ValueError('Invalid parameter for full_dating_regression: %r' % k)

            if k == 'to_decorr':
                added_to_decorr = self.added_to_params[k]
                for decorrfunc, decorr_args in v.items():
                    current_decorr_args = current.setdefault(decorrfunc, Args())
                    generic_update(current_decorr_args, decorr_args)
                    # Register the additions:
                    added_decorr_args = added_to_decorr.setdefault(decorrfunc, Args())
                    generic_update(added_decorr_args, decorr_args)
            else:
                generic_update(current, v)
                # Register the additions
                try:
                    generic_update(self.added_to_params[k], v)
                except KeyError:
                    self.added_to_params[k] = v  # copy?


    def remove_from_param(self, **kwargs):
        for k,v in kwargs.items():
            try:
                current = self.full_dating_kwargs[k]
            except KeyError:
                raise ValueError('Invalid parameter for full_dating_regression: %r' % k)
            if k == 'to_decorr':
                removed_from_decorr = self.removed_from_params[k]
                for decorrfunc, decorr_args in v.items():
                    current_decorr_args = current.setdefault(decorrfunc, Args())
                    for decorr_item in decorr_args.items():
                        try:
                            current_decorr_args.remove_item(decorr_item)
                        except ValueError as err:
                            err.args += (str(current_decorr_args), current.keys())
                            raise
                    # Register the removals
                    removed_decorr_args = removed_from_decorr.setdefault(decorrfunc, Args())
                    generic_update(removed_decorr_args, decorr_args)
            else:
                # Remove keys from dict, values from set/list, items from Args.
                generic_remove_items(current, v)
                # Register the removals
                try:
                    generic_update(self.removed_from_params[k], v)
                except KeyError:
                    self.removed_from_params[k] = v  # copy?


    def print_parametrisation(self, file=None):
        def fmt_dict(d, valfunc=type):
            return '{\n  ' + ',\n  '.join(
                    '%r: %s' %(k,valfunc(v)) for k,v in d.items()) + '\n}'
        print('# Parametrisation\n\n## Main attributes\n', file=file)
        print('Description: %s\nout_template: %r\ndataset: %r\n' % (
                self.description, self.out_template, self.dataset), file=file)
        print('ancestors = [%s]\n' % ' '.join(self.ancestors), file=file)
        print('age_analyses = %s\n' % fmt_dict(self.age_analyses, fmt_dict), file=file)
        print('## full_dating_regression parameters\n', file=file)
        for p in self.param_keys:
            print('%25s:\t%s\n' % (p, self.parametrisation[p]), file=file)
        for p in sorted(set(self.parametrisation).difference(self.param_keys)):
            print('[Extra] %25s:\t%s\n' % (p, self.parametrisation[p]), file=file)
        print('## Added to full_dating_regression parameters', file=file)
        for p, v in self.added_to_params.items():
            print('%s: %s' % (p,v), file=file)
        print('## Removed from full_dating_regression parameters', file=file)
        for p, v in self.removed_from_params.items():
            print('%s: %s' % (p,v), file=file)

        print('## Global variables\n', file=file)
        print('loggers (name) = %s\n' % [lgr.name for lgr in loggers], file=file)
        #print('myslideshow_css, myslideshow_js', file=file)
        #print('same_alls_fsahmmc = %s\n' % same_alls_fsahmmc.info(), file=file)
        print('features = %s\n' % ' '.join(features), file=file)
        print('dataset_params_dS = %s\n' % ' '.join(dataset_params_dS), file=file)
        print('renames = %s\n' % type(renames), file=file)

    def do(self, refresh=False):
        #global loggers, myslideshow_css, myslideshow_js,
        # age_analyses, same_alls_fsahmmc

        ref_anc = self.parametrisation['ref_suggested_transform_anc']
        ref_anc_i = self.ancestors.index(ref_anc)
        other_ancs = self.ancestors[:ref_anc_i] + self.ancestors[(ref_anc_i+1):]
        self.ref_anc, self.ref_anc_i, self.other_ancs = ref_anc, ref_anc_i, other_ancs

        self.ref_suggested_transform = None
        self.all_reports = dict()
        self.all_reports.update({(anc, filename): t for (anc, filename, t) in
                                 getattr(self, 'reports', [])})
        self.reports = []  # (anc, <html filename>, <timestamp>).
        
        try:
            for anc_i, anc in enumerate([ref_anc] + other_ancs):
                print('# %d. %s #\n------------------' % (anc_i, anc))
                reportfile = self.out_template.format(anc=anc, desc=self.description)
                self.reports.append((anc, reportfile, dt.now()))
                if anc in self.all_coeffs and not refresh:
                    print('Already done. SKIP.')
                    continue
                try:
                    with reroute_loggers(
                            HtmlReport(reportfile,
                                       css=[myslideshow_css], scripts_embedded=[myslideshow_js]),
                            loggers) as hr:

                        # Add a header information with link to previous/next file
                        try:
                            prev_anc, prev_report, _ = self.reports[anc_i-1]
                        except IndexError:
                            prev_anc = 'allspeciations'
                            prev_report = self.out_template.format(anc=prev_anc,
                                                                desc=self.description)
                        try:
                            next_anc = other_ancs[anc_i]
                        except IndexError:
                            prev_anc = 'allspeciations'
                        next_report = self.out_template.format(anc=next_anc,
                                                            desc=self.description)

                        hr.begin('header')
                        #hr.cols()
                        hr.raw('<a href="%s">Previous: %s</a><br />' %
                                (op.basename(prev_report), prev_anc))
                        hr.raw('<a href="%s">Next: %s</a><br />' %
                                (op.basename(next_report), next_anc))
                        hr.end('header')

                        self.all_regressions[anc] = \
                        reg = full_dating_regression(
                                        self.age_analyses[anc][self.dataset],
                                        out=hr, widget=slideshow_generator(hr),
                                        **self.full_dating_kwargs)
                        reg.do_rates(**self.parametrisation['rate_settings'])
                        old_rate_shape = reg.cs_rates.shape
                        reg.cs_rates = reg.cs_rates.assign(**self.parametrisation['rate_reprocess'])
                        hr.print('cs_rates.shape: %s -> %s' % (old_rate_shape, reg.cs_rates.shape))
                        reg.do_undecorred_fit()
                        reg.do_decorr()
                        self.parametrisation['fix_decorr'](reg)
                        reg.do_decorred_fit()
                        coeffs = reg.do_bestfit()
                        reg.do_extreme_points()
                        if anc == ref_anc:
                            self.ref_suggested_transform = reg.suggested_transform

                    self.all_coeffs[anc] = coeffs
                except BaseException:
                    print('ERROR at', anc)
                    raise

            self.summarize()
        finally:
            display_markdown(('Created %d reports:\n\n- ' % len(self.reports)) +
                    '\n- '.join('[%s](%s) (%s)' % (a,f,t.strftime('%Y-%m-%d %H:%M:%S'))
                                for (a,f,t) in self.reports),
                             raw=True)


    def summarize(self, ncoefs=8,  # No more than 8 in the main text.
                  startcoef=0, ncols=2, fixed_order=None, save=False):
        """fixed_order: ordered list of coefficients (strictly included in sorted_coefs)"""
        reportfile = self.out_template.format(anc='allspeciations', desc=self.description)
        self.reports.append(('allspeciations', reportfile, dt.now()))
        self.all_reports[('allspeciations', reportfile)] = self.reports[-1][-1]
        all_coeffs = self.all_coeffs
        with reroute_loggers(
                HtmlReport(reportfile,
                           css=[myslideshow_css], scripts=[myslideshow_js]),
                loggers) as hr:

            prev_anc = self.other_ancs[-1]
            prev_report = op.basename(self.out_template.format(anc=prev_anc, desc=self.description))
            hr.begin('header')
            #hr.cols()
            hr.raw('<a href="%s">Previous: %s</a><br />' % (prev_report, prev_anc))
            hr.raw('<a href="%s">Next: %s</a><br />' % (
                   op.basename(self.out_template.format(anc=self.ref_anc,
                                                        desc=self.description)),
                   self.ref_anc))
            hr.end('header')

            hr.print('\n# adjusted R²\n')
            self.all_adjR2 = pd.Series([self.all_regressions[anc].refitlasso2.rsquared_adj
                                        for anc in self.ancestors],
                                       index=self.ancestors, name='all_adjR2')
            hr.html(self.all_adjR2.to_frame().style.background_gradient(cmap='YlGn'))

            hr.print('\n# Top coefficients\n')
            self.spec_coefs = \
            spec_coefs = pd.concat([all_coeffs[anc][['coef']] for anc in self.ancestors],
                                keys=self.ancestors, axis=1, join='outer', sort=False)\
                      .set_axis(self.ancestors, axis=1, inplace=False)\
                      .assign(absmean=lambda df: df.abs().mean(axis=1))\
                      .sort_values('absmean', ascending=False)

            spec_top5_coefs = [all_coeffs[anc]['coef'].index[1:6].tolist() for anc in self.ancestors]
            union_top5_coefs = set.union(*(set(s) for s in spec_top5_coefs))
            self.union_top5_coefs = union_top5_coefs
            spec_top3_coefs = [all_coeffs[anc]['coef'].index[1:4].tolist() for anc in self.ancestors]
            union_top3_coefs = set.union(*(set(s) for s in spec_top3_coefs))
            self.union_top3_coefs = union_top3_coefs
            
            hr.print('%d in union top5; %d in union top3.' % (
                        len(union_top5_coefs), len(union_top3_coefs)))
            hr.print('Specific in tops VS global ranking:\ntop5: %s\ntop3: %s\n' % (
                union_top5_coefs.difference(spec_coefs.drop('const').index[:5]),
                union_top3_coefs.difference(spec_coefs.drop('const').index[:3])
                ))
            #spec_coefs.drop('const').index[:5].difference(union_top5_coefs)
            hr.html(spec_coefs.style.bar(align='zero')\
                .set_caption('Top coefficients (sorted by absmean)')\
                .set_table_styles([dict(selector='caption',
                                       props=[('caption-side', 'top')])
                                 ])
                )
            #.set_table_attributes(title=))

            # Here, because of .loc[], we need to re-sort the dataframe '--
            self.sorted_union_top5_coefs = \
            sorted_union_top5_coefs = spec_coefs.loc[union_top5_coefs].sort_values('absmean', ascending=False).index.tolist()

            if fixed_order is not None:
                bad_feature_names = set(fixed_order).difference(sorted_union_top5_coefs)
                assert not bad_feature_names, bad_feature_names
                sorted_union_top5_coefs = fixed_order + [ft for ft in sorted_union_top5_coefs if ft not in fixed_order]

            # Also check the simple regression coefs, for inconsistencies
            self.spec_simple_coefs = \
            spec_simple_coefs = pd.concat([all_coeffs[anc][['Simple regression coef']] for anc in self.ancestors],
                                    keys=self.ancestors, axis=1, join='outer', sort=False)\
                          .set_axis(self.ancestors, axis=1, inplace=False)\
                          .assign(absmean=lambda df: df.abs().mean(axis=1))\
                          .sort_values('absmean', ascending=False)
            #sorted_features = spec_coefs.index.tolist()
            hr.html(spec_simple_coefs.style.bar(align='zero')\
                        .set_caption('Top coefficients from simple regressions'
                                     '(sorted by absmean)')\
                        .set_table_styles([dict(selector='caption',
                                               props=[('caption-side', 'top')])
                                         ])
                        )
                        #.set_table_attributes(title=))

            # Paper figures: top coefficients
            hr.print('## Top 5 coefficients, sorted by their mean absolute value across ancestors/imposed order')

            if ncoefs is None:
                ncoefs = len(sorted_union_top5_coefs)
            plotdata = spec_coefs.loc[union_top5_coefs].drop('absmean', 1) #.rename(renames)
            # The error bars for the barplot:
            spec_conf_int = pd.concat(
                                [self.all_coeffs[anc][['[0.025', '0.975]']]
                                 for anc in self.ancestors],
                                axis=1, sort=False, keys=self.ancestors, verify_integrity=True)
            # -> index: coefficients
            # -> columns: MultiIndex(( anc, inf/sup ))

            nrows = ncoefs // ncols + (ncoefs % ncols >0)
            fig, axes = plt.subplots(nrows, ncols,
                                     figsize=(ncols*5, 2+nrows*1.5), sharex=True, sharey=True,
                                     gridspec_kw={'hspace': 0}, squeeze=False)
            n_anc = len(self.ancestors)
            cmap = plt.get_cmap('tab20b', n_anc)  # 'set3', 'viridis'
            barcolors = [cmap(i) for i in range(n_anc)]

            letters = 'abcdefghijklmnopqrstuvwxyz'
            x = np.arange(n_anc)
            for lett, ax, coef in zip(letters, axes.flat, sorted_union_top5_coefs[startcoef:]):
                heights = plotdata.loc[coef, self.ancestors]
                yerr = np.abs(spec_conf_int.loc[coef].unstack(level=-1).values.T - heights.values)
                ax.bar(x, heights, color=barcolors, width=0.9, yerr=yerr,
                       capsize=2, error_kw=dict(elinewidth=1, alpha=0.5))
                ax.set_ylabel('\n'.join(textwrap.wrap(renames.get(coef, coef), width=25, break_long_words=False))
                             , rotation='horizontal', ha='right', va='center')
                ax.yaxis.set_ticks_position('right')
                ax.grid(True, axis='y')
                ax.grid(False, axis='x')
                ax.annotate(lett+')', (-0.1, 1.), weight='bold', xycoords='axes fraction', va='top')
            ax.set_xticks(x)

            #x_shared_group = axes[-1, -1].get_shared_x_axes()
            for i, ax in enumerate(axes[-1]):
                if i >= (ncoefs % ncols) > 0:
                    # Nothing was plotted: remove the axis
                    ax.axis('off')
                    #x_shared_group.remove(ax)
                    ax = axes[-2, i]
                    ax.tick_params('x', labelbottom=True)
                ax.set_xticklabels(self.ancestors, rotation=45, va='top', ha='right');

            fig.tight_layout()
            #return fig

            if save:
                for ext in ('pdf', 'png'):
                    fig.savefig('../fig/speciations_coefs-top5_%s_%s_ylab-h-%d-%d.%s'
                                % (dataset_to_pipeline[self.dataset].replace(',','+'),
                                   self.description, startcoef+1, startcoef+ncoefs, ext),
                                bbox_inches='tight')
            hr.show(fig)
            plt.close(fig)

            #fig.savefig('../fig/speciations_coefs-top5_fsahmmc+mpl_ylab-h-8first.pdf', bbox_inches='tight')
            #fig.savefig('../fig/speciations_coefs-top5_fsahmmc+mpl_ylab-h-8first.png', bbox_inches='tight')

            hr.print('Top 5 coefficients, from the multiple VS simple regression')
        #def plot_anc_coeffs(plotdata, ncoefs=8, startcoef=0, ncols=2,
                             #plotdata_simple=None):
            
            plotdata2 = spec_simple_coefs.loc[union_top5_coefs].drop('absmean', 1) #.rename(renames)
            spec_conf_int2 = pd.concat(
                                [self.all_coeffs[anc][
                                    ['Simple regression [0.025',
                                     'Simple regression 0.975]']]
                                 for anc in self.ancestors],
                                axis=1, sort=False, keys=self.ancestors, verify_integrity=True)
            # -> index: coefficients
            # -> columns: MultiIndex(( anc, inf/sup ))

            nrows = ncoefs // ncols + (ncoefs % ncols >0)
            fig, axes = plt.subplots(nrows, ncols,
                                     figsize=(ncols*5, 2+nrows*1.5), sharex=True, sharey=True,
                                     gridspec_kw={'hspace': 0}, squeeze=False)
            n_anc = len(self.ancestors)
            cmap = plt.get_cmap('tab20b', n_anc)  # 'set3', 'viridis'
            barcolors = [cmap(i) for i in range(n_anc)]

            x = np.arange(n_anc)
            for ax, coef in zip(axes.flat, sorted_union_top5_coefs[startcoef:]):
                heights = plotdata.loc[coef, self.ancestors]
                yerr = np.abs(spec_conf_int.loc[coef].unstack(level=-1).values.T - heights.values)
                ax.bar(x, heights, color=barcolors, width=-0.4, align='edge',
                       yerr=yerr, capsize=2, error_kw=dict(elinewidth=1, alpha=0.5))
                heights2 = plotdata2.loc[coef, self.ancestors]
                yerr2 = np.abs(spec_conf_int2.loc[coef].unstack(level=-1).values.T - heights2.values)
                ax.bar(x, heights2, color='gray', width=0.4, align='edge',
                       yerr=yerr2, capsize=2, error_kw=dict(elinewidth=1, alpha=0.3))
                
                ax.set_ylabel('\n'.join(textwrap.wrap(renames.get(coef, coef), width=25, break_long_words=False))
                             , rotation='horizontal', ha='right', va='center')
                ax.yaxis.set_ticks_position('right')
                ax.grid(True, axis='y')
                ax.grid(False, axis='x')
            ax.set_xticks(x)

            #x_shared_group = axes[-1, -1].get_shared_x_axes()
            for i, ax in enumerate(axes[-1]):
                if i >= (ncoefs % ncols) > 0:
                    # Nothing was plotted: remove the axis
                    ax.axis('off')
                    #x_shared_group.remove(ax)
                    ax = axes[-2, i]
                    ax.tick_params('x', labelbottom=True)
                ax.set_xticklabels(self.ancestors, rotation=45, va='top', ha='right');
                    
            fig.tight_layout()
            if save:
                for ext in ('pdf', 'png'):
                    fig.savefig('../fig/speciations_simple+multiple_coefs-top5_%s_%s_ylab-h-%d-%d.%s'
                                % (dataset_to_pipeline[self.dataset].replace(',','+'),
                                   self.description, startcoef+1, startcoef+ncoefs, ext),
                                bbox_inches='tight')
            hr.show(fig)
            plt.close(fig)



from siphon import dependency_func, dependency, auto_internalmethod
from functools import wraps


class analysis(object):
    """Analysis pipeline and environment"""
    
    @property
    def measures(self):
        return self._measures
    
    @property
    def dist_measures(self):
        return self._dist_measures

    @property
    def age_measures(self):
        return self._age_measures

    @property
    def rate_measures(self):
        return self._rate_measures

    @measures.setter
    def set_measures(self, value):
        self._measures = value
        self._dist_measures = ['branch_%s' %m for m in self._measures]
        self._rate_measures = ['rate_%s' %m for m in self._measures]
        self._age_measures =  ['age_%s' %m for m in self._measures]


    def __init__(self, age_file, stats_template,
                 stattypes=('al', 'tree', 'codeml'),
                 measures = ['dist', 'dS', 'dN', 't'],
                 control_condition='really_robust && aberrant_dists == 0',
                 phyltreefile=None,
                 control_ages_CI=None,
                 common_info=None,
                 al_params=None,
                 tree_params=None,
                 cl_params_restricted=None):
        # Add a settings attribute that stores the keys
        self.settings = set()
        for k,v in locals():
            if k != 'self':
                setattr(self, k, v)
                self.settings.add(k)

    @dependency
    def phyltree(self):
        return myPhylTree.PhylogeneticTree(phyltreefile)
    
    @dependency
    def ts(self):
        return load_stats_tree(stats_template.format(stattype=stattypes[1]))

    @dependency
    def aS(self):
        return load_stats_al(stats_template.format(stattype=stattypes[0]))

    @dependency
    def cs(self):
        return load_stats_codeml(stats_template.format(stattype=stattypes[2]))

    @dependency
    def clS(self):
        return load_any_stats(stats_template.format(stattype=stattypes[3]))

    @dependency_func
    @wraps(load_prepare_ages)
    def load_prepare_ages(self):
        return load_prepare_ages(self.ages_file, self.ts)

    @dependency
    def ages_treestats(self):
        return self.load_prepare_ages()[0]

    @dependency
    def ns(self):
        return self.load_prepare_ages()[1]

    @auto_internalmethod
    @wraps(add_control_dates_lengths)
    def get_control_dates_length(self, *args, **kwargs): 
        return add_control_dates_lengths(*args, **kwargs)

    @dependency_func
    def getcheck_control_dates_lengths(self, root):
        a_ctrled, ctrl_ages, ctrl_brlen = self.get_control_dates_length()
        check_control_dates_lengths(ctrl_brlen, self.phyltree, root)
        return a_ctrled, ctrl_ages, ctrl_brlen

    @property
    def ages_controled(self):
        return self.get_control_dates_length[0]

    @property
    def control_ages(self):
        return self.get_control_dates_length[1]
    
    @property
    def control_brlen(self):
        return self.get_control_dates_length[2]

    @dependency_func
    def mean_errors(self, control='median'):
        return compute_dating_errors(self.ages_controled, control=control)

    @dependency_func
    def mean_speciation_errors(self, control='median'):
        pass

    def compare_error_with(other_analysis):
        pass

    compute_branchrate_std = dependency(auto_internalmethod(compute_branchrate_std))

    @dependency
    def csrates(self):
        return self.compute_branchrate_std

    @dependency
    def alls(self):
        return pd.concat((self.mean_errors,
                          self.aS[self.common_info + self.al_params],
                          self.ts[self.tree_params],
                          self.cs[self.cl_params_restricted],
                          self.cs_rates),
                         axis=1, join='inner')


#for k in _must_decorr.keys():
#    print('** 0x%x %s' % (id(k), k))
