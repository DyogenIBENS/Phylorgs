#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Display a multiple sequence alignment (nucleotide or codon), and the 
conservation scores along it."""


from sys import stdin
#from time import clock # to display execution time

import os.path
import re

from itertools import product, combinations

import numpy as np
import matplotlib as mpl
try:
    mpl.use('Agg', warn=False)
except TypeError:
    # keyword warn disappeared in recent matplotlib
    mpl.use('Agg')
import matplotlib.pyplot as plt
from Bio import AlignIO
try:
    import argparse_custom as argparse
except ImportError:
    import argparse

from seqtools.symbols import GAPS, NUCLEOTIDES, AA
from seqtools.arrayal import *
from phylorg import parsimony_score
from datasci.graphs import stackedbar, plottree

import logging
logger = logging.getLogger(__name__)

#try:
#    mpl.style.use('softer')
#except OSError:
# From Seaborn
mpl.rcParams['axes.prop_cycle'] = mpl.cycler('color',
                ['4C72B0', '55A868', 'C44E52', '8172B2', 'CCB974', '64B5CD'])
mpl.rcParams['grid.alpha'] = 0.5
mpl.rcParams['grid.linestyle'] = '--'
#    pass

mpl.rcParams['axes.grid'] = True
mpl.rcParams['axes.grid.axis'] = 'x'

# Change all black to dark grey
grey10 = '#1a1a1a'
for param, paramval in mpl.rcParamsDefault.items():
    if paramval == 'k':
        mpl.rcParams[param] = grey10


EXT2FMT = {'.fa':    'fasta',
           '.faa':   'fasta',
           '.fasta': 'fasta',
           '.mfa':   'fasta',
           '.phy':   'phylip-relaxed',
           '.sto':   'stockholm',
           '.nx':    'nexus',
           '.nex':   'nexus'}


def filename2format(filename):
    _, ext = os.path.splitext(filename)
    return EXT2FMT[ext]


def pairs_score(vint, indexpairs, dist_mat):
    nsites = vint.shape[1]
    seqpairs = [np.stack((vint[i,:], vint[j,:])) for i, j in indexpairs]
    pair_scores = np.array([[dist_mat[tuple(p[:,site])] for site in range(nsites)]
                            for p in seqpairs])
    return pair_scores.sum(axis=0)

def sp_score(vint, dist_mat):
    """Compute the sum of all pairs score, column-wise"""
    if dist_mat.shape[0] - 1 < vint.max():
        raise ValueError('dist_mat size (current: %d) must be larger than the max residue code (current: %d); requires 1(gap) + alphabet size + 1(unknown residue)' % (dist_mat.shape[0], vint.max()))
    return pairs_score(vint, combinations(range(vint.shape[0]), 2), dist_mat)

def part_sp_score(vint, split, dist_mat):
    assert len(split) == 2
    return pairs_score(vint, product(*split), dist_mat)

# ~~> datasci.stats?
def cov(X,Y, axis=0):
    mx = X.mean(axis=axis)
    my = Y.mean(axis=axis)
    return ((X - mx) * (Y - my)).mean(axis=axis)


# ~~> datasci.stats?
def pearson_coeff(X, Y, axis=0):
    """Column-wise pearson coefficient"""
    coeff = cov(X, Y, axis)
    coeff /= (X.std(axis=axis) * Y.std(axis=axis))
    return coeff


#from dendro.any import biophylo
def get_items_biophylo(tree, nodedist):
    return [(child, child.branch_length) for child in nodedist[0].clades]

def get_label_biophylo(tree, node):
    return node.name


def annotate_summary(ax, values):
    fmt = '%g' if len(np.asarray(values).shape) <= 1 else '%s'
    with np.printoptions(precision=4):
        summary = 'Mean: {0}\nMedian: {0}\nStd-dev: {0}'.format(fmt) % (
                    np.nanmean(values, axis=-1),  # axis=-1
                    np.nanmedian(values, axis=-1),
                    np.nanstd(values, axis=-1))
        #FIXME: no effect of precision when shape of arrays have dimension 1.
    ax.text(0.995, 0.99, summary, transform=ax.transAxes,
            verticalalignment='top', horizontalalignment='right')


PLOTS = ['gap_prop', 'entropy', 'gap_entropy', 'al', 'sp_score']
COMP_PLOTS = ['manhattan', 'pearson', 'split_pairs_score']

def al_colorbar(ax, cax=None, altype='nucl', orientation='vertical', **kwargs):
    kwargs = {'shrink':0.05, **kwargs}
    im = ax.images[-1]
    cbar = plt.colorbar(im, ax=(ax if cax is None else None),
                        cax=cax,
                        orientation=orientation, extend='both', **kwargs)
    if cax is not None and orientation == 'horizontal':
        cax_pos = cax.get_position()
        cax.set_position(cax_pos.shrunk_to_aspect(1/20))
    cbar.set_label('Residues')
    #ticks = cbar.get_ticks()
    #cbar.set_ticks()
    #nticks = len(ticks)
    logger.debug('cbar lim= %s; ylim= %s; xlim= %s', cbar.get_clim(),
                 cbar.ax.get_ylim(), cbar.ax.get_xlim())
    logger.debug('norm.vmax = %g; norm.N = %g; cmap.N = %g',
                 im.norm.vmax, im.norm.N, im.cmap.N)
    gap = GAPS[0]
    if altype=='nucl':
        residues = NUCLEOTIDES
        rescode = NUCL2INT
    elif altype=='aa':
        residues = AA
        rescode = AA2INT
    elif altype=='codon':
        residues = CODONS + list(CODONS_STOP)
        rescode = CODON2INT
        gap *= 3

    ticks = np.arange(1, len(residues)+1)
    ticklabels = residues
    # Put multiple residues in one label if redundant colors.
    if isinstance(im.cmap, mpl.colors.ListedColormap):
        sep = '\n' if orientation=='horizontal' else ' '
        ticklabels = [[]]
        prev_resid_col = tuple(im.cmap(im.norm(1)))
        ticks = [1]
        logger.debug('residue color 0: %s', prev_resid_col)
        for i, resid in enumerate(residues, start=1):
            if rescode[resid] != i:
                raise RuntimeError('Residue integer code is mismatched! Wrong color legend.')

            resid_col = tuple(im.cmap(im.norm(i)))
            #logger.debug('resid %d %r col %s; != %s', i, resid, resid_col,
            #             resid_col != prev_resid_col)
            if resid_col == prev_resid_col:
                ticklabels[-1].append(resid)
            else:
                ticks.append(i)
                #try:
                ticklabels[-1] = sep.join(ticklabels[-1])
                #except IndexError:
                #    pass
                ticklabels.append([resid])
                prev_resid_col = resid_col
        #ticks.append(i)
        #if not ticklabels[-1]:
        #    ticklabels[-1].append(resid)
        ticklabels[-1] = sep.join(ticklabels[-1])
        logger.debug('Pooled %d ticks %s, %d labels %s', len(ticks), ticks,
                     len(ticklabels), ticklabels)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(ticklabels)
    logger.debug('cbar lim= %s; ylim= %s; xlim= %s', cbar.get_clim(),
                 cbar.ax.get_ylim(), cbar.ax.get_xlim())
    #TODO
    cbar.ax.text(1, -0.05, gap, fontsize='x-small', clip_on=False) # under
    #cbar.ax.annotate(gap, (1, -0.05), xytext=(2, 0),
    #                 textcoords='offset points',
    #                 fontsize='x-small', annotation_clip=False) # under
    cbar.ax.tick_params(labelsize='x-small')
    return cbar


class AlignPlotter(object):

    #colorcycle = plt.rcParams['axes.prop_cycle']

    # FIXME: these should be instance attributes. Otherwise, when modified in 
    # an instance, they are changed for each newly created instance...
    @classmethod
    def set_plots(cls):
        plot_properties = {'al': {'title': 'Global scoring (all sequences)'},
                                  #'ylabel': 'Alignment'},
                           'gap': {'ylabel': 'Proportion of gaps\n(G)'},
                                   #'ylim': (-0.05,1.05)},
                           'pars': {'ylabel': 'Parsimony score\n(min number of changes per branch)'},
                           'entropy': {'ylabel': 'Entropy\n(H)'},
                           'gap_entropy': {'ylabel': 'Gap-entropy score:\n(1-H)*(1-G)'},
                           'sp': {'ylabel': 'SP score\n(Sum-of-pair differences)'},
                           'manh': {'title': 'Difference scoring between parts',
                                    'ylabel': 'manhattan distance'},
                           'pearson': {'ylabel': "Pearson's correlation coefficient"},
                           'part_sp': {'ylabel': 'sum of pair differences'},
                           'part_pars': {'ylabel': 'parsimony score within parts'},
                           'tree': {}}

        default_step = ('step', {'where': 'mid', 'alpha': 0.65})
        default_bar = ('bar', {'width': 1})
        default_stacked = (stackedbar, {'width': 1, 'edgecolor': 'none'})

        #from collections import namedtuple
        #funcArgs = namedtuple('funcArgs', 'func args')

        # { plot: [list of tuples(function, kwargs)] }
        plot_funcs = {'al':          [('imshow', {'aspect': 'auto'})], #(al_colorbar, {})],
                      #'al':          [('pcolormesh',   {'edgecolors': 'None'}),
                      #                ('invert_yaxis', {})],
                      #'gap':         default_step, #default_bar, "stackplot"
                      'gap':         [default_stacked],
                      'pars':        [default_bar, (annotate_summary,)],
                      'entropy':     [default_bar, (annotate_summary,)],
                      'gap_entropy': [default_bar, (annotate_summary,)],
                      'sp':          [default_bar],
                      'manh':        [default_step],
                      'pearson':     [default_step],
                      'part_sp':     [default_step],
                      'part_pars':   [default_stacked,
                                      (annotate_summary,),
                                      ('legend', {'bbox_to_anchor': (0,0,0.5,1)})],
                      'tree':        [(plottree, {'get_items': get_items_biophylo,
                                                  'get_label': get_label_biophylo,
                                                  'label_params': {'fontsize': 'x-small'}}),
                                      ('tick_params', {'labelright': False}),
                                      ('grid', {'b': False})]}
        return plot_properties, plot_funcs


    def __init__(self, alint, seqlabels=None, valid_range=(1,64), tree=None,
                 topology_only=False):
        """Initialize the Plotter instance from a matrix of integers.
        
        - alint: alignment provided as a numpy 2D array of integers.
        - tree: tree in Bio.Phylo format.
        """
        self.plot_properties, self.plot_funcs = self.set_plots()
        self.alint = alint
        self.x = np.arange(self.alint.shape[1]) # X values for plotting
        logger.debug('Alignment shape: %s', alint.shape)
        
        self.malint = np.ma.array(alint, mask=(alint==0))
        self.seqlabels = seqlabels
        self.plotlist = ['al']
        self.plotdata = {'al': (self.malint,)}
        #print(np.unique(self.malint))
        #print(np.unique(self.malint - self.malint.max()))
        self.valid_range = valid_range
        self.tree = tree
        if tree is not None:
            self.plotlist.insert(0, 'tree')
            self.plotdata['tree'] = (tree,)
            self.plot_funcs['tree'][0][1].update(topology_only=topology_only)

        #cmap_size = valid_range[1] - valid_range[0] # If len(valid_range) > 1
        # Dark2, tab20b, tab20, Set1
        #FIXME: valid_range for codons should exclude the stops: (1,61)
        cmap = plt.get_cmap('tab20b', valid_range[1] - valid_range[0] + 1) #, cmap_size)
        cmap.set_over((0.2, 0.2, 0.2))
        cmap.set_under('w')
        # a norm is needed to set bounds!!!
        
        # Directly use integer values to get colors in cmap: Same speed but less clear.
        #self.plotdata = {'al': (self.malint - valid_range[0],)} # With NoNorm
        #norm = mpl.colors.NoNorm(0, valid_range[1] - valid_range[0])

        bounds = np.arange(valid_range[0]-0.5, valid_range[1] + 0.6)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        # Maybe a ListedColormap is better.
        self.plot_funcs['al'][0][1].update(cmap=cmap, norm=norm)
        self.plot_properties['al']['yticks'] = np.arange(alint.shape[0])
        #TODO: Add the *last* xtick showing alignment length.
        #self.plot_properties['al']['xlim'] = (0, alint.shape[1]) #pcolormesh
        self.plot_properties['al']['xlim'] = (-0.5, alint.shape[1]-0.5) #imshow
        #TODO: Add plt.colorbar for residue colors.

        if seqlabels is not None:
            self.plot_funcs['al'].append(('set_yticklabels',
                                          {'labels': seqlabels,
                                           'fontsize': 'x-small'}))

    @classmethod
    def fromfile(cls, infile, format=None, altype='nucl', allow_N=False,
                 ungap=True, slice=None, records=None, recordsfile=None,
                 treefile=None, topology_only=False):
        align = AlignIO.read(infile,
                             format=(format or
                                     filename2format(getattr(infile, 'name',
                                                             infile)))
                            )

        al_len = align.get_alignment_length()
        if altype == 'codon' and al_len % 3:
            logger.error("Not a codon alignment!")

        tree = None

        if records or recordsfile or treefile:
            record_order = None
            if records:
                record_order = [int(r) for r in records.split(',')]
                records = None
            elif recordsfile:
                with open(recordsfile) as recf:
                    records = [line.rstrip() for line in recf]
            elif treefile:
                #try:
                #    import ete3
                #    tree = ete3.Tree(treefile, format=1)
                #    get_leaves = tree.get_leaves
                #except ImportError:
                from Bio.Phylo import read as phyloread
                tree = phyloread(treefile, 'newick')
                get_leaves = tree.get_terminals
                records = [leaf.name for leaf in get_leaves()]

            align = reorder_al(align, records, record_order)

        if slice:
            slstart, slend = [int(pos) for pos in slice.split(':')]
            if not nucl:
                slstart *= 3
                slend   *= 3
            if not (0 <= slstart < slend <= al_len):
                raise IndexError("Slice indices (%s,%s) out of bounds." \
                                 % (slstart, slend))
            align = align[:,slstart:slend]
        else:
            slstart, slend = None, None

        seqlabels = [record.name for record in align]
        # Convert matrix of codon/nucleotide strings to matrix of integers.
        alint = al2int(align, altype, allow_N)

        # ungap: remove columns being 100% gaps
        if ungap:
            alint = alint[:, alint.sum(axis=0) > 0]

        # Init and setup the class instance
        valid_range = (1, len(NUCLEOTIDES)) if altype=='nucl' else (1, 64) if altype=='codon' else (1,len(AA))
        instance = cls(alint, seqlabels, valid_range, tree, topology_only)
        instance.fromfile_params = {'infile': infile, 'format': format,
                                    'slice': (slstart, slend),
                                    'records': records}
        instance.align = align
        instance.al_len = al_len
        instance.altype = altype
        instance.residue_distmat = make_unif_codon_dist() if altype=='codon' else (1 - np.eye(valid_range[1]+2))  #TODO: distance matrices
        instance.residue_distmat[-1,:] = instance.residue_distmat[:,-1] = np.NaN
        instance.slstart = slstart
        return instance


    def measure(self):
        """Compute column-wise global measures: residue frequencies, entropy, 
        gap proportion, entropy-gap score, SP-score"""

        max_valid = self.valid_range[1]
        
        al_freqs = freq_matrix(self.alint, minlength=(max_valid+2))
        gap_prop = al_freqs[0,:]
        invalid_prop = al_freqs[(max_valid+1):, :].sum(axis=0)
        if invalid_prop.any():
            #print(al_freqs[(max_valid+1):, :])
            self.gap_prop = np.stack((gap_prop,
                                      invalid_prop), axis=0)
                                      #gap_prop + invalid_prop), axis=1)
            self.plot_funcs['gap'].append(('legend',
                                           {'labels': ('Gaps only',
                                                       #'Gaps + 
                                                       'Invalid values (N)')}))
        else:
            # change from 1D to 2D. Here, same as array.reshape((array.shape[0], 1))
            self.gap_prop = gap_prop[np.newaxis]

        # Convert zeros to ones so that x * log(x) = 0
        self.freqs = al_freqs[1:,:]
        # scipy.stats.entropy
        self.entropy = freqs2entropy(self.freqs)
        #self.gap_entropy = (1 - self.entropy) * (1 - self.gap_prop[:, 1])
        self.gap_entropy = (1 - self.entropy) * (1 - self.gap_prop.sum(axis=0))
        self.sp_score = sp_score(self.alint, self.residue_distmat)

        self.plotlist.extend(('gap', 'entropy', 'gap_entropy', 'sp'))
        self.plotdata.update(gap=(self.x, self.gap_prop,),
                             entropy=(self.x, self.entropy,),
                             gap_entropy=(self.x, self.gap_entropy,),
                             sp=(self.x, self.sp_score,))
        self.plot_funcs['entropy'][1] = self.plot_funcs['entropy'][1] \
                                     + ({'values': self.entropy},) 
        self.plot_funcs['gap_entropy'][1] = self.plot_funcs['gap_entropy'][1] \
                                     + ({'values': self.gap_entropy},) 

        if self.tree is not None:
            if hasattr(self.tree, 'children'):
                get_children = lambda tree, node: node.children
            else:
                get_children = lambda tree, node: node.clades
            self.parsimony_score = parsimony_score(self.alint, self.tree,
                                                   seqlabels=self.seqlabels,
                                                   minlength=(max_valid+2),
                                                   get_children=get_children)
            self.plotlist.insert(self.plotlist.index('al')+1, 'pars')
            self.plotdata['pars'] = (self.x, self.parsimony_score)
            # annotate_summary:
            self.plot_funcs['pars'][1] += ({'values': self.parsimony_score},)


    def annotate_parts(self, ax):
        """Draw arrows with annotation to show 'compare_parts' on the
        alignment plot."""
        #x0, x1 = ax.get_xlim()
        #xtext = x1 + (x1-x0) * 1.05
        x1, xtext = 1, 1.02
        hlines = []
        for part_number, part in enumerate(self.parts, start=1):
            ytext = np.mean(part) / self.alint.shape[0]
            ax.annotate(str(part_number), xy=(x1, ytext), xytext=(xtext, ytext),
                        xycoords='axes fraction', alpha=0.8)
            for row in part:
                ax.annotate("", xy=(x1, (row+0.5)/self.alint.shape[0]),
                            xytext=(xtext, ytext),
                            xycoords='axes fraction',
                            arrowprops={'arrowstyle': '->', 'alpha':0.8},
                            alpha=0.8)
            # if part is contiguous:
            row0, row_end = min(part), max(part)+1
            if set(part) == set(range(row0, row_end)):
                if row0 not in hlines:
                    ax.axhline(row0-0.5, linestyle='--', color='#191919', alpha=0.5)
                    hlines.append(row0)
                if row_end not in hlines:
                    ax.axhline(row_end-0.5, linestyle='--', color='#191919', alpha=0.5)
                    hlines.append(row_end)


    def comp_parts(self, compare_parts):
        """Compute column-wise metrics of difference between two parts (groups
        of records) of the dataset:
            - manhattan distance between frequencies
            - pearson coefficient between counts
            - sum of pair differences between any possible pair of distinct parts
        """
        ### TODO: add a parsimony comparison: 1 if null intersection.
        if isinstance(compare_parts, str):
            compare_parts = [part.rstrip(')').lstrip('(').split(',')
                             for part in compare_parts.split(';')]
        # Convert to the proper type
        self.parts = []
        for part in compare_parts:
            parti = []
            for subpart in part:
                try:
                    parti.append(int(subpart))
                except ValueError:
                    try:
                        parti.extend(range(*(int(i) for i in subpart.split(':'))))
                    except ValueError:
                        label_reg = re.compile(subpart)
                        parti.extend((i for i,lab in enumerate(self.seqlabels)
                                      if label_reg.match(lab)))
            self.parts.append(parti)

        # The complementary sequences.
        out_part = list(set(range(self.alint.shape[0]))
                        .difference(set().union(*(self.parts))))
        all_parts = self.parts + [out_part]

        #assert len(self.parts) == 2

        self.part_al = [self.alint[part,:] for part in all_parts[:2]]

        self.part_freqs = fmat1, fmat2  = [freq_matrix(alpart) for alpart in self.part_al]
        self.part_counts = cmat1, cmat2 = [count_matrix(alpart) for alpart in self.part_al]

        self.part_manh_dist = np.abs(fmat1 - fmat2).sum(axis=0)
        self.part_pearson_c = pearson_coeff(cmat1, cmat2, axis=0)
        self.part_sp_score = part_sp_score(self.alint, all_parts[:2], self.residue_distmat)

        self.plotlist.extend(('manh', 'pearson', 'part_sp'))
        self.plotdata.update(manh=(self.x, self.part_manh_dist,),
                             pearson=(self.x, self.part_pearson_c,),
                             part_sp=(self.x, self.part_sp_score,))

        self.plot_funcs['al'].append((self.annotate_parts, {}))

        if self.tree is not None: # and 'part_pars' in self.plotlist:
            if hasattr(self.tree, 'children'):
                get_children = lambda tree, node: node.children
            else:
                get_children = lambda tree, node: node.clades
            self.part_pars_scores, part_branch_nbs = parsimony_score(self.alint, self.tree,
                                                seqlabels=self.seqlabels,
                                                minlength=(self.valid_range[1]+2),
                                                get_children=get_children,
                                                parts=self.parts)
            logger.debug('Finished partial parsimony scores')
            self.plotlist.append('part_pars')
            #self.plotdata['part_pars'] = tuple((arr for p_sc in self.part_pars_scores
            #                                    for arr in (self.x, p_sc)))
            self.plotdata['part_pars'] = (self.x, self.part_pars_scores)

            #self.plot_funcs['part_pars'][0][1].update(alpha=0.6)
            self.plot_funcs['part_pars'][2][1].update(labels=['outgroup']+
                                                      list(range(1, 1+len(self.parts))))
            self.plot_funcs['part_pars'][1] += ({'values': self.part_pars_scores},)


    def makefig(self, plotlist=None, figwidth=16, plotheight=4):

        plotlist = plotlist if plotlist is not None else self.plotlist
        nplots = len(plotlist)
        fig = plt.figure(figsize=(figwidth, plotheight*nplots))

        logger.debug('Nplots: %s, %s', nplots, plotlist)
        #fig, axes = plt.subplots(nplots, sharex=True,
        #                         figsize=(figwidth, 4*nplots), squeeze=False)
        
        if 'tree' in plotlist:
            rows = nplots - 1
            cols = 5  # tree plot spans 1 col, alignment spans 4.
            colspan = cols - 1
            gridspec_kw = {'width_ratios': [1, colspan]}
            # And hide alignment labels
            #self.plot_funcs['al'].append(('tick_params', {'label_left'}))
        else:
            rows = nplots
            cols = 1
            colspan = 1
            gridspec_kw = {}
        logger.debug("subplots: %d, %d" % (rows, cols))
        pos = (0,0)
        #ax = fig.add_subplot(rows, cols, pos, gridspec_kw=gridspec_kw)
        ax = plt.subplot2grid((rows,cols), pos, colspan=1, fig=fig)
        axes = [ax]

        for plot in plotlist:
            logger.debug('Plot=%r; pos=%s', plot, pos)
            logger.debug('len(axes): %d', len(axes))

            plotfuncname, plot_kwargs = self.plot_funcs[plot][0]
            try:
                plotfunc = getattr(ax, plotfuncname)
            except TypeError:
                plotfunc = plotfuncname
                # ax argument required:
                plot_kwargs.update(ax=ax)

            plotdata = self.plotdata[plot]

            plotfunc(*plotdata, **plot_kwargs)

            for i, (extrafunc, extra_kwargs) in enumerate(self.plot_funcs[plot][1:],
                                                          start=1):
                try:
                    if isinstance(extrafunc, str):
                        extrafunc = getattr(ax, extrafunc)
                    else:
                        extra_kwargs.update(ax=ax)
                        #extrafunc = lambda **kw: extrafuncname(ax, **kw)
                    extrafunc(**extra_kwargs)
                except BaseException as err:
                    err.args += ('extra func %d %r' % (i, extrafunc),)
                    raise

            plot_prop = self.plot_properties[plot]
            ax.set(**plot_prop)

            # Next ax:
            if len(axes) < nplots:
                if plot == 'tree':
                    #pos += 1
                    #ax = fig.add_subplot(rows, cols, pos, sharey=ax)
                    pos = (pos[0], pos[1]+1)
                    ax = plt.subplot2grid((rows,cols), pos, colspan=colspan, fig=fig, sharey=ax, autoscalex_on=False)
                else:
                    #pos += 2 if 'tree' in plotlist else 1
                    #ax = fig.add_subplot(rows, cols, pos, sharex=ax)
                    if plot == 'al' and rows>1 and cols>1:
                        # Plot the colorbar below the tree
                        cax = plt.subplot2grid((rows,cols),(pos[0]+1, 0), fig=fig,
                                               aspect=20)
                        al_colorbar(ax, cax, self.altype)#, orientation='horizontal')
                    pos = (pos[0]+1, pos[1])
                    logger.debug('rows,cols: (%d,%d); pos: %s; colspan: %d',
                                 rows, cols, pos, colspan)
                    ax = plt.subplot2grid((rows,cols), pos, colspan=colspan, fig=fig, sharex=ax, autoscalex_on=False)
                axes.append(ax)

        # On the last plot:
        if self.slstart:
            # Update the displayed xticklabels if start position > 0
            slstart = self.slstart // 3 if self.altype=='codon' else self.slstart
                
            ax.set_xticklabels([(slstart + int(xt)) for xt in \
                                ax.get_xticks()])

        fig.tight_layout()
        self.fig, self.axes = fig, axes


    def display(self, outfile=None, **kwargs):
        logger.debug('Displaying')
        if outfile is None:
            plt.show(block=True)
            #print(clock())
        else:
            self.fig.savefig(outfile, bbox_inches='tight', **kwargs)


def plot_al_conservation(infile, format=None, altype='nucl', allow_N=False,
         ungap=True, records=None, recordsfile=None, treefile=None, slice=None,
         topology_only=False, compare_parts=None,
         compare_only=False, plotlist=None, figwidth=16, plotheight=4):

    align_plot = AlignPlotter.fromfile(infile, format, altype, allow_N, ungap,
                                       slice, records, recordsfile, treefile,
                                       topology_only)
    if logger.getEffectiveLevel() <= logging.DEBUG:
        del align_plot.plot_funcs['tree'][1]  # tickparams(labelright=False)
    if not compare_only:
        align_plot.measure()
    if compare_parts:
        align_plot.comp_parts(compare_parts)
    if plotlist is not None:
        plotlist = plotlist.split(',')
        if set(plotlist).intersection(('pars', 'tree', 'part_pars')):
            assert treefile is not None, \
                    "A treefile must be given to compute parsimony score or plot a tree."
            # plotlist.remobe('pars') ## 'tree'
            if 'tree' not in plotlist and 'al' in plotlist:
                plotlist.insert(plotlist.index('al'), 'tree')
    align_plot.makefig(plotlist, figwidth, plotheight)
    return align_plot


def main():
    logging.basicConfig(format="%(levelname)s:%(funcName)s:%(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile', nargs='?', default=stdin,
                        type=argparse.FileType('r'))
    parser.add_argument('-o', '--outfile')
    parser.add_argument('-f', '--format', help='Force format usage.' \
                        ' Can be any format accepted by Bio.alignIO')
    parser.add_argument('-a', '--altype', choices=['nucl', 'aa', 'codon'],
                        help='Process per nucleotide position (instead of codon)')
    parser.add_argument('-N', '--allow-N', action='store_true',
                        help='Tolerate "N" nucleotides as NA value.')
    parser.add_argument('-g', '--no-ungap', dest='ungap', action='store_false',
                        help='Do not remove columns of gaps.')
    rec_g = parser.add_mutually_exclusive_group()
    rec_g.add_argument('-r', '--records',
                        help='select records (comma-sep list of 0-based integers)')
    rec_g.add_argument('-R', '--recordsfile',
                        help='select records from a file (one record name per line)')
    rec_g.add_argument('-t', '--treefile',
                        help='select records from a newick tree (reorder the sequences)')
    parser.add_argument('-T', '--topology-only', action='store_true',
                        help='Plot unit branch lengths instead of the real ones.')
    parser.add_argument('-s', '--slice',
                        help='select positions (start:end). '\
                             '0-based, end excluded, in number of '\
                             'codons/nucleotides.')
    parser.add_argument('-c', '--compare-parts',
                        help='Plot the per-column correlation between two ' \
                             'groups of data, e.g "(0,1);(2,3,4)". The last '\
                             'part is calculated automatically. slices and '\
                             'regex (matching sequence names) allowed')
    parser.add_argument('-C', '--compare-only', action='store_true',
                        help='Do not display global scores')
    parser.add_argument('-W', '--figwidth', default=16, type=float,
                        help='Figure width (inches) [%(default)s]')
    parser.add_argument('-H', '--plotheight', default=1.8, type=float,
                        help='Individual plot height (inches) [%(default)s]')
    parser.add_argument('-p', '--plotlist', default='al,gap,pars',
                        help='comma-sep list of plots. Valid values are: '\
                             'al,gap,pars,entropy,gap_entropy,sp,manh,pearson'
                             'part_sp,part_pars. [%(default)s]')
    parser.add_argument('-v', '--verbose', action='count', default=0)
    
    args = parser.parse_args()

    if args.verbose > 1:
        logger.setLevel(logging.DEBUG)
    elif args.verbose > 0:
        logger.setLevel(logging.INFO)
    delattr(args, 'verbose')

    outfile = args.outfile
    delattr(args, 'outfile')
    if not outfile:
        plt.switch_backend('TkAgg')

    plot_al_conservation(**vars(args)).display(outfile)


if __name__ == '__main__':
    main()
